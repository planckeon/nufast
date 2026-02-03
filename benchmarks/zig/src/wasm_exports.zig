//! WASM exports for NuFast
//!
//! Provides a C-ABI interface for calling NuFast from JavaScript/WebAssembly.
//! All functions use simple scalar types for easy interop.

const nufast = @import("nufast.zig");

// =============================================================================
// Exported Parameter Setters (build params in WASM memory)
// =============================================================================

/// Global vacuum parameters (set via exported functions)
var g_vacuum_params: nufast.VacuumParams = nufast.VacuumParams.default;

/// Global matter parameters
var g_matter_params: nufast.MatterParams = nufast.MatterParams.default;

/// Result buffer for probability matrix (9 f64s, row-major)
var g_result: [9]f64 = undefined;

// =============================================================================
// Parameter Setters
// =============================================================================

/// Set vacuum oscillation parameters
export fn set_vacuum_params(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    antineutrino: bool,
) void {
    g_vacuum_params = .{
        .s12sq = s12sq,
        .s13sq = s13sq,
        .s23sq = s23sq,
        .delta = delta,
        .Dmsq21 = Dmsq21,
        .Dmsq31 = Dmsq31,
        .antineutrino = antineutrino,
    };
}

/// Set to default NuFIT 5.2 parameters
export fn set_default_params() void {
    g_vacuum_params = nufast.VacuumParams.default;
    g_matter_params = nufast.MatterParams.default;
}

/// Set matter-specific parameters (vacuum params should be set first)
export fn set_matter_params(rho: f64, Ye: f64, n_newton: u8, antineutrino: bool) void {
    g_matter_params = .{
        .vacuum = g_vacuum_params,
        .rho = rho,
        .Ye = Ye,
        .n_newton = n_newton,
        .antineutrino = antineutrino,
    };
}

// =============================================================================
// Probability Calculations
// =============================================================================

/// Calculate vacuum oscillation probability and store in result buffer
/// Returns pointer to result buffer (9 f64s, row-major: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt])
export fn vacuum_probability(L: f64, E: f64) *const [9]f64 {
    const probs = nufast.vacuumProbability(g_vacuum_params, L, E);

    // Flatten to row-major
    g_result[0] = probs[0][0]; // Pee
    g_result[1] = probs[0][1]; // Pem
    g_result[2] = probs[0][2]; // Pet
    g_result[3] = probs[1][0]; // Pme
    g_result[4] = probs[1][1]; // Pmm
    g_result[5] = probs[1][2]; // Pmt
    g_result[6] = probs[2][0]; // Pte
    g_result[7] = probs[2][1]; // Ptm
    g_result[8] = probs[2][2]; // Ptt

    return &g_result;
}

/// Calculate matter oscillation probability and store in result buffer
export fn matter_probability(L: f64, E: f64) *const [9]f64 {
    const probs = nufast.matterProbability(g_matter_params, L, E);

    g_result[0] = probs[0][0];
    g_result[1] = probs[0][1];
    g_result[2] = probs[0][2];
    g_result[3] = probs[1][0];
    g_result[4] = probs[1][1];
    g_result[5] = probs[1][2];
    g_result[6] = probs[2][0];
    g_result[7] = probs[2][1];
    g_result[8] = probs[2][2];

    return &g_result;
}

// =============================================================================
// Direct Calculations (no global state, for advanced use)
// =============================================================================

/// Direct vacuum probability calculation - returns P(mu->e)
/// This is the most common use case for neutrino experiments
export fn vacuum_Pme(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    L: f64,
    E: f64,
) f64 {
    const params = nufast.VacuumParams{
        .s12sq = s12sq,
        .s13sq = s13sq,
        .s23sq = s23sq,
        .delta = delta,
        .Dmsq21 = Dmsq21,
        .Dmsq31 = Dmsq31,
        .antineutrino = false,
    };
    const probs = nufast.vacuumProbability(params, L, E);
    return probs[1][0]; // P(mu -> e)
}

/// Direct vacuum probability with default params - returns P(mu->e)
export fn vacuum_Pme_default(L: f64, E: f64) f64 {
    const probs = nufast.vacuumProbability(nufast.VacuumParams.default, L, E);
    return probs[1][0];
}

/// Direct matter probability with default params - returns P(mu->e)
export fn matter_Pme_default(L: f64, E: f64, rho: f64) f64 {
    const params = nufast.MatterParams{
        .vacuum = nufast.VacuumParams.default,
        .rho = rho,
        .Ye = 0.5,
        .n_newton = 0,
        .antineutrino = false,
    };
    const probs = nufast.matterProbability(params, L, E);
    return probs[1][0];
}

// =============================================================================
// Memory Access (for reading result buffer from JS)
// =============================================================================

/// Get pointer to result buffer
export fn get_result_ptr() *const [9]f64 {
    return &g_result;
}

/// Get a specific probability from the result buffer
/// idx: 0=Pee, 1=Pem, 2=Pet, 3=Pme, 4=Pmm, 5=Pmt, 6=Pte, 7=Ptm, 8=Ptt
export fn get_result(idx: usize) f64 {
    if (idx < 9) {
        return g_result[idx];
    }
    return 0.0;
}

// =============================================================================
// Batch Processing (for maximum throughput)
// =============================================================================

/// Maximum batch size (determines buffer allocation)
const MAX_BATCH_SIZE: usize = 1024;

/// Pre-computed vacuum batch for repeated calculations
var g_vacuum_batch: nufast.VacuumBatch = nufast.VacuumBatch.init(nufast.VacuumParams.default);

/// Pre-computed matter batch for repeated calculations
var g_matter_batch: nufast.MatterBatch = nufast.MatterBatch.init(nufast.MatterParams.default);

/// Energy input buffer (written by JS)
var g_energies: [MAX_BATCH_SIZE]f64 = undefined;

/// Output buffer for batch Pme results
var g_batch_output: [MAX_BATCH_SIZE]f64 = undefined;

/// Initialize vacuum batch from current global params
/// Call this after set_vacuum_params for best batch performance
export fn init_vacuum_batch() void {
    g_vacuum_batch = nufast.VacuumBatch.init(g_vacuum_params);
}

/// Initialize matter batch from current global params
/// Call this after set_matter_params for best batch performance
export fn init_matter_batch() void {
    g_matter_batch = nufast.MatterBatch.init(g_matter_params);
}

/// Get pointer to energy input buffer (JS writes energies here)
export fn get_energies_ptr() *[MAX_BATCH_SIZE]f64 {
    return &g_energies;
}

/// Get pointer to batch output buffer (JS reads results here)
export fn get_batch_output_ptr() *const [MAX_BATCH_SIZE]f64 {
    return &g_batch_output;
}

/// Get maximum batch size
export fn get_max_batch_size() usize {
    return MAX_BATCH_SIZE;
}

/// Calculate vacuum Pme for a batch of energies
/// Energies should be written to get_energies_ptr() first
/// Results are written to get_batch_output_ptr()
/// Returns pointer to output buffer
export fn vacuum_batch_Pme(L: f64, count: usize) *const [MAX_BATCH_SIZE]f64 {
    const n = @min(count, MAX_BATCH_SIZE);

    for (0..n) |i| {
        const probs = g_vacuum_batch.probabilityAt(L, g_energies[i]);
        g_batch_output[i] = probs[1][0]; // Pme
    }

    return &g_batch_output;
}

/// Calculate matter Pme for a batch of energies
/// Uses pre-computed matter batch (call init_matter_batch first)
export fn matter_batch_Pme(L: f64, count: usize) *const [MAX_BATCH_SIZE]f64 {
    const n = @min(count, MAX_BATCH_SIZE);

    for (0..n) |i| {
        const probs = g_matter_batch.probabilityAt(L, g_energies[i]);
        g_batch_output[i] = probs[1][0]; // Pme
    }

    return &g_batch_output;
}

/// Calculate full probability matrices for a batch of energies (vacuum)
/// Output: 9 * count f64s in row-major order per matrix
/// Buffer must be pre-allocated by caller at get_batch_matrix_ptr()
var g_batch_matrix: [MAX_BATCH_SIZE * 9]f64 = undefined;

export fn get_batch_matrix_ptr() *const [MAX_BATCH_SIZE * 9]f64 {
    return &g_batch_matrix;
}

export fn vacuum_batch_full(L: f64, count: usize) *const [MAX_BATCH_SIZE * 9]f64 {
    const n = @min(count, MAX_BATCH_SIZE);

    for (0..n) |i| {
        const probs = g_vacuum_batch.probabilityAt(L, g_energies[i]);
        const base = i * 9;
        g_batch_matrix[base + 0] = probs[0][0];
        g_batch_matrix[base + 1] = probs[0][1];
        g_batch_matrix[base + 2] = probs[0][2];
        g_batch_matrix[base + 3] = probs[1][0];
        g_batch_matrix[base + 4] = probs[1][1];
        g_batch_matrix[base + 5] = probs[1][2];
        g_batch_matrix[base + 6] = probs[2][0];
        g_batch_matrix[base + 7] = probs[2][1];
        g_batch_matrix[base + 8] = probs[2][2];
    }

    return &g_batch_matrix;
}

export fn matter_batch_full(L: f64, count: usize) *const [MAX_BATCH_SIZE * 9]f64 {
    const n = @min(count, MAX_BATCH_SIZE);

    for (0..n) |i| {
        const probs = g_matter_batch.probabilityAt(L, g_energies[i]);
        const base = i * 9;
        g_batch_matrix[base + 0] = probs[0][0];
        g_batch_matrix[base + 1] = probs[0][1];
        g_batch_matrix[base + 2] = probs[0][2];
        g_batch_matrix[base + 3] = probs[1][0];
        g_batch_matrix[base + 4] = probs[1][1];
        g_batch_matrix[base + 5] = probs[1][2];
        g_batch_matrix[base + 6] = probs[2][0];
        g_batch_matrix[base + 7] = probs[2][1];
        g_batch_matrix[base + 8] = probs[2][2];
    }

    return &g_batch_matrix;
}
