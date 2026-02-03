//! C-ABI exports for NuFast
//!
//! Provides a C-compatible interface for calling NuFast from Python/ctypes.
//! All functions use simple scalar types and pointers for easy FFI.

const nufast = @import("nufast.zig");

// =============================================================================
// Global State (for stateful API)
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
export fn nufast_set_vacuum_params(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    antineutrino: bool,
) callconv(.c) void {
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
export fn nufast_set_default_params() callconv(.c) void {
    g_vacuum_params = nufast.VacuumParams.default;
    g_matter_params = nufast.MatterParams.default;
}

/// Set matter-specific parameters (vacuum params should be set first)
export fn nufast_set_matter_params(rho: f64, Ye: f64, n_newton: u8, antineutrino: bool) callconv(.c) void {
    g_matter_params = .{
        .vacuum = g_vacuum_params,
        .rho = rho,
        .Ye = Ye,
        .n_newton = n_newton,
        .antineutrino = antineutrino,
    };
}

// =============================================================================
// Probability Calculations (stateful)
// =============================================================================

/// Calculate vacuum oscillation probability and store in result buffer
/// Returns pointer to result buffer (9 f64s, row-major: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt])
export fn nufast_vacuum_probability(L: f64, E: f64) callconv(.c) [*]const f64 {
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
export fn nufast_matter_probability(L: f64, E: f64) callconv(.c) [*]const f64 {
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
// Direct Calculations (stateless, for advanced use)
// =============================================================================

/// Direct vacuum probability calculation - writes full 3x3 matrix to output buffer
/// Output buffer must be at least 9 f64s (row-major)
export fn nufast_vacuum_prob_direct(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    L: f64,
    E: f64,
    antineutrino: bool,
    out: [*]f64,
) callconv(.c) void {
    const params = nufast.VacuumParams{
        .s12sq = s12sq,
        .s13sq = s13sq,
        .s23sq = s23sq,
        .delta = delta,
        .Dmsq21 = Dmsq21,
        .Dmsq31 = Dmsq31,
        .antineutrino = antineutrino,
    };
    const probs = nufast.vacuumProbability(params, L, E);

    out[0] = probs[0][0];
    out[1] = probs[0][1];
    out[2] = probs[0][2];
    out[3] = probs[1][0];
    out[4] = probs[1][1];
    out[5] = probs[1][2];
    out[6] = probs[2][0];
    out[7] = probs[2][1];
    out[8] = probs[2][2];
}

/// Direct matter probability calculation - writes full 3x3 matrix to output buffer
export fn nufast_matter_prob_direct(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    L: f64,
    E: f64,
    rho: f64,
    Ye: f64,
    n_newton: u8,
    antineutrino: bool,
    out: [*]f64,
) callconv(.c) void {
    const vacuum = nufast.VacuumParams{
        .s12sq = s12sq,
        .s13sq = s13sq,
        .s23sq = s23sq,
        .delta = delta,
        .Dmsq21 = Dmsq21,
        .Dmsq31 = Dmsq31,
        .antineutrino = antineutrino,
    };
    const params = nufast.MatterParams{
        .vacuum = vacuum,
        .rho = rho,
        .Ye = Ye,
        .n_newton = n_newton,
        .antineutrino = antineutrino,
    };
    const probs = nufast.matterProbability(params, L, E);

    out[0] = probs[0][0];
    out[1] = probs[0][1];
    out[2] = probs[0][2];
    out[3] = probs[1][0];
    out[4] = probs[1][1];
    out[5] = probs[1][2];
    out[6] = probs[2][0];
    out[7] = probs[2][1];
    out[8] = probs[2][2];
}

/// Direct vacuum probability - returns P(mu->e) (most common use case)
export fn nufast_vacuum_Pme(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    L: f64,
    E: f64,
) callconv(.c) f64 {
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
export fn nufast_vacuum_Pme_default(L: f64, E: f64) callconv(.c) f64 {
    const probs = nufast.vacuumProbability(nufast.VacuumParams.default, L, E);
    return probs[1][0];
}

/// Direct matter probability with default params - returns P(mu->e)
export fn nufast_matter_Pme_default(L: f64, E: f64, rho: f64) callconv(.c) f64 {
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
// Batch Processing
// =============================================================================

/// Maximum batch size for stateful batch API
const MAX_BATCH_SIZE: usize = 1024;

/// Pre-computed vacuum batch for repeated calculations
var g_vacuum_batch: nufast.VacuumBatch = nufast.VacuumBatch.init(nufast.VacuumParams.default);

/// Pre-computed matter batch for repeated calculations
var g_matter_batch: nufast.MatterBatch = nufast.MatterBatch.init(nufast.MatterParams.default);

/// Energy input buffer (written by caller)
var g_energies: [MAX_BATCH_SIZE]f64 = undefined;

/// Output buffer for batch Pme results
var g_batch_output: [MAX_BATCH_SIZE]f64 = undefined;

/// Initialize vacuum batch from current global params
export fn nufast_init_vacuum_batch() callconv(.c) void {
    g_vacuum_batch = nufast.VacuumBatch.init(g_vacuum_params);
}

/// Initialize matter batch from current global params
export fn nufast_init_matter_batch() callconv(.c) void {
    g_matter_batch = nufast.MatterBatch.init(g_matter_params);
}

/// Get pointer to energy input buffer
export fn nufast_get_energies_ptr() callconv(.c) [*]f64 {
    return &g_energies;
}

/// Get pointer to batch output buffer
export fn nufast_get_batch_output_ptr() callconv(.c) [*]const f64 {
    return &g_batch_output;
}

/// Get maximum batch size
export fn nufast_get_max_batch_size() callconv(.c) usize {
    return MAX_BATCH_SIZE;
}

/// Calculate vacuum Pme for a batch of energies
export fn nufast_vacuum_batch_Pme(L: f64, count: usize) callconv(.c) [*]const f64 {
    const n = @min(count, MAX_BATCH_SIZE);

    for (0..n) |i| {
        const probs = g_vacuum_batch.probabilityAt(L, g_energies[i]);
        g_batch_output[i] = probs[1][0]; // Pme
    }

    return &g_batch_output;
}

/// Calculate matter Pme for a batch of energies
export fn nufast_matter_batch_Pme(L: f64, count: usize) callconv(.c) [*]const f64 {
    const n = @min(count, MAX_BATCH_SIZE);

    for (0..n) |i| {
        const probs = g_matter_batch.probabilityAt(L, g_energies[i]);
        g_batch_output[i] = probs[1][0]; // Pme
    }

    return &g_batch_output;
}

// =============================================================================
// Stateless Batch Processing (caller provides buffers)
// =============================================================================

/// Calculate vacuum Pme for a batch of energies (stateless)
/// Caller provides input (energies) and output (results) buffers
export fn nufast_vacuum_batch_Pme_direct(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    L: f64,
    energies: [*]const f64,
    count: usize,
    antineutrino: bool,
    out: [*]f64,
) callconv(.c) void {
    const params = nufast.VacuumParams{
        .s12sq = s12sq,
        .s13sq = s13sq,
        .s23sq = s23sq,
        .delta = delta,
        .Dmsq21 = Dmsq21,
        .Dmsq31 = Dmsq31,
        .antineutrino = antineutrino,
    };
    const batch = nufast.VacuumBatch.init(params);

    for (0..count) |i| {
        const probs = batch.probabilityAt(L, energies[i]);
        out[i] = probs[1][0];
    }
}

/// Calculate matter Pme for a batch of energies (stateless)
export fn nufast_matter_batch_Pme_direct(
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    delta: f64,
    Dmsq21: f64,
    Dmsq31: f64,
    L: f64,
    energies: [*]const f64,
    count: usize,
    rho: f64,
    Ye: f64,
    n_newton: u8,
    antineutrino: bool,
    out: [*]f64,
) callconv(.c) void {
    const vacuum = nufast.VacuumParams{
        .s12sq = s12sq,
        .s13sq = s13sq,
        .s23sq = s23sq,
        .delta = delta,
        .Dmsq21 = Dmsq21,
        .Dmsq31 = Dmsq31,
        .antineutrino = antineutrino,
    };
    const params = nufast.MatterParams{
        .vacuum = vacuum,
        .rho = rho,
        .Ye = Ye,
        .n_newton = n_newton,
        .antineutrino = antineutrino,
    };
    const batch = nufast.MatterBatch.init(params);

    for (0..count) |i| {
        const probs = batch.probabilityAt(L, energies[i]);
        out[i] = probs[1][0];
    }
}

// =============================================================================
// Default Parameters Access (for Python to query defaults)
// =============================================================================

/// Get default s12sq
export fn nufast_default_s12sq() callconv(.c) f64 {
    return nufast.VacuumParams.default.s12sq;
}

/// Get default s13sq
export fn nufast_default_s13sq() callconv(.c) f64 {
    return nufast.VacuumParams.default.s13sq;
}

/// Get default s23sq
export fn nufast_default_s23sq() callconv(.c) f64 {
    return nufast.VacuumParams.default.s23sq;
}

/// Get default delta
export fn nufast_default_delta() callconv(.c) f64 {
    return nufast.VacuumParams.default.delta;
}

/// Get default Dmsq21
export fn nufast_default_Dmsq21() callconv(.c) f64 {
    return nufast.VacuumParams.default.Dmsq21;
}

/// Get default Dmsq31
export fn nufast_default_Dmsq31() callconv(.c) f64 {
    return nufast.VacuumParams.default.Dmsq31;
}

/// Get default rho
export fn nufast_default_rho() callconv(.c) f64 {
    return nufast.MatterParams.default.rho;
}

/// Get default Ye
export fn nufast_default_Ye() callconv(.c) f64 {
    return nufast.MatterParams.default.Ye;
}
