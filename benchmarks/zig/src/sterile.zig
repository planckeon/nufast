//! Sterile Neutrino Oscillations (3+1 Model)
//!
//! Extension of NuFast for 4-flavor neutrino oscillations including one sterile neutrino.
//! This module provides EXACT vacuum oscillation probabilities for the 3+1 scenario.
//!
//! ## Physics Background
//!
//! The 3+1 model extends the standard 3-flavor framework with one sterile neutrino (ν_s).
//! The sterile neutrino does not participate in weak interactions but mixes with active
//! neutrinos through a 4×4 PMNS matrix.
//!
//! Key parameters:
//! - θ14, θ24, θ34: mixing angles between active and sterile states
//! - δ14, δ24: additional CP-violating phases (δ34 can be absorbed)
//! - Δm²41: mass-squared difference ~1 eV² (typical for short-baseline anomalies)
//!
//! ## Limitations
//!
//! **The NuFast approximation does NOT apply to 4-flavor oscillations.**
//!
//! The NuFast algorithm (arXiv:2405.02400) is specifically optimized for 3-flavor
//! matter oscillations using:
//! 1. DMP approximation for the largest eigenvalue
//! 2. Newton-Raphson refinement based on 3×3 characteristic polynomial
//!
//! For 4-flavor, the characteristic polynomial is 4th-degree (quartic), requiring:
//! - Exact diagonalization (O(n³) per evaluation), or
//! - Perturbative expansions (assuming hierarchy of mass splittings)
//!
//! This module implements EXACT VACUUM oscillations using the general formula.
//! Matter effects for 4-flavor would require numerical diagonalization.
//!
//! ## References
//!
//! - arXiv:1303.3011 - "Sterile Neutrino Oscillations: The Global Picture"
//! - arXiv:1609.01211 - "Anomalies in neutrino experiments"
//! - PDG 2023 - Review of Neutrino Mixing

const std = @import("std");
const math = std.math;
const nufast = @import("nufast.zig");

// =============================================================================
// Physical Constants
// =============================================================================

/// Conversion factor: eV² km to GeV / 4 (same as 3-flavor)
const eVsqkm_to_GeV_over4: f64 = 1e-9 / 1.97327e-7 * 1e3 / 4.0;

// =============================================================================
// Parameter Structures
// =============================================================================

/// Oscillation parameters for 3+1 sterile neutrino vacuum calculations.
///
/// Extends the standard 3-flavor parameters with:
/// - Three additional mixing angles (θ14, θ24, θ34)
/// - Two additional CP phases (δ14, δ24)
/// - One additional mass-squared difference (Δm²41)
pub const SterileParams = struct {
    // Standard 3-flavor parameters
    /// sin²θ12 (solar)
    s12sq: f64,
    /// sin²θ13 (reactor)
    s13sq: f64,
    /// sin²θ23 (atmospheric)
    s23sq: f64,
    /// CP phase δ13 (radians), often just called δ
    delta13: f64,
    /// Δm²21 (eV²) - solar mass splitting
    Dmsq21: f64,
    /// Δm²31 (eV²) - atmospheric mass splitting
    Dmsq31: f64,

    // Sterile sector parameters
    /// sin²θ14 (e-s mixing)
    s14sq: f64,
    /// sin²θ24 (μ-s mixing)
    s24sq: f64,
    /// sin²θ34 (τ-s mixing)
    s34sq: f64,
    /// CP phase δ14 (radians)
    delta14: f64,
    /// CP phase δ24 (radians)
    delta24: f64,
    /// Δm²41 (eV²) - sterile mass splitting (~1 eV² for LSND/MiniBooNE)
    Dmsq41: f64,

    /// Anti-neutrino mode (flips sign of all CP phases)
    antineutrino: bool = false,

    /// Default NuFIT 5.2 + typical sterile parameters
    /// Sterile mixing set to values motivated by short-baseline anomalies
    pub const default = SterileParams{
        // Standard 3-flavor (NuFIT 5.2 normal ordering)
        .s12sq = 0.307,
        .s13sq = 0.0220,
        .s23sq = 0.546,
        .delta13 = -0.7 * math.pi,
        .Dmsq21 = 7.53e-5,
        .Dmsq31 = 2.453e-3,
        // Sterile sector (motivated by LSND/MiniBooNE hints)
        // Note: These are illustrative; actual values are controversial
        .s14sq = 0.025, // |Ue4|² ~ 0.025
        .s24sq = 0.015, // |Uμ4|² ~ 0.015
        .s34sq = 0.0, // Usually unconstrained
        .delta14 = 0.0,
        .delta24 = 0.0,
        .Dmsq41 = 1.0, // ~1 eV² for short-baseline
    };

    /// Create from standard VacuumParams with additional sterile parameters
    pub fn fromVacuumParams(
        vacuum: nufast.VacuumParams,
        s14sq: f64,
        s24sq: f64,
        s34sq: f64,
        Dmsq41: f64,
    ) SterileParams {
        return .{
            .s12sq = vacuum.s12sq,
            .s13sq = vacuum.s13sq,
            .s23sq = vacuum.s23sq,
            .delta13 = vacuum.delta,
            .Dmsq21 = vacuum.Dmsq21,
            .Dmsq31 = vacuum.Dmsq31,
            .s14sq = s14sq,
            .s24sq = s24sq,
            .s34sq = s34sq,
            .delta14 = 0.0,
            .delta24 = 0.0,
            .Dmsq41 = Dmsq41,
            .antineutrino = vacuum.antineutrino,
        };
    }

    /// Convert to 3-flavor VacuumParams (ignoring sterile sector)
    pub fn toVacuumParams(self: SterileParams) nufast.VacuumParams {
        return .{
            .s12sq = self.s12sq,
            .s13sq = self.s13sq,
            .s23sq = self.s23sq,
            .delta = self.delta13,
            .Dmsq21 = self.Dmsq21,
            .Dmsq31 = self.Dmsq31,
            .antineutrino = self.antineutrino,
        };
    }

    /// Create parameters with zero sterile mixing (reduces to 3-flavor)
    pub fn noSterileMixing() SterileParams {
        var p = default;
        p.s14sq = 0.0;
        p.s24sq = 0.0;
        p.s34sq = 0.0;
        return p;
    }
};

/// 4×4 probability matrix [from][to]
/// Indices: 0=e, 1=μ, 2=τ, 3=s
pub const ProbabilityMatrix4 = [4][4]f64;

/// Complex number for PMNS matrix calculations
const Complex = std.math.Complex(f64);

/// 4×4 complex PMNS matrix
pub const PMNS4 = [4][4]Complex;

// =============================================================================
// PMNS Matrix Construction
// =============================================================================

/// Construct the 4×4 PMNS matrix from sterile parameters.
///
/// The 4×4 PMNS matrix is parameterized as a product of rotation matrices.
/// Following the convention commonly used in sterile neutrino literature:
///
/// U = R34 × W24 × W14 × R23 × W13 × R12
///
/// where Wij = Rij with CP phase (e^{-iδ} in the ij element).
///
/// The rotation matrices act on mass eigenstates to produce flavor eigenstates:
/// |να⟩ = Σᵢ Uαᵢ |νᵢ⟩
///
/// Indices: 0=ν1, 1=ν2, 2=ν3, 3=ν4 for mass eigenstates
///          0=νe, 1=νμ, 2=ντ, 3=νs for flavor eigenstates
pub fn constructPMNS4(params: SterileParams) PMNS4 {
    const antinu_sign: f64 = if (params.antineutrino) -1.0 else 1.0;

    // Extract sin and cos of mixing angles
    const s12 = @sqrt(params.s12sq);
    const c12 = @sqrt(1.0 - params.s12sq);
    const s13 = @sqrt(params.s13sq);
    const c13 = @sqrt(1.0 - params.s13sq);
    const s23 = @sqrt(params.s23sq);
    const c23 = @sqrt(1.0 - params.s23sq);
    const s14 = @sqrt(params.s14sq);
    const c14 = @sqrt(1.0 - params.s14sq);
    const s24 = @sqrt(params.s24sq);
    const c24 = @sqrt(1.0 - params.s24sq);
    const s34 = @sqrt(params.s34sq);
    const c34 = @sqrt(1.0 - params.s34sq);

    // CP phases (sign-flipped for antineutrinos)
    const d13 = antinu_sign * params.delta13;
    const d14 = antinu_sign * params.delta14;
    const d24 = antinu_sign * params.delta24;

    // Helper to create real complex number
    const real = struct {
        fn f(x: f64) Complex {
            return Complex.init(x, 0.0);
        }
    }.f;

    // Build rotation matrices
    // R12: rotation in 1-2 plane (solar angle θ12)
    var R12 = identityPMNS4();
    R12[0][0] = real(c12);
    R12[0][1] = real(s12);
    R12[1][0] = real(-s12);
    R12[1][1] = real(c12);

    // W13: rotation in 1-3 plane with CP phase δ13 (reactor angle θ13)
    // The phase convention: U_e3 = s13 * e^{-iδ}
    var W13 = identityPMNS4();
    W13[0][0] = real(c13);
    W13[0][2] = Complex.init(s13 * @cos(-d13), s13 * @sin(-d13)); // s13 * e^{-iδ13}
    W13[2][0] = Complex.init(-s13 * @cos(-d13), s13 * @sin(-d13)); // -s13 * e^{iδ13}
    W13[2][2] = real(c13);

    // R23: rotation in 2-3 plane (atmospheric angle θ23)
    var R23 = identityPMNS4();
    R23[1][1] = real(c23);
    R23[1][2] = real(s23);
    R23[2][1] = real(-s23);
    R23[2][2] = real(c23);

    // W14: rotation in 1-4 plane with CP phase δ14
    var W14 = identityPMNS4();
    W14[0][0] = real(c14);
    W14[0][3] = Complex.init(s14 * @cos(-d14), s14 * @sin(-d14));
    W14[3][0] = Complex.init(-s14 * @cos(-d14), s14 * @sin(-d14));
    W14[3][3] = real(c14);

    // W24: rotation in 2-4 plane with CP phase δ24
    var W24 = identityPMNS4();
    W24[1][1] = real(c24);
    W24[1][3] = Complex.init(s24 * @cos(-d24), s24 * @sin(-d24));
    W24[3][1] = Complex.init(-s24 * @cos(-d24), s24 * @sin(-d24));
    W24[3][3] = real(c24);

    // R34: rotation in 3-4 plane (no additional CP phase needed)
    var R34 = identityPMNS4();
    R34[2][2] = real(c34);
    R34[2][3] = real(s34);
    R34[3][2] = real(-s34);
    R34[3][3] = real(c34);

    // Construct U = R34 × W24 × W14 × R23 × W13 × R12
    // Matrix multiplication from right to left (rightmost applied first)
    var U = R12;
    U = matmul4(W13, U);
    U = matmul4(R23, U);
    U = matmul4(W14, U);
    U = matmul4(W24, U);
    U = matmul4(R34, U);

    return U;
}

/// Create a 4×4 identity matrix
fn identityPMNS4() PMNS4 {
    var I: PMNS4 = undefined;
    for (0..4) |i| {
        for (0..4) |j| {
            I[i][j] = if (i == j) Complex.init(1.0, 0.0) else Complex.init(0.0, 0.0);
        }
    }
    return I;
}

/// 4×4 complex matrix multiplication
fn matmul4(A: PMNS4, B: PMNS4) PMNS4 {
    var C: PMNS4 = undefined;
    for (0..4) |i| {
        for (0..4) |j| {
            C[i][j] = Complex.init(0.0, 0.0);
            for (0..4) |k| {
                C[i][j] = C[i][j].add(A[i][k].mul(B[k][j]));
            }
        }
    }
    return C;
}

// =============================================================================
// Vacuum Oscillation Probabilities (Exact 4-Flavor)
// =============================================================================

/// Calculate exact 4-flavor vacuum oscillation probability.
///
/// Uses the general formula:
/// P(α→β) = δ_αβ - 4 Σ_{j>k} Re(U*αj Uβj Uαk U*βk) sin²(Δjk)
///                + 2 Σ_{j>k} Im(U*αj Uβj Uαk U*βk) sin(2Δjk)
///
/// where Δjk = Δm²jk × L / (4E)
///
/// Arguments:
///   params: Sterile oscillation parameters
///   L: baseline distance (km)
///   E: neutrino energy (GeV)
///
/// Returns: 4×4 probability matrix P[α][β] = P(να → νβ)
pub fn sterileProbabilityVacuum(params: SterileParams, L: f64, E: f64) ProbabilityMatrix4 {
    const U = constructPMNS4(params);

    // Kinematic factor
    const Lover4E = eVsqkm_to_GeV_over4 * L / E;

    // Mass-squared differences relative to m1
    // Δm²21 = m2² - m1²
    // Δm²31 = m3² - m1²
    // Δm²41 = m4² - m1²
    // Δm²32 = Δm²31 - Δm²21
    // Δm²42 = Δm²41 - Δm²21
    // Δm²43 = Δm²41 - Δm²31
    const Dmsq21 = params.Dmsq21;
    const Dmsq31 = params.Dmsq31;
    const Dmsq41 = params.Dmsq41;
    const Dmsq32 = Dmsq31 - Dmsq21;
    const Dmsq42 = Dmsq41 - Dmsq21;
    const Dmsq43 = Dmsq41 - Dmsq31;

    // Phase arguments Δjk = Δm²jk × L / (4E)
    const D21 = Dmsq21 * Lover4E;
    const D31 = Dmsq31 * Lover4E;
    const D41 = Dmsq41 * Lover4E;
    const D32 = Dmsq32 * Lover4E;
    const D42 = Dmsq42 * Lover4E;
    const D43 = Dmsq43 * Lover4E;

    var probs: ProbabilityMatrix4 = undefined;

    // Calculate each transition probability
    for (0..4) |alpha| {
        for (0..4) |beta| {
            // Start with δ_αβ
            var P: f64 = if (alpha == beta) 1.0 else 0.0;

            // Sum over pairs (j > k) for j, k ∈ {1, 2, 3, 4} (using 0-indexing: 0, 1, 2, 3)
            // Pairs: (1,0), (2,0), (2,1), (3,0), (3,1), (3,2)

            // Helper function to compute the contribution from pair (j, k)
            const addContribution = struct {
                fn f(
                    P_ptr: *f64,
                    U_mat: PMNS4,
                    a: usize,
                    b: usize,
                    j: usize,
                    k: usize,
                    Djk: f64,
                ) void {
                    // U*αj Uβj Uαk U*βk
                    const Ua_j_conj = U_mat[a][j].conjugate();
                    const Ub_j = U_mat[b][j];
                    const Ua_k = U_mat[a][k];
                    const Ub_k_conj = U_mat[b][k].conjugate();

                    const term = Ua_j_conj.mul(Ub_j).mul(Ua_k).mul(Ub_k_conj);

                    const sin_Djk = @sin(Djk);
                    const sin_2Djk = @sin(2.0 * Djk);

                    // -4 Re(...) sin²(Δjk) + 2 Im(...) sin(2Δjk)
                    P_ptr.* -= 4.0 * term.re * sin_Djk * sin_Djk;
                    P_ptr.* += 2.0 * term.im * sin_2Djk;
                }
            }.f;

            // (j=1, k=0) corresponds to Δm²21
            addContribution(&P, U, alpha, beta, 1, 0, D21);
            // (j=2, k=0) corresponds to Δm²31
            addContribution(&P, U, alpha, beta, 2, 0, D31);
            // (j=2, k=1) corresponds to Δm²32
            addContribution(&P, U, alpha, beta, 2, 1, D32);
            // (j=3, k=0) corresponds to Δm²41
            addContribution(&P, U, alpha, beta, 3, 0, D41);
            // (j=3, k=1) corresponds to Δm²42
            addContribution(&P, U, alpha, beta, 3, 1, D42);
            // (j=3, k=2) corresponds to Δm²43
            addContribution(&P, U, alpha, beta, 3, 2, D43);

            probs[alpha][beta] = P;
        }
    }

    return probs;
}

/// Get just the active-flavor 3×3 submatrix from 4-flavor probabilities.
/// Useful for comparing with 3-flavor NuFast results.
pub fn getActiveProbabilities(probs4: ProbabilityMatrix4) nufast.ProbabilityMatrix {
    var probs3: nufast.ProbabilityMatrix = undefined;
    for (0..3) |i| {
        for (0..3) |j| {
            probs3[i][j] = probs4[i][j];
        }
    }
    return probs3;
}

/// Calculate the probability of active-to-sterile transitions.
/// P(νe → νs), P(νμ → νs), P(ντ → νs)
pub fn getSterileAppearance(probs4: ProbabilityMatrix4) [3]f64 {
    return .{ probs4[0][3], probs4[1][3], probs4[2][3] };
}

/// Calculate the effective disappearance due to sterile mixing.
/// This is the total probability lost to the sterile sector.
pub fn getTotalSterileDisappearance(probs4: ProbabilityMatrix4) [3]f64 {
    // For each active flavor, the sterile disappearance is P(α→s)
    // which equals 1 - Σ_β P(α→β) for active β
    return .{
        1.0 - (probs4[0][0] + probs4[0][1] + probs4[0][2]),
        1.0 - (probs4[1][0] + probs4[1][1] + probs4[1][2]),
        1.0 - (probs4[2][0] + probs4[2][1] + probs4[2][2]),
    };
}

// =============================================================================
// Convenience Functions for Short-Baseline Experiments
// =============================================================================

/// Short-baseline (SBL) effective 2-flavor approximation for νe disappearance.
///
/// For L/E << 1/Δm²21, the oscillations are dominated by Δm²41, giving:
/// P(νe → νe) ≈ 1 - sin²(2θ14) sin²(Δm²41 L / 4E)
///
/// This is the formula often used for reactor anomaly analyses.
pub fn sblElectronDisappearance(params: SterileParams, L: f64, E: f64) f64 {
    const Lover4E = eVsqkm_to_GeV_over4 * L / E;
    const D41 = params.Dmsq41 * Lover4E;
    const sin_2theta14_sq = 4.0 * params.s14sq * (1.0 - params.s14sq);
    const sin_D41_sq = @sin(D41) * @sin(D41);
    return 1.0 - sin_2theta14_sq * sin_D41_sq;
}

/// Short-baseline effective 2-flavor approximation for νμ → νe appearance.
///
/// For L/E << 1/Δm²21, the appearance probability is approximately:
/// P(νμ → νe) ≈ sin²(2θμe) sin²(Δm²41 L / 4E)
///
/// where sin²(2θμe) = 4|Ue4|²|Uμ4|² = 4 s14² s24²
///
/// This is the formula often used for LSND/MiniBooNE analyses.
pub fn sblMuonToElectronAppearance(params: SterileParams, L: f64, E: f64) f64 {
    const Lover4E = eVsqkm_to_GeV_over4 * L / E;
    const D41 = params.Dmsq41 * Lover4E;
    // Approximate effective mixing angle for appearance
    // sin²(2θμe) ≈ 4 |Ue4|² |Uμ4|² for small mixing
    const sin_2theta_mue_sq = 4.0 * params.s14sq * params.s24sq;
    const sin_D41_sq = @sin(D41) * @sin(D41);
    return sin_2theta_mue_sq * sin_D41_sq;
}

/// Short-baseline effective 2-flavor approximation for νμ disappearance.
///
/// P(νμ → νμ) ≈ 1 - sin²(2θ24) sin²(Δm²41 L / 4E)
pub fn sblMuonDisappearance(params: SterileParams, L: f64, E: f64) f64 {
    const Lover4E = eVsqkm_to_GeV_over4 * L / E;
    const D41 = params.Dmsq41 * Lover4E;
    const sin_2theta24_sq = 4.0 * params.s24sq * (1.0 - params.s24sq);
    const sin_D41_sq = @sin(D41) * @sin(D41);
    return 1.0 - sin_2theta24_sq * sin_D41_sq;
}

// =============================================================================
// Pre-computed Batch Structure
// =============================================================================

/// Pre-computed structure for batch 4-flavor vacuum calculations.
/// Stores the PMNS matrix to avoid recomputation.
pub const SterileBatch = struct {
    U: PMNS4,
    Dmsq21: f64,
    Dmsq31: f64,
    Dmsq41: f64,

    pub fn init(params: SterileParams) SterileBatch {
        return .{
            .U = constructPMNS4(params),
            .Dmsq21 = params.Dmsq21,
            .Dmsq31 = params.Dmsq31,
            .Dmsq41 = params.Dmsq41,
        };
    }

    pub fn probabilityAt(self: SterileBatch, L: f64, E: f64) ProbabilityMatrix4 {
        const Lover4E = eVsqkm_to_GeV_over4 * L / E;

        const D21 = self.Dmsq21 * Lover4E;
        const D31 = self.Dmsq31 * Lover4E;
        const D41 = self.Dmsq41 * Lover4E;
        const D32 = D31 - D21;
        const D42 = D41 - D21;
        const D43 = D41 - D31;

        var probs: ProbabilityMatrix4 = undefined;

        for (0..4) |alpha| {
            for (0..4) |beta| {
                var P: f64 = if (alpha == beta) 1.0 else 0.0;

                inline for (.{ .{ 1, 0, D21 }, .{ 2, 0, D31 }, .{ 2, 1, D32 }, .{ 3, 0, D41 }, .{ 3, 1, D42 }, .{ 3, 2, D43 } }) |pair| {
                    const j = pair[0];
                    const k = pair[1];
                    const Djk = pair[2];

                    const Ua_j_conj = self.U[alpha][j].conjugate();
                    const Ub_j = self.U[beta][j];
                    const Ua_k = self.U[alpha][k];
                    const Ub_k_conj = self.U[beta][k].conjugate();

                    const term = Ua_j_conj.mul(Ub_j).mul(Ua_k).mul(Ub_k_conj);

                    const sin_Djk = @sin(Djk);
                    const sin_2Djk = @sin(2.0 * Djk);

                    P -= 4.0 * term.re * sin_Djk * sin_Djk;
                    P += 2.0 * term.im * sin_2Djk;
                }

                probs[alpha][beta] = P;
            }
        }

        return probs;
    }
};

// =============================================================================
// Tests
// =============================================================================

test "4-flavor reduces to 3-flavor when sterile mixing is zero" {
    // When s14 = s24 = s34 = 0, the 4-flavor should match 3-flavor
    var sterile_params = SterileParams.noSterileMixing();
    sterile_params.s14sq = 0.0;
    sterile_params.s24sq = 0.0;
    sterile_params.s34sq = 0.0;

    const vacuum_params = sterile_params.toVacuumParams();

    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const probs4 = sterileProbabilityVacuum(sterile_params, L, E);
    const probs3 = nufast.vacuumProbability(vacuum_params, L, E);

    // Active sector should match
    for (0..3) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(probs3[i][j], probs4[i][j], 1e-10);
        }
    }

    // Sterile column and row should be zero/one appropriately
    // P(α→s) = 0 for all active α
    for (0..3) |alpha| {
        try std.testing.expectApproxEqAbs(probs4[alpha][3], 0.0, 1e-10);
    }
    // P(s→s) = 1
    try std.testing.expectApproxEqAbs(probs4[3][3], 1.0, 1e-10);
}

test "4-flavor probability conservation (unitarity)" {
    const params = SterileParams.default;
    const probs = sterileProbabilityVacuum(params, 1300.0, 2.5);

    // Rows should sum to 1
    for (0..4) |alpha| {
        var sum: f64 = 0;
        for (0..4) |beta| sum += probs[alpha][beta];
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }

    // Columns should sum to 1
    for (0..4) |beta| {
        var sum: f64 = 0;
        for (0..4) |alpha| sum += probs[alpha][beta];
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }
}

test "4-flavor batch matches single calculation" {
    const params = SterileParams.default;
    const batch = SterileBatch.init(params);

    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const single = sterileProbabilityVacuum(params, L, E);
    const batched = batch.probabilityAt(L, E);

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(single[i][j], batched[i][j], 1e-14);
        }
    }
}

test "short-baseline approximations" {
    var params = SterileParams.default;
    params.Dmsq41 = 1.0; // 1 eV²

    // At very short baseline, should see oscillation pattern
    const L: f64 = 0.5; // 500 m = 0.5 km
    const E: f64 = 0.030; // 30 MeV

    const Pee_approx = sblElectronDisappearance(params, L, E);
    const Pme_approx = sblMuonToElectronAppearance(params, L, E);
    const Pmm_approx = sblMuonDisappearance(params, L, E);

    // Probabilities should be in valid range
    try std.testing.expect(Pee_approx >= 0.0 and Pee_approx <= 1.0);
    try std.testing.expect(Pme_approx >= 0.0 and Pme_approx <= 1.0);
    try std.testing.expect(Pmm_approx >= 0.0 and Pmm_approx <= 1.0);

    // At the peak of oscillation, should see significant effect
    // L/E = 0.5/0.03 ≈ 16.7 km/GeV
    // Δ = 1.27 × 1.0 × 16.7 ≈ 21 rad (multiple oscillations)
}

test "PMNS matrix is unitary" {
    const params = SterileParams.default;
    const U = constructPMNS4(params);

    // Check U × U† = I
    for (0..4) |i| {
        for (0..4) |j| {
            var sum = Complex.init(0.0, 0.0);
            for (0..4) |k| {
                sum = sum.add(U[i][k].mul(U[j][k].conjugate()));
            }
            const expected: f64 = if (i == j) 1.0 else 0.0;
            try std.testing.expectApproxEqAbs(sum.re, expected, 1e-10);
            try std.testing.expectApproxEqAbs(sum.im, 0.0, 1e-10);
        }
    }
}

test "sterile appearance is significant with non-zero mixing" {
    const params = SterileParams.default;

    // Short baseline where sterile oscillations are maximal
    // L/E ≈ π/(2.54 × Δm²41) for first maximum
    const Dmsq41 = params.Dmsq41;
    const L_over_E_max = math.pi / (2.54 * Dmsq41);
    const E: f64 = 0.030; // 30 MeV
    const L = L_over_E_max * E;

    const probs = sterileProbabilityVacuum(params, L, E);

    // Should see non-zero sterile appearance
    const sterile_app = getSterileAppearance(probs);
    try std.testing.expect(sterile_app[0] > 0.001); // P(νe → νs) > 0
    try std.testing.expect(sterile_app[1] > 0.001); // P(νμ → νs) > 0
}

test "antineutrino mode changes CP-violating terms" {
    var params_nu = SterileParams.default;
    params_nu.delta13 = 0.5 * math.pi;
    params_nu.delta14 = 0.3 * math.pi;
    params_nu.antineutrino = false;

    var params_nubar = params_nu;
    params_nubar.antineutrino = true;

    const L: f64 = 500.0;
    const E: f64 = 1.0;

    const probs_nu = sterileProbabilityVacuum(params_nu, L, E);
    const probs_nubar = sterileProbabilityVacuum(params_nubar, L, E);

    // Survival probabilities should be the same (CPT)
    try std.testing.expectApproxEqAbs(probs_nu[0][0], probs_nubar[0][0], 1e-10);
    try std.testing.expectApproxEqAbs(probs_nu[1][1], probs_nubar[1][1], 1e-10);

    // Transition probabilities: P(α→β) for ν should equal P(β→α) for ν̄
    try std.testing.expectApproxEqAbs(probs_nu[1][0], probs_nubar[0][1], 1e-10);
    try std.testing.expectApproxEqAbs(probs_nu[0][1], probs_nubar[1][0], 1e-10);
}

test "near-zero baseline gives identity" {
    const params = SterileParams.default;
    const probs = sterileProbabilityVacuum(params, 0.0001, 2.5);

    // At L→0, P should approach identity matrix
    // Due to Δm²41 ~ 1 eV², we need very small L to get close to identity
    for (0..4) |i| {
        for (0..4) |j| {
            const expected: f64 = if (i == j) 1.0 else 0.0;
            try std.testing.expectApproxEqAbs(probs[i][j], expected, 1e-6);
        }
    }
}
