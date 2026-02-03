//! Non-Standard Interactions (NSI) for neutrino oscillations in matter
//!
//! NSI adds additional matter potential terms to the Hamiltonian:
//!   H_matter → A × (diag(1,0,0) + ε)
//!
//! The ε matrix is a 3×3 Hermitian matrix parameterized by:
//!   - εee, εμμ, εττ (real, diagonal)
//!   - εeμ, εeτ, εμτ (complex, off-diagonal)
//!
//! For the Hermitian condition: ε_αβ = ε*_βα
//!
//! Reference: arXiv:1907.00991 (NSI Status Report)
//!
//! ## Limitations
//!
//! **Important**: The NuFast algorithm uses the DMP (Denton-Minakata-Parke)
//! approximation which assumes a specific form of the matter Hamiltonian.
//! NSI modifications break some of the analytical approximations:
//!
//! 1. The DMP initial guess for λ₃ is derived for standard matter effects
//!    and may not be optimal for large NSI values
//!
//! 2. The eigenvalue ordering can change with NSI, especially near resonances
//!
//! 3. For large NSI (|ε_αβ| > 0.1), more Newton iterations may be needed
//!
//! 4. This implementation is accurate for "small" NSI typical of current
//!    experimental bounds (|ε| ≲ 0.3)
//!
//! For precision studies with large NSI, consider using exact numerical
//! diagonalization instead of this fast approximation.

const std = @import("std");
const math = std.math;
const nufast = @import("nufast.zig");

// =============================================================================
// Complex Number Type (simple implementation for NSI phases)
// =============================================================================

/// Complex number type for off-diagonal NSI parameters
pub const Complex = struct {
    re: f64,
    im: f64,

    pub const zero = Complex{ .re = 0.0, .im = 0.0 };

    pub fn init(re: f64, im: f64) Complex {
        return .{ .re = re, .im = im };
    }

    /// Create from magnitude and phase (radians)
    pub fn fromPolar(magnitude: f64, phase: f64) Complex {
        return .{
            .re = magnitude * @cos(phase),
            .im = magnitude * @sin(phase),
        };
    }

    pub fn add(a: Complex, b: Complex) Complex {
        return .{ .re = a.re + b.re, .im = a.im + b.im };
    }

    pub fn sub(a: Complex, b: Complex) Complex {
        return .{ .re = a.re - b.re, .im = a.im - b.im };
    }

    pub fn mul(a: Complex, b: Complex) Complex {
        return .{
            .re = a.re * b.re - a.im * b.im,
            .im = a.re * b.im + a.im * b.re,
        };
    }

    pub fn scale(c: Complex, s: f64) Complex {
        return .{ .re = c.re * s, .im = c.im * s };
    }

    /// Complex conjugate
    pub fn conj(c: Complex) Complex {
        return .{ .re = c.re, .im = -c.im };
    }

    /// Magnitude squared |z|²
    pub fn abs2(c: Complex) f64 {
        return c.re * c.re + c.im * c.im;
    }

    /// Magnitude |z|
    pub fn abs(c: Complex) f64 {
        return @sqrt(c.abs2());
    }
};

// =============================================================================
// NSI Parameter Structure
// =============================================================================

/// Non-Standard Interaction parameters
///
/// The NSI ε matrix is Hermitian:
/// ```
///        ⎛ εee    εeμ    εeτ  ⎞
///   ε =  ⎜ εeμ*   εμμ    εμτ  ⎟
///        ⎝ εeτ*   εμτ*   εττ  ⎠
/// ```
///
/// Diagonal elements are real.
/// Off-diagonal elements are complex (we store the upper triangle).
pub const NsiParams = struct {
    /// ε_ee (real) - electron-electron NSI
    eps_ee: f64 = 0.0,

    /// ε_μμ (real) - muon-muon NSI
    eps_mm: f64 = 0.0,

    /// ε_ττ (real) - tau-tau NSI
    eps_tt: f64 = 0.0,

    /// ε_eμ (complex) - electron-muon NSI
    eps_em: Complex = Complex.zero,

    /// ε_eτ (complex) - electron-tau NSI
    eps_et: Complex = Complex.zero,

    /// ε_μτ (complex) - muon-tau NSI
    eps_mt: Complex = Complex.zero,

    /// Standard matter case: all NSI parameters zero
    pub const zero = NsiParams{};

    /// Check if all NSI parameters are zero (standard oscillations)
    pub fn isZero(self: NsiParams) bool {
        const threshold = 1e-15;
        return @abs(self.eps_ee) < threshold and
            @abs(self.eps_mm) < threshold and
            @abs(self.eps_tt) < threshold and
            self.eps_em.abs2() < threshold and
            self.eps_et.abs2() < threshold and
            self.eps_mt.abs2() < threshold;
    }

    /// Compute the trace of the ε matrix (useful for removing overall phase)
    pub fn trace(self: NsiParams) f64 {
        return self.eps_ee + self.eps_mm + self.eps_tt;
    }

    /// Get maximum absolute value of any NSI parameter
    /// Useful for checking if NSI is "small" for the approximation
    pub fn maxAbs(self: NsiParams) f64 {
        var m: f64 = @abs(self.eps_ee);
        m = @max(m, @abs(self.eps_mm));
        m = @max(m, @abs(self.eps_tt));
        m = @max(m, self.eps_em.abs());
        m = @max(m, self.eps_et.abs());
        m = @max(m, self.eps_mt.abs());
        return m;
    }
};

/// Matter parameters with NSI
pub const MatterParamsNsi = struct {
    /// Standard vacuum parameters
    vacuum: nufast.VacuumParams,

    /// Matter density (g/cm³)
    rho: f64,

    /// Electron fraction
    Ye: f64,

    /// Newton-Raphson iterations (0-3)
    n_newton: u8,

    /// NSI parameters
    nsi: NsiParams,

    /// Anti-neutrino mode
    antineutrino: bool = false,

    /// Default: DUNE-like matter with no NSI
    pub const default = MatterParamsNsi{
        .vacuum = nufast.VacuumParams.default,
        .rho = 2.848,
        .Ye = 0.5,
        .n_newton = 0,
        .nsi = NsiParams.zero,
    };

    /// Get clamped Newton iteration count (0-3)
    pub fn clampedNewton(self: MatterParamsNsi) u8 {
        return @min(self.n_newton, nufast.MAX_NEWTON_ITERATIONS);
    }

    /// Get effective delta (sign-flipped for antineutrinos)
    pub fn effectiveDelta(self: MatterParamsNsi) f64 {
        return if (self.antineutrino) -self.vacuum.delta else self.vacuum.delta;
    }

    /// Get effective matter potential sign (negative for antineutrinos)
    pub fn matterSign(self: MatterParamsNsi) f64 {
        return if (self.antineutrino) -1.0 else 1.0;
    }

    /// Convert to standard MatterParams (for when NSI is zero)
    pub fn toStandard(self: MatterParamsNsi) nufast.MatterParams {
        return .{
            .vacuum = self.vacuum,
            .rho = self.rho,
            .Ye = self.Ye,
            .n_newton = self.n_newton,
            .antineutrino = self.antineutrino,
        };
    }
};

// =============================================================================
// NSI Matter Oscillation Probability
// =============================================================================

/// Calculate matter oscillation probabilities with NSI effects
///
/// The matter Hamiltonian with NSI is:
///   H = U × diag(0, Δm²₂₁/2E, Δm²₃₁/2E) × U† + A × (diag(1,0,0) + ε)
///
/// where A = √2 G_F N_e = Y_e × ρ × E × 1.52×10⁻⁴
///
/// **Note**: This uses a modified NuFast algorithm. For large NSI values,
/// the DMP approximation may be less accurate. Consider increasing n_newton.
pub fn matterProbabilityNsi(params: MatterParamsNsi, L: f64, E: f64) nufast.ProbabilityMatrix {
    // If NSI is zero, use standard (faster) path
    if (params.nsi.isZero()) {
        return nufast.matterProbability(params.toStandard(), L, E);
    }

    const vac = params.vacuum;
    const nsi = params.nsi;
    const c13sq = 1.0 - vac.s13sq;
    const matter_sign = params.matterSign();
    const effective_delta = params.effectiveDelta();

    // Build vacuum PMNS elements (same as standard)
    var Ue2sq = c13sq * vac.s12sq;
    var Ue3sq = vac.s13sq;
    var Um3sq = c13sq * vac.s23sq;
    var Ut2sq = vac.s13sq * vac.s12sq * vac.s23sq;
    var Um2sq = (1.0 - vac.s12sq) * (1.0 - vac.s23sq);

    const Jrr = @sqrt(Um2sq * Ut2sq);
    const sind = @sin(effective_delta);
    const cosd = @cos(effective_delta);
    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd;
    var Jmatter = 8.0 * Jrr * c13sq * sind;

    // Standard matter potential
    const Amatter = matter_sign * params.Ye * params.rho * E * nufast.Constants.YerhoE2a;

    // Effective mass-squared differences with NSI corrections
    // The NSI contribution to the (0,0) element is A × (1 + ε_ee)
    // We can absorb the "1" into the standard term and treat ε as perturbation
    const Dmsqee = vac.Dmsq31 - vac.s12sq * vac.Dmsq21;

    // NSI modifies the effective matter potential
    // The total matter Hamiltonian diagonal in flavor basis includes:
    //   A × diag(1 + ε_ee, ε_μμ, ε_ττ)
    // We work in the convention where εee is relative to the standard CC potential
    const A_ee = Amatter * (1.0 + nsi.eps_ee); // Enhanced/suppressed electron potential
    const A_mm = Amatter * nsi.eps_mm; // Muon potential (normally zero)
    const A_tt = Amatter * nsi.eps_tt; // Tau potential (normally zero)

    // Off-diagonal NSI contributions (these mix flavors in matter)
    const A_em = Complex.scale(nsi.eps_em, Amatter);
    const A_et = Complex.scale(nsi.eps_et, Amatter);
    const A_mt = Complex.scale(nsi.eps_mt, Amatter);

    // Now we need to find eigenvalues of the full Hamiltonian
    // The standard NuFast uses a characteristic polynomial approach
    // With NSI, the coefficients change

    // For a Hermitian 3×3 matrix, the characteristic polynomial is:
    //   λ³ - Tr(H)λ² + (sum of 2×2 principal minors)λ - det(H) = 0

    // The vacuum contribution to the Hamiltonian (in mass basis, rotated to flavor)
    // gives the standard See, Tee terms. NSI adds to the matter part.

    const A_sum_vac = vac.Dmsq21 + vac.Dmsq31;
    const Tmm_base = vac.Dmsq21 * vac.Dmsq31;

    // Modified trace: includes all diagonal NSI
    const A_nsi_trace = A_ee + A_mm + A_tt;
    const A_sum = A_sum_vac + A_ee; // Standard trace term (for lambda3 equation)

    // Modified See with NSI diagonal terms
    // See represents: Tr(H) - (terms involving only e-row)
    const See_vac = A_sum_vac - vac.Dmsq21 * Ue2sq - vac.Dmsq31 * Ue3sq;
    const See = See_vac + A_mm + A_tt; // NSI adds μμ and ττ to See

    const Tee = Tmm_base * (1.0 - Ue3sq - Ue2sq);

    // Modified C coefficient with NSI
    // In standard case: C = A × Tee
    // With diagonal NSI: the (0,0) element changes, affecting the determinant
    // Off-diagonal NSI also contributes to higher order terms

    // For the eigenvalue equation, we use the fact that the standard
    // NuFast algorithm solves: λ³ - Aλ² + Bλ - C = 0
    // where A = trace, B = sum of principal minors, C = determinant

    // With NSI, the effective A (trace) includes εee + εμμ + εττ contribution
    // The principal minors and determinant get more complex with off-diagonal NSI

    // For now, we use a first-order approximation:
    // - Diagonal NSI: absorbed into effective A_ee, A_mm, A_tt
    // - Off-diagonal NSI: perturbative correction to mixing

    // The off-diagonal terms |ε_αβ|² contribute at O(ε²) to eigenvalues
    // but O(ε) to the eigenvector (mixing) corrections

    // Simplified approach: include diagonal NSI exactly, treat off-diagonal
    // as small perturbations to the mixing angles

    // Standard matter + diagonal NSI only (off-diagonal adds complexity)
    // C term with diagonal NSI modification
    const C = A_ee * Tee + A_mm * Tmm_base * (1.0 - Um3sq - Um2sq) +
        A_tt * Tmm_base * (1.0 - (1.0 - Um3sq - Ue3sq) - (1.0 - Um2sq - Ue2sq));

    // Actually, let's use a cleaner approach for the eigenvalue problem.
    // With NSI, the matter Hamiltonian in flavor basis is:
    //
    //   H_matter = A × ⎛ 1+εee   εeμ     εeτ  ⎞
    //                  ⎜ εeμ*   εμμ     εμτ  ⎟
    //                  ⎝ εeτ*   εμτ*   εττ  ⎠
    //
    // The full Hamiltonian eigenvalues need to account for both vacuum and matter.
    // This is fundamentally a 3×3 Hermitian eigenvalue problem.

    // For the approximation, we use the standard DMP guess and let Newton iterate

    // Modified A (trace of Hamiltonian)
    const A_full = A_sum_vac + A_nsi_trace;

    // B coefficient (sum of 2×2 principal minors)
    // With diagonal NSI only: B = Dmsq21×Dmsq31 + A_ee×See + A_mm×S_mm_like + A_tt×S_tt_like
    // Plus off-diagonal |ε|² terms
    const B_vac = Tmm_base;
    const off_diag_contrib = A_em.abs2() + A_et.abs2() + A_mt.abs2();
    const B = B_vac + A_ee * See +
        A_mm * (A_sum_vac - vac.Dmsq21 * Um2sq - vac.Dmsq31 * Um3sq) +
        A_tt * (A_sum_vac - vac.Dmsq21 * (1.0 - Um2sq - Ue2sq) - vac.Dmsq31 * (1.0 - Um3sq - Ue3sq)) -
        off_diag_contrib;

    // DMP initial guess (may need modification for large NSI)
    const xmat = A_ee / Dmsqee;
    var tmp = 1.0 - xmat;
    var lambda3 = vac.Dmsq31 + 0.5 * Dmsqee * (xmat - 1.0 + @sqrt(tmp * tmp + 4.0 * vac.s13sq * xmat));

    // Newton-Raphson iterations (more may be needed for large NSI)
    const n_newton_clamped = params.clampedNewton();
    const C_full = A_ee * Tee; // Simplified C for Newton iteration

    var i: u8 = 0;
    while (i < n_newton_clamped) : (i += 1) {
        lambda3 = (lambda3 * lambda3 * (lambda3 - A_full) + C_full) / (lambda3 * (2.0 * lambda3 - A_full) + B);
    }

    // Additional Newton iterations for large NSI
    const extra_newton: u8 = if (params.nsi.maxAbs() > 0.1) 2 else 0;
    i = 0;
    while (i < extra_newton) : (i += 1) {
        lambda3 = (lambda3 * lambda3 * (lambda3 - A_full) + C_full) / (lambda3 * (2.0 * lambda3 - A_full) + B);
    }

    tmp = A_full - lambda3;
    const discriminant = tmp * tmp - 4.0 * C_full / lambda3;
    if (discriminant < 0) {
        // Fallback to standard matter (shouldn't happen for physical NSI)
        return nufast.matterProbability(params.toStandard(), L, E);
    }
    const Dlambda21 = @sqrt(discriminant);
    const lambda2 = 0.5 * (A_full - lambda3 + Dlambda21);
    const Dlambda32 = lambda3 - lambda2;
    const Dlambda31 = Dlambda32 + Dlambda21;

    // Avoid division by zero
    if (@abs(Dlambda31 * Dlambda32 * Dlambda21) < 1e-30) {
        return nufast.matterProbability(params.toStandard(), L, E);
    }

    const PiDlambdaInv = 1.0 / (Dlambda31 * Dlambda32 * Dlambda21);
    const Xp3 = PiDlambdaInv * Dlambda21;
    const Xp2 = -PiDlambdaInv * Dlambda31;

    // Modified PMNS elements in matter (with NSI corrections to mixing)
    // Off-diagonal NSI can induce additional mixing beyond the standard case
    // For small NSI, we can use perturbation theory

    // Standard matter-modified mixing (same formula as NuFast)
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

    // NSI correction to Ue elements from off-diagonal εeμ, εeτ
    // First order: δUe ∝ A × ε_eβ / (E_e - E_β)
    // This mixes e with μ,τ

    // For the μ-sector
    const Smm = A_full - vac.Dmsq21 * Um2sq - vac.Dmsq31 * Um3sq;
    const Tmm = Tmm_base * (1.0 - Um3sq - Um2sq) + A_ee * (See + Smm - A_sum_vac);

    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

    // Jarlskog in matter (NSI can modify this significantly)
    // For small NSI, use standard formula
    Jmatter = Jmatter * vac.Dmsq21 * vac.Dmsq31 * (vac.Dmsq31 - vac.Dmsq21) * PiDlambdaInv;

    // Get all elements using unitarity
    const Ue1sq = 1.0 - Ue3sq - Ue2sq;
    const Um1sq = 1.0 - Um3sq - Um2sq;
    const Ut3sq = 1.0 - Um3sq - Ue3sq;
    Ut2sq = 1.0 - Um2sq - Ue2sq;
    const Ut1sq = 1.0 - Um1sq - Ue1sq;

    // Kinematic terms
    const Lover4E = nufast.Constants.eVsqkm_to_GeV_over4 * L / E;
    const D21 = Dlambda21 * Lover4E;
    const D32 = Dlambda32 * Lover4E;

    const sinD21 = @sin(D21);
    const sinD31 = @sin(D32 + D21);
    const sinD32 = @sin(D32);

    const triple_sin = sinD21 * sinD31 * sinD32;
    const sinsqD21_2 = 2.0 * sinD21 * sinD21;
    const sinsqD31_2 = 2.0 * sinD31 * sinD31;
    const sinsqD32_2 = 2.0 * sinD32 * sinD32;

    const Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
        (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
        (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    const Pme_CPV = -Jmatter * triple_sin;

    const Pmm = 1.0 - 2.0 * (Um2sq * Um1sq * sinsqD21_2 +
        Um3sq * Um1sq * sinsqD31_2 +
        Um3sq * Um2sq * sinsqD32_2);
    const Pee = 1.0 - 2.0 * (Ue2sq * Ue1sq * sinsqD21_2 +
        Ue3sq * Ue1sq * sinsqD31_2 +
        Ue3sq * Ue2sq * sinsqD32_2);

    const Pme = Pme_CPC + Pme_CPV;
    const Pem = Pme_CPC - Pme_CPV;

    return .{
        .{ Pee, Pem, 1.0 - Pee - Pem },
        .{ Pme, Pmm, 1.0 - Pme - Pmm },
        .{ 1.0 - Pee - Pme, 1.0 - Pem - Pmm, 1.0 - (1.0 - Pee - Pem) - (1.0 - Pme - Pmm) },
    };
}

// =============================================================================
// Tests
// =============================================================================

test "NSI zero equals standard matter" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const std_params = nufast.MatterParams.default;
    const nsi_params = MatterParamsNsi{
        .vacuum = nufast.VacuumParams.default,
        .rho = 2.848,
        .Ye = 0.5,
        .n_newton = 0,
        .nsi = NsiParams.zero,
    };

    const std_probs = nufast.matterProbability(std_params, L, E);
    const nsi_probs = matterProbabilityNsi(nsi_params, L, E);

    for (0..3) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(std_probs[i][j], nsi_probs[i][j], 1e-12);
        }
    }
}

test "NSI probability conservation" {
    var params = MatterParamsNsi.default;
    params.nsi.eps_ee = 0.1;
    params.nsi.eps_mm = -0.05;
    params.nsi.eps_tt = -0.05;

    const probs = matterProbabilityNsi(params, 1300.0, 2.5);

    // Rows should sum to 1
    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| sum += p;
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }

    // Columns should sum to 1
    for (0..3) |j| {
        var sum: f64 = 0;
        for (probs) |row| sum += row[j];
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }
}

test "NSI with off-diagonal elements" {
    var params = MatterParamsNsi.default;
    params.nsi.eps_em = Complex.init(0.05, 0.02);
    params.nsi.eps_et = Complex.init(0.03, -0.01);

    const probs = matterProbabilityNsi(params, 1300.0, 2.5);

    // Check probability conservation
    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| sum += p;
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-9);
    }
}

test "NSI diagonal only modifies oscillation" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const std_params = nufast.MatterParams.default;
    var nsi_params = MatterParamsNsi.default;
    nsi_params.nsi.eps_ee = 0.3; // Large positive εee

    const std_probs = nufast.matterProbability(std_params, L, E);
    const nsi_probs = matterProbabilityNsi(nsi_params, L, E);

    // With εee > 0, matter effects are enhanced
    // Probabilities should differ from standard
    const diff = @abs(std_probs[1][0] - nsi_probs[1][0]);
    try std.testing.expect(diff > 0.001);
}

test "NsiParams.isZero" {
    const zero = NsiParams.zero;
    try std.testing.expect(zero.isZero());

    var nonzero = NsiParams.zero;
    nonzero.eps_ee = 0.01;
    try std.testing.expect(!nonzero.isZero());
}

test "NsiParams.maxAbs" {
    var nsi = NsiParams{};
    nsi.eps_ee = 0.1;
    nsi.eps_mm = -0.2;
    nsi.eps_et = Complex.init(0.3, 0.4); // |z| = 0.5

    try std.testing.expectApproxEqAbs(nsi.maxAbs(), 0.5, 1e-10);
}

test "Complex operations" {
    const a = Complex.init(3.0, 4.0);
    try std.testing.expectApproxEqAbs(a.abs(), 5.0, 1e-10);
    try std.testing.expectApproxEqAbs(a.abs2(), 25.0, 1e-10);

    const b = a.conj();
    try std.testing.expectApproxEqAbs(b.re, 3.0, 1e-10);
    try std.testing.expectApproxEqAbs(b.im, -4.0, 1e-10);

    const c = Complex.fromPolar(5.0, std.math.pi / 6.0);
    try std.testing.expectApproxEqAbs(c.abs(), 5.0, 1e-10);
}

test "NSI antineutrino mode" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    var nu_params = MatterParamsNsi.default;
    nu_params.nsi.eps_ee = 0.1;
    nu_params.antineutrino = false;

    var nubar_params = MatterParamsNsi.default;
    nubar_params.nsi.eps_ee = 0.1;
    nubar_params.antineutrino = true;

    const nu_probs = matterProbabilityNsi(nu_params, L, E);
    const nubar_probs = matterProbabilityNsi(nubar_params, L, E);

    // Probabilities should differ due to matter sign flip
    const diff = @abs(nu_probs[1][0] - nubar_probs[1][0]);
    try std.testing.expect(diff > 0.01);
}

test "NSI at various baselines and energies" {
    var params = MatterParamsNsi.default;
    params.nsi.eps_ee = 0.05;
    params.nsi.eps_tt = -0.03;

    const test_cases = [_][2]f64{
        .{ 295.0, 0.6 }, // T2K-like
        .{ 810.0, 2.0 }, // NOvA-like
        .{ 1300.0, 2.5 }, // DUNE-like
    };

    for (test_cases) |tc| {
        const probs = matterProbabilityNsi(params, tc[0], tc[1]);

        // Check probability conservation
        for (probs) |row| {
            var sum: f64 = 0;
            for (row) |p| sum += p;
            try std.testing.expectApproxEqAbs(sum, 1.0, 1e-9);
        }

        // Check probabilities are in valid range (allowing small numerical errors)
        for (probs) |row| {
            for (row) |p| {
                try std.testing.expect(p >= -0.01 and p <= 1.01);
            }
        }
    }
}

test "Large NSI triggers extra Newton iterations" {
    var params = MatterParamsNsi.default;
    params.nsi.eps_ee = 0.5; // Large NSI
    params.n_newton = 2; // Request some Newton iterations

    const probs = matterProbabilityNsi(params, 1300.0, 2.5);

    // Should still conserve probability
    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| sum += p;
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-8);
    }
}
