//! NuFast - Fast three-flavor neutrino oscillation probabilities
//!
//! Zig implementation of the NuFast algorithm by Denton & Parke.
//! Provides high-performance oscillation probability calculations
//! for vacuum and constant-density matter.
//!
//! Features:
//! - SIMD-optimized batch calculations (vacuum and matter)
//! - Compile-time generic f32/f64 modes (f32 = 2× SIMD lanes)
//! - Pre-computed batch structures for repeated calculations
//! - Anti-neutrino mode (flips sign of matter potential and δCP)
//!
//! Reference: arXiv:2405.02400

const std = @import("std");
const math = std.math;

// =============================================================================
// Physical Constants (compile-time generic)
// =============================================================================

/// Physical constants parameterized by float type
pub fn ConstantsFor(comptime T: type) type {
    return struct {
        /// Conversion factor: Y_e * rho * E to matter potential
        pub const YerhoE2a: T = 1.52e-4;
        /// Conversion: eV² km to GeV / 4
        pub const eVsqkm_to_GeV_over4: T = 1e-9 / 1.97327e-7 * 1e3 / 4.0;
    };
}

/// f64 constants (default)
pub const Constants = ConstantsFor(f64);

// =============================================================================
// Parameter Structures
// =============================================================================

/// Oscillation parameters for vacuum calculations
pub const VacuumParams = struct {
    /// sin²θ₁₂
    s12sq: f64,
    /// sin²θ₁₃
    s13sq: f64,
    /// sin²θ₂₃
    s23sq: f64,
    /// CP phase δ (radians)
    delta: f64,
    /// Δm²₂₁ (eV²)
    Dmsq21: f64,
    /// Δm²₃₁ (eV²)
    Dmsq31: f64,
    /// Anti-neutrino mode (flips sign of δCP)
    antineutrino: bool = false,

    /// Default NuFIT 5.2 parameters
    pub const default = VacuumParams{
        .s12sq = 0.307,
        .s13sq = 0.0220,
        .s23sq = 0.546,
        .delta = -0.7 * math.pi,
        .Dmsq21 = 7.53e-5,
        .Dmsq31 = 2.453e-3,
    };

    /// Get effective delta (sign-flipped for antineutrinos)
    pub fn effectiveDelta(self: VacuumParams) f64 {
        return if (self.antineutrino) -self.delta else self.delta;
    }
};

/// Maximum allowed Newton iterations (for numerical stability)
pub const MAX_NEWTON_ITERATIONS: u8 = 3;

/// Oscillation parameters for matter calculations
pub const MatterParams = struct {
    /// Vacuum parameters
    vacuum: VacuumParams,
    /// Matter density (g/cm³)
    rho: f64,
    /// Electron fraction
    Ye: f64,
    /// Newton-Raphson iterations (0-3, clamped to MAX_NEWTON_ITERATIONS)
    n_newton: u8,
    /// Anti-neutrino mode (flips sign of matter potential and δCP)
    antineutrino: bool = false,

    /// Default: DUNE-like matter
    pub const default = MatterParams{
        .vacuum = VacuumParams.default,
        .rho = 2.848,
        .Ye = 0.5,
        .n_newton = 0,
    };
    
    /// Get clamped Newton iteration count (0-3)
    pub fn clampedNewton(self: MatterParams) u8 {
        return @min(self.n_newton, MAX_NEWTON_ITERATIONS);
    }

    /// Get effective delta (sign-flipped for antineutrinos)
    pub fn effectiveDelta(self: MatterParams) f64 {
        return if (self.antineutrino) -self.vacuum.delta else self.vacuum.delta;
    }

    /// Get effective matter potential sign (negative for antineutrinos)
    pub fn matterSign(self: MatterParams) f64 {
        return if (self.antineutrino) -1.0 else 1.0;
    }
};

/// 3×3 probability matrix [from][to]
pub const ProbabilityMatrix = [3][3]f64;

// =============================================================================
// Pre-computed Batch Structures
// =============================================================================

/// Pre-computed mixing matrix elements for batch vacuum calculations
pub const VacuumBatch = struct {
    // PMNS matrix elements squared
    Ue1sq: f64,
    Ue2sq: f64,
    Ue3sq: f64,
    Um1sq: f64,
    Um2sq: f64,
    Um3sq: f64,
    Ut1sq: f64,
    Ut2sq: f64,
    Ut3sq: f64,
    // Jarlskog invariant for vacuum
    Jvac: f64,
    // Mass squared differences
    Dmsq21: f64,
    Dmsq31: f64,
    // Anti-neutrino mode
    antineutrino: bool,

    /// Create batch calculator from vacuum parameters
    pub fn init(params: VacuumParams) VacuumBatch {
        const c13sq = 1.0 - params.s13sq;
        const effective_delta = params.effectiveDelta();

        const Ue3sq = params.s13sq;
        const Ue2sq = c13sq * params.s12sq;
        const Um3sq = c13sq * params.s23sq;
        var Ut2sq = params.s13sq * params.s12sq * params.s23sq;
        var Um2sq = (1.0 - params.s12sq) * (1.0 - params.s23sq);

        const Jrr = @sqrt(Um2sq * Ut2sq);
        const sind = @sin(effective_delta);
        const cosd = @cos(effective_delta);
        Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd;
        const Jvac = 8.0 * Jrr * c13sq * sind;

        const Ue1sq = 1.0 - Ue3sq - Ue2sq;
        const Um1sq = 1.0 - Um3sq - Um2sq;
        const Ut3sq = 1.0 - Um3sq - Ue3sq;
        Ut2sq = 1.0 - Um2sq - Ue2sq;
        const Ut1sq = 1.0 - Um1sq - Ue1sq;

        return .{
            .Ue1sq = Ue1sq,
            .Ue2sq = Ue2sq,
            .Ue3sq = Ue3sq,
            .Um1sq = Um1sq,
            .Um2sq = Um2sq,
            .Um3sq = Um3sq,
            .Ut1sq = Ut1sq,
            .Ut2sq = Ut2sq,
            .Ut3sq = Ut3sq,
            .Jvac = Jvac,
            .Dmsq21 = params.Dmsq21,
            .Dmsq31 = params.Dmsq31,
            .antineutrino = params.antineutrino,
        };
    }

    /// Calculate probability at given L and E
    pub fn probabilityAt(self: VacuumBatch, L: f64, E: f64) ProbabilityMatrix {
        const Lover4E = Constants.eVsqkm_to_GeV_over4 * L / E;
        const D21 = self.Dmsq21 * Lover4E;
        const D31 = self.Dmsq31 * Lover4E;

        const sinD21 = @sin(D21);
        const sinD31 = @sin(D31);
        const sinD32 = @sin(D31 - D21);

        const triple_sin = sinD21 * sinD31 * sinD32;
        const sinsqD21_2 = 2.0 * sinD21 * sinD21;
        const sinsqD31_2 = 2.0 * sinD31 * sinD31;
        const sinsqD32_2 = 2.0 * sinD32 * sinD32;

        const Pme_CPC = (self.Ut3sq - self.Um2sq * self.Ue1sq - self.Um1sq * self.Ue2sq) * sinsqD21_2 +
            (self.Ut2sq - self.Um3sq * self.Ue1sq - self.Um1sq * self.Ue3sq) * sinsqD31_2 +
            (self.Ut1sq - self.Um3sq * self.Ue2sq - self.Um2sq * self.Ue3sq) * sinsqD32_2;
        const Pme_CPV = -self.Jvac * triple_sin;

        const Pmm = 1.0 - 2.0 * (self.Um2sq * self.Um1sq * sinsqD21_2 +
            self.Um3sq * self.Um1sq * sinsqD31_2 +
            self.Um3sq * self.Um2sq * sinsqD32_2);
        const Pee = 1.0 - 2.0 * (self.Ue2sq * self.Ue1sq * sinsqD21_2 +
            self.Ue3sq * self.Ue1sq * sinsqD31_2 +
            self.Ue3sq * self.Ue2sq * sinsqD32_2);

        const Pme = Pme_CPC + Pme_CPV;
        const Pem = Pme_CPC - Pme_CPV;

        return .{
            .{ Pee, Pem, 1.0 - Pee - Pem },
            .{ Pme, Pmm, 1.0 - Pme - Pmm },
            .{ 1.0 - Pee - Pme, 1.0 - Pem - Pmm, 1.0 - (1.0 - Pee - Pem) - (1.0 - Pme - Pmm) },
        };
    }
};

/// Pre-computed structure for batch matter calculations at constant density.
/// 
/// Pre-computes all terms that don't depend on energy:
/// - Vacuum mixing elements and Jarlskog
/// - Matter density terms
/// - Vacuum mass-squared products
///
/// This is ~30-40% faster than repeated matterProbability calls.
pub const MatterBatch = struct {
    // Original vacuum mixing (before matter modification)
    s12sq: f64,
    s13sq: f64,
    s23sq: f64,
    c13sq: f64,

    // Pre-computed vacuum PMNS-like terms
    Ue2sq_vac: f64,
    Ue3sq_vac: f64,
    Um3sq_vac: f64,
    Um2sq_vac: f64,
    Jrr: f64,
    sind: f64,
    cosd: f64,

    // Mass-squared differences
    Dmsq21: f64,
    Dmsq31: f64,
    Dmsqee: f64,

    // Pre-computed products
    Tmm_base: f64, // Dmsq21 * Dmsq31
    A_sum: f64, // Dmsq21 + Dmsq31

    // Matter parameters
    Ye_rho: f64, // Ye * rho (pre-multiplied)
    n_newton: u8,

    // Anti-neutrino sign for matter potential
    matter_sign: f64,

    /// Initialize from matter parameters (pre-compute everything possible)
    pub fn init(params: MatterParams) MatterBatch {
        const vac = params.vacuum;
        const c13sq = 1.0 - vac.s13sq;
        const matter_sign: f64 = params.matterSign();
        const effective_delta = params.effectiveDelta();

        const Ue2sq_vac = c13sq * vac.s12sq;
        const Ue3sq_vac = vac.s13sq;
        const Um3sq_vac = c13sq * vac.s23sq;
        const Ut2sq_temp = vac.s13sq * vac.s12sq * vac.s23sq;
        const Um2sq_temp = (1.0 - vac.s12sq) * (1.0 - vac.s23sq);

        const Jrr = @sqrt(Um2sq_temp * Ut2sq_temp);
        const sind = @sin(effective_delta);
        const cosd = @cos(effective_delta);
        const Um2sq_vac = Um2sq_temp + Ut2sq_temp - 2.0 * Jrr * cosd;

        const Dmsqee = vac.Dmsq31 - vac.s12sq * vac.Dmsq21;
        const Tmm_base = vac.Dmsq21 * vac.Dmsq31;
        const A_sum = vac.Dmsq21 + vac.Dmsq31;

        return .{
            .s12sq = vac.s12sq,
            .s13sq = vac.s13sq,
            .s23sq = vac.s23sq,
            .c13sq = c13sq,
            .Ue2sq_vac = Ue2sq_vac,
            .Ue3sq_vac = Ue3sq_vac,
            .Um3sq_vac = Um3sq_vac,
            .Um2sq_vac = Um2sq_vac,
            .Jrr = Jrr,
            .sind = sind,
            .cosd = cosd,
            .Dmsq21 = vac.Dmsq21,
            .Dmsq31 = vac.Dmsq31,
            .Dmsqee = Dmsqee,
            .Tmm_base = Tmm_base,
            .A_sum = A_sum,
            .Ye_rho = params.Ye * params.rho,
            .n_newton = params.clampedNewton(),
            .matter_sign = matter_sign,
        };
    }

    /// Calculate probability at given L and E
    pub fn probabilityAt(self: MatterBatch, L: f64, E: f64) ProbabilityMatrix {
        // Energy-dependent matter potential
        const Amatter = self.matter_sign * self.Ye_rho * E * Constants.YerhoE2a;

        // Compute See and Tee
        const See = self.A_sum - self.Dmsq21 * self.Ue2sq_vac - self.Dmsq31 * self.Ue3sq_vac;
        const Tee = self.Tmm_base * (1.0 - self.Ue3sq_vac - self.Ue2sq_vac);
        const C = Amatter * Tee;
        const A = self.A_sum + Amatter;

        // Get lambda3 from DMP approximation
        const xmat = Amatter / self.Dmsqee;
        var tmp = 1.0 - xmat;
        var lambda3 = self.Dmsq31 + 0.5 * self.Dmsqee * (xmat - 1.0 + @sqrt(tmp * tmp + 4.0 * self.s13sq * xmat));

        const B = self.Tmm_base + Amatter * See;

        // Newton-Raphson iterations
        var i: u8 = 0;
        while (i < self.n_newton) : (i += 1) {
            lambda3 = (lambda3 * lambda3 * (lambda3 - A) + C) / (lambda3 * (2.0 * lambda3 - A) + B);
        }

        tmp = A - lambda3;
        const Dlambda21 = @sqrt(tmp * tmp - 4.0 * C / lambda3);
        const lambda2 = 0.5 * (A - lambda3 + Dlambda21);
        const Dlambda32 = lambda3 - lambda2;
        const Dlambda31 = Dlambda32 + Dlambda21;

        const PiDlambdaInv = 1.0 / (Dlambda31 * Dlambda32 * Dlambda21);
        const Xp3 = PiDlambdaInv * Dlambda21;
        const Xp2 = -PiDlambdaInv * Dlambda31;

        const Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
        const Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

        const Smm = A - self.Dmsq21 * self.Um2sq_vac - self.Dmsq31 * self.Um3sq_vac;
        const Tmm = self.Tmm_base * (1.0 - self.Um3sq_vac - self.Um2sq_vac) + Amatter * (See + Smm - self.A_sum);

        const Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
        const Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

        // Jarlskog in matter
        var Jmatter = 8.0 * self.Jrr * self.c13sq * self.sind;
        Jmatter = Jmatter * self.Dmsq21 * self.Dmsq31 * (self.Dmsq31 - self.Dmsq21) * PiDlambdaInv;

        // Get all elements using unitarity
        const Ue1sq = 1.0 - Ue3sq - Ue2sq;
        const Um1sq = 1.0 - Um3sq - Um2sq;
        const Ut3sq = 1.0 - Um3sq - Ue3sq;
        const Ut2sq = 1.0 - Um2sq - Ue2sq;
        const Ut1sq = 1.0 - Um1sq - Ue1sq;

        // Kinematic terms
        const Lover4E = Constants.eVsqkm_to_GeV_over4 * L / E;
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
};

// =============================================================================
// Scalar Functions
// =============================================================================

/// Calculate vacuum oscillation probabilities
pub fn vacuumProbability(params: VacuumParams, L: f64, E: f64) ProbabilityMatrix {
    const batch = VacuumBatch.init(params);
    return batch.probabilityAt(L, E);
}

/// Calculate matter oscillation probabilities
pub fn matterProbability(params: MatterParams, L: f64, E: f64) ProbabilityMatrix {
    const vac = params.vacuum;
    const c13sq = 1.0 - vac.s13sq;
    const matter_sign = params.matterSign();
    const effective_delta = params.effectiveDelta();

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

    const Amatter = matter_sign * params.Ye * params.rho * E * Constants.YerhoE2a;
    const Dmsqee = vac.Dmsq31 - vac.s12sq * vac.Dmsq21;

    const A_sum = vac.Dmsq21 + vac.Dmsq31;
    const See = A_sum - vac.Dmsq21 * Ue2sq - vac.Dmsq31 * Ue3sq;
    const Tmm_base = vac.Dmsq21 * vac.Dmsq31;
    const Tee = Tmm_base * (1.0 - Ue3sq - Ue2sq);
    const C = Amatter * Tee;
    const A = A_sum + Amatter;

    const xmat = Amatter / Dmsqee;
    var tmp = 1.0 - xmat;
    var lambda3 = vac.Dmsq31 + 0.5 * Dmsqee * (xmat - 1.0 + @sqrt(tmp * tmp + 4.0 * vac.s13sq * xmat));

    const B = Tmm_base + Amatter * See;

    // Newton-Raphson iterations (clamped)
    const n_newton_clamped = params.clampedNewton();
    var i: u8 = 0;
    while (i < n_newton_clamped) : (i += 1) {
        lambda3 = (lambda3 * lambda3 * (lambda3 - A) + C) / (lambda3 * (2.0 * lambda3 - A) + B);
    }

    tmp = A - lambda3;
    const Dlambda21 = @sqrt(tmp * tmp - 4.0 * C / lambda3);
    const lambda2 = 0.5 * (A - lambda3 + Dlambda21);
    const Dlambda32 = lambda3 - lambda2;
    const Dlambda31 = Dlambda32 + Dlambda21;

    const PiDlambdaInv = 1.0 / (Dlambda31 * Dlambda32 * Dlambda21);
    const Xp3 = PiDlambdaInv * Dlambda21;
    const Xp2 = -PiDlambdaInv * Dlambda31;

    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

    const Smm = A - vac.Dmsq21 * Um2sq - vac.Dmsq31 * Um3sq;
    const Tmm = Tmm_base * (1.0 - Um3sq - Um2sq) + Amatter * (See + Smm - A_sum);

    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

    Jmatter = Jmatter * vac.Dmsq21 * vac.Dmsq31 * (vac.Dmsq31 - vac.Dmsq21) * PiDlambdaInv;

    const Ue1sq = 1.0 - Ue3sq - Ue2sq;
    const Um1sq = 1.0 - Um3sq - Um2sq;
    const Ut3sq = 1.0 - Um3sq - Ue3sq;
    Ut2sq = 1.0 - Um2sq - Ue2sq;
    const Ut1sq = 1.0 - Um1sq - Ue1sq;

    const Lover4E = Constants.eVsqkm_to_GeV_over4 * L / E;
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
// SIMD-optimized batch calculations (f64)
// =============================================================================

/// SIMD vector size for f64
pub const simd_len = std.simd.suggestVectorLength(f64) orelse 4;

/// SIMD vector of f64
pub const F64Vec = @Vector(simd_len, f64);

/// SIMD batch vacuum calculation
/// Calculates `simd_len` energy points simultaneously
pub fn vacuumProbabilitySimd(batch: VacuumBatch, L: f64, energies: F64Vec) [3][3]F64Vec {
    const L_vec: F64Vec = @splat(L);
    const conv: F64Vec = @splat(Constants.eVsqkm_to_GeV_over4);
    const Dmsq21_vec: F64Vec = @splat(batch.Dmsq21);
    const Dmsq31_vec: F64Vec = @splat(batch.Dmsq31);

    const Lover4E = conv * L_vec / energies;
    const D21 = Dmsq21_vec * Lover4E;
    const D31 = Dmsq31_vec * Lover4E;

    const sinD21 = @sin(D21);
    const sinD31 = @sin(D31);
    const sinD32 = @sin(D31 - D21);

    const triple_sin = sinD21 * sinD31 * sinD32;
    const two: F64Vec = @splat(2.0);
    const sinsqD21_2 = two * sinD21 * sinD21;
    const sinsqD31_2 = two * sinD31 * sinD31;
    const sinsqD32_2 = two * sinD32 * sinD32;

    const Ut3sq: F64Vec = @splat(batch.Ut3sq);
    const Ut2sq: F64Vec = @splat(batch.Ut2sq);
    const Ut1sq: F64Vec = @splat(batch.Ut1sq);
    const Um3sq: F64Vec = @splat(batch.Um3sq);
    const Um2sq: F64Vec = @splat(batch.Um2sq);
    const Um1sq: F64Vec = @splat(batch.Um1sq);
    const Ue3sq: F64Vec = @splat(batch.Ue3sq);
    const Ue2sq: F64Vec = @splat(batch.Ue2sq);
    const Ue1sq: F64Vec = @splat(batch.Ue1sq);
    const Jvac: F64Vec = @splat(batch.Jvac);

    const Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
        (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
        (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    const Pme_CPV = -Jvac * triple_sin;

    const one: F64Vec = @splat(1.0);
    const Pmm = one - two * (Um2sq * Um1sq * sinsqD21_2 +
        Um3sq * Um1sq * sinsqD31_2 +
        Um3sq * Um2sq * sinsqD32_2);
    const Pee = one - two * (Ue2sq * Ue1sq * sinsqD21_2 +
        Ue3sq * Ue1sq * sinsqD31_2 +
        Ue3sq * Ue2sq * sinsqD32_2);

    const Pme = Pme_CPC + Pme_CPV;
    const Pem = Pme_CPC - Pme_CPV;

    return .{
        .{ Pee, Pem, one - Pee - Pem },
        .{ Pme, Pmm, one - Pme - Pmm },
        .{ one - Pee - Pme, one - Pem - Pmm, one - (one - Pee - Pem) - (one - Pme - Pmm) },
    };
}

/// SIMD batch matter calculation
/// Calculates `simd_len` energy points simultaneously
pub fn matterProbabilitySimd(batch: MatterBatch, L: f64, energies: F64Vec) [3][3]F64Vec {
    const conv: F64Vec = @splat(Constants.eVsqkm_to_GeV_over4);
    const YerhoE2a: F64Vec = @splat(Constants.YerhoE2a);
    const L_vec: F64Vec = @splat(L);
    const Ye_rho: F64Vec = @splat(batch.Ye_rho);
    const matter_sign: F64Vec = @splat(batch.matter_sign);

    // Energy-dependent matter potentials
    const Amatter = matter_sign * Ye_rho * energies * YerhoE2a;

    // Splat scalar precomputed values
    const Dmsq21: F64Vec = @splat(batch.Dmsq21);
    const Dmsq31: F64Vec = @splat(batch.Dmsq31);
    const Dmsqee: F64Vec = @splat(batch.Dmsqee);
    const A_sum: F64Vec = @splat(batch.A_sum);
    const Tmm_base: F64Vec = @splat(batch.Tmm_base);
    const s13sq: F64Vec = @splat(batch.s13sq);
    const c13sq: F64Vec = @splat(batch.c13sq);
    const Ue2sq_vac: F64Vec = @splat(batch.Ue2sq_vac);
    const Ue3sq_vac: F64Vec = @splat(batch.Ue3sq_vac);
    const Um3sq_vac: F64Vec = @splat(batch.Um3sq_vac);
    const Um2sq_vac: F64Vec = @splat(batch.Um2sq_vac);
    const Jrr: F64Vec = @splat(batch.Jrr);
    const sind: F64Vec = @splat(batch.sind);

    // Compute See and Tee
    const See = A_sum - Dmsq21 * Ue2sq_vac - Dmsq31 * Ue3sq_vac;
    const one: F64Vec = @splat(1.0);
    const Tee = Tmm_base * (one - Ue3sq_vac - Ue2sq_vac);
    const C = Amatter * Tee;
    const A = A_sum + Amatter;

    // DMP approximation
    const xmat = Amatter / Dmsqee;
    var tmp = one - xmat;
    const four: F64Vec = @splat(4.0);
    const half: F64Vec = @splat(0.5);
    var lambda3 = Dmsq31 + half * Dmsqee * (xmat - one + @sqrt(tmp * tmp + four * s13sq * xmat));

    const B = Tmm_base + Amatter * See;

    // Newton-Raphson iterations (unrolled for SIMD)
    const two: F64Vec = @splat(2.0);
    comptime var iter: u8 = 0;
    inline while (iter < 3) : (iter += 1) {
        if (iter < batch.n_newton) {
            lambda3 = (lambda3 * lambda3 * (lambda3 - A) + C) / (lambda3 * (two * lambda3 - A) + B);
        }
    }

    tmp = A - lambda3;
    const Dlambda21 = @sqrt(tmp * tmp - four * C / lambda3);
    const lambda2 = half * (A - lambda3 + Dlambda21);
    const Dlambda32 = lambda3 - lambda2;
    const Dlambda31 = Dlambda32 + Dlambda21;

    const PiDlambdaInv = one / (Dlambda31 * Dlambda32 * Dlambda21);
    const Xp3 = PiDlambdaInv * Dlambda21;
    const Xp2 = -PiDlambdaInv * Dlambda31;

    const Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    const Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

    const Smm = A - Dmsq21 * Um2sq_vac - Dmsq31 * Um3sq_vac;
    const Tmm = Tmm_base * (one - Um3sq_vac - Um2sq_vac) + Amatter * (See + Smm - A_sum);

    const Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    const Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

    // Jarlskog in matter
    const eight: F64Vec = @splat(8.0);
    var Jmatter = eight * Jrr * c13sq * sind;
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;

    // Get all elements using unitarity
    const Ue1sq = one - Ue3sq - Ue2sq;
    const Um1sq = one - Um3sq - Um2sq;
    const Ut3sq = one - Um3sq - Ue3sq;
    const Ut2sq = one - Um2sq - Ue2sq;
    const Ut1sq = one - Um1sq - Ue1sq;

    // Kinematic terms
    const Lover4E = conv * L_vec / energies;
    const D21 = Dlambda21 * Lover4E;
    const D32 = Dlambda32 * Lover4E;

    const sinD21 = @sin(D21);
    const sinD31 = @sin(D32 + D21);
    const sinD32 = @sin(D32);

    const triple_sin = sinD21 * sinD31 * sinD32;
    const sinsqD21_2 = two * sinD21 * sinD21;
    const sinsqD31_2 = two * sinD31 * sinD31;
    const sinsqD32_2 = two * sinD32 * sinD32;

    const Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
        (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
        (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    const Pme_CPV = -Jmatter * triple_sin;

    const Pmm = one - two * (Um2sq * Um1sq * sinsqD21_2 +
        Um3sq * Um1sq * sinsqD31_2 +
        Um3sq * Um2sq * sinsqD32_2);
    const Pee = one - two * (Ue2sq * Ue1sq * sinsqD21_2 +
        Ue3sq * Ue1sq * sinsqD31_2 +
        Ue3sq * Ue2sq * sinsqD32_2);

    const Pme = Pme_CPC + Pme_CPV;
    const Pem = Pme_CPC - Pme_CPV;

    return .{
        .{ Pee, Pem, one - Pee - Pem },
        .{ Pme, Pmm, one - Pme - Pmm },
        .{ one - Pee - Pme, one - Pem - Pmm, one - (one - Pee - Pem) - (one - Pme - Pmm) },
    };
}

// =============================================================================
// f32 SIMD Mode (2× lanes: 8×f32 vs 4×f64)
// =============================================================================

/// SIMD vector size for f32 (typically 2× f64 lanes)
pub const simd_len_f32 = std.simd.suggestVectorLength(f32) orelse 8;

/// SIMD vector of f32
pub const F32Vec = @Vector(simd_len_f32, f32);

/// f32 probability matrix for SIMD results
pub const ProbabilityMatrixF32 = [3][3]F32Vec;

/// f32 version of vacuum batch (for faster computation with 2× lanes)
pub const VacuumBatchF32 = struct {
    Ue1sq: f32,
    Ue2sq: f32,
    Ue3sq: f32,
    Um1sq: f32,
    Um2sq: f32,
    Um3sq: f32,
    Ut1sq: f32,
    Ut2sq: f32,
    Ut3sq: f32,
    Jvac: f32,
    Dmsq21: f32,
    Dmsq31: f32,

    /// Create from f64 VacuumBatch
    pub fn fromF64(batch: VacuumBatch) VacuumBatchF32 {
        return .{
            .Ue1sq = @floatCast(batch.Ue1sq),
            .Ue2sq = @floatCast(batch.Ue2sq),
            .Ue3sq = @floatCast(batch.Ue3sq),
            .Um1sq = @floatCast(batch.Um1sq),
            .Um2sq = @floatCast(batch.Um2sq),
            .Um3sq = @floatCast(batch.Um3sq),
            .Ut1sq = @floatCast(batch.Ut1sq),
            .Ut2sq = @floatCast(batch.Ut2sq),
            .Ut3sq = @floatCast(batch.Ut3sq),
            .Jvac = @floatCast(batch.Jvac),
            .Dmsq21 = @floatCast(batch.Dmsq21),
            .Dmsq31 = @floatCast(batch.Dmsq31),
        };
    }

    /// Create from vacuum parameters
    pub fn init(params: VacuumParams) VacuumBatchF32 {
        return fromF64(VacuumBatch.init(params));
    }
};

/// f32 version of matter batch (for faster computation with 2× lanes)
pub const MatterBatchF32 = struct {
    s12sq: f32,
    s13sq: f32,
    s23sq: f32,
    c13sq: f32,
    Ue2sq_vac: f32,
    Ue3sq_vac: f32,
    Um3sq_vac: f32,
    Um2sq_vac: f32,
    Jrr: f32,
    sind: f32,
    cosd: f32,
    Dmsq21: f32,
    Dmsq31: f32,
    Dmsqee: f32,
    Tmm_base: f32,
    A_sum: f32,
    Ye_rho: f32,
    n_newton: u8,
    matter_sign: f32,

    /// Create from f64 MatterBatch
    pub fn fromF64(batch: MatterBatch) MatterBatchF32 {
        return .{
            .s12sq = @floatCast(batch.s12sq),
            .s13sq = @floatCast(batch.s13sq),
            .s23sq = @floatCast(batch.s23sq),
            .c13sq = @floatCast(batch.c13sq),
            .Ue2sq_vac = @floatCast(batch.Ue2sq_vac),
            .Ue3sq_vac = @floatCast(batch.Ue3sq_vac),
            .Um3sq_vac = @floatCast(batch.Um3sq_vac),
            .Um2sq_vac = @floatCast(batch.Um2sq_vac),
            .Jrr = @floatCast(batch.Jrr),
            .sind = @floatCast(batch.sind),
            .cosd = @floatCast(batch.cosd),
            .Dmsq21 = @floatCast(batch.Dmsq21),
            .Dmsq31 = @floatCast(batch.Dmsq31),
            .Dmsqee = @floatCast(batch.Dmsqee),
            .Tmm_base = @floatCast(batch.Tmm_base),
            .A_sum = @floatCast(batch.A_sum),
            .Ye_rho = @floatCast(batch.Ye_rho),
            .n_newton = batch.n_newton,
            .matter_sign = @floatCast(batch.matter_sign),
        };
    }

    /// Create from matter parameters
    pub fn init(params: MatterParams) MatterBatchF32 {
        return fromF64(MatterBatch.init(params));
    }
};

/// f32 SIMD vacuum calculation (2× lanes vs f64)
pub fn vacuumProbabilitySimdF32(batch: VacuumBatchF32, L: f32, energies: F32Vec) ProbabilityMatrixF32 {
    const conv = ConstantsFor(f32).eVsqkm_to_GeV_over4;
    const L_vec: F32Vec = @splat(L);
    const conv_vec: F32Vec = @splat(conv);
    const Dmsq21_vec: F32Vec = @splat(batch.Dmsq21);
    const Dmsq31_vec: F32Vec = @splat(batch.Dmsq31);

    const Lover4E = conv_vec * L_vec / energies;
    const D21 = Dmsq21_vec * Lover4E;
    const D31 = Dmsq31_vec * Lover4E;

    const sinD21 = @sin(D21);
    const sinD31 = @sin(D31);
    const sinD32 = @sin(D31 - D21);

    const triple_sin = sinD21 * sinD31 * sinD32;
    const two: F32Vec = @splat(2.0);
    const sinsqD21_2 = two * sinD21 * sinD21;
    const sinsqD31_2 = two * sinD31 * sinD31;
    const sinsqD32_2 = two * sinD32 * sinD32;

    const Ut3sq: F32Vec = @splat(batch.Ut3sq);
    const Ut2sq: F32Vec = @splat(batch.Ut2sq);
    const Ut1sq: F32Vec = @splat(batch.Ut1sq);
    const Um3sq: F32Vec = @splat(batch.Um3sq);
    const Um2sq: F32Vec = @splat(batch.Um2sq);
    const Um1sq: F32Vec = @splat(batch.Um1sq);
    const Ue3sq: F32Vec = @splat(batch.Ue3sq);
    const Ue2sq: F32Vec = @splat(batch.Ue2sq);
    const Ue1sq: F32Vec = @splat(batch.Ue1sq);
    const Jvac: F32Vec = @splat(batch.Jvac);

    const Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
        (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
        (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    const Pme_CPV = -Jvac * triple_sin;

    const one: F32Vec = @splat(1.0);
    const Pmm = one - two * (Um2sq * Um1sq * sinsqD21_2 +
        Um3sq * Um1sq * sinsqD31_2 +
        Um3sq * Um2sq * sinsqD32_2);
    const Pee = one - two * (Ue2sq * Ue1sq * sinsqD21_2 +
        Ue3sq * Ue1sq * sinsqD31_2 +
        Ue3sq * Ue2sq * sinsqD32_2);

    const Pme = Pme_CPC + Pme_CPV;
    const Pem = Pme_CPC - Pme_CPV;

    return .{
        .{ Pee, Pem, one - Pee - Pem },
        .{ Pme, Pmm, one - Pme - Pmm },
        .{ one - Pee - Pme, one - Pem - Pmm, one - (one - Pee - Pem) - (one - Pme - Pmm) },
    };
}

/// f32 SIMD matter calculation (2× lanes vs f64)
pub fn matterProbabilitySimdF32(batch: MatterBatchF32, L: f32, energies: F32Vec) ProbabilityMatrixF32 {
    const conv = ConstantsFor(f32).eVsqkm_to_GeV_over4;
    const YerhoE2a = ConstantsFor(f32).YerhoE2a;
    const L_vec: F32Vec = @splat(L);
    const conv_vec: F32Vec = @splat(conv);
    const YerhoE2a_vec: F32Vec = @splat(YerhoE2a);
    const Ye_rho: F32Vec = @splat(batch.Ye_rho);
    const matter_sign: F32Vec = @splat(batch.matter_sign);

    const Amatter = matter_sign * Ye_rho * energies * YerhoE2a_vec;

    const Dmsq21: F32Vec = @splat(batch.Dmsq21);
    const Dmsq31: F32Vec = @splat(batch.Dmsq31);
    const Dmsqee: F32Vec = @splat(batch.Dmsqee);
    const A_sum: F32Vec = @splat(batch.A_sum);
    const Tmm_base: F32Vec = @splat(batch.Tmm_base);
    const s13sq: F32Vec = @splat(batch.s13sq);
    const c13sq: F32Vec = @splat(batch.c13sq);
    const Ue2sq_vac: F32Vec = @splat(batch.Ue2sq_vac);
    const Ue3sq_vac: F32Vec = @splat(batch.Ue3sq_vac);
    const Um3sq_vac: F32Vec = @splat(batch.Um3sq_vac);
    const Um2sq_vac: F32Vec = @splat(batch.Um2sq_vac);
    const Jrr: F32Vec = @splat(batch.Jrr);
    const sind: F32Vec = @splat(batch.sind);

    const See = A_sum - Dmsq21 * Ue2sq_vac - Dmsq31 * Ue3sq_vac;
    const one: F32Vec = @splat(1.0);
    const Tee = Tmm_base * (one - Ue3sq_vac - Ue2sq_vac);
    const C = Amatter * Tee;
    const A = A_sum + Amatter;

    const xmat = Amatter / Dmsqee;
    var tmp = one - xmat;
    const four: F32Vec = @splat(4.0);
    const half: F32Vec = @splat(0.5);
    var lambda3 = Dmsq31 + half * Dmsqee * (xmat - one + @sqrt(tmp * tmp + four * s13sq * xmat));

    const B = Tmm_base + Amatter * See;

    const two: F32Vec = @splat(2.0);
    comptime var iter: u8 = 0;
    inline while (iter < 3) : (iter += 1) {
        if (iter < batch.n_newton) {
            lambda3 = (lambda3 * lambda3 * (lambda3 - A) + C) / (lambda3 * (two * lambda3 - A) + B);
        }
    }

    tmp = A - lambda3;
    const Dlambda21 = @sqrt(tmp * tmp - four * C / lambda3);
    const lambda2 = half * (A - lambda3 + Dlambda21);
    const Dlambda32 = lambda3 - lambda2;
    const Dlambda31 = Dlambda32 + Dlambda21;

    const PiDlambdaInv = one / (Dlambda31 * Dlambda32 * Dlambda21);
    const Xp3 = PiDlambdaInv * Dlambda21;
    const Xp2 = -PiDlambdaInv * Dlambda31;

    const Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    const Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

    const Smm = A - Dmsq21 * Um2sq_vac - Dmsq31 * Um3sq_vac;
    const Tmm = Tmm_base * (one - Um3sq_vac - Um2sq_vac) + Amatter * (See + Smm - A_sum);

    const Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    const Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

    const eight: F32Vec = @splat(8.0);
    var Jmatter = eight * Jrr * c13sq * sind;
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;

    const Ue1sq = one - Ue3sq - Ue2sq;
    const Um1sq = one - Um3sq - Um2sq;
    const Ut3sq = one - Um3sq - Ue3sq;
    const Ut2sq = one - Um2sq - Ue2sq;
    const Ut1sq = one - Um1sq - Ue1sq;

    const Lover4E = conv_vec * L_vec / energies;
    const D21 = Dlambda21 * Lover4E;
    const D32 = Dlambda32 * Lover4E;

    const sinD21 = @sin(D21);
    const sinD31 = @sin(D32 + D21);
    const sinD32 = @sin(D32);

    const triple_sin = sinD21 * sinD31 * sinD32;
    const sinsqD21_2 = two * sinD21 * sinD21;
    const sinsqD31_2 = two * sinD31 * sinD31;
    const sinsqD32_2 = two * sinD32 * sinD32;

    const Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
        (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
        (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    const Pme_CPV = -Jmatter * triple_sin;

    const Pmm = one - two * (Um2sq * Um1sq * sinsqD21_2 +
        Um3sq * Um1sq * sinsqD31_2 +
        Um3sq * Um2sq * sinsqD32_2);
    const Pee = one - two * (Ue2sq * Ue1sq * sinsqD21_2 +
        Ue3sq * Ue1sq * sinsqD31_2 +
        Ue3sq * Ue2sq * sinsqD32_2);

    const Pme = Pme_CPC + Pme_CPV;
    const Pem = Pme_CPC - Pme_CPV;

    return .{
        .{ Pee, Pem, one - Pee - Pem },
        .{ Pme, Pmm, one - Pme - Pmm },
        .{ one - Pee - Pme, one - Pem - Pmm, one - (one - Pee - Pem) - (one - Pme - Pmm) },
    };
}

// =============================================================================
// Tests
// =============================================================================

test "vacuum probability conservation" {
    const params = VacuumParams.default;
    const probs = vacuumProbability(params, 1300.0, 2.5);

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

test "matter probability conservation" {
    const params = MatterParams.default;
    const probs = matterProbability(params, 1300.0, 2.5);

    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| sum += p;
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }
}

test "vacuum batch matches single" {
    const params = VacuumParams.default;
    const batch = VacuumBatch.init(params);

    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const single = vacuumProbability(params, L, E);
    const batched = batch.probabilityAt(L, E);

    for (0..3) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(single[i][j], batched[i][j], 1e-15);
        }
    }
}

test "matter batch matches single" {
    const params = MatterParams.default;
    const batch = MatterBatch.init(params);

    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const single = matterProbability(params, L, E);
    const batched = batch.probabilityAt(L, E);

    for (0..3) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(single[i][j], batched[i][j], 1e-14);
        }
    }
}

test "SIMD vacuum matches scalar" {
    const params = VacuumParams.default;
    const batch = VacuumBatch.init(params);
    const L: f64 = 1300.0;

    var energies: F64Vec = undefined;
    for (0..simd_len) |i| {
        energies[i] = 0.5 + @as(f64, @floatFromInt(i)) * 0.5;
    }

    const simd_probs = vacuumProbabilitySimd(batch, L, energies);

    for (0..simd_len) |i| {
        const E = energies[i];
        const scalar = batch.probabilityAt(L, E);

        for (0..3) |r| {
            for (0..3) |c| {
                try std.testing.expectApproxEqAbs(scalar[r][c], simd_probs[r][c][i], 1e-14);
            }
        }
    }
}

test "SIMD matter matches scalar" {
    const params = MatterParams.default;
    const batch = MatterBatch.init(params);
    const L: f64 = 1300.0;

    var energies: F64Vec = undefined;
    for (0..simd_len) |i| {
        energies[i] = 0.5 + @as(f64, @floatFromInt(i)) * 0.5;
    }

    const simd_probs = matterProbabilitySimd(batch, L, energies);

    for (0..simd_len) |i| {
        const E = energies[i];
        const scalar = batch.probabilityAt(L, E);

        for (0..3) |r| {
            for (0..3) |c| {
                try std.testing.expectApproxEqAbs(scalar[r][c], simd_probs[r][c][i], 1e-12);
            }
        }
    }
}

test "f32 SIMD vacuum matches f64 scalar" {
    const params = VacuumParams.default;
    const batch_f64 = VacuumBatch.init(params);
    const batch_f32 = VacuumBatchF32.fromF64(batch_f64);
    const L: f32 = 1300.0;

    var energies: F32Vec = undefined;
    for (0..simd_len_f32) |i| {
        energies[i] = 0.5 + @as(f32, @floatFromInt(i)) * 0.3;
    }

    const simd_probs = vacuumProbabilitySimdF32(batch_f32, L, energies);

    for (0..simd_len_f32) |i| {
        const E = energies[i];
        const scalar = batch_f64.probabilityAt(L, E);

        for (0..3) |r| {
            for (0..3) |c| {
                // f32 has ~7 digits of precision, so we allow more tolerance
                try std.testing.expectApproxEqAbs(@as(f32, @floatCast(scalar[r][c])), simd_probs[r][c][i], 1e-5);
            }
        }
    }
}

test "antineutrino mode flips probabilities" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    var nu_params = VacuumParams.default;
    nu_params.antineutrino = false;

    var nubar_params = VacuumParams.default;
    nubar_params.antineutrino = true;

    const nu_probs = vacuumProbability(nu_params, L, E);
    const nubar_probs = vacuumProbability(nubar_params, L, E);

    // Pee should be the same (no CP violation in ee channel)
    try std.testing.expectApproxEqAbs(nu_probs[0][0], nubar_probs[0][0], 1e-10);

    // Pmm should be the same
    try std.testing.expectApproxEqAbs(nu_probs[1][1], nubar_probs[1][1], 1e-10);

    // Pme and Pem should swap (CPT theorem)
    try std.testing.expectApproxEqAbs(nu_probs[1][0], nubar_probs[0][1], 1e-10);
    try std.testing.expectApproxEqAbs(nu_probs[0][1], nubar_probs[1][0], 1e-10);
}

test "antineutrino matter mode" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    var nu_params = MatterParams.default;
    nu_params.antineutrino = false;

    var nubar_params = MatterParams.default;
    nubar_params.antineutrino = true;

    const nu_probs = matterProbability(nu_params, L, E);
    const nubar_probs = matterProbability(nubar_params, L, E);

    // Probabilities should differ due to matter effects
    // (matter has opposite sign for antineutrinos)
    const diff = @abs(nu_probs[1][0] - nubar_probs[1][0]);
    try std.testing.expect(diff > 0.01); // Should be noticeably different
}

// =============================================================================
// Cross-Validation Tests (against original NuFast Python implementation)
// =============================================================================

test "cross-validation: T2K-like vacuum (L=295 km, E=0.6 GeV)" {
    // Reference values from original NuFast Python implementation
    const params = VacuumParams.default;
    const probs = vacuumProbability(params, 295.0, 0.6);

    // Expected: Pee = 0.912467..., Pme = 0.057131..., Pmm = 0.010845...
    // These are from the original NuFast algorithm
    try std.testing.expectApproxEqAbs(probs[0][0], 0.9124677609000021, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 0.05713195504421340, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 0.01084557487519722, 1e-10);
}

test "cross-validation: DUNE-like vacuum (L=1300 km, E=2.5 GeV)" {
    const params = VacuumParams.default;
    const probs = vacuumProbability(params, 1300.0, 2.5);

    // Expected from original NuFast
    try std.testing.expectApproxEqAbs(probs[0][0], 0.9120570461953496, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 0.05859330228592972, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 0.004726368319321272, 1e-10);
}

test "cross-validation: inverted ordering (L=1300 km, E=2.5 GeV)" {
    var params = VacuumParams.default;
    params.Dmsq31 = -2.498e-3; // Inverted ordering

    const probs = vacuumProbability(params, 1300.0, 2.5);

    // Expected from original NuFast
    try std.testing.expectApproxEqAbs(probs[0][0], 0.9126759218475518, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 0.05722581249769201, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 0.01733314218311000, 1e-10);
}

// =============================================================================
// Edge Case Tests
// =============================================================================

test "edge case: near-zero baseline (L ≈ 0)" {
    const params = VacuumParams.default;
    const probs = vacuumProbability(params, 0.001, 2.5);

    // At L→0, should approach identity: P_αα → 1, P_αβ → 0
    try std.testing.expectApproxEqAbs(probs[0][0], 1.0, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 1.0, 1e-10);
    try std.testing.expectApproxEqAbs(probs[2][2], 1.0, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 0.0, 1e-10);
    try std.testing.expectApproxEqAbs(probs[0][1], 0.0, 1e-10);
}

test "edge case: very high energy (E = 100 GeV)" {
    const params = VacuumParams.default;
    const probs = vacuumProbability(params, 1300.0, 100.0);

    // At high energy, oscillations are suppressed (L/E small)
    // Should approach identity
    try std.testing.expectApproxEqAbs(probs[0][0], 0.9998609501566286, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 6.845956204204167e-05, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 0.9984513776819962, 1e-10);
}

test "edge case: zero theta13 (s13sq = 0)" {
    var params = VacuumParams.default;
    params.s13sq = 0.0;

    const probs = vacuumProbability(params, 1300.0, 2.5);

    // With θ₁₃ = 0, there's no νμ → νe appearance at atmospheric Δm²
    // (only solar oscillations contribute)
    // All probabilities should still be in [0, 1]
    for (probs) |row| {
        for (row) |p| {
            try std.testing.expect(p >= -1e-10 and p <= 1.0 + 1e-10);
        }
    }
    // Pee and Pmm should sum rows/cols to 1
    var row_sum: f64 = 0;
    for (probs[0]) |p| row_sum += p;
    try std.testing.expectApproxEqAbs(row_sum, 1.0, 1e-10);
}

test "edge case: maximal mixing (s23sq = 0.5)" {
    var params = VacuumParams.default;
    params.s23sq = 0.5; // Exactly maximal

    const probs = vacuumProbability(params, 1300.0, 2.5);

    // Should still have valid probabilities
    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| {
            try std.testing.expect(p >= -1e-10 and p <= 1.0 + 1e-10);
            sum += p;
        }
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }
}

test "edge case: delta = 0 (no CP violation)" {
    var params = VacuumParams.default;
    params.delta = 0.0;

    const probs = vacuumProbability(params, 1300.0, 2.5);

    // With δ = 0, Pme should equal Pem (T symmetry)
    try std.testing.expectApproxEqAbs(probs[1][0], probs[0][1], 1e-10);
}

test "edge case: delta = pi (maximal CP violation for this phase)" {
    var params = VacuumParams.default;
    params.delta = math.pi;

    const probs = vacuumProbability(params, 1300.0, 2.5);

    // Should still have valid probabilities
    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| sum += p;
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }
}

test "n_newton clamping: values > 3 are clamped" {
    var params = MatterParams.default;
    params.n_newton = 100; // Way over the limit

    // Should be clamped to MAX_NEWTON_ITERATIONS (3)
    try std.testing.expectEqual(params.clampedNewton(), MAX_NEWTON_ITERATIONS);

    // Should still compute valid probabilities
    const probs = matterProbability(params, 1300.0, 2.5);
    for (probs) |row| {
        var sum: f64 = 0;
        for (row) |p| sum += p;
        try std.testing.expectApproxEqAbs(sum, 1.0, 1e-10);
    }
}

test "matter: Newton iterations improve accuracy" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    // Get results for different Newton iteration counts
    var params_n0 = MatterParams.default;
    params_n0.n_newton = 0;
    const probs_n0 = matterProbability(params_n0, L, E);

    var params_n1 = MatterParams.default;
    params_n1.n_newton = 1;
    const probs_n1 = matterProbability(params_n1, L, E);

    var params_n2 = MatterParams.default;
    params_n2.n_newton = 2;
    const probs_n2 = matterProbability(params_n2, L, E);

    // N=1 and N=2 should be closer to each other than N=0 and N=1
    // (Newton iterations converge)
    const diff_01 = @abs(probs_n0[1][0] - probs_n1[1][0]);
    const diff_12 = @abs(probs_n1[1][0] - probs_n2[1][0]);

    try std.testing.expect(diff_12 < diff_01);
}

test "matter: zero density equals vacuum" {
    const L: f64 = 1300.0;
    const E: f64 = 2.5;

    const vac_params = VacuumParams.default;
    const vac_probs = vacuumProbability(vac_params, L, E);

    var mat_params = MatterParams.default;
    mat_params.rho = 0.0; // Zero density
    const mat_probs = matterProbability(mat_params, L, E);

    // With ρ = 0, matter should equal vacuum
    for (0..3) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(vac_probs[i][j], mat_probs[i][j], 1e-10);
        }
    }
}

// =============================================================================
// Cross-Validation Tests: MATTER (against original NuFast Python)
// =============================================================================

test "cross-validation: DUNE-like matter (L=1300 km, E=2.5 GeV, rho=2.848)" {
    // Reference values from original NuFast Python implementation
    // Parameters: NuFit 5.2 (s12sq=0.307, s13sq=0.0220, s23sq=0.546, delta=-0.7*pi)
    // Dmsq21=7.53e-5, Dmsq31=2.453e-3, rho=2.848 g/cm³, Ye=0.5, N_Newton=0
    const params = MatterParams.default;
    const probs = matterProbability(params, 1300.0, 2.5);

    // Python reference:
    //   P(e->e)  = 8.727627013058699e-01
    //   P(e->mu) = 5.731530983046839e-02
    //   P(mu->e) = 8.249293182552475e-02
    //   P(mu->mu)= 5.260049583365356e-02
    try std.testing.expectApproxEqAbs(probs[0][0], 8.727627013058699e-01, 1e-10);
    try std.testing.expectApproxEqAbs(probs[0][1], 5.731530983046839e-02, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 8.249293182552475e-02, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 5.260049583365356e-02, 1e-10);
}

test "cross-validation: T2K-like matter (L=295 km, E=0.6 GeV, rho=2.6)" {
    // Reference values from original NuFast Python implementation
    // Parameters: NuFit 5.2, rho=2.6 g/cm³, Ye=0.5, N_Newton=0
    var params = MatterParams.default;
    params.rho = 2.6;
    const probs = matterProbability(params, 295.0, 0.6);

    // Python reference:
    //   P(e->e)  = 9.051251298911471e-01
    //   P(e->mu) = 4.191271845601963e-02
    //   P(mu->e) = 6.284592128957496e-02
    //   P(mu->mu)= 1.004304479827633e-02
    try std.testing.expectApproxEqAbs(probs[0][0], 9.051251298911471e-01, 1e-10);
    try std.testing.expectApproxEqAbs(probs[0][1], 4.191271845601963e-02, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][0], 6.284592128957496e-02, 1e-10);
    try std.testing.expectApproxEqAbs(probs[1][1], 1.004304479827633e-02, 1e-10);
}

