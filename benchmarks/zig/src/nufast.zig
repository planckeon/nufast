//! NuFast - Fast three-flavor neutrino oscillation probabilities
//!
//! Zig implementation of the NuFast algorithm by Denton & Parke.
//! Provides high-performance oscillation probability calculations
//! for vacuum and constant-density matter.
//!
//! Reference: arXiv:2405.02400

const std = @import("std");
const math = std.math;

/// Physical constants
pub const Constants = struct {
    /// Conversion factor: Y_e * rho * E to matter potential
    pub const YerhoE2a: f64 = 1.52e-4;
    /// Conversion: eV² km to GeV / 4
    pub const eVsqkm_to_GeV_over4: f64 = 1e-9 / 1.97327e-7 * 1e3 / 4.0;
};

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

    /// Default NuFIT 5.2 parameters
    pub const default = VacuumParams{
        .s12sq = 0.307,
        .s13sq = 0.0220,
        .s23sq = 0.546,
        .delta = -0.7 * math.pi,
        .Dmsq21 = 7.53e-5,
        .Dmsq31 = 2.453e-3,
    };
};

/// Oscillation parameters for matter calculations
pub const MatterParams = struct {
    /// Vacuum parameters
    vacuum: VacuumParams,
    /// Matter density (g/cm³)
    rho: f64,
    /// Electron fraction
    Ye: f64,
    /// Newton-Raphson iterations (0-3)
    n_newton: u8,

    /// Default: DUNE-like matter
    pub const default = MatterParams{
        .vacuum = VacuumParams.default,
        .rho = 2.848,
        .Ye = 0.5,
        .n_newton = 0,
    };
};

/// 3×3 probability matrix [from][to]
pub const ProbabilityMatrix = [3][3]f64;

/// Pre-computed mixing matrix elements for batch calculations
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

    /// Create batch calculator from vacuum parameters
    pub fn init(params: VacuumParams) VacuumBatch {
        const c13sq = 1.0 - params.s13sq;

        const Ue3sq = params.s13sq;
        const Ue2sq = c13sq * params.s12sq;
        const Um3sq = c13sq * params.s23sq;
        var Ut2sq = params.s13sq * params.s12sq * params.s23sq;
        var Um2sq = (1.0 - params.s12sq) * (1.0 - params.s23sq);

        const Jrr = @sqrt(Um2sq * Ut2sq);
        const sind = @sin(params.delta);
        const cosd = @cos(params.delta);
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

/// Calculate vacuum oscillation probabilities
pub fn vacuumProbability(params: VacuumParams, L: f64, E: f64) ProbabilityMatrix {
    const batch = VacuumBatch.init(params);
    return batch.probabilityAt(L, E);
}

/// Calculate matter oscillation probabilities
pub fn matterProbability(params: MatterParams, L: f64, E: f64) ProbabilityMatrix {
    const vac = params.vacuum;
    const c13sq = 1.0 - vac.s13sq;

    var Ue2sq = c13sq * vac.s12sq;
    var Ue3sq = vac.s13sq;
    var Um3sq = c13sq * vac.s23sq;
    var Ut2sq = vac.s13sq * vac.s12sq * vac.s23sq;
    var Um2sq = (1.0 - vac.s12sq) * (1.0 - vac.s23sq);

    const Jrr = @sqrt(Um2sq * Ut2sq);
    const sind = @sin(vac.delta);
    const cosd = @cos(vac.delta);
    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd;
    var Jmatter = 8.0 * Jrr * c13sq * sind;

    const Amatter = params.Ye * params.rho * E * Constants.YerhoE2a;
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

    // Newton-Raphson iterations
    var i: u8 = 0;
    while (i < params.n_newton) : (i += 1) {
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
// SIMD-optimized batch calculations
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

test "SIMD matches scalar" {
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
