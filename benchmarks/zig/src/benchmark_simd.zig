//! NuFast SIMD Benchmark
//!
//! Benchmarks SIMD batch calculations vs scalar.
//! Run with: zig build simd

const std = @import("std");
const nufast = @import("nufast");
const time = std.time;

const simd_len = nufast.simd_len;
const F64Vec = nufast.F64Vec;

// Total calculations (must be divisible by simd_len)
const n_scalar: usize = 10_000_000;
const n_simd: usize = n_scalar / simd_len;

const Emin: f64 = 0.5;
const Emax: f64 = 5.0;
const L: f64 = 1300.0;

var sink: f64 = 0;

pub fn main() !void {
    std.debug.print("NuFast SIMD Benchmark\n", .{});
    std.debug.print("=====================\n\n", .{});
    std.debug.print("SIMD vector length: {} × f64 ({} bytes)\n", .{ simd_len, simd_len * 8 });
    std.debug.print("Total calculations: {}\n\n", .{n_scalar});

    const params = nufast.VacuumParams.default;
    const batch = nufast.VacuumBatch.init(params);

    // Scalar benchmark
    {
        var total: f64 = 0;
        const start = time.nanoTimestamp();

        for (0..n_scalar) |i| {
            const E = Emin + (Emax - Emin) * @as(f64, @floatFromInt(i % 1000)) / 1000.0;
            const probs = batch.probabilityAt(L, E);
            total += probs[1][0];
        }

        const end = time.nanoTimestamp();
        sink += total;

        const elapsed_ns: f64 = @floatFromInt(end - start);
        const ns_per_call = elapsed_ns / @as(f64, @floatFromInt(n_scalar));
        const throughput = 1e9 / ns_per_call / 1e6;

        std.debug.print("Scalar:\n", .{});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_call});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // SIMD benchmark
    {
        var total: f64 = 0;
        const start = time.nanoTimestamp();

        for (0..n_simd) |batch_idx| {
            // Create energy vector
            var energies: F64Vec = undefined;
            for (0..simd_len) |i| {
                const idx = batch_idx * simd_len + i;
                energies[i] = Emin + (Emax - Emin) * @as(f64, @floatFromInt(idx % 1000)) / 1000.0;
            }

            const probs = nufast.vacuumProbabilitySimd(batch, L, energies);

            // Sum all Pme values
            var lane_sum: f64 = 0;
            for (0..simd_len) |i| {
                lane_sum += probs[1][0][i];
            }
            total += lane_sum;
        }

        const end = time.nanoTimestamp();
        sink += total;

        const elapsed_ns: f64 = @floatFromInt(end - start);
        const ns_per_calc = elapsed_ns / @as(f64, @floatFromInt(n_scalar));
        const ns_per_batch = elapsed_ns / @as(f64, @floatFromInt(n_simd));
        const throughput = 1e9 / ns_per_calc / 1e6;

        std.debug.print("SIMD ({} lanes):\n", .{simd_len});
        std.debug.print("  Time per batch: {d:8.2} ns\n", .{ns_per_batch});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_calc});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // Speedup
    std.debug.print("Expected speedup: {}×\n", .{simd_len});
    std.debug.print("(Actual speedup depends on memory bandwidth and instruction mix)\n", .{});

    if (sink == 0) std.debug.print(".", .{});
}
