//! NuFast SIMD Benchmark
//!
//! Benchmarks SIMD batch calculations vs scalar for both vacuum and matter.
//! Tests f64 and f32 modes.
//! Run with: zig build simd

const std = @import("std");
const nufast = @import("nufast");
const time = std.time;

const simd_len = nufast.simd_len;
const simd_len_f32 = nufast.simd_len_f32;
const F64Vec = nufast.F64Vec;
const F32Vec = nufast.F32Vec;

// Total calculations (must be divisible by simd_len)
const n_scalar: usize = 10_000_000;
const n_simd_f64: usize = n_scalar / simd_len;
const n_simd_f32: usize = n_scalar / simd_len_f32;

const Emin: f64 = 0.5;
const Emax: f64 = 5.0;
const L: f64 = 1300.0;

var sink: f64 = 0;
var sink_f32: f32 = 0;

pub fn main() !void {
    std.debug.print("NuFast SIMD Benchmark\n", .{});
    std.debug.print("=====================\n\n", .{});
    std.debug.print("SIMD vector lengths:\n", .{});
    std.debug.print("  f64: {} lanes ({} bytes)\n", .{ simd_len, simd_len * 8 });
    std.debug.print("  f32: {} lanes ({} bytes)\n", .{ simd_len_f32, simd_len_f32 * 4 });
    std.debug.print("Total calculations: {}\n\n", .{n_scalar});

    const params = nufast.VacuumParams.default;
    const batch = nufast.VacuumBatch.init(params);
    const batch_f32 = nufast.VacuumBatchF32.fromF64(batch);

    const matter_params = nufast.MatterParams.default;
    const matter_batch = nufast.MatterBatch.init(matter_params);
    const matter_batch_f32 = nufast.MatterBatchF32.fromF64(matter_batch);

    // ==========================================================================
    // VACUUM BENCHMARKS
    // ==========================================================================
    std.debug.print("=== VACUUM OSCILLATIONS ===\n\n", .{});

    // Scalar vacuum benchmark
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

        std.debug.print("Scalar (f64):\n", .{});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_call});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // SIMD f64 vacuum benchmark
    {
        var total: f64 = 0;
        const start = time.nanoTimestamp();

        for (0..n_simd_f64) |batch_idx| {
            var energies: F64Vec = undefined;
            for (0..simd_len) |i| {
                const idx = batch_idx * simd_len + i;
                energies[i] = Emin + (Emax - Emin) * @as(f64, @floatFromInt(idx % 1000)) / 1000.0;
            }

            const probs = nufast.vacuumProbabilitySimd(batch, L, energies);

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
        const ns_per_batch = elapsed_ns / @as(f64, @floatFromInt(n_simd_f64));
        const throughput = 1e9 / ns_per_calc / 1e6;

        std.debug.print("SIMD f64 ({} lanes):\n", .{simd_len});
        std.debug.print("  Time per batch: {d:8.2} ns\n", .{ns_per_batch});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_calc});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // SIMD f32 vacuum benchmark
    {
        var total: f32 = 0;
        const start = time.nanoTimestamp();
        const L_f32: f32 = @floatCast(L);

        for (0..n_simd_f32) |batch_idx| {
            var energies: F32Vec = undefined;
            for (0..simd_len_f32) |i| {
                const idx = batch_idx * simd_len_f32 + i;
                energies[i] = @as(f32, @floatCast(Emin)) + @as(f32, @floatCast(Emax - Emin)) * @as(f32, @floatFromInt(idx % 1000)) / 1000.0;
            }

            const probs = nufast.vacuumProbabilitySimdF32(batch_f32, L_f32, energies);

            var lane_sum: f32 = 0;
            for (0..simd_len_f32) |i| {
                lane_sum += probs[1][0][i];
            }
            total += lane_sum;
        }

        const end = time.nanoTimestamp();
        sink_f32 += total;

        const elapsed_ns: f64 = @floatFromInt(end - start);
        const ns_per_calc = elapsed_ns / @as(f64, @floatFromInt(n_scalar));
        const ns_per_batch = elapsed_ns / @as(f64, @floatFromInt(n_simd_f32));
        const throughput = 1e9 / ns_per_calc / 1e6;

        std.debug.print("SIMD f32 ({} lanes):\n", .{simd_len_f32});
        std.debug.print("  Time per batch: {d:8.2} ns\n", .{ns_per_batch});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_calc});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // ==========================================================================
    // MATTER BENCHMARKS
    // ==========================================================================
    std.debug.print("=== MATTER OSCILLATIONS (N=0) ===\n\n", .{});

    // Scalar matter benchmark
    {
        var total: f64 = 0;
        const start = time.nanoTimestamp();

        for (0..n_scalar) |i| {
            const E = Emin + (Emax - Emin) * @as(f64, @floatFromInt(i % 1000)) / 1000.0;
            const probs = matter_batch.probabilityAt(L, E);
            total += probs[1][0];
        }

        const end = time.nanoTimestamp();
        sink += total;

        const elapsed_ns: f64 = @floatFromInt(end - start);
        const ns_per_call = elapsed_ns / @as(f64, @floatFromInt(n_scalar));
        const throughput = 1e9 / ns_per_call / 1e6;

        std.debug.print("Scalar (f64):\n", .{});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_call});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // SIMD f64 matter benchmark
    {
        var total: f64 = 0;
        const start = time.nanoTimestamp();

        for (0..n_simd_f64) |batch_idx| {
            var energies: F64Vec = undefined;
            for (0..simd_len) |i| {
                const idx = batch_idx * simd_len + i;
                energies[i] = Emin + (Emax - Emin) * @as(f64, @floatFromInt(idx % 1000)) / 1000.0;
            }

            const probs = nufast.matterProbabilitySimd(matter_batch, L, energies);

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
        const ns_per_batch = elapsed_ns / @as(f64, @floatFromInt(n_simd_f64));
        const throughput = 1e9 / ns_per_calc / 1e6;

        std.debug.print("SIMD f64 ({} lanes):\n", .{simd_len});
        std.debug.print("  Time per batch: {d:8.2} ns\n", .{ns_per_batch});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_calc});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // SIMD f32 matter benchmark
    {
        var total: f32 = 0;
        const start = time.nanoTimestamp();
        const L_f32: f32 = @floatCast(L);

        for (0..n_simd_f32) |batch_idx| {
            var energies: F32Vec = undefined;
            for (0..simd_len_f32) |i| {
                const idx = batch_idx * simd_len_f32 + i;
                energies[i] = @as(f32, @floatCast(Emin)) + @as(f32, @floatCast(Emax - Emin)) * @as(f32, @floatFromInt(idx % 1000)) / 1000.0;
            }

            const probs = nufast.matterProbabilitySimdF32(matter_batch_f32, L_f32, energies);

            var lane_sum: f32 = 0;
            for (0..simd_len_f32) |i| {
                lane_sum += probs[1][0][i];
            }
            total += lane_sum;
        }

        const end = time.nanoTimestamp();
        sink_f32 += total;

        const elapsed_ns: f64 = @floatFromInt(end - start);
        const ns_per_calc = elapsed_ns / @as(f64, @floatFromInt(n_scalar));
        const ns_per_batch = elapsed_ns / @as(f64, @floatFromInt(n_simd_f32));
        const throughput = 1e9 / ns_per_calc / 1e6;

        std.debug.print("SIMD f32 ({} lanes):\n", .{simd_len_f32});
        std.debug.print("  Time per batch: {d:8.2} ns\n", .{ns_per_batch});
        std.debug.print("  Time per calc:  {d:8.2} ns\n", .{ns_per_calc});
        std.debug.print("  Throughput:     {d:8.2} M/s\n\n", .{throughput});
    }

    // Summary
    std.debug.print("=== SUMMARY ===\n\n", .{});
    std.debug.print("SIMD f64 provides {}× parallelism\n", .{simd_len});
    std.debug.print("SIMD f32 provides {}× parallelism (2× more lanes)\n", .{simd_len_f32});
    std.debug.print("f32 mode trades precision for throughput\n", .{});

    if (sink == 0) std.debug.print(".", .{});
    if (sink_f32 == 0) std.debug.print(".", .{});
}
