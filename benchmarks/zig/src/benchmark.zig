//! NuFast Zig Benchmark
//!
//! Benchmarks the Zig implementation against other languages.
//! Run with: zig build bench

const std = @import("std");
const nufast = @import("nufast");
const time = std.time;

const n: usize = 10_000_000;
const Emin: f64 = 0.5;
const Emax: f64 = 5.0;

// Prevent optimizer from removing computations
var sink: f64 = 0;

fn getStdout() std.fs.File.Writer {
    return std.io.getStdOut().writer();
}

fn benchmark(comptime name: []const u8, comptime func: anytype, args: anytype) void {
    // Warmup
    var warmup_total: f64 = 0;
    for (0..10_000) |i| {
        const E = Emin + (Emax - Emin) * @as(f64, @floatFromInt(i % 1000)) / 1000.0;
        const probs = @call(.auto, func, args ++ .{E});
        warmup_total += probs[1][0];
    }
    sink += warmup_total;

    // Timed run
    var total: f64 = 0;
    const start = time.nanoTimestamp();

    for (0..n) |i| {
        const E = Emin + (Emax - Emin) * @as(f64, @floatFromInt(i % 1000)) / 1000.0;
        const probs = @call(.auto, func, args ++ .{E});
        total += probs[1][0];
    }

    const end = time.nanoTimestamp();
    sink += total;

    const elapsed_ns: f64 = @floatFromInt(end - start);
    const ns_per_call = elapsed_ns / @as(f64, @floatFromInt(n));

    std.debug.print("  {s}: {d:10.2} ns\n", .{ name, ns_per_call });
}

fn vacuumWrapper(params: nufast.VacuumParams, L: f64, E: f64) nufast.ProbabilityMatrix {
    return nufast.vacuumProbability(params, L, E);
}

fn matterWrapper(params: nufast.MatterParams, L: f64, E: f64) nufast.ProbabilityMatrix {
    return nufast.matterProbability(params, L, E);
}

fn batchWrapper(batch: nufast.VacuumBatch, L: f64, E: f64) nufast.ProbabilityMatrix {
    return batch.probabilityAt(L, E);
}

pub fn main() !void {
    std.debug.print("NuFast Zig Benchmark (n={} iterations)\n", .{n});
    std.debug.print("============================================\n\n", .{});
    std.debug.print("SIMD vector length: {} × f64\n\n", .{nufast.simd_len});

    const L: f64 = 1300.0;
    const vacuum_params = nufast.VacuumParams.default;
    const batch = nufast.VacuumBatch.init(vacuum_params);

    std.debug.print("Single-point calculations:\n", .{});

    // Vacuum (direct)
    benchmark("Vacuum (direct)", vacuumWrapper, .{ vacuum_params, L });

    // Vacuum (batch precomputed)
    benchmark("Vacuum (batch) ", batchWrapper, .{ batch, L });

    // Matter with different Newton iterations
    inline for (0..4) |N| {
        const matter_params = nufast.MatterParams{
            .vacuum = vacuum_params,
            .rho = 3.0,
            .Ye = 0.5,
            .n_newton = N,
        };
        const name = std.fmt.comptimePrint("Matter N={d}      ", .{N});
        benchmark(name, matterWrapper, .{ matter_params, L });
    }

    std.debug.print("\n", .{});

    // Final sink to prevent DCE
    if (sink == 0) std.debug.print(".", .{});
}
