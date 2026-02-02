const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Create the main nufast module
    const nufast_mod = b.addModule("nufast", .{
        .root_source_file = b.path("src/nufast.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Unit tests
    const lib_unit_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/nufast.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    // Benchmark executable
    const bench_mod = b.createModule(.{
        .root_source_file = b.path("src/benchmark.zig"),
        .target = target,
        .optimize = .ReleaseFast,
    });
    bench_mod.addImport("nufast", nufast_mod);

    const bench = b.addExecutable(.{
        .name = "benchmark",
        .root_module = bench_mod,
    });
    b.installArtifact(bench);

    const run_bench = b.addRunArtifact(bench);
    const bench_step = b.step("bench", "Run benchmarks");
    bench_step.dependOn(&run_bench.step);

    // SIMD benchmark
    const simd_mod = b.createModule(.{
        .root_source_file = b.path("src/benchmark_simd.zig"),
        .target = target,
        .optimize = .ReleaseFast,
    });
    simd_mod.addImport("nufast", nufast_mod);

    const simd_bench = b.addExecutable(.{
        .name = "benchmark_simd",
        .root_module = simd_mod,
    });
    b.installArtifact(simd_bench);

    const run_simd = b.addRunArtifact(simd_bench);
    const simd_step = b.step("simd", "Run SIMD benchmarks");
    simd_step.dependOn(&run_simd.step);
}
