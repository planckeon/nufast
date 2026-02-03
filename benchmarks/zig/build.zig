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

    // =============================================================================
    // WASM Target
    // =============================================================================
    
    // WASM library (exports C-ABI functions for JS interop)
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });
    
    const wasm_mod = b.createModule(.{
        .root_source_file = b.path("src/wasm_exports.zig"),
        .target = wasm_target,
        .optimize = .ReleaseSmall,
    });

    const wasm = b.addExecutable(.{
        .name = "nufast",
        .root_module = wasm_mod,
    });
    
    // Important: export all symbols and don't use entry point
    wasm.entry = .disabled;
    wasm.rdynamic = true;
    
    b.installArtifact(wasm);

    const wasm_step = b.step("wasm", "Build WASM library");
    wasm_step.dependOn(&wasm.step);
    
    // WASM with SIMD support (for browsers with WASM SIMD)
    const wasm_simd_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
        .cpu_features_add = std.Target.wasm.featureSet(&.{.simd128}),
    });
    
    const wasm_simd_mod = b.createModule(.{
        .root_source_file = b.path("src/wasm_exports.zig"),
        .target = wasm_simd_target,
        .optimize = .ReleaseSmall,
    });

    const wasm_simd = b.addExecutable(.{
        .name = "nufast-simd",
        .root_module = wasm_simd_mod,
    });
    
    wasm_simd.entry = .disabled;
    wasm_simd.rdynamic = true;
    
    b.installArtifact(wasm_simd);

    const wasm_simd_step = b.step("wasm-simd", "Build WASM library with SIMD");
    wasm_simd_step.dependOn(&wasm_simd.step);
}
