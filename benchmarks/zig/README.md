# NuFast Zig Implementation

High-performance three-flavor neutrino oscillation probabilities in Zig.

## Features

- **Zero dependencies** - Pure Zig implementation
- **SIMD vectorization** - Batch calculations using native SIMD
- **VacuumBatch API** - Pre-computed mixing matrix for repeated calculations
- **Allocator-free** - No heap allocations in hot paths
- **Cross-platform** - Works on any Zig 0.15.2+ target

## Usage

### As a Dependency

Add to your `build.zig.zon`:

```zig
.dependencies = .{
    .nufast = .{
        .url = "https://github.com/planckeon/nufast/archive/refs/heads/main.tar.gz",
        .hash = "...", // Get hash from: zig fetch <url>
    },
},
```

In `build.zig`:

```zig
const nufast = b.dependency("nufast", .{
    .target = target,
    .optimize = optimize,
});
exe.root_module.addImport("nufast", nufast.module("nufast"));
```

### API

```zig
const nufast = @import("nufast");

// Vacuum oscillation
const params = nufast.VacuumParams.default;
const probs = nufast.vacuumProbability(params, 1300.0, 2.5);
// probs[1][0] = P(νμ → νe)

// Batch calculations (pre-computed mixing matrix)
const batch = nufast.VacuumBatch.init(params);
for (energies) |E| {
    const p = batch.probabilityAt(1300.0, E);
}

// SIMD batch (4-8 energies simultaneously)
var E_vec: nufast.F64Vec = .{ 1.0, 2.0, 3.0, 4.0 };
const p_vec = nufast.vacuumProbabilitySimd(batch, 1300.0, E_vec);

// Matter oscillation
const matter = nufast.MatterParams{
    .vacuum = params,
    .rho = 2.848, // g/cm³
    .Ye = 0.5,
    .n_newton = 0,
};
const p_matter = nufast.matterProbability(matter, 1300.0, 2.5);
```

## Building

```bash
# Run tests
zig build test

# Run scalar benchmark
zig build bench && ./zig-out/bin/benchmark

# Run SIMD benchmark
zig build simd && ./zig-out/bin/benchmark_simd
```

## Benchmark Results (AMD Ryzen, WSL2, 10M iterations)

| Implementation | Vacuum | Matter (N=0) | Matter (N=1) |
|---------------|--------|--------------|--------------|
| **Zig scalar**| **25.6 ns** | **73.3 ns** | **77.7 ns** |
| Zig SIMD (4×) | 23.3 ns | N/A          | N/A          |
| Rust          | 61 ns  | 95 ns        | 106 ns       |
| C++           | 49 ns  | 130 ns       | 143 ns       |
| Fortran       | 51 ns  | 107 ns       | 123 ns       |

### Key Findings

1. **Zig is 2.4× faster than Rust** for vacuum oscillations
2. **Zig is 1.3× faster than Rust** for matter oscillations
3. **Zig is 1.8× faster than C++** for matter (the more complex case)
4. SIMD provides ~10% additional speedup for vacuum

### Throughput

| Implementation | Vacuum | Matter (N=0) |
|---------------|--------|--------------|
| **Zig**       | **39.2 M/s** | **13.6 M/s** |
| Rust          | 17.5 M/s | 10.5 M/s     |
| C++           | 20.3 M/s | 7.7 M/s      |

## Algorithm

This implements the NuFast algorithm by Denton & Parke ([arXiv:2405.02400](https://arxiv.org/abs/2405.02400)).

Key optimizations:
1. **Eigenvalue-Eigenvector Identity (EEI)** - Avoids cubic equation solving
2. **DMP approximation** - Fast initial eigenvalue estimate
3. **Newton-Raphson refinement** - Optional precision iterations
4. **SIMD vectorization** - Process multiple energies in parallel

## License

MIT
