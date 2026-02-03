# NuFast Zig Implementation

Fast three-flavor neutrino oscillation probabilities in Zig.

## Features

- **Vacuum oscillations**: Full 3-flavor PMNS probabilities
- **Matter effects (MSW)**: Constant-density matter with arbitrary electron fraction
- **SIMD acceleration**: Vectorized calculations for multiple energies
  - f64: 4 lanes (AVX2) 
  - f32: 8 lanes (2× throughput, reduced precision)
- **Batch APIs**: Pre-computed mixing for repeated calculations
- **Anti-neutrino mode**: Sign-flipped matter potential and δCP
- **Zero allocations**: All computation on the stack

## Quick Start

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

// SIMD: 4 energies simultaneously (f64)
var E_vec: nufast.F64Vec = .{ 1.0, 2.0, 3.0, 4.0 };
const p_vec = nufast.vacuumProbabilitySimd(batch, 1300.0, E_vec);

// f32 SIMD: 8 energies simultaneously (2× throughput)
const batch_f32 = nufast.VacuumBatchF32.fromF64(batch);
var E_vec_f32: nufast.F32Vec = .{ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 };
const p_vec_f32 = nufast.vacuumProbabilitySimdF32(batch_f32, 1300.0, E_vec_f32);
```

## Anti-Neutrino Mode

```zig
// Set antineutrino = true to flip sign of matter potential and δCP
var matter_params = nufast.MatterParams.default;
matter_params.antineutrino = true;
const nubar_probs = nufast.matterProbability(matter_params, 1300.0, 2.5);
```

## Matter Batch API

For repeated calculations at constant density (e.g., energy spectra):

```zig
const matter_params = nufast.MatterParams{
    .vacuum = nufast.VacuumParams.default,
    .rho = 2.848,
    .Ye = 0.5,
    .n_newton = 0,
};
const matter_batch = nufast.MatterBatch.init(matter_params);

// 30-40% faster than repeated matterProbability calls
for (energies) |E| {
    const p = matter_batch.probabilityAt(1300.0, E);
}

// SIMD matter calculation
const p_simd = nufast.matterProbabilitySimd(matter_batch, 1300.0, E_vec);
```

## Performance

Benchmarks on AMD Ryzen (WSL2), 10M iterations:

### Scalar Performance

| Mode | Time/calc | Throughput |
|------|-----------|------------|
| Vacuum (f64) | 31 ns | 32 M/s |
| Matter N=0 (f64) | 85 ns | 12 M/s |
| Matter N=1 (f64) | 98 ns | 10 M/s |

### SIMD Performance (per-calculation, amortized)

| Mode | f64 (4 lanes) | f32 (8 lanes) |
|------|---------------|---------------|
| Vacuum | 30 ns → 33 M/s | **14 ns → 70 M/s** |
| Matter N=0 | 38 ns → 26 M/s | **26 ns → 39 M/s** |

**Key findings:**
- f32 SIMD achieves **70 M/s** vacuum throughput (2× lanes, 2.3× speedup)
- f32 SIMD achieves **39 M/s** matter throughput (3× faster than scalar)
- SIMD provides near-linear speedup for embarrassingly parallel energy spectra

## Building

Requires Zig 0.15.2 or later.

```bash
# Run tests
zig build test

# Run scalar benchmark
zig build bench

# Run SIMD benchmark
zig build simd
```

## Using as a Dependency

Add to your `build.zig.zon`:

```zig
.dependencies = .{
    .nufast = .{
        .url = "https://github.com/planckeon/nufast/archive/refs/tags/v0.4.0.tar.gz",
        .hash = "...",
    },
},
```

Then in `build.zig`:

```zig
const nufast = b.dependency("nufast", .{
    .target = target,
    .optimize = optimize,
});
exe.root_module.addImport("nufast", nufast.module("nufast"));
```

## API Reference

### Types

- `VacuumParams` - Oscillation parameters for vacuum
- `MatterParams` - Parameters including matter density
- `VacuumBatch` - Pre-computed vacuum mixing (f64)
- `MatterBatch` - Pre-computed matter mixing (f64)
- `VacuumBatchF32` - Pre-computed vacuum mixing (f32)
- `MatterBatchF32` - Pre-computed matter mixing (f32)
- `F64Vec` - SIMD vector of f64 (4 lanes on AVX2)
- `F32Vec` - SIMD vector of f32 (8 lanes on AVX2)

### Functions

- `vacuumProbability(params, L, E)` - Single vacuum calculation
- `matterProbability(params, L, E)` - Single matter calculation
- `vacuumProbabilitySimd(batch, L, energies)` - SIMD vacuum (f64)
- `matterProbabilitySimd(batch, L, energies)` - SIMD matter (f64)
- `vacuumProbabilitySimdF32(batch, L, energies)` - SIMD vacuum (f32)
- `matterProbabilitySimdF32(batch, L, energies)` - SIMD matter (f32)

## License

MIT
