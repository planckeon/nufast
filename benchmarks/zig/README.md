# nufast (Zig)

Three-flavor neutrino oscillation probabilities in vacuum and constant-density matter.

Zig implementation of the NuFast algorithm by P.B. Denton ([arXiv:2405.02400](https://arxiv.org/abs/2405.02400)).

## Features

- Vacuum oscillations (exact analytic)
- Matter effects via DMP approximation with Newton refinement
- Anti-neutrino mode (sign-flipped δCP and matter potential)
- SIMD vectorization (4×f64 or 8×f32 on AVX2)
- f32 mode for increased throughput
- Batch APIs for pre-computed mixing matrices
- Zero dependencies (std only)
- Zero heap allocations in hot paths

## Requirements

- Zig 0.15.0 or later

## Installation

### Using `zig fetch` (recommended)

```bash
zig fetch --save git+https://github.com/planckeon/nufast.git#main
```

Or with a specific version:

```bash
zig fetch --save git+https://github.com/planckeon/nufast.git#v0.4.0
```

To save with a custom dependency name:

```bash
zig fetch --save=nufast git+https://github.com/planckeon/nufast.git
```

### Manual configuration

Add to `build.zig.zon`:

```zig
.dependencies = .{
    .nufast = .{
        .url = "git+https://github.com/planckeon/nufast.git",
        .hash = "...", // Run `zig build` to get the correct hash
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

## Usage

### Vacuum oscillations

```zig
const nufast = @import("nufast");

const params = nufast.VacuumParams.default;
const probs = nufast.vacuumProbability(params, 1300.0, 2.5);
// probs[1][0] = P(νμ → νe)
```

### Matter oscillations

```zig
const matter = nufast.MatterParams{
    .vacuum = nufast.VacuumParams.default,
    .rho = 2.848,  // g/cm³
    .Ye = 0.5,
    .n_newton = 0, // 0-3, higher = more accurate
};
const probs = nufast.matterProbability(matter, 1300.0, 2.5);
```

### Batch API (pre-computed mixing)

```zig
// Pre-compute mixing matrix for repeated calculations
const batch = nufast.VacuumBatch.init(params);
for (energies) |E| {
    const p = batch.probabilityAt(1300.0, E);
}

// Matter batch
const matter_batch = nufast.MatterBatch.init(matter_params);
```

### SIMD (vectorized)

```zig
// f64 SIMD: 4 energies at once
var E_vec: nufast.F64Vec = .{ 1.0, 2.0, 3.0, 4.0 };
const p_vec = nufast.vacuumProbabilitySimd(batch, 1300.0, E_vec);

// f32 SIMD: 8 energies at once (2× throughput)
const batch_f32 = nufast.VacuumBatchF32.fromF64(batch);
var E_f32: nufast.F32Vec = .{ 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 };
const p_f32 = nufast.vacuumProbabilitySimdF32(batch_f32, 1300.0, E_f32);

// SIMD matter
const mat_probs = nufast.matterProbabilitySimd(matter_batch, 1300.0, E_vec);
```

### Anti-neutrino mode

```zig
var params = nufast.MatterParams.default;
params.antineutrino = true;
const probs = nufast.matterProbability(params, 1300.0, 2.5);
```

## Benchmarks

AMD Ryzen, WSL2, 10M iterations.

### Scalar

| Mode | Time | Throughput |
|------|------|------------|
| Vacuum | 42 ns | 24 M/s |
| Matter N=0 | 108 ns | 9 M/s |

### SIMD

| Mode | f64 (4 lanes) | f32 (8 lanes) |
|------|---------------|---------------|
| Vacuum | 44 ns/calc | 21 ns/calc, 48 M/s |
| Matter | 56 ns/calc | 37 ns/calc, 27 M/s |

## Physics

Default parameters: NuFit 5.2 (2022) normal ordering.

The probability matrix `probs[α][β]` gives P(ν_α → ν_β):

```
        e       μ       τ
    ┌───────────────────────┐
e   │ P_ee    P_eμ    P_eτ  │
μ   │ P_μe    P_μμ    P_μτ  │
τ   │ P_τe    P_τμ    P_ττ  │
    └───────────────────────┘
```

`n_newton` controls matter eigenvalue accuracy:
- 0: DMP approximation (~0.1%)
- 1: One Newton iteration (~0.01%)
- 2-3: Machine precision

Values > 3 are clamped to 3.

## Tests

```bash
zig build test
```

## Documentation

Generate API documentation:

```bash
zig build docs
```

The documentation will be generated in the `docs/` folder. Open `docs/index.html` in a browser to view.

22 tests covering:
- Probability conservation (unitarity)
- Batch/SIMD consistency with scalar
- CPT theorem for anti-neutrinos
- Cross-validation against original NuFast
- Edge cases (L→0, high E, zero θ₁₃, etc.)

## License

MIT
