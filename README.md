# nufast

[![Crates.io](https://img.shields.io/crates/v/nufast.svg)](https://crates.io/crates/nufast)
[![docs.rs](https://docs.rs/nufast/badge.svg)](https://docs.rs/nufast)
[![Zig](https://img.shields.io/badge/zig-0.15.2-f7a41d?logo=zig)](https://ziglang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Three-flavor neutrino oscillation probabilities in vacuum and constant-density matter.

Implementations in **Rust**, **Zig**, and **WebAssembly**, based on NuFast by P.B. Denton ([arXiv:2405.02400](https://arxiv.org/abs/2405.02400)).

## Overview

This repository provides implementations of the NuFast algorithm for computing neutrino oscillation probabilities. The algorithm achieves sub-microsecond performance while maintaining numerical accuracy suitable for experimental analyses.

The original NuFast is used in production by T2K (via MaCH3) and JUNO.

## Implementations

| Platform | Location | Features |
|----------|----------|----------|
| **Rust** | `src/lib.rs` | crates.io, no_std, VacuumBatch |
| **Zig** | `benchmarks/zig/` | SIMD, f32 mode, MatterBatch |
| **WASM** | `benchmarks/zig/wasm/` | Browser, Node.js, TypeScript, batch API |

All implementations support:
- Vacuum oscillations (exact analytic)
- Matter effects via the DMP approximation with optional Newton refinement
- Anti-neutrino mode (sign-flipped δCP and matter potential)
- Batch APIs for pre-computed mixing matrices

## Usage

### Rust

```toml
[dependencies]
nufast = "0.5"
```

```rust
use nufast::{VacuumParameters, probability_vacuum_lbl};
use std::f64::consts::PI;

let params = VacuumParameters {
    s12sq: 0.307,
    s13sq: 0.0218,
    s23sq: 0.545,
    delta: 1.36 * PI,
    Dmsq21: 7.42e-5,
    Dmsq31: 2.517e-3,
    L: 295.0,
    E: 0.6,
    antineutrino: false,
};

let probs = probability_vacuum_lbl(&params);
// probs[1][0] = P(νμ → νe)
```

### Zig

```zig
const nufast = @import("nufast");

const params = nufast.VacuumParams.default;
const probs = nufast.vacuumProbability(params, 1300.0, 2.5);
// probs[1][0] = P(νμ → νe)

// SIMD: compute 4 energies simultaneously
const batch = nufast.VacuumBatch.init(params);
var energies: nufast.F64Vec = .{ 1.0, 2.0, 3.0, 4.0 };
const p_vec = nufast.vacuumProbabilitySimd(batch, 1300.0, energies);
```

### WebAssembly (TypeScript)

```typescript
import { loadNuFast } from '@nufast/wasm';

const nufast = await loadNuFast();
const Pme = nufast.vacuumPmeDefault(1300, 2.5);
console.log(`P(νμ → νe) = ${(Pme * 100).toFixed(2)}%`);

// Batch mode (2× faster)
const energies = new Float64Array([0.5, 1.0, 1.5, 2.0, 2.5]);
nufast.initVacuumBatch();
const results = nufast.vacuumBatchPme(1300, energies);
```

## Benchmarks

AMD Ryzen, WSL2, 10M iterations per measurement.

### Scalar Performance

| Implementation | Vacuum | Matter (N=0) |
|----------------|--------|--------------|
| Zig            | 42 ns  | 108 ns       |
| Rust           | 61 ns  | 95 ns        |
| C++ (original) | 49 ns  | 130 ns       |
| Fortran (original) | 51 ns | 107 ns    |
| Python (original)  | 14.7 µs | 21.9 µs |

### Zig SIMD

| Mode | f64 (4 lanes) | f32 (8 lanes) |
|------|---------------|---------------|
| Vacuum | 44 ns/calc | 21 ns/calc |
| Matter | 56 ns/calc | 37 ns/calc |

Throughput with f32 SIMD: **~48M vacuum/sec**, **~27M matter/sec**.

### WebAssembly

| Mode | Single-point | Batch (1000) |
|------|-------------|--------------|
| Vacuum | ~100 ns | ~50 ns/point |
| Matter | ~150 ns | ~110 ns/point |

Batch mode gives **2× speedup** by amortizing JS↔WASM overhead.

## Physics

The probability matrix `probs[α][β]` gives P(ν_α → ν_β):

```
        e       μ       τ
    ┌───────────────────────┐
e   │ P_ee    P_eμ    P_eτ  │
μ   │ P_μe    P_μμ    P_μτ  │
τ   │ P_τe    P_τμ    P_ττ  │
    └───────────────────────┘
```

Default parameters use NuFit 5.2 (2022) best-fit values for normal ordering.

For anti-neutrinos, set `antineutrino = true`. This flips the sign of δCP and the matter potential.

The `N_Newton` parameter controls matter eigenvalue accuracy:
- 0: DMP approximation (~0.1% accuracy)
- 1: One Newton iteration (~0.01%)
- ≥2: Machine precision

## Repository Structure

```
src/              Rust implementation
benchmarks/
  zig/            Zig implementation + WASM
    src/          Core physics + WASM exports
    wasm/         WASM binaries, TypeScript, npm package
  cpp/            C++ (original NuFast)
  fortran/        Fortran (original NuFast)
  python/         Python (original NuFast)
paper/            Benchmark methodology and results
```

## Building WASM

```bash
cd benchmarks/zig

# Build WASM (baseline + SIMD)
zig build wasm wasm-simd

# Copy to wasm/ directory
cp .zig-cache/o/*/nufast*.wasm wasm/

# Test
cd wasm && bun test-ts.ts
```

## References

1. P.B. Denton, [arXiv:2405.02400](https://arxiv.org/abs/2405.02400) (2024)
2. [NuFast](https://github.com/PeterDenton/NuFast) — original implementations
3. [NuFit 5.2](http://www.nu-fit.org) (2022)

## License

MIT
