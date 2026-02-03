# nufast

[![Crates.io](https://img.shields.io/crates/v/nufast.svg)](https://crates.io/crates/nufast)
[![Documentation](https://docs.rs/nufast/badge.svg)](https://docs.rs/nufast)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Fast and accurate three-flavor neutrino oscillation probabilities in vacuum and matter.

**Implementations: Rust + Zig** — both first-class, production-ready.

Based on [NuFast](https://github.com/PeterDenton/NuFast) by Peter Denton ([arXiv:2405.02400](https://arxiv.org/abs/2405.02400)). The original algorithm is used in production by major neutrino experiments including T2K (via MaCH3) and JUNO.

## Features

| Feature | Rust | Zig |
|---------|------|-----|
| Vacuum oscillations | ✅ | ✅ |
| Matter effects (MSW) | ✅ | ✅ |
| Anti-neutrino mode | ✅ | ✅ |
| Batch API (pre-computed mixing) | ✅ | ✅ |
| SIMD vectorization | — | ✅ (f64 & f32) |
| f32 mode (2× SIMD lanes) | — | ✅ |
| Zero dependencies | ✅ | ✅ |

## Performance

Benchmarks on AMD Ryzen (WSL2), 10M iterations:

### Scalar Performance

| Language | Vacuum | Matter N=0 |
|----------|--------|------------|
| **Zig** | **42 ns** | **108 ns** |
| Rust | 61 ns | 95 ns |
| C++ | 49 ns | 130 ns |
| Fortran | 51 ns | 107 ns |
| Python | 14,700 ns | 21,900 ns |

### Zig SIMD Performance

| Mode | Scalar | SIMD f64 (4×) | SIMD f32 (8×) |
|------|--------|---------------|---------------|
| Vacuum | 42 ns | 44 ns | **21 ns (48 M/s)** |
| Matter | 108 ns | 56 ns | **37 ns (27 M/s)** |

**Key findings:**
- Zig SIMD f32 achieves **48 M/s** vacuum throughput
- Zig SIMD f32 is **3× faster** than scalar for matter
- Rust is **27% faster than C++** for matter calculations

---

## Rust

### Installation

```toml
[dependencies]
nufast = "0.3"
```

### Quick Start

```rust
use nufast::{VacuumParameters, MatterParameters, probability_vacuum_lbl, probability_matter_lbl};
use std::f64::consts::PI;

// Vacuum oscillation (T2K-like)
let vacuum = VacuumParameters {
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

let probs = probability_vacuum_lbl(&vacuum);
println!("P(νμ → νe) = {:.4}", probs[1][0]);

// Matter oscillation (DUNE-like)
let matter = MatterParameters {
    s12sq: 0.307,
    s13sq: 0.0218,
    s23sq: 0.545,
    delta: 1.36 * PI,
    Dmsq21: 7.42e-5,
    Dmsq31: 2.517e-3,
    L: 1300.0,
    E: 2.5,
    rho: 2.848,
    Ye: 0.5,
    N_Newton: 0,
    antineutrino: false,
};

let probs = probability_matter_lbl(&matter);
println!("P(νμ → νe) in matter = {:.4}", probs[1][0]);
```

### Batch API

Pre-compute mixing matrix for repeated calculations:

```rust
use nufast::VacuumBatch;

let batch = VacuumBatch::new(0.307, 0.0218, 0.545, 1.36 * PI, 7.42e-5, 2.517e-3, false);

// ~30% faster for energy spectra
for e in [0.5, 1.0, 2.0, 3.0, 5.0] {
    let probs = batch.probability_at(1300.0, e);
    println!("P(νμ→νe) at {} GeV: {:.4}", e, probs[1][0]);
}
```

### Anti-Neutrino Mode

```rust
let mut params = MatterParameters::nufit52_no(1300.0, 2.5);
params.antineutrino = true;
let probs = probability_matter_lbl(&params);
```

---

## Zig

### Installation

Add to `build.zig.zon`:

```zig
.dependencies = .{
    .nufast = .{
        .url = "https://github.com/planckeon/nufast/archive/refs/heads/main.tar.gz",
        .hash = "...",
    },
},
```

Then in `build.zig`:

```zig
const nufast = b.dependency("nufast", .{ .target = target, .optimize = optimize });
exe.root_module.addImport("nufast", nufast.module("nufast"));
```

### Quick Start

```zig
const nufast = @import("nufast");

// Vacuum oscillation
const params = nufast.VacuumParams.default;
const probs = nufast.vacuumProbability(params, 1300.0, 2.5);
// probs[1][0] = P(νμ → νe)

// Matter oscillation
const matter = nufast.MatterParams{
    .vacuum = params,
    .rho = 2.848,
    .Ye = 0.5,
    .n_newton = 0,
};
const matter_probs = nufast.matterProbability(matter, 1300.0, 2.5);
```

### Batch API

```zig
// Pre-computed mixing (30-40% faster for spectra)
const batch = nufast.VacuumBatch.init(params);
for (energies) |E| {
    const p = batch.probabilityAt(1300.0, E);
}

// Matter batch
const matter_batch = nufast.MatterBatch.init(matter_params);
```

### SIMD

```zig
// f64 SIMD: 4 energies at once
var E_vec: nufast.F64Vec = .{ 1.0, 2.0, 3.0, 4.0 };
const p_vec = nufast.vacuumProbabilitySimd(batch, 1300.0, E_vec);

// f32 SIMD: 8 energies at once (2× throughput)
const batch_f32 = nufast.VacuumBatchF32.fromF64(batch);
var E_f32: nufast.F32Vec = .{ 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 };
const p_f32 = nufast.vacuumProbabilitySimdF32(batch_f32, 1300.0, E_f32);

// SIMD matter
const matter_probs = nufast.matterProbabilitySimd(matter_batch, 1300.0, E_vec);
```

### Anti-Neutrino Mode

```zig
var params = nufast.MatterParams.default;
params.antineutrino = true;
const probs = nufast.matterProbability(params, 1300.0, 2.5);
```

---

## Probability Matrix

The returned `[[f64; 3]; 3]` matrix is indexed as `probs[α][β]` = P(ν_α → ν_β):

|       | e (0) | μ (1) | τ (2) |
|-------|-------|-------|-------|
| **e** | P_ee  | P_eμ  | P_eτ  |
| **μ** | P_μe  | P_μμ  | P_μτ  |
| **τ** | P_τe  | P_τμ  | P_ττ  |

## Physics Notes

### NuFit 5.2 Best-Fit Values (2022)

| Parameter | Normal Ordering | Inverted Ordering |
|-----------|-----------------|-------------------|
| sin²θ₁₂ | 0.307 | 0.307 |
| sin²θ₁₃ | 0.02203 | 0.02219 |
| sin²θ₂₃ | 0.546 | 0.539 |
| δ_CP | 1.36π | 1.56π |
| Δm²₂₁ | 7.42×10⁻⁵ eV² | 7.42×10⁻⁵ eV² |
| Δm²₃₁ | 2.517×10⁻³ eV² | −2.498×10⁻³ eV² |

### Matter Effect (MSW)

The `N_Newton` parameter controls accuracy:
- `N_Newton = 0`: DMP approximation (fast, ~0.1% accuracy)
- `N_Newton = 1`: One Newton iteration (~0.01% accuracy)
- `N_Newton ≥ 2`: Machine precision

### Anti-Neutrino Physics

For anti-neutrinos, two sign flips are applied:
1. **CP phase**: δ → −δ
2. **Matter potential**: A → −A

This enables direct CP asymmetry calculation: A_CP = P(νμ→νe) − P(ν̄μ→ν̄e)

## Project Structure

```
nufast/
├── src/                    # Rust implementation
│   ├── lib.rs
│   └── main.rs
├── benchmarks/
│   ├── zig/                # Zig implementation (SIMD, f32)
│   │   ├── src/nufast.zig
│   │   └── README.md
│   ├── cpp/                # C++ (original NuFast)
│   ├── fortran/            # Fortran (original NuFast)
│   └── python/             # Python (original NuFast)
├── paper/
│   └── nufast-benchmark.pdf
└── README.md
```

## Running Benchmarks

```bash
# Rust
cargo bench

# Zig
cd benchmarks/zig
zig build bench   # Scalar
zig build simd    # SIMD (f64 & f32)

# C++
cd benchmarks/cpp && g++ -O3 -march=native -o bench benchmark.cpp && ./bench
```

See [`benchmarks/README.md`](benchmarks/README.md) for full details.

## References

1. P.B. Denton, "Neutrino Oscillation Probabilities: A Compact Multi-algorithm Approach," [arXiv:2405.02400](https://arxiv.org/abs/2405.02400) (2024)
2. [NuFast](https://github.com/PeterDenton/NuFast) — Original C++/Fortran/Python implementations
3. [NuFit 5.2](http://www.nu-fit.org) (2022)

## Acknowledgments

This implementation was inspired by correspondence with Dr. Peter B. Denton, who recommended NuFast as the optimal algorithm for efficient neutrino oscillation probability calculations.

## License

MIT
