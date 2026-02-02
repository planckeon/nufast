# nufast

[![Crates.io](https://img.shields.io/crates/v/nufast.svg)](https://crates.io/crates/nufast)
[![Documentation](https://docs.rs/nufast/badge.svg)](https://docs.rs/nufast)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Fast and accurate three-flavor neutrino oscillation probabilities in vacuum and matter.

Rust port of [NuFast](https://github.com/PeterDenton/NuFast) by Peter Denton.

## Features

- **Vacuum oscillations**: Full 3-flavor PMNS oscillation probabilities
- **Matter effects (MSW)**: Constant-density matter with arbitrary electron fraction
- **CP violation**: Full δ_CP phase support
- **Newton refinement**: Optional iterative improvement for matter eigenvalues
- **Zero dependencies**: Pure Rust, no external crates required
- **Fast**: Optimized for batch calculations

## Installation

```toml
[dependencies]
nufast = "0.2"
```

## Quick Start

```rust
use nufast::{VacuumParameters, MatterParameters, probability_vacuum_lbl, probability_matter_lbl};
use std::f64::consts::PI;

// Vacuum oscillation example (T2K-like)
let vacuum = VacuumParameters {
    s12sq: 0.307,      // sin²θ₁₂
    s13sq: 0.0218,     // sin²θ₁₃
    s23sq: 0.545,      // sin²θ₂₃
    delta: 1.36 * PI,  // δ_CP in radians
    Dmsq21: 7.42e-5,   // Δm²₂₁ in eV²
    Dmsq31: 2.517e-3,  // Δm²₃₁ in eV² (normal ordering)
    L: 295.0,          // baseline in km
    E: 0.6,            // energy in GeV
};

let probs = probability_vacuum_lbl(&vacuum);
println!("P(νμ → νe) = {:.4}", probs[1][0]);

// Matter oscillation example (DUNE-like)
let matter = MatterParameters {
    s12sq: 0.307,
    s13sq: 0.0218,
    s23sq: 0.545,
    delta: 1.36 * PI,
    Dmsq21: 7.42e-5,
    Dmsq31: 2.517e-3,
    L: 1300.0,         // km
    E: 2.5,            // GeV
    rho: 2.848,        // matter density in g/cm³
    Ye: 0.5,           // electron fraction
    N_Newton: 0,       // Newton iterations (0 = DMP approximation)
};

let probs = probability_matter_lbl(&matter);
println!("P(νμ → νe) in matter = {:.4}", probs[1][0]);
```

## Probability Matrix

The returned `[[f64; 3]; 3]` matrix is indexed as `probs[α][β]` = P(ν_α → ν_β):

|       | e (0) | μ (1) | τ (2) |
|-------|-------|-------|-------|
| **e** | P_ee  | P_eμ  | P_eτ  |
| **μ** | P_μe  | P_μμ  | P_μτ  |
| **τ** | P_τe  | P_τμ  | P_ττ  |

## Physics Notes

### NuFit 5.2 Best-Fit Values (2022)

```rust
// Normal Ordering
let s12sq = 0.307;      // sin²θ₁₂
let s13sq = 0.0218;     // sin²θ₁₃  
let s23sq = 0.545;      // sin²θ₂₃
let delta = 1.36 * PI;  // δ_CP
let Dmsq21 = 7.42e-5;   // eV²
let Dmsq31 = 2.517e-3;  // eV² (positive for NO)

// Inverted Ordering: use Dmsq31 = -2.498e-3 eV²
```

### Matter Effect

The MSW (Mikheyev-Smirnov-Wolfenstein) effect modifies oscillation probabilities in matter. The `N_Newton` parameter controls accuracy:
- `N_Newton = 0`: DMP approximation (fast, ~0.1% accuracy)
- `N_Newton = 1`: One Newton iteration (~0.01% accuracy)
- `N_Newton ≥ 2`: Machine precision

### Antineutrinos

For antineutrinos, flip the sign of δ_CP:
```rust
let delta_antinu = -delta;
```

## Performance

Benchmarked on AMD Ryzen 9 5900X:
- Vacuum: ~15 ns per calculation
- Matter (N_Newton=0): ~25 ns per calculation
- Matter (N_Newton=1): ~30 ns per calculation

## References

1. P.B. Denton, [NuFast: Fast and Accurate Neutrino Oscillation Probabilities](https://github.com/PeterDenton/NuFast)
2. P.B. Denton, H. Minakata, S.J. Parke, "Compact Perturbative Expressions for Neutrino Oscillations in Matter" (2016)
3. NuFit 5.2 (2022): [http://www.nu-fit.org](http://www.nu-fit.org)

## License

MIT License - see [LICENSE](LICENSE) for details.
