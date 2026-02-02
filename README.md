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
- **Fast**: ~61 ns vacuum, ~95 ns matter — **27% faster than C++ for matter!**

## Installation

```toml
[dependencies]
nufast = "0.3"
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

## Performance

Comprehensive benchmarks on AMD Ryzen (WSL2), 10M iterations each:

| Language | Vacuum | Matter N=0 | Matter N=1 | Matter N=2 | Matter N=3 |
|----------|--------|------------|------------|------------|------------|
| **Rust** | 61 ns  | **95 ns**  | 106 ns     | 113 ns     | 117 ns     |
| C++      | 49 ns  | 130 ns     | 143 ns     | 154 ns     | 164 ns     |
| Fortran  | 51 ns  | 107 ns     | 123 ns     | 146 ns     | 167 ns     |
| Python   | 14,700 ns | 21,900 ns | 21,200 ns | 18,500 ns | 16,300 ns |

**Key finding**: Rust is **~27% faster than C++** for matter calculations, likely due to LLVM's better optimization of the Newton iteration loop and Rust's stricter aliasing rules.

### Throughput

| Language | Vacuum | Matter (N=0) |
|----------|--------|--------------|
| Rust     | 17.5 M/s | 10.5 M/s   |
| C++      | 20.3 M/s | 7.7 M/s    |
| Fortran  | 19.7 M/s | 9.4 M/s    |
| Python   | 0.07 M/s | 0.05 M/s   |

### Run Benchmarks

```bash
# Rust (Criterion)
cargo bench

# All languages (see benchmarks/ directory)
cd benchmarks/cpp && g++ -O3 -march=native -o benchmark benchmark.cpp && ./benchmark
cd benchmarks/fortran && gfortran -O3 -march=native -o benchmark benchmark.f90 && ./benchmark
cd benchmarks/python && python benchmark.py
```

See [`benchmarks/README.md`](benchmarks/README.md) for details and [`paper/nufast-benchmark.pdf`](paper/nufast-benchmark.pdf) for the full research paper.

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

Or use the convenience constructors:
```rust
let params = VacuumParameters::nufit52_no(1300.0, 2.5);  // Normal Ordering
let params = VacuumParameters::nufit52_io(1300.0, 2.5);  // Inverted Ordering
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

## Project Structure

```
nufast/
├── src/
│   ├── lib.rs          # Core library
│   └── main.rs         # CLI example
├── benches/
│   └── oscillation.rs  # Criterion benchmarks
├── benchmarks/
│   ├── README.md       # Benchmark documentation
│   ├── rust/           # Rust benchmarks (Criterion)
│   ├── cpp/            # C++ implementation
│   ├── fortran/        # Fortran implementation
│   └── python/         # Python implementation
├── paper/
│   ├── nufast-benchmark.typ  # Typst source
│   └── nufast-benchmark.pdf  # Research paper
├── CHANGELOG.md
└── README.md
```

## References

1. P.B. Denton, "Neutrino Oscillation Probabilities: A Compact Multi-algorithm Approach," [arXiv:2405.02400](https://arxiv.org/abs/2405.02400) (2024)
2. P.B. Denton, [NuFast: Fast and Accurate Neutrino Oscillation Probabilities](https://github.com/PeterDenton/NuFast) - Original C++, Fortran, and Python implementations
3. NuFit 5.2 (2022): [http://www.nu-fit.org](http://www.nu-fit.org)

## Acknowledgments

This Rust implementation was inspired by correspondence with Dr. Peter B. Denton, who recommended the NuFast algorithm as the optimal approach for efficient neutrino oscillation probability calculations.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.
