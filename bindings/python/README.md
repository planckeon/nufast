# NuFast Python Bindings

Python bindings for the NuFast neutrino oscillation library.

NuFast provides fast and accurate three-flavor neutrino oscillation probabilities
in vacuum and constant-density matter, implementing the algorithm from 
[Denton & Parke (arXiv:2405.02400)](https://arxiv.org/abs/2405.02400).

## Installation

### 1. Build the shared library

From the repository root:

```bash
cd benchmarks/zig
zig build lib -Doptimize=ReleaseFast
```

This creates:
- Linux: `zig-out/lib/libnufast.so`
- macOS: `zig-out/lib/libnufast.dylib`
- Windows: `zig-out/lib/nufast.dll`

### 2. Set up Python

Copy `nufast.py` to your project, or add this directory to your `PYTHONPATH`.

The library will be found automatically if:
- It's in the same directory as `nufast.py`
- It's in `benchmarks/zig/zig-out/lib/` (relative to `nufast.py`)
- You set `NUFAST_LIB` environment variable to the full path

## Quick Start

```python
import nufast

# Vacuum oscillation at DUNE baseline
probs = nufast.vacuum_probability(L=1300, E=2.5)
print(f"P(νμ → νe) = {probs.Pme:.4f}")
print(f"P(νμ → νμ) = {probs.Pmm:.4f}")

# Matter effects (constant density)
probs = nufast.matter_probability(L=1300, E=2.5, rho=2.848)
print(f"P(νμ → νe) with matter = {probs.Pme:.4f}")

# Quick Pme calculation (most common use case)
pme = nufast.vacuum_Pme(L=1300, E=2.5)
print(f"P(νμ → νe) = {pme:.4f}")
```

## API Reference

### Main Functions

#### `vacuum_probability(L, E, *, params=None, antineutrino=False)`

Calculate vacuum oscillation probabilities.

- `L`: Baseline distance (km)
- `E`: Neutrino energy (GeV)
- `params`: `OscillationParams` (default: NuFIT 5.2)
- `antineutrino`: Calculate for ν̄ instead of ν

Returns: `ProbabilityMatrix` with Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt

#### `matter_probability(L, E, *, rho=None, Ye=0.5, n_newton=0, params=None, antineutrino=False)`

Calculate matter-affected oscillation probabilities.

- `rho`: Matter density (g/cm³), default: 2.848 (DUNE-like)
- `Ye`: Electron fraction (default: 0.5)
- `n_newton`: Newton-Raphson iterations (0-3, default: 0)

#### `vacuum_Pme(L, E, *, params=None)` / `matter_Pme(...)`

Quick calculation of just P(νμ → νe).

### Batch Processing

For computing many energies at the same baseline:

```python
energies = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
pme_values = nufast.vacuum_Pme_batch(L=1300, energies=energies)
pme_values = nufast.matter_Pme_batch(L=1300, energies=energies, rho=2.848)
```

### Custom Parameters

```python
params = nufast.OscillationParams(
    s12sq=0.307,      # sin²θ₁₂
    s13sq=0.0220,     # sin²θ₁₃
    s23sq=0.546,      # sin²θ₂₃
    delta=-0.7*π,     # δCP (radians)
    Dmsq21=7.53e-5,   # Δm²₂₁ (eV²)
    Dmsq31=2.453e-3,  # Δm²₃₁ (eV², positive=NO, negative=IO)
)

probs = nufast.vacuum_probability(1300, 2.5, params=params)
```

### Experiment Presets

```python
# Use preset experiment configurations
probs = nufast.experiment_probability("DUNE")
probs = nufast.experiment_probability("T2K")
probs = nufast.experiment_probability("NOvA")

# Available presets:
# - T2K:     L=295 km,  E=0.6 GeV,  ρ=2.6 g/cm³
# - NOvA:    L=810 km,  E=2.0 GeV,  ρ=2.84 g/cm³
# - DUNE:    L=1300 km, E=2.5 GeV,  ρ=2.848 g/cm³
# - Hyper-K: L=295 km,  E=0.6 GeV,  ρ=2.6 g/cm³
# - JUNO:    L=52.5 km, E=4 MeV,    ρ=2.6 g/cm³
```

### Antineutrinos

```python
# Neutrino
probs_nu = nufast.vacuum_probability(1300, 2.5, antineutrino=False)

# Antineutrino (flips sign of δCP and matter potential)
probs_nubar = nufast.vacuum_probability(1300, 2.5, antineutrino=True)
```

## Default Parameters

NuFIT 5.2 (normal ordering):

| Parameter | Value |
|-----------|-------|
| sin²θ₁₂ | 0.307 |
| sin²θ₁₃ | 0.0220 |
| sin²θ₂₃ | 0.546 |
| δCP | -0.7π |
| Δm²₂₁ | 7.53 × 10⁻⁵ eV² |
| Δm²₃₁ | 2.453 × 10⁻³ eV² |

## Testing

```bash
python test_nufast.py
```

## Performance

The Zig library is optimized for high performance:
- SIMD-vectorized batch calculations
- Pre-computed mixing matrices for repeated calculations
- Minimal memory allocations

For best performance with many energy points, use the batch functions.

## License

MIT License - see the main repository.
