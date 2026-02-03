# Cross-Language Benchmarks

This directory contains NuFast implementations and benchmarks in multiple languages:

- **zig/** - Zig implementation (fastest, with SIMD + WASM support)
- **rust/** - Rust implementation (this crate, uses Criterion)
- **cpp/** - C++ implementation (original NuFast)
- **fortran/** - Fortran implementation (original NuFast)
- **python/** - Python implementation (original NuFast)

## Running Benchmarks

### Zig

```bash
cd zig
zig build bench  # Scalar benchmark
zig build simd   # SIMD benchmark (f64 and f32 modes)
```

### Zig WASM

```bash
cd zig
zig build wasm wasm-simd
cd wasm && bun test-ts.ts  # Run TypeScript test suite
```

### Rust (Criterion)

```bash
cd .. && cargo bench
```

### C++

```bash
cd cpp
g++ -O3 -march=native -o benchmark benchmark.cpp
./benchmark
```

### Fortran

```bash
cd fortran
gfortran -O3 -march=native -o benchmark benchmark.f90
./benchmark
```

### Python

```bash
cd python
python benchmark.py
```

## Results (AMD Ryzen, WSL2, 10M iterations)

### Scalar Performance

| Language | Vacuum | Matter N=0 | Matter N=1 | Matter N=2 | Matter N=3 |
|----------|--------|------------|------------|------------|------------|
| **Zig**  | **31 ns** | **85 ns** | **98 ns** | **114 ns** | **121 ns** |
| Rust     | 61 ns  | 95 ns      | 106 ns     | 113 ns     | 117 ns     |
| C++      | 49 ns  | 130 ns     | 143 ns     | 154 ns     | 164 ns     |
| Fortran  | 51 ns  | 107 ns     | 123 ns     | 146 ns     | 167 ns     |
| Python   | 14,700 ns | 21,900 ns | 21,200 ns | 18,500 ns | 16,300 ns |

### Zig SIMD Performance

| Mode | Scalar | SIMD f64 (4×) | SIMD f32 (8×) |
|------|--------|---------------|---------------|
| Vacuum | 33 ns, 30 M/s | 30 ns, 33 M/s | **14 ns, 70 M/s** |
| Matter N=0 | 77 ns, 13 M/s | 38 ns, 26 M/s | **26 ns, 39 M/s** |

### Zig WASM Performance

| Mode | Single-point | Batch (1000) | Speedup |
|------|--------------|--------------|---------|
| Vacuum | ~100 ns | ~50 ns/point | **2×** |
| Matter | ~150 ns | ~110 ns/point | **1.4×** |

WASM binaries:
- `nufast.wasm`: 13.6 KB (baseline)
- `nufast-simd.wasm`: 13.4 KB (with SIMD128)

### Key Findings

1. **Zig is 2× faster than Rust** for vacuum calculations (31 ns vs 61 ns)
2. **Zig SIMD f32 achieves 70 M/s** vacuum throughput (2.3× scalar speedup)
3. **Zig SIMD matter achieves 3× speedup** over scalar with f32 mode
4. **Rust is 27% faster than C++** for matter calculations
5. **WASM batch mode achieves 2× speedup** over single-point calls
6. **Python is ~600× slower than Zig**

### Throughput Summary

| Language | Vacuum | Matter (N=0) |
|----------|--------|--------------|
| **Zig SIMD f32** | **70 M/s** | **39 M/s** |
| **Zig SIMD f64** | 33 M/s | 26 M/s |
| Zig scalar | 32 M/s | 12 M/s |
| **WASM batch** | **20 M/s** | **9 M/s** |
| WASM single | 10 M/s | 6 M/s |
| Rust     | 17.5 M/s | 10.5 M/s |
| C++      | 20.3 M/s | 7.7 M/s |
| Fortran  | 19.7 M/s | 9.4 M/s |
| Python   | 0.07 M/s | 0.05 M/s |

## New Features (v0.5.0)

### Experiment Presets (Zig)

Pre-configured parameters for common neutrino experiments:

```zig
const nufast = @import("nufast");

// Use a preset directly
const dune = nufast.experiments.dune;
const probs = nufast.matterProbability(dune.toMatterParams(), dune.L, dune.E);

// Available presets:
// - experiments.t2k      — T2K (295 km, 0.6 GeV)
// - experiments.nova     — NOvA (810 km, 2.0 GeV)  
// - experiments.dune     — DUNE (1300 km, 2.5 GeV)
// - experiments.hyper_k  — Hyper-K (295 km, 0.6 GeV)
// - experiments.juno     — JUNO (52.5 km, 4 MeV)
```

### PREM Earth Model (Zig)

Variable density calculations using the Preliminary Reference Earth Model:

```zig
const nufast = @import("nufast");

// Automatic path integration through Earth layers
const probs = nufast.matterProbabilityPrem(params, 1300.0, 2.5);

// Get average density for a baseline
const avg = nufast.getAverageDensityAlongPath(1300.0);
// avg.rho ≈ 2.8 g/cm³, avg.Ye ≈ 0.49

// Query density at specific depth/radius
const surface = nufast.premDensityAtDepth(0);      // Crust: ~2.6 g/cm³
const deep = nufast.premDensityAtDepth(3000.0);    // Lower mantle: ~4.9 g/cm³
```

Supports baselines up to Earth diameter (~12,742 km) for atmospheric neutrinos.

### C FFI / Python Bindings

Build a shared library for Python/ctypes:

```bash
cd zig
zig build lib  # Produces libnufast.so or nufast.dll
```

Use from Python:

```python
# Option 1: Low-level ctypes (see above)

# Option 2: High-level Python wrapper
from nufast import vacuum_probability, matter_probability, experiment_probability

# Full 3×3 probability matrix
probs = vacuum_probability(L=1300, E=2.5)
print(f"P(νμ → νe) = {probs.Pme:.4f}")

# Matter effects
probs = matter_probability(L=1300, E=2.5, rho=2.848)

# Use experiment presets
probs = experiment_probability("DUNE")
probs = experiment_probability("T2K", antineutrino=True)

# Batch processing (optimized)
energies = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
pme_values = vacuum_Pme_batch(1300, energies)
```

Note: Set `NUFAST_LIB` environment variable to library path if not auto-detected.

### WebAssembly Support

Build WASM with:

```bash
cd zig
zig build wasm        # Baseline (max compatibility)
zig build wasm-simd   # With SIMD128 (faster, modern browsers)
```

Use from TypeScript:

```typescript
import { loadNuFast } from '@nufast/wasm';

const nufast = await loadNuFast();
nufast.setDefaultParams();

// Single point
const Pme = nufast.vacuumPmeDefault(1300, 2.5);

// Batch (2× faster)
const energies = new Float64Array(1000);
nufast.initVacuumBatch();
const results = nufast.vacuumBatchPme(1300, energies);
```

### Anti-Neutrino Mode

Both Rust and Zig now support anti-neutrino calculations:

```rust
// Rust
let mut params = MatterParameters::nufit52_no(1300.0, 2.5);
params.antineutrino = true;
let probs = probability_matter_lbl(&params);
```

```zig
// Zig
var params = nufast.MatterParams.default;
params.antineutrino = true;
const probs = nufast.matterProbability(params, 1300.0, 2.5);
```

For anti-neutrinos:
- δCP sign is flipped
- Matter potential sign is flipped (A → -A)

### SIMD Matter Calculations (Zig)

```zig
const matter_batch = nufast.MatterBatch.init(params);
const probs = nufast.matterProbabilitySimd(matter_batch, L, energies);
```

### f32 Mode (Zig)

For maximum throughput when precision isn't critical:

```zig
const batch_f32 = nufast.VacuumBatchF32.fromF64(batch);
const probs = nufast.vacuumProbabilitySimdF32(batch_f32, L, energies);
// 8 energies at once (vs 4 for f64)
```

### Batch APIs

Pre-compute mixing matrix elements for repeated calculations:

```zig
// Vacuum batch (30-40% faster for energy spectra)
const batch = nufast.VacuumBatch.init(params);

// Matter batch (30-40% faster for constant density)
const matter_batch = nufast.MatterBatch.init(matter_params);
```

## Notes

- Zig requires version 0.15.2+
- All benchmarks use NuFIT 5.2 parameters
- SIMD lane count depends on CPU (4×f64 / 8×f32 on AVX2)
- f32 mode provides ~7 digits of precision (sufficient for most applications)
- WASM batch mode pre-allocates 1024-point buffers
