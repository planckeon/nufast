# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2026-02-03

### Added (Zig WebAssembly)
- **WASM build target** — `zig build wasm` produces 13.6 KB WASM binary
- **WASM SIMD build** — `zig build wasm-simd` with SIMD128 support
- **wasm_exports.zig** — C-ABI wrapper for JS interop
- **Batch processing API** — `vacuum_batch_Pme()`, `matter_batch_Pme()` for 2× throughput
- **TypeScript bindings** — Full type definitions in `nufast.ts`
- **NPM package structure** — Ready for `@nufast/wasm` publication

### WASM Performance
- Single-point vacuum: ~80-100 ns/call, 10-12 M/sec
- Single-point matter: ~120-170 ns/call, 6-8 M/sec
- **Batch vacuum: ~50-60 ns/point, 17-20 M/sec** (2× speedup)
- **Batch matter: ~100-120 ns/point, 8-10 M/sec**

### WASM Features
- Zero dependencies (pure math, no libc)
- Pre-allocated buffers (1024 points max batch)
- Global state for simplified JS↔WASM interop
- `nufast.ts` class with `NuFast.vacuumBatchPme()` and friends

### Files Added
- `benchmarks/zig/WASM.md` — Full documentation
- `benchmarks/zig/src/wasm_exports.zig` — WASM exports
- `benchmarks/zig/wasm/nufast.wasm` — Baseline WASM
- `benchmarks/zig/wasm/nufast-simd.wasm` — SIMD WASM
- `benchmarks/zig/wasm/nufast.ts` — TypeScript API
- `benchmarks/zig/wasm/nufast.js` — Bundled JavaScript
- `benchmarks/zig/wasm/nufast.d.ts` — Type definitions
- `benchmarks/zig/wasm/package.json` — NPM package
- `benchmarks/zig/wasm/test.html` — Interactive browser demo
- `benchmarks/zig/wasm/test-ts.ts` — TypeScript test suite

## [0.4.0] - 2026-02-03

### Added (Rust)
- **Anti-neutrino mode** — `antineutrino: bool` field in `VacuumParameters` and `MatterParameters`
  - Flips sign of δCP and matter potential (A → −A)
  - `VacuumBatch::new()` now accepts `antineutrino` parameter
  - `VacuumBatch::from_params()` for construction from `VacuumParameters`
- New tests for CPT theorem validation and matter asymmetry

### Added (Zig)
- **SIMD matter calculations** — `matterProbabilitySimd()` for vectorized matter oscillations
- **f32 mode** — 8 SIMD lanes (vs 4 for f64), 2× throughput
  - `VacuumBatchF32`, `MatterBatchF32`, `F32Vec`
  - `vacuumProbabilitySimdF32()`, `matterProbabilitySimdF32()`
- **MatterBatch API** — pre-computed mixing for constant-density calculations
- **Anti-neutrino mode** — same physics as Rust implementation

### Performance (Zig SIMD)
- Vacuum f32 SIMD: 21 ns/calc, **48 M/s** (2× faster than scalar)
- Matter f32 SIMD: 37 ns/calc, **27 M/s** (3× faster than scalar)
- Matter f64 SIMD: 56 ns/calc, 18 M/s (2× faster than scalar)

### Documentation
- README completely rewritten with Rust and Zig at parity
- Updated benchmarks/README.md with SIMD results
- Paper updated with anti-neutrino and SIMD sections

## [0.3.1] - 2026-02-03

### Added
- **VacuumBatch API** — 45% faster batch calculations via pre-computed mixing matrix elements
  - `VacuumBatch::new()` pre-computes all 9 Usq elements and Jvac
  - `probability_at(L, E)` computes single point using cached values
  - `spectrum(L, &energies, flavor)` batch computes P(E) at fixed L
  - `baseline_scan(E, &baselines, flavor)` batch computes P(L) at fixed E
- Factory constructors `VacuumBatch::nufit52_no()` and `nufit52_io()`
- New benchmark: `vacuum_batch_spectrum_1000`

### Performance
- Energy spectrum (1000 points): 72µs → **72ns/point** (was 131ns/point)
- Single vacuum probability: ~46ns (unchanged)
- Throughput: 21M calculations/second for vacuum

## [0.3.0] - 2026-02-03

### Added
- Complete cross-language benchmark suite in `benchmarks/` directory
  - C++ implementation (`benchmarks/cpp/benchmark.cpp`)
  - Fortran implementation (`benchmarks/fortran/benchmark.f90`)
  - Python implementation (`benchmarks/python/benchmark.py`)
- Research paper documenting benchmarks (`paper/nufast-benchmark.pdf`)
- Improved CLI example with formatted output

### Changed
- Added `#![allow(non_snake_case)]` to suppress physics naming warnings
- Enhanced lib.rs documentation with performance metrics
- Reorganized project structure for multi-language benchmarks

### Documentation
- Comprehensive README with benchmark tables and project structure
- Added `benchmarks/README.md` with instructions for all languages
- Updated lib.rs rustdoc with features and performance section

## [0.2.1] - 2026-02-03

### Added
- Comprehensive Criterion benchmark suite (`benches/oscillation.rs`)
- Cross-language benchmark comparison (Rust vs C++, Fortran, Python)

### Documentation
- Updated README with detailed performance comparison table
- Added throughput metrics (17.5M vacuum, 10.5M matter calculations/sec)
- Documented key finding: Rust is ~27% faster than C++ for matter calculations

## [0.2.0] - 2026-02-02

### Added
- Initial public release on crates.io
- `VacuumParameters` and `MatterParameters` structs
- `probability_vacuum_lbl()` for vacuum oscillations
- `probability_matter_lbl()` for matter effects (MSW)
- `normalize_probabilities()` utility function
- NuFit 5.2 best-fit constructors for Normal and Inverted ordering
- `no_std` feature flag support
- Comprehensive unit tests

### Performance
- Vacuum calculations: ~61 ns per call
- Matter calculations (N=0): ~95 ns per call

## [0.1.0] - 2026-02-01

### Added
- Initial development version (unpublished)
- Core algorithm port from NuFast by Peter Denton
