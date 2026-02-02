# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
