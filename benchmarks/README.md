# Cross-Language Benchmarks

This directory contains NuFast implementations and benchmarks in multiple languages:

- **zig/** - Zig implementation (fastest!)
- **rust/** - Rust implementation (this crate, uses Criterion)
- **cpp/** - C++ implementation (original NuFast)
- **fortran/** - Fortran implementation (original NuFast)
- **python/** - Python implementation (original NuFast)

## Running Benchmarks

### Zig
```bash
cd zig
zig build bench  # Scalar benchmark
zig build simd   # SIMD benchmark
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

| Language | Vacuum | Matter N=0 | Matter N=1 | Matter N=2 | Matter N=3 |
|----------|--------|------------|------------|------------|------------|
| **Zig**  | **25.6 ns** | **73.3 ns** | **77.7 ns** | **81.1 ns** | **82.4 ns** |
| Rust     | 61 ns  | 95 ns      | 106 ns     | 113 ns     | 117 ns     |
| C++      | 49 ns  | 130 ns     | 143 ns     | 154 ns     | 164 ns     |
| Fortran  | 51 ns  | 107 ns     | 123 ns     | 146 ns     | 167 ns     |
| Python   | 14,700 ns | 21,900 ns | 21,200 ns | 18,500 ns | 16,300 ns |

### Key Findings

1. **Zig is 2.4× faster than Rust** for vacuum calculations
2. **Zig is 1.3× faster than Rust** for matter calculations
3. **Zig is 1.8× faster than C++** for matter (the complex case)
4. Zig SIMD provides ~10% additional speedup
5. Python is ~600× slower than Zig

### Throughput

| Language | Vacuum | Matter (N=0) |
|----------|--------|--------------|
| **Zig**  | **39.2 M/s** | **13.6 M/s** |
| Rust     | 17.5 M/s | 10.5 M/s   |
| C++      | 20.3 M/s | 7.7 M/s    |
| Fortran  | 19.7 M/s | 9.4 M/s    |
| Python   | 0.07 M/s | 0.05 M/s   |
