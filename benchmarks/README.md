# Cross-Language Benchmarks

This directory contains NuFast implementations and benchmarks in multiple languages:

- **rust/** - Rust implementation (this crate, uses Criterion)
- **cpp/** - C++ implementation (original NuFast)
- **fortran/** - Fortran implementation (original NuFast)
- **python/** - Python implementation (original NuFast)

## Running Benchmarks

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
| **Rust** | 61 ns  | **95 ns**  | 106 ns     | 113 ns     | 117 ns     |
| C++      | 49 ns  | 130 ns     | 143 ns     | 154 ns     | 164 ns     |
| Fortran  | 51 ns  | 107 ns     | 123 ns     | 146 ns     | 167 ns     |
| Python   | 14,700 ns | 21,900 ns | 21,200 ns | 18,500 ns | 16,300 ns |

### Key Findings

1. **Rust is 27% faster than C++** for matter calculations (the more complex case)
2. C++ and Fortran are ~20% faster for vacuum (simpler calculation)
3. Python is ~240× slower than compiled languages

### Throughput

| Language | Vacuum | Matter (N=0) |
|----------|--------|--------------|
| Rust     | 17.5 M/s | 10.5 M/s   |
| C++      | 20.3 M/s | 7.7 M/s    |
| Fortran  | 19.7 M/s | 9.4 M/s    |
| Python   | 0.07 M/s | 0.05 M/s   |
