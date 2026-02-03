# NuFast Zig → WebAssembly

This directory contains the WebAssembly port of NuFast.

## Status: ✅ Production Ready

The core NuFast Zig implementation compiles cleanly to WASM with:
- Full vacuum and matter oscillation support
- Batch processing for 2× throughput
- TypeScript type definitions
- NPM-ready package structure

## Quick Start

```typescript
import { loadNuFast } from '@nufast/wasm';

const nufast = await loadNuFast();
const Pme = nufast.vacuumPmeDefault(1300, 2.5);
console.log(`P(νμ → νe) = ${(Pme * 100).toFixed(2)}%`);
```

## Build

```bash
# Build both variants
zig build wasm wasm-simd

# Copy to wasm/ directory
cp .zig-cache/o/*/nufast*.wasm wasm/
```

## File Sizes

| Variant | Size | Notes |
|---------|------|-------|
| `nufast.wasm` | ~13.6 KB | Baseline, works everywhere |
| `nufast-simd.wasm` | ~13.4 KB | Uses WASM SIMD128 |

## Performance

Benchmark results (Bun on WSL2, AMD Threadripper):

| Mode | Speed | Throughput |
|------|-------|------------|
| Single-point vacuum | ~80-100 ns | 10-12 M/sec |
| Single-point matter | ~120-170 ns | 6-8 M/sec |
| **Batch vacuum** | ~50-60 ns | **17-20 M/sec** |
| **Batch matter** | ~100-120 ns | **8-10 M/sec** |

Batch mode gives **2× speedup** by amortizing JS↔WASM overhead.

## API Overview

### Loading

```typescript
// Browser
import { loadNuFast } from '@nufast/wasm';
const nufast = await loadNuFast('/path/to/nufast.wasm');

// Node.js / Bun  
import { loadNuFastFromBytes } from '@nufast/wasm';
import { readFileSync } from 'fs';
const nufast = await loadNuFastFromBytes(readFileSync('nufast.wasm'));
```

### Single-Point Calculations

```typescript
// Quick P(νμ → νe)
const Pme = nufast.vacuumPmeDefault(L, E);
const PmeMatter = nufast.matterPmeDefault(L, E, rho);

// Full 3×3 matrix
const matrix = nufast.vacuumProbability(L, E);
// matrix[1][0] = Pme, matrix[0][0] = Pee, etc.
```

### Batch Processing

```typescript
const energies = new Float64Array(1000);
for (let i = 0; i < 1000; i++) {
  energies[i] = 0.5 + i * 0.005;
}

// Initialize batch (once after setting params)
nufast.initVacuumBatch();

// Process 1000 energies at once
const results = nufast.vacuumBatchPme(1300, energies);
```

### Custom Parameters

```typescript
import { DEFAULT_PARAMS } from '@nufast/wasm';

// Inverted ordering
nufast.setVacuumParams({
  ...DEFAULT_PARAMS,
  Dmsq31: -2.498e-3,
});

// Anti-neutrino
nufast.setVacuumParams({
  ...DEFAULT_PARAMS,
  antineutrino: true,
});

// Matter
nufast.setMatterParams({
  rho: 2.848,
  Ye: 0.5,
  nNewton: 0,
});
nufast.initMatterBatch();
```

## Files

| File | Purpose |
|------|---------|
| `nufast.wasm` | Baseline WASM binary |
| `nufast-simd.wasm` | WASM with SIMD128 |
| `nufast.ts` | TypeScript source |
| `nufast.js` | Bundled JavaScript |
| `nufast.d.ts` | Type definitions |
| `package.json` | NPM package config |
| `README.md` | Package documentation |
| `test.html` | Browser test page |
| `test.js` | Node.js test script |
| `test-ts.ts` | TypeScript test suite |

## Architecture

### wasm_exports.zig

Thin wrapper providing C-ABI exports:
- Global state for parameters (simpler JS interop)
- Pre-allocated buffers for batch processing
- Direct single-value accessors

### Batch Buffers

```
g_energies[1024]      - Input energies (written by JS)
g_batch_output[1024]  - Output Pme values
g_batch_matrix[1024×9] - Full matrices (optional)
```

### No libc

Uses `wasm32-freestanding` - pure math, no runtime dependencies.

## Exported Functions

### Parameter Setters
- `set_vacuum_params(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, antineutrino)`
- `set_default_params()`
- `set_matter_params(rho, Ye, n_newton, antineutrino)`

### Single-Point
- `vacuum_probability(L, E) → ptr`
- `matter_probability(L, E) → ptr`
- `vacuum_Pme_default(L, E) → f64`
- `matter_Pme_default(L, E, rho) → f64`

### Batch
- `init_vacuum_batch()`
- `init_matter_batch()`
- `get_energies_ptr() → ptr` (write energies here)
- `get_max_batch_size() → 1024`
- `vacuum_batch_Pme(L, count) → ptr`
- `matter_batch_Pme(L, count) → ptr`
- `vacuum_batch_full(L, count) → ptr`
- `matter_batch_full(L, count) → ptr`

### Memory Access
- `get_result_ptr() → ptr`
- `get_result(idx) → f64`
- `get_batch_output_ptr() → ptr`
- `get_batch_matrix_ptr() → ptr`

## Running Tests

```bash
cd wasm

# Node.js/Bun basic test
bun test.js

# TypeScript test suite with benchmarks
bun test-ts.ts

# Browser (start server, open test.html)
python3 -m http.server 8765
```

## NPM Publishing

```bash
cd wasm
npm publish --access public
```

## Future Work

- [ ] Streaming WASM compilation (`instantiateStreaming`)
- [ ] Web Worker wrapper for off-main-thread
- [ ] SIMD batch functions (process vector lanes in parallel)
- [ ] Deno/Cloudflare Workers support
