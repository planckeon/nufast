# nufast

Fast three-flavor neutrino oscillation probabilities running in WebAssembly.

## Installation

```bash
npm install nufast
# or
bun add nufast
```

## Quick Start

```typescript
import { loadNuFast } from 'nufast';

const nufast = await loadNuFast();

// DUNE-like calculation (L=1300 km, E=2.5 GeV)
const Pme = nufast.vacuumPmeDefault(1300, 2.5);
console.log(`P(νμ → νe) = ${(Pme * 100).toFixed(2)}%`);
// Output: P(νμ → νe) = 5.86%
```

## Features

- **Tiny**: ~13 KB WASM binary
- **Fast**: 5-10 million calculations per second
- **Accurate**: Full double precision, matches original NuFast
- **Complete**: Vacuum + constant-density matter oscillations
- **Batch mode**: 2× faster for processing multiple energies

## API

### Loading

```typescript
// Browser
import { loadNuFast } from 'nufast';
const nufast = await loadNuFast('path/to/nufast.wasm');

// Node.js / Bun
import { loadNuFastFromBytes } from 'nufast';
import { readFileSync } from 'fs';
const bytes = readFileSync('nufast.wasm');
const nufast = await loadNuFastFromBytes(bytes);
```

### Single-Point Calculations

```typescript
// Quick P(νμ → νe) with default parameters
const Pme = nufast.vacuumPmeDefault(L, E);

// Full 3×3 probability matrix
const matrix = nufast.vacuumProbability(L, E);
// matrix[1][0] = P(νμ → νe)
// matrix[0][0] = P(νe → νe)
// etc.

// Matter effects
const PmeMatter = nufast.matterPmeDefault(L, E, rho);
```

### Custom Parameters

```typescript
import { DEFAULT_PARAMS } from 'nufast';

// Inverted mass ordering
nufast.setVacuumParams({
  ...DEFAULT_PARAMS,
  Dmsq31: -2.498e-3,
});

// Anti-neutrino mode
nufast.setVacuumParams({
  ...DEFAULT_PARAMS,
  antineutrino: true,
});

// Matter parameters
nufast.setMatterParams({
  rho: 2.848,    // g/cm³
  Ye: 0.5,       // electron fraction
  nNewton: 0,    // Newton iterations (0-3)
});
```

### Batch Processing

For maximum throughput, use batch mode:

```typescript
const energies = new Float64Array(1000);
for (let i = 0; i < 1000; i++) {
  energies[i] = 0.5 + i * 0.005; // 0.5 to 5.5 GeV
}

// Initialize batch (call once after setting parameters)
nufast.initVacuumBatch();

// Calculate 1000 probabilities at once
const results = nufast.vacuumBatchPme(1300, energies);
// results[i] = P(νμ → νe) at energies[i]
```

## Performance

| Mode | Speed | Throughput |
|------|-------|------------|
| Single-point vacuum | ~80-100 ns | 10-12 M/sec |
| Single-point matter | ~120-170 ns | 6-8 M/sec |
| Batch vacuum | ~50-60 ns | 17-20 M/sec |
| Batch matter | ~100-120 ns | 8-10 M/sec |

## Physics

This is a WebAssembly port of the NuFast algorithm (arXiv:2405.02400).

Parameters use NuFIT 5.2 defaults:
- sin²θ₁₂ = 0.307
- sin²θ₁₃ = 0.0220  
- sin²θ₂₃ = 0.546
- δCP = -0.7π
- Δm²₂₁ = 7.53×10⁻⁵ eV²
- Δm²₃₁ = 2.453×10⁻³ eV² (normal ordering)

## License

MIT
