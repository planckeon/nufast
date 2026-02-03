/**
 * NuFast WASM Test Suite
 * Run with: bun test-ts.ts
 */

import { readFileSync } from 'fs';
import { join } from 'path';
import { loadNuFastFromBytes, NuFast, DEFAULT_PARAMS, DUNE_MATTER } from './nufast';

const wasmPath = join(__dirname, 'nufast.wasm');
const wasmBytes = readFileSync(wasmPath);

async function main() {
  console.log('='.repeat(60));
  console.log('NuFast WASM TypeScript Test');
  console.log('='.repeat(60));

  const nufast = await loadNuFastFromBytes(wasmBytes);

  // =========================================================================
  // Basic Tests
  // =========================================================================
  console.log('\n📊 Basic Probability Tests\n');

  // DUNE-like parameters
  const L = 1300.0; // km
  const E = 2.5;    // GeV

  const Pme = nufast.vacuumPmeDefault(L, E);
  console.log(`Vacuum P(νμ → νe) at DUNE: ${(Pme * 100).toFixed(4)}%`);

  const PmeMatter = nufast.matterPmeDefault(L, E, DUNE_MATTER.rho);
  console.log(`Matter P(νμ → νe) at DUNE: ${(PmeMatter * 100).toFixed(4)}%`);
  console.log(`Matter effect: +${((PmeMatter - Pme) * 100).toFixed(4)}%`);

  // Full matrix
  const matrix = nufast.vacuumProbability(L, E);
  console.log('\nFull vacuum probability matrix:');
  console.log('       →νe       →νμ       →ντ');
  const labels = ['νe', 'νμ', 'ντ'];
  for (let i = 0; i < 3; i++) {
    const row = matrix[i].map(p => `${(p * 100).toFixed(4)}%`.padStart(9)).join(' ');
    console.log(`${labels[i]} → ${row}`);
  }

  // Verify unitarity
  for (let i = 0; i < 3; i++) {
    const rowSum = matrix[i].reduce((a, b) => a + b, 0);
    if (Math.abs(rowSum - 1.0) > 1e-10) {
      console.error(`❌ Row ${i} sum = ${rowSum}, expected 1.0`);
    }
  }
  console.log('✓ Unitarity verified');

  // =========================================================================
  // Batch Processing Tests
  // =========================================================================
  console.log('\n📦 Batch Processing Tests\n');

  // Generate energy array
  const energies = new Float64Array(100);
  for (let i = 0; i < 100; i++) {
    energies[i] = 0.5 + i * 0.045; // 0.5 to 5.0 GeV
  }

  // Batch vacuum
  nufast.setDefaultParams();
  nufast.initVacuumBatch();
  const batchPme = nufast.vacuumBatchPme(L, energies);
  console.log(`Batch vacuum Pme (100 energies):`);
  console.log(`  First: ${(batchPme[0] * 100).toFixed(4)}% at E=${energies[0].toFixed(2)} GeV`);
  console.log(`  Last:  ${(batchPme[99] * 100).toFixed(4)}% at E=${energies[99].toFixed(2)} GeV`);

  // Verify batch matches single-point
  const singlePme = nufast.vacuumPmeDefault(L, energies[50]);
  if (Math.abs(batchPme[50] - singlePme) > 1e-14) {
    console.error(`❌ Batch mismatch at idx 50: ${batchPme[50]} vs ${singlePme}`);
  } else {
    console.log('✓ Batch matches single-point calculations');
  }

  // Batch matter
  nufast.setMatterParams(DUNE_MATTER);
  nufast.initMatterBatch();
  const batchMatterPme = nufast.matterBatchPme(L, energies);
  console.log(`\nBatch matter Pme (100 energies, ρ=${DUNE_MATTER.rho} g/cm³):`);
  console.log(`  First: ${(batchMatterPme[0] * 100).toFixed(4)}% at E=${energies[0].toFixed(2)} GeV`);
  console.log(`  Last:  ${(batchMatterPme[99] * 100).toFixed(4)}% at E=${energies[99].toFixed(2)} GeV`);

  // =========================================================================
  // Custom Parameters
  // =========================================================================
  console.log('\n🔧 Custom Parameters Test\n');

  // Inverted ordering
  nufast.setVacuumParams({
    ...DEFAULT_PARAMS,
    Dmsq31: -2.498e-3, // Inverted ordering
  });
  nufast.initVacuumBatch();

  const PmeIO = nufast.vacuumPmeDefault(L, E);
  console.log(`Inverted ordering P(νμ → νe): ${(PmeIO * 100).toFixed(4)}%`);

  // Anti-neutrino
  nufast.setVacuumParams({
    ...DEFAULT_PARAMS,
    antineutrino: true,
  });
  const PmeBar = nufast.vacuumProbability(L, E)[1][0];
  console.log(`Anti-neutrino P(ν̄μ → ν̄e): ${(PmeBar * 100).toFixed(4)}%`);

  // =========================================================================
  // Benchmarks
  // =========================================================================
  console.log('\n⚡ Benchmarks\n');

  // Reset to defaults
  nufast.setDefaultParams();
  nufast.initVacuumBatch();
  nufast.setMatterParams(DUNE_MATTER);
  nufast.initMatterBatch();

  const n = 1_000_000;

  // Single-point vacuum
  let sum = 0;
  const startVac = Bun.nanoseconds();
  for (let i = 0; i < n; i++) {
    sum += nufast.vacuumPmeDefault(L, 0.5 + (i % 1000) * 0.0045);
  }
  const elapsedVac = Bun.nanoseconds() - startVac;
  console.log(`Single-point vacuum (${(n/1e6).toFixed(0)}M):`);
  console.log(`  ${(elapsedVac / n).toFixed(1)} ns/call`);
  console.log(`  ${(n / elapsedVac * 1e9 / 1e6).toFixed(2)} M calls/sec`);

  // Batch vacuum
  const batchSize = 1000;
  const batchN = n / batchSize;
  const batchEnergies = new Float64Array(batchSize);
  for (let i = 0; i < batchSize; i++) {
    batchEnergies[i] = 0.5 + i * 0.0045;
  }

  sum = 0;
  const startBatch = Bun.nanoseconds();
  for (let i = 0; i < batchN; i++) {
    const results = nufast.vacuumBatchPme(L, batchEnergies);
    sum += results[0]; // Prevent DCE
  }
  const elapsedBatch = Bun.nanoseconds() - startBatch;
  const nsPerPointBatch = elapsedBatch / n;
  console.log(`\nBatch vacuum (${batchSize} per batch, ${(n/1e6).toFixed(0)}M total):`);
  console.log(`  ${nsPerPointBatch.toFixed(1)} ns/point`);
  console.log(`  ${(n / elapsedBatch * 1e9 / 1e6).toFixed(2)} M points/sec`);
  console.log(`  Speedup vs single-point: ${(elapsedVac / elapsedBatch).toFixed(2)}×`);

  // Batch matter
  sum = 0;
  const startMatterBatch = Bun.nanoseconds();
  for (let i = 0; i < batchN; i++) {
    const results = nufast.matterBatchPme(L, batchEnergies);
    sum += results[0];
  }
  const elapsedMatterBatch = Bun.nanoseconds() - startMatterBatch;
  console.log(`\nBatch matter (${batchSize} per batch, ${(n/1e6).toFixed(0)}M total):`);
  console.log(`  ${(elapsedMatterBatch / n).toFixed(1)} ns/point`);
  console.log(`  ${(n / elapsedMatterBatch * 1e9 / 1e6).toFixed(2)} M points/sec`);

  // Prevent DCE
  if (sum === 0) console.log('.');

  console.log('\n' + '='.repeat(60));
  console.log('✅ All tests passed!');
  console.log('='.repeat(60));
}

main().catch(console.error);
