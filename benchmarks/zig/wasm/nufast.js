// nufast.ts
var ProbIndex;
((ProbIndex2) => {
  ProbIndex2[ProbIndex2["Pee"] = 0] = "Pee";
  ProbIndex2[ProbIndex2["Pem"] = 1] = "Pem";
  ProbIndex2[ProbIndex2["Pet"] = 2] = "Pet";
  ProbIndex2[ProbIndex2["Pme"] = 3] = "Pme";
  ProbIndex2[ProbIndex2["Pmm"] = 4] = "Pmm";
  ProbIndex2[ProbIndex2["Pmt"] = 5] = "Pmt";
  ProbIndex2[ProbIndex2["Pte"] = 6] = "Pte";
  ProbIndex2[ProbIndex2["Ptm"] = 7] = "Ptm";
  ProbIndex2[ProbIndex2["Ptt"] = 8] = "Ptt";
})(ProbIndex ||= {});
var DEFAULT_PARAMS = {
  s12sq: 0.307,
  s13sq: 0.022,
  s23sq: 0.546,
  delta: -0.7 * Math.PI,
  Dmsq21: 0.0000753,
  Dmsq31: 0.002453,
  antineutrino: false
};
var DUNE_MATTER = {
  rho: 2.848,
  Ye: 0.5,
  nNewton: 0
};

class NuFast {
  wasm;
  maxBatchSize;
  constructor(wasmExports) {
    this.wasm = wasmExports;
    this.maxBatchSize = wasmExports.get_max_batch_size();
  }
  setVacuumParams(params) {
    this.wasm.set_vacuum_params(params.s12sq, params.s13sq, params.s23sq, params.delta, params.Dmsq21, params.Dmsq31, params.antineutrino ?? false);
  }
  setDefaultParams() {
    this.wasm.set_default_params();
  }
  setMatterParams(params) {
    this.wasm.set_matter_params(params.rho, params.Ye ?? 0.5, params.nNewton ?? 0, params.antineutrino ?? false);
  }
  vacuumProbability(L, E) {
    const ptr = this.wasm.vacuum_probability(L, E);
    return this.readMatrix(ptr);
  }
  matterProbability(L, E) {
    const ptr = this.wasm.matter_probability(L, E);
    return this.readMatrix(ptr);
  }
  vacuumPmeDefault(L, E) {
    return this.wasm.vacuum_Pme_default(L, E);
  }
  matterPmeDefault(L, E, rho) {
    return this.wasm.matter_Pme_default(L, E, rho);
  }
  initVacuumBatch() {
    this.wasm.init_vacuum_batch();
  }
  initMatterBatch() {
    this.wasm.init_matter_batch();
  }
  vacuumBatchPme(L, energies) {
    const count = Math.min(energies.length, this.maxBatchSize);
    const energiesPtr = this.wasm.get_energies_ptr();
    const energiesView = new Float64Array(this.wasm.memory.buffer, energiesPtr, count);
    for (let i = 0;i < count; i++) {
      energiesView[i] = energies[i];
    }
    const outputPtr = this.wasm.vacuum_batch_Pme(L, count);
    const output = new Float64Array(this.wasm.memory.buffer, outputPtr, count);
    return new Float64Array(output);
  }
  matterBatchPme(L, energies) {
    const count = Math.min(energies.length, this.maxBatchSize);
    const energiesPtr = this.wasm.get_energies_ptr();
    const energiesView = new Float64Array(this.wasm.memory.buffer, energiesPtr, count);
    for (let i = 0;i < count; i++) {
      energiesView[i] = energies[i];
    }
    const outputPtr = this.wasm.matter_batch_Pme(L, count);
    const output = new Float64Array(this.wasm.memory.buffer, outputPtr, count);
    return new Float64Array(output);
  }
  vacuumBatchFull(L, energies) {
    const count = Math.min(energies.length, this.maxBatchSize);
    const energiesPtr = this.wasm.get_energies_ptr();
    const energiesView = new Float64Array(this.wasm.memory.buffer, energiesPtr, count);
    for (let i = 0;i < count; i++) {
      energiesView[i] = energies[i];
    }
    const outputPtr = this.wasm.vacuum_batch_full(L, count);
    const output = new Float64Array(this.wasm.memory.buffer, outputPtr, count * 9);
    const matrices = [];
    for (let i = 0;i < count; i++) {
      const base = i * 9;
      matrices.push([
        [output[base + 0], output[base + 1], output[base + 2]],
        [output[base + 3], output[base + 4], output[base + 5]],
        [output[base + 6], output[base + 7], output[base + 8]]
      ]);
    }
    return matrices;
  }
  getMaxBatchSize() {
    return this.maxBatchSize;
  }
  readMatrix(ptr) {
    const view = new Float64Array(this.wasm.memory.buffer, ptr, 9);
    return [
      [view[0], view[1], view[2]],
      [view[3], view[4], view[5]],
      [view[6], view[7], view[8]]
    ];
  }
}
async function loadNuFast(wasmPath = "nufast.wasm") {
  const response = await fetch(wasmPath);
  const bytes = await response.arrayBuffer();
  const { instance } = await WebAssembly.instantiate(bytes, { env: {} });
  const wasm = instance.exports;
  const nufast = new NuFast(wasm);
  nufast.setDefaultParams();
  nufast.initVacuumBatch();
  return nufast;
}
async function loadNuFastFromBytes(wasmBytes) {
  const { instance } = await WebAssembly.instantiate(wasmBytes, { env: {} });
  const wasm = instance.exports;
  const nufast = new NuFast(wasm);
  nufast.setDefaultParams();
  nufast.initVacuumBatch();
  return nufast;
}
export {
  loadNuFastFromBytes,
  loadNuFast,
  ProbIndex,
  NuFast,
  DUNE_MATTER,
  DEFAULT_PARAMS
};
