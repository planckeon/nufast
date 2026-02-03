/**
 * NuFast WASM - Three-flavor neutrino oscillation probabilities
 * 
 * TypeScript type definitions for the NuFast WebAssembly module.
 * 
 * @example
 * ```typescript
 * import { loadNuFast, NuFast } from './nufast';
 * 
 * const nufast = await loadNuFast();
 * nufast.setDefaultParams();
 * const Pme = nufast.vacuumPmeDefault(1300, 2.5);
 * console.log(`P(νμ → νe) = ${(Pme * 100).toFixed(4)}%`);
 * ```
 */

/** Probability matrix indices */
export const enum ProbIndex {
  Pee = 0, Pem = 1, Pet = 2,
  Pme = 3, Pmm = 4, Pmt = 5,
  Pte = 6, Ptm = 7, Ptt = 8,
}

/** Raw WASM exports (low-level) */
export interface NuFastWasmExports {
  memory: WebAssembly.Memory;

  // Parameter setters
  set_vacuum_params(
    s12sq: number, s13sq: number, s23sq: number,
    delta: number, Dmsq21: number, Dmsq31: number,
    antineutrino: boolean
  ): void;
  set_default_params(): void;
  set_matter_params(rho: number, Ye: number, n_newton: number, antineutrino: boolean): void;

  // Single-point calculations
  vacuum_probability(L: number, E: number): number;
  matter_probability(L: number, E: number): number;
  vacuum_Pme(
    s12sq: number, s13sq: number, s23sq: number,
    delta: number, Dmsq21: number, Dmsq31: number,
    L: number, E: number
  ): number;
  vacuum_Pme_default(L: number, E: number): number;
  matter_Pme_default(L: number, E: number, rho: number): number;

  // Result access
  get_result_ptr(): number;
  get_result(idx: number): number;

  // Batch processing
  init_vacuum_batch(): void;
  init_matter_batch(): void;
  get_energies_ptr(): number;
  get_batch_output_ptr(): number;
  get_batch_matrix_ptr(): number;
  get_max_batch_size(): number;
  vacuum_batch_Pme(L: number, count: number): number;
  matter_batch_Pme(L: number, count: number): number;
  vacuum_batch_full(L: number, count: number): number;
  matter_batch_full(L: number, count: number): number;
}

/** Vacuum oscillation parameters */
export interface VacuumParams {
  /** sin²θ₁₂ (default: 0.307) */
  s12sq: number;
  /** sin²θ₁₃ (default: 0.0220) */
  s13sq: number;
  /** sin²θ₂₃ (default: 0.546) */
  s23sq: number;
  /** CP phase δ in radians (default: -0.7π) */
  delta: number;
  /** Δm²₂₁ in eV² (default: 7.53e-5) */
  Dmsq21: number;
  /** Δm²₃₁ in eV² (default: 2.453e-3 for NO, negative for IO) */
  Dmsq31: number;
  /** Anti-neutrino mode (default: false) */
  antineutrino?: boolean;
}

/** Matter oscillation parameters */
export interface MatterParams {
  /** Matter density in g/cm³ */
  rho: number;
  /** Electron fraction (default: 0.5) */
  Ye?: number;
  /** Newton-Raphson iterations 0-3 (default: 0) */
  nNewton?: number;
  /** Anti-neutrino mode (default: false) */
  antineutrino?: boolean;
}

/** 3×3 probability matrix [from][to] */
export type ProbabilityMatrix = [
  [number, number, number],  // νe → [νe, νμ, ντ]
  [number, number, number],  // νμ → [νe, νμ, ντ]
  [number, number, number],  // ντ → [νe, νμ, ντ]
];

/** NuFIT 5.2 default parameters */
export const DEFAULT_PARAMS: VacuumParams = {
  s12sq: 0.307,
  s13sq: 0.0220,
  s23sq: 0.546,
  delta: -0.7 * Math.PI,
  Dmsq21: 7.53e-5,
  Dmsq31: 2.453e-3,
  antineutrino: false,
};

/** DUNE-like matter parameters */
export const DUNE_MATTER: MatterParams = {
  rho: 2.848,
  Ye: 0.5,
  nNewton: 0,
};

/**
 * High-level NuFast wrapper with ergonomic TypeScript API
 */
export class NuFast {
  private readonly wasm: NuFastWasmExports;
  private readonly maxBatchSize: number;

  constructor(wasmExports: NuFastWasmExports) {
    this.wasm = wasmExports;
    this.maxBatchSize = wasmExports.get_max_batch_size();
  }

  // =========================================================================
  // Parameter Configuration
  // =========================================================================

  /** Set vacuum oscillation parameters */
  setVacuumParams(params: VacuumParams): void {
    this.wasm.set_vacuum_params(
      params.s12sq,
      params.s13sq,
      params.s23sq,
      params.delta,
      params.Dmsq21,
      params.Dmsq31,
      params.antineutrino ?? false
    );
  }

  /** Set to NuFIT 5.2 default parameters */
  setDefaultParams(): void {
    this.wasm.set_default_params();
  }

  /** Set matter parameters (call setVacuumParams first) */
  setMatterParams(params: MatterParams): void {
    this.wasm.set_matter_params(
      params.rho,
      params.Ye ?? 0.5,
      params.nNewton ?? 0,
      params.antineutrino ?? false
    );
  }

  // =========================================================================
  // Single-Point Calculations
  // =========================================================================

  /**
   * Calculate vacuum oscillation probability matrix
   * @param L Baseline in km
   * @param E Energy in GeV
   * @returns 3×3 probability matrix
   */
  vacuumProbability(L: number, E: number): ProbabilityMatrix {
    const ptr = this.wasm.vacuum_probability(L, E);
    return this.readMatrix(ptr);
  }

  /**
   * Calculate matter oscillation probability matrix
   * @param L Baseline in km
   * @param E Energy in GeV
   * @returns 3×3 probability matrix
   */
  matterProbability(L: number, E: number): ProbabilityMatrix {
    const ptr = this.wasm.matter_probability(L, E);
    return this.readMatrix(ptr);
  }

  /**
   * Quick P(νμ → νe) with default parameters
   * @param L Baseline in km
   * @param E Energy in GeV
   */
  vacuumPmeDefault(L: number, E: number): number {
    return this.wasm.vacuum_Pme_default(L, E);
  }

  /**
   * Quick P(νμ → νe) in matter with default parameters
   * @param L Baseline in km
   * @param E Energy in GeV
   * @param rho Matter density in g/cm³
   */
  matterPmeDefault(L: number, E: number, rho: number): number {
    return this.wasm.matter_Pme_default(L, E, rho);
  }

  // =========================================================================
  // Batch Processing (for maximum throughput)
  // =========================================================================

  /**
   * Initialize vacuum batch calculator
   * Call after setVacuumParams for best batch performance
   */
  initVacuumBatch(): void {
    this.wasm.init_vacuum_batch();
  }

  /**
   * Initialize matter batch calculator
   * Call after setMatterParams for best batch performance
   */
  initMatterBatch(): void {
    this.wasm.init_matter_batch();
  }

  /**
   * Calculate P(νμ → νe) for multiple energies (vacuum)
   * @param L Baseline in km
   * @param energies Array of energies in GeV
   * @returns Array of P(νμ → νe) values
   */
  vacuumBatchPme(L: number, energies: Float64Array | number[]): Float64Array {
    const count = Math.min(energies.length, this.maxBatchSize);
    
    // Write energies to WASM memory
    const energiesPtr = this.wasm.get_energies_ptr();
    const energiesView = new Float64Array(this.wasm.memory.buffer, energiesPtr, count);
    for (let i = 0; i < count; i++) {
      energiesView[i] = energies[i];
    }
    
    // Calculate
    const outputPtr = this.wasm.vacuum_batch_Pme(L, count);
    
    // Read results (copy to avoid issues with buffer detachment)
    const output = new Float64Array(this.wasm.memory.buffer, outputPtr, count);
    return new Float64Array(output);
  }

  /**
   * Calculate P(νμ → νe) for multiple energies (matter)
   * @param L Baseline in km
   * @param energies Array of energies in GeV
   * @returns Array of P(νμ → νe) values
   */
  matterBatchPme(L: number, energies: Float64Array | number[]): Float64Array {
    const count = Math.min(energies.length, this.maxBatchSize);
    
    const energiesPtr = this.wasm.get_energies_ptr();
    const energiesView = new Float64Array(this.wasm.memory.buffer, energiesPtr, count);
    for (let i = 0; i < count; i++) {
      energiesView[i] = energies[i];
    }
    
    const outputPtr = this.wasm.matter_batch_Pme(L, count);
    const output = new Float64Array(this.wasm.memory.buffer, outputPtr, count);
    return new Float64Array(output);
  }

  /**
   * Calculate full probability matrices for multiple energies (vacuum)
   * @param L Baseline in km
   * @param energies Array of energies in GeV
   * @returns Array of 3×3 probability matrices
   */
  vacuumBatchFull(L: number, energies: Float64Array | number[]): ProbabilityMatrix[] {
    const count = Math.min(energies.length, this.maxBatchSize);
    
    const energiesPtr = this.wasm.get_energies_ptr();
    const energiesView = new Float64Array(this.wasm.memory.buffer, energiesPtr, count);
    for (let i = 0; i < count; i++) {
      energiesView[i] = energies[i];
    }
    
    const outputPtr = this.wasm.vacuum_batch_full(L, count);
    const output = new Float64Array(this.wasm.memory.buffer, outputPtr, count * 9);
    
    const matrices: ProbabilityMatrix[] = [];
    for (let i = 0; i < count; i++) {
      const base = i * 9;
      matrices.push([
        [output[base + 0], output[base + 1], output[base + 2]],
        [output[base + 3], output[base + 4], output[base + 5]],
        [output[base + 6], output[base + 7], output[base + 8]],
      ]);
    }
    return matrices;
  }

  // =========================================================================
  // Utilities
  // =========================================================================

  /** Maximum batch size (1024) */
  getMaxBatchSize(): number {
    return this.maxBatchSize;
  }

  private readMatrix(ptr: number): ProbabilityMatrix {
    const view = new Float64Array(this.wasm.memory.buffer, ptr, 9);
    return [
      [view[0], view[1], view[2]],
      [view[3], view[4], view[5]],
      [view[6], view[7], view[8]],
    ];
  }
}

/**
 * Load NuFast WASM module
 * @param wasmPath Path to nufast.wasm file (default: 'nufast.wasm')
 * @returns Initialized NuFast instance
 */
export async function loadNuFast(wasmPath = 'nufast.wasm'): Promise<NuFast> {
  const response = await fetch(wasmPath);
  const bytes = await response.arrayBuffer();
  const { instance } = await WebAssembly.instantiate(bytes, { env: {} });
  const wasm = instance.exports as unknown as NuFastWasmExports;
  
  const nufast = new NuFast(wasm);
  nufast.setDefaultParams();
  nufast.initVacuumBatch();
  
  return nufast;
}

/**
 * Load NuFast WASM module from bytes (Node.js / Bun)
 * @param wasmBytes WASM binary data
 * @returns Initialized NuFast instance
 */
export async function loadNuFastFromBytes(wasmBytes: ArrayBuffer): Promise<NuFast> {
  const { instance } = await WebAssembly.instantiate(wasmBytes, { env: {} });
  const wasm = instance.exports as unknown as NuFastWasmExports;
  
  const nufast = new NuFast(wasm);
  nufast.setDefaultParams();
  nufast.initVacuumBatch();
  
  return nufast;
}
