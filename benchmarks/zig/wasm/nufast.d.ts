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
export declare const enum ProbIndex {
    Pee = 0,
    Pem = 1,
    Pet = 2,
    Pme = 3,
    Pmm = 4,
    Pmt = 5,
    Pte = 6,
    Ptm = 7,
    Ptt = 8
}
/** Raw WASM exports (low-level) */
export interface NuFastWasmExports {
    memory: WebAssembly.Memory;
    set_vacuum_params(s12sq: number, s13sq: number, s23sq: number, delta: number, Dmsq21: number, Dmsq31: number, antineutrino: boolean): void;
    set_default_params(): void;
    set_matter_params(rho: number, Ye: number, n_newton: number, antineutrino: boolean): void;
    vacuum_probability(L: number, E: number): number;
    matter_probability(L: number, E: number): number;
    vacuum_Pme(s12sq: number, s13sq: number, s23sq: number, delta: number, Dmsq21: number, Dmsq31: number, L: number, E: number): number;
    vacuum_Pme_default(L: number, E: number): number;
    matter_Pme_default(L: number, E: number, rho: number): number;
    get_result_ptr(): number;
    get_result(idx: number): number;
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
    [
        number,
        number,
        number
    ],
    [
        number,
        number,
        number
    ],
    [
        number,
        number,
        number
    ]
];
/** NuFIT 5.2 default parameters */
export declare const DEFAULT_PARAMS: VacuumParams;
/** DUNE-like matter parameters */
export declare const DUNE_MATTER: MatterParams;
/**
 * High-level NuFast wrapper with ergonomic TypeScript API
 */
export declare class NuFast {
    private readonly wasm;
    private readonly maxBatchSize;
    constructor(wasmExports: NuFastWasmExports);
    /** Set vacuum oscillation parameters */
    setVacuumParams(params: VacuumParams): void;
    /** Set to NuFIT 5.2 default parameters */
    setDefaultParams(): void;
    /** Set matter parameters (call setVacuumParams first) */
    setMatterParams(params: MatterParams): void;
    /**
     * Calculate vacuum oscillation probability matrix
     * @param L Baseline in km
     * @param E Energy in GeV
     * @returns 3×3 probability matrix
     */
    vacuumProbability(L: number, E: number): ProbabilityMatrix;
    /**
     * Calculate matter oscillation probability matrix
     * @param L Baseline in km
     * @param E Energy in GeV
     * @returns 3×3 probability matrix
     */
    matterProbability(L: number, E: number): ProbabilityMatrix;
    /**
     * Quick P(νμ → νe) with default parameters
     * @param L Baseline in km
     * @param E Energy in GeV
     */
    vacuumPmeDefault(L: number, E: number): number;
    /**
     * Quick P(νμ → νe) in matter with default parameters
     * @param L Baseline in km
     * @param E Energy in GeV
     * @param rho Matter density in g/cm³
     */
    matterPmeDefault(L: number, E: number, rho: number): number;
    /**
     * Initialize vacuum batch calculator
     * Call after setVacuumParams for best batch performance
     */
    initVacuumBatch(): void;
    /**
     * Initialize matter batch calculator
     * Call after setMatterParams for best batch performance
     */
    initMatterBatch(): void;
    /**
     * Calculate P(νμ → νe) for multiple energies (vacuum)
     * @param L Baseline in km
     * @param energies Array of energies in GeV
     * @returns Array of P(νμ → νe) values
     */
    vacuumBatchPme(L: number, energies: Float64Array | number[]): Float64Array;
    /**
     * Calculate P(νμ → νe) for multiple energies (matter)
     * @param L Baseline in km
     * @param energies Array of energies in GeV
     * @returns Array of P(νμ → νe) values
     */
    matterBatchPme(L: number, energies: Float64Array | number[]): Float64Array;
    /**
     * Calculate full probability matrices for multiple energies (vacuum)
     * @param L Baseline in km
     * @param energies Array of energies in GeV
     * @returns Array of 3×3 probability matrices
     */
    vacuumBatchFull(L: number, energies: Float64Array | number[]): ProbabilityMatrix[];
    /** Maximum batch size (1024) */
    getMaxBatchSize(): number;
    private readMatrix;
}
/**
 * Load NuFast WASM module
 * @param wasmPath Path to nufast.wasm file (default: 'nufast.wasm')
 * @returns Initialized NuFast instance
 */
export declare function loadNuFast(wasmPath?: string): Promise<NuFast>;
/**
 * Load NuFast WASM module from bytes (Node.js / Bun)
 * @param wasmBytes WASM binary data
 * @returns Initialized NuFast instance
 */
export declare function loadNuFastFromBytes(wasmBytes: ArrayBuffer): Promise<NuFast>;
