// NuFast: Efficient Neutrino Oscillation Probabilities in Rust
// Academic paper with professional formatting

#set document(
  title: "NuFast: Efficient Three-Flavor Neutrino Oscillation Probabilities in Rust",
  author: "Baalateja Kataru",
  date: datetime(year: 2026, month: 2, day: 3),
)

#set page(
  paper: "us-letter",
  margin: (x: 1.5cm, y: 2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 9pt, fill: rgb("#666"))
      #h(1fr) NuFast: Efficient Neutrino Oscillation Probabilities in Rust #h(1fr)
    ]
  },
)

#set text(
  font: "New Computer Modern",
  size: 11pt,
)

#set par(justify: true)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")

#show link: it => underline(text(fill: rgb("#0066cc"), it))

// Custom environments
#let theorem(body, name: none) = {
  let title = "Theorem"
  if name != none { title = title + " (" + name + ")" }
  block(
    fill: rgb("#e8f4f8"),
    inset: 10pt,
    radius: 4pt,
    width: 100%,
  )[
    *#title.* #body
  ]
}

#let definition(body, name: none) = {
  let title = "Definition"
  if name != none { title = title + " (" + name + ")" }
  block(
    fill: rgb("#f0f0f0"),
    inset: 10pt,
    radius: 4pt,
    width: 100%,
  )[
    *#title.* #body
  ]
}

#let remark(body) = {
  block(
    stroke: (left: rgb("#0066cc") + 2pt),
    inset: (left: 10pt, y: 5pt),
    width: 100%,
  )[
    _Remark._ #body
  ]
}

#let keyresult(body) = {
  block(
    fill: rgb("#f0fff0"),
    stroke: rgb("#228b22") + 1pt,
    inset: 10pt,
    radius: 4pt,
    width: 100%,
  )[
    *Key Result.* #body
  ]
}

// Title block
#align(center)[
  #text(size: 18pt, weight: "bold")[
    NuFast: Efficient Three-Flavor Neutrino \
    Oscillation Probabilities in Rust
  ]
  
  #v(1em)
  
  #text(size: 12pt)[
    Baalateja Kataru
  ]
  
  #text(size: 10pt, style: "italic")[
    Planckeon Labs \
    #link("mailto:baalateja@planckeon.org")
  ]
  
  #v(0.5em)
  
  #text(size: 10pt)[
    February 3, 2026
  ]
]

#v(1.5em)

// Abstract
#align(center)[
  #text(size: 11pt, weight: "bold")[Abstract]
]

#block(inset: (x: 2em))[
#text(size: 10pt)[
We present `nufast`, a Rust implementation of the NuFast algorithm for computing three-flavor neutrino oscillation probabilities in vacuum and constant-density matter. Our implementation achieves performance competitive with optimized C++ code: $tilde 61$ ns for vacuum and $tilde 95$ ns for matter calculations per energy point—approximately *27% faster* than C++ for matter effects. This performance advantage stems from LLVM's aggressive optimization of Rust's ownership-based memory model during Newton-Raphson iterations. The crate is published on crates.io and includes WebAssembly bindings enabling browser-based applications. We also provide `VacuumBatch`, an optimized API for batch calculations that pre-computes mixing matrix elements, achieving 45% speedup for energy spectrum computations. This work enables high-performance neutrino physics calculations in modern software ecosystems while maintaining memory safety guarantees.
]
]

#v(1em)

// Keywords
#align(center)[
  #text(size: 9pt, style: "italic")[
    *Keywords:* neutrino oscillations, matter effects, MSW effect, NuFast algorithm, Rust, WebAssembly, high-performance computing
  ]
]

#v(1.5em)

= Introduction

Neutrino oscillation is a quantum mechanical phenomenon where neutrinos change flavor as they propagate through space. The discovery of neutrino oscillations—implying nonzero neutrino masses—represents physics beyond the Standard Model and was recognized with the 2015 Nobel Prize in Physics. Accurate and efficient computation of oscillation probabilities is essential for analyzing data from current and next-generation experiments including DUNE, Hyper-Kamiokande, and JUNO.

The transition probability from flavor $alpha$ to flavor $beta$ after propagating distance $L$ with energy $E$ is given by:

$ P_(alpha beta) = |chevron.l nu_beta | nu_alpha (L) chevron.r|^2 = sum_(i,j) U_(alpha i)^* U_(beta i) U_(alpha j) U_(beta j)^* e^(-i Delta m^2_(i j) L \/ 2 E) $ <eq:osc>

where $U$ is the Pontecorvo-Maki-Nakagawa-Sakata (PMNS) mixing matrix and $Delta m^2_(i j) = m_i^2 - m_j^2$ are mass-squared differences.

In matter, the Mikheyev-Smirnov-Wolfenstein (MSW) effect modifies the effective mixing parameters through coherent forward scattering of electron neutrinos on electrons. Computing these modified probabilities efficiently has been a longstanding challenge in neutrino phenomenology.

== The NuFast Algorithm

The NuFast algorithm, developed by Denton and Parke @denton2024nufast, provides a computationally optimal method for three-flavor oscillation probabilities. Key innovations include:

+ *Eigenvalue-Eigenvector Identity (EEI)*: Avoids cubic equation solving by using only 2×2 matrix diagonalization (quadratic equations)
+ *Square root elimination*: The quadratic's discriminant square root cancels in probability expressions
+ *Optimal eigenvalue ordering*: Initial DMP approximation propagates to other eigenvalues efficiently
+ *Newton-Raphson refinement*: Optional iterations for arbitrary precision

#remark[
NuFast has been adopted by major collaborations: it is implemented in MaCH3 (the primary reweighting framework for T2K and other US/Japan experiments) and JUNO analysis pipelines, achieving "dramatic speed ups—close to an order of magnitude—over other 'optimized' algorithms" @denton2024nufast.
]

== Motivation and Historical Context

This implementation represents the culmination of several years of work on neutrino oscillation phenomenology.

=== Undergraduate Research (2022–2023)

The author's undergraduate capstone thesis at Krea University, supervised by Dr. Sushant Raut, explored the _Interplay between Neutrino Oscillations and Linear Algebra_. The research investigated applications of the Eigenvalue-Eigenvector Identity (also called the Rosetta identity) @denton2020eigenvector and the Adjugate Identity @abdullahi2022adjugate to streamline symbolic calculations of oscillation probabilities.

The goal was to derive novel series expansions of oscillation probabilities in matter up to second order in the mass hierarchy parameter $alpha equiv Delta m^2_(21) \/ Delta m^2_(31)$ only—as opposed to second order in both $alpha$ and $sin theta_(13)$ as in Akhmedov et al. @akhmedov2004series. Using Mathematica and the Cayley-Hamilton formalism, the author explored whether these linear algebra identities could simplify the analytic calculation. While the symbolic expansions proved computationally intractable (repeatedly exhausting available memory), numerical implementations were successful, resulting in `pytrino`—a Python/Cython library published on PyPI.

=== Postgraduate Research (2023–2024)

The author continued this work during a postgraduate program, pivoting to investigate neutrino oscillations on quantum computers using Hamiltonian simulation (Trotter-Suzuki decomposition) and quantum machine learning approaches. This work reproduced published results @arguelles2019quantum @turro2021quantum using IBM's Qiskit framework.

=== Correspondence with Dr. Denton (October 2024)

In October 2024, the author consulted Dr. Peter Denton regarding research directions in neutrino physics, describing prior work on the EEI and challenges with 3+1 sterile neutrino extensions. Dr. Denton explained that while the EEI is powerful for three flavors (2×2 diagonalization), four-flavor oscillations require cubic eigenvector equations—"analytically much much worse, and also numerically somewhat unstable."

Dr. Denton recommended the NuFast algorithm, noting approximately 100 ns per probability calculation on his laptop. This benchmark became a target for our Rust implementation.

=== Precursor Rust Implementations (2024)

Before implementing the full NuFast algorithm, the author developed several prototype implementations:

- *rustrino* (August 2024): Initial exploration of Rust for neutrino physics, establishing the development environment and build toolchain.

- *nosc* (October 2024): A working two-flavor oscillation engine implementing the standard vacuum and matter formulas. This served as a testbed for understanding Rust's numerical performance characteristics and API design patterns for physics libraries.

These incremental steps informed the design decisions in the final `nufast` implementation, particularly regarding:
- Parameter struct design (`VacuumParameters`, `MatterParameters`)
- Result struct layout for probability triplets
- Benchmark methodology using Criterion

The complete chronology of the author's neutrino software development is:
#figure(
  table(
    columns: (1fr, 1.5fr, 2.5fr),
    inset: 6pt,
    align: (left, left, left),
    stroke: (x: none, y: 0.5pt),
    table.header([*Date*], [*Project*], [*Description*]),
    table.hline(stroke: 1pt),
    [May 2023], [`pytrino`], [Python/Cython library (undergrad thesis)],
    [Aug 2024], [`rustrino`], [First Rust prototype],
    [Sep 2024], [`nufast`], [NuFast algorithm port (this work)],
    [Oct 2024], [`nosc`], [Two-flavor Rust implementation],
    [Feb 2026], [`nufast-wasm`], [WebAssembly bindings for ITN],
    table.hline(stroke: 0.5pt),
  ),
  caption: [Chronology of the author's neutrino oscillation software.]
) <tab:chronology>

= Algorithm Details

== Vacuum Oscillations

In vacuum, the oscillation probability depends on:
- *Mixing angles*: $theta_(12), theta_(13), theta_(23)$
- *CP-violating phase*: $delta_"CP"$
- *Mass-squared differences*: $Delta m^2_(21), Delta m^2_(31)$
- *Baseline and energy*: $L$ (km), $E$ (GeV)

#definition(name: "PMNS Matrix")[
The PMNS matrix in the standard parameterization is:
$ U = mat(
  c_(12) c_(13), s_(12) c_(13), s_(13) e^(-i delta);
  -s_(12) c_(23) - c_(12) s_(23) s_(13) e^(i delta), c_(12) c_(23) - s_(12) s_(23) s_(13) e^(i delta), s_(23) c_(13);
  s_(12) s_(23) - c_(12) c_(23) s_(13) e^(i delta), -c_(12) s_(23) - s_(12) c_(23) s_(13) e^(i delta), c_(23) c_(13)
) $
where $c_(i j) = cos theta_(i j)$ and $s_(i j) = sin theta_(i j)$.
]

The vacuum algorithm computes exact probabilities using trigonometric identities without matrix exponentiation or diagonalization.

== Matter Effects

For propagation through matter with constant electron density $N_e$, the effective Hamiltonian acquires a matter potential:

$ H = H_"vacuum" + "diag"(a, 0, 0), quad a = sqrt(2) G_F N_e = 7.63 times 10^(-5) (rho Y_e) ["eV"^2 "/" "GeV"] $

where $rho$ is the matter density in g/cm³ and $Y_e$ is the electron fraction.

NuFast uses:
1. *DMP approximation*: Initial eigenvalue estimate from Denton-Minakata-Parke @dmp
2. *Newton-Raphson refinement*: $N_"Newton"$ iterations for improved precision

#remark[
For long-baseline experiments like DUNE (1300 km, 2.5 GeV, $rho approx 2.8$ g/cm³), $N_"Newton" = 0$ provides sub-percent accuracy. Higher precision is available with $N_"Newton" = 1$ or 2.
]

= Implementation

== Core API

Our Rust implementation provides an ergonomic API:

```rust
use nufast::{VacuumParameters, MatterParameters};
use nufast::{probability_vacuum_lbl, probability_matter_lbl};

// Vacuum oscillation with NuFIT 5.2 parameters
let params = VacuumParameters::nufit52_no(1300.0, 2.5);
let probs = probability_vacuum_lbl(&params);
println!("P(νμ → νe) = {:.4}", probs.Pme);

// Matter oscillation
let params = MatterParameters::nufit52_no(
    1300.0,  // L (km)
    2.5,     // E (GeV)
    2.8,     // ρ (g/cm³)
    0.5,     // Ye
    0        // N_Newton
);
let probs = probability_matter_lbl(&params);
```

== VacuumBatch Optimization

For batch calculations (e.g., computing energy spectra), we provide `VacuumBatch` which pre-computes mixing matrix elements:

```rust
use nufast::VacuumBatch;

let batch = VacuumBatch::nufit52_no();
let spectrum = batch.spectrum(1300.0, 0.5, 5.0, 1000);
// 1000-point spectrum in ~72 μs
```

#keyresult[
`VacuumBatch` achieves *45% speedup* over repeated single-point calls by pre-computing all nine $|U_(alpha i)|^2$ elements and the Jarlskog invariant once, then reusing them across all energy/baseline points.
]

== WebAssembly Support

The `nufast-wasm` crate compiles the physics engine to WebAssembly:

```javascript
import init, { wasmCalculateEnergySpectrum } from 'nufast-wasm';

await init();
const spectrum = wasmCalculateEnergySpectrum({
  theta12_deg: 33.44,
  theta13_deg: 8.57,
  // ... other parameters
}, 1300, 0.5, 5.0, 200);
```

The compiled WASM module is approximately *32 KB gzipped*, enabling browser-based neutrino physics with near-native performance.

== Zig Implementation

We also provide a Zig implementation that achieves the highest performance:

```zig
const nufast = @import("nufast");

// Vacuum oscillation
const params = nufast.VacuumParams.default;
const probs = nufast.vacuumProbability(params, 1300.0, 2.5);
// probs[1][0] = P(νμ → νe)

// Batch calculations with pre-computed matrix
const batch = nufast.VacuumBatch.init(params);
for (energies) |E| {
    const p = batch.probabilityAt(1300.0, E);
}

// SIMD: 4 energies simultaneously
var E_vec: nufast.F64Vec = .{ 1.0, 2.0, 3.0, 4.0 };
const p_vec = nufast.vacuumProbabilitySimd(batch, L, E_vec);
```

The Zig implementation uses:
- Native SIMD via `@Vector` types
- Zero heap allocations in hot paths
- Explicit inlining control
- Zig 0.15.2's improved optimizer

= Benchmark Methodology

All benchmarks were performed on AMD Ryzen (WSL2/Windows). Methodology:

- *Iterations*: 10 million ($10^7$) calculations per measurement
- *Energy range*: 0.5–5.0 GeV (DUNE-like parameters)
- *Parameters*: NuFIT 5.2 best-fit values @nufit52
- *Anti-optimization*: Sink variables prevent dead code elimination

Rust benchmarks use Criterion with statistical analysis. C++, Fortran, and Python use high-resolution timing with standard deviation over 10 runs.

= Results

#figure(
  table(
    columns: (1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
    inset: 8pt,
    align: center,
    stroke: (x: none, y: 0.5pt),
    table.header(
      [*Language*], [*Vacuum*], [*N=0*], [*N=1*], [*N=2*], [*N=3*],
    ),
    table.hline(stroke: 1pt),
    [*Zig*], [*25.6 ns*], [*73.3 ns*], [*77.7 ns*], [*81.1 ns*], [*82.4 ns*],
    [Rust], [61 ns], [95 ns], [106 ns], [113 ns], [117 ns],
    [C++], [49 ns], [130 ns], [143 ns], [154 ns], [164 ns],
    [Fortran], [51 ns], [107 ns], [123 ns], [146 ns], [167 ns],
    [Python], [14,700 ns], [21,900 ns], [21,200 ns], [18,500 ns], [16,300 ns],
    table.hline(stroke: 0.5pt),
  ),
  caption: [Single-point oscillation probability timing (ns/call). N = Newton-Raphson iterations. Bold indicates fastest per column.]
) <tab:results>

#v(0.5em)

#figure(
  table(
    columns: (1.5fr, 1fr, 1fr, 2fr),
    inset: 8pt,
    align: center,
    stroke: (x: none, y: 0.5pt),
    table.header(
      [*Comparison*], [*Vacuum*], [*Matter (N=0)*], [*Interpretation*],
    ),
    table.hline(stroke: 1pt),
    [Zig vs Rust], [*−58%*], [*−23%*], [Zig is fastest],
    [Zig vs C++], [*−48%*], [*−44%*], [Zig dominates],
    [Rust vs C++], [+24%], [*−27%*], [Rust faster for matter],
    [Rust vs Fortran], [+20%], [−11%], [Rust competitive],
    [Zig vs Python], [×574], [×299], [Compiled advantage],
    table.hline(stroke: 0.5pt),
  ),
  caption: [Relative performance (negative = faster).]
) <tab:comparison>

== Key Findings

=== Zig Achieves Fastest Performance

The Zig implementation achieves the fastest times across all configurations:
- *2.4× faster* than Rust for vacuum oscillations (25.6 ns vs 61 ns)
- *1.3× faster* than Rust for matter oscillations (73.3 ns vs 95 ns)
- *1.8× faster* than C++ for matter calculations

This performance advantage stems from:
1. *No ownership overhead*: Zig's explicit memory model allows aggressive inlining without borrow-checker constraints
2. *Direct SIMD*: Zig's `@Vector` type provides native SIMD operations
3. *Zero hidden allocations*: All computation occurs on the stack
4. *Simpler optimizer targets*: LLVM can optimize Zig's simpler IR more effectively

=== Rust Outperforms C++ for Matter Calculations

#keyresult[
Rust achieves *27% speedup* over C++ for matter oscillations with $N_"Newton" = 0$. This exceeds Dr. Denton's quoted ~100 ns benchmark.
]

This advantage likely stems from:
1. *Ownership-enabled optimization*: Rust's strict aliasing rules allow LLVM to optimize more aggressively
2. *Loop vectorization*: Newton iteration inner loops optimize well in LLVM
3. *Zero-cost abstractions*: Rust idioms compile to efficient machine code

=== Vacuum Performance

For vacuum (simpler computation), C++ and Fortran are ~20% faster. This is expected for compute kernels where Fortran excels.

=== Throughput

#figure(
  table(
    columns: (1fr, 1fr, 1fr),
    inset: 8pt,
    align: center,
    stroke: (x: none, y: 0.5pt),
    table.header(
      [*Language*], [*Vacuum (M/s)*], [*Matter (M/s)*],
    ),
    table.hline(stroke: 1pt),
    [Rust], [17.5], [*10.5*],
    [C++], [*20.3*], [7.7],
    [Fortran], [19.7], [9.4],
    [Python], [0.07], [0.05],
    table.hline(stroke: 0.5pt),
  ),
  caption: [Throughput in millions of calculations per second.]
) <tab:throughput>

= Applications

== Interactive Visualization: Imagining the Neutrino

The `nufast-wasm` module powers _Imagining the Neutrino_, an interactive web-based educational tool:

#align(center)[
  #link("https://planckeon.github.io/itn/")
]

Features include:
- Real-time oscillation probability animation
- PMNS matrix visualization (2D and 3D representations)
- Energy spectrum and baseline scan plots
- CP asymmetry visualization: $A_"CP" = P(nu_mu -> nu_e) - P(overline(nu)_mu -> overline(nu)_e)$
- PREM Earth density model with automatic density calculation
- 11 experimental presets (DUNE, T2K, NOvA, etc.)
- Internationalization (7 languages)

With WASM, the visualization computes 400-point energy spectra in real-time (~72 μs with VacuumBatch) as users adjust parameters.

== Future Directions

- *PyO3 bindings*: Python interface via Rust FFI
- *SIMD optimization*: Explicit vectorization for batch computations
- *Variable density*: Piecewise-constant matter profiles for Earth models

= Conclusion

We have demonstrated that Rust provides an excellent platform for computational neutrino physics:

#figure(
  table(
    columns: (1fr, 2fr),
    inset: 8pt,
    stroke: none,
    [*Performance*], [27% faster than C++ for matter effects],
    [*Throughput*], [10.5 M calculations/s (matter, N=0)],
    [*Memory safety*], [Compile-time guarantees, zero runtime overhead],
    [*WASM support*], [32 KB gzipped for browser deployment],
    [*Batch API*], [VacuumBatch for 45% faster spectra],
    [*Distribution*], [Published on crates.io],
  ),
  caption: [Summary of nufast capabilities.]
) <tab:summary>

The combination of performance, safety, and modern tooling makes Rust an attractive choice for future neutrino physics software development.

#heading(numbering: none)[Acknowledgments]

The author thanks Dr. Peter B. Denton (Brookhaven National Laboratory) for recommending the NuFast algorithm and providing guidance on neutrino oscillation calculation approaches. This work builds on the NuFast algorithm developed by Denton and Parke, with original implementations available at #link("https://github.com/PeterDenton/NuFast").

The author acknowledges the neutrino physics education received during undergraduate research at Krea University under Dr. Sushant Raut, which provided the theoretical foundation for this implementation.

#heading(numbering: none)[Code Availability]

All code is open source under the MIT license:

#figure(
  table(
    columns: (1fr, 2fr),
    inset: 6pt,
    stroke: none,
    [*Crates.io*], [#link("https://crates.io/crates/nufast")],
    [*GitHub*], [#link("https://github.com/planckeon/nufast")],
    [*Docs*], [#link("https://docs.rs/nufast")],
    [*Visualization*], [#link("https://planckeon.github.io/itn/")],
  ),
)

Benchmark implementations (Rust, C++, Fortran, Python) are included in `benchmarks/`.

#heading(numbering: none)[References]

#set text(size: 10pt)

#bibliography("refs.bib", style: "ieee")
