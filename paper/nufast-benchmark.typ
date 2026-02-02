#set document(
  title: "NuFast: Efficient Neutrino Oscillation Probabilities in Rust",
  author: "Baalateja Kataru",
  date: datetime.today(),
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 2.5cm),
)

#set text(
  font: "New Computer Modern",
  size: 11pt,
)

#set heading(numbering: "1.1")
#set par(justify: true)

#show link: it => underline(text(fill: blue, it))

#align(center)[
  #text(size: 20pt, weight: "bold")[
    NuFast: Efficient Three-Flavor Neutrino Oscillation Probabilities in Rust
  ]
  
  #v(0.5em)
  
  #text(size: 12pt)[Baalateja Kataru]
  
  #text(size: 10pt, style: "italic")[
    Planckeon Labs \
    #link("mailto:baalateja@planckeon.org")
  ]
  
  #v(0.5em)
  
  #text(size: 10pt)[February 2026]
]

#v(1em)

#align(center)[
  #box(width: 85%)[
    *Abstract.* We present `nufast`, a Rust implementation of the NuFast algorithm
    for computing three-flavor neutrino oscillation probabilities in vacuum and
    matter. Our implementation achieves performance competitive with optimized
    C++ code, with ~61 ns for vacuum and ~95 ns for matter calculations per
    energy point—approximately 27% faster than C++ for matter effects. The crate
    is published on crates.io and includes WebAssembly bindings for browser-based
    applications. This work enables high-performance neutrino physics calculations
    in modern software ecosystems.
  ]
]

#v(1em)

= Introduction

Neutrino oscillation is a quantum mechanical phenomenon where neutrinos change
flavor as they propagate through space. Accurate and efficient computation of
oscillation probabilities is essential for analyzing data from neutrino
experiments such as DUNE, Hyper-Kamiokande, and JUNO.

The NuFast algorithm, developed by Denton and Parke [1], provides a
computationally efficient method for calculating three-flavor oscillation
probabilities including matter effects. The original implementation was
provided in C++, Python, and Fortran.

We present `nufast`, a Rust port of this algorithm, designed to leverage Rust's
memory safety guarantees and zero-cost abstractions while achieving performance
comparable to or better than the original implementations.

== Motivation

This project originated from correspondence with Dr. Peter Denton in 2024
regarding research opportunities in neutrino physics. During our discussion of
computational methods for oscillation probability calculations, Dr. Denton
recommended the NuFast algorithm as "the optimal algorithm" for this purpose,
pointing to the preprint arXiv:2405.02400 and the reference implementations
at #link("https://github.com/PeterDenton/NuFast").

The author's background in computational physics—specifically neutrino
oscillation phenomenology during undergraduate research—provided the foundation
for understanding and implementing this algorithm. The decision to port NuFast
to Rust was motivated by several factors:

1. *Performance exploration*: Testing whether Rust's ownership model and LLVM
   backend could match or exceed traditional numerical computing languages
2. *Modern tooling*: Enabling integration with Rust-based physics tools and
   WebAssembly deployment for browser applications
3. *Educational visualization*: Creating "Imagining the Neutrino," an
   interactive web-based tool for teaching neutrino oscillation concepts

= Algorithm Overview

The NuFast algorithm computes the full 3×3 oscillation probability matrix
$P_(alpha beta)$ for neutrino flavor transitions $nu_alpha -> nu_beta$ where
$alpha, beta in {e, mu, tau}$.

== Vacuum Oscillations

In vacuum, the oscillation probability depends on:
- Mixing angles: $theta_(12), theta_(13), theta_(23)$
- CP-violating phase: $delta$
- Mass-squared differences: $Delta m^2_(21), Delta m^2_(31)$
- Baseline $L$ and energy $E$

The algorithm computes exact probabilities using trigonometric identities
without matrix diagonalization.

== Matter Effects

For propagation through matter with constant density $rho$ and electron
fraction $Y_e$, the Mikheyev-Smirnov-Wolfenstein (MSW) effect modifies the
effective mixing parameters. NuFast uses:

1. An initial estimate from the DMP (Denton-Minakata-Parke) approximation
2. Optional Newton-Raphson iterations to improve precision

The number of Newton iterations $N_"Newton"$ controls the trade-off between
speed and precision. For most long-baseline experiments, $N_"Newton" = 0$ or $1$
is sufficient.

= Implementation

Our Rust implementation provides a simple API:

```rust
use nufast::{VacuumParameters, MatterParameters};
use nufast::{probability_vacuum_lbl, probability_matter_lbl};

// Vacuum oscillation
let params = VacuumParameters::nufit52_no(1300.0, 2.5);
let probs = probability_vacuum_lbl(&params);
println!("P(μ→e) = {}", probs.Pme);

// Matter oscillation
let params = MatterParameters::nufit52_no(1300.0, 2.5, 3.0, 0.5, 0);
let probs = probability_matter_lbl(&params);
```

The crate is published on crates.io at #link("https://crates.io/crates/nufast").

== WebAssembly Support

For web applications, we provide `nufast-wasm`, which compiles the core physics
engine to WebAssembly. This enables browser-based neutrino physics tools with
near-native performance. The compiled WASM module is approximately 32 KB gzipped.

= Benchmark Methodology

All benchmarks were performed on a system with AMD Ryzen CPU running WSL2 on
Windows. Each benchmark:

1. Iterates over 10 million (10⁷) calculations
2. Varies energy from 0.5–5.0 GeV (DUNE-like range)
3. Uses standard oscillation parameters (NuFIT 5.2)
4. Includes a sink variable to prevent compiler optimization

The Rust benchmarks used the Criterion library with statistical analysis.
C++, Fortran, and Python benchmarks used high-resolution timing.

= Results

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    inset: 8pt,
    align: horizon,
    table.header(
      [*Language*], [*Vacuum*], [*N=0*], [*N=1*], [*N=2*], [*N=3*],
    ),
    [*Rust*], [61 ns], [95 ns], [106 ns], [113 ns], [117 ns],
    [C++], [49 ns], [130 ns], [143 ns], [154 ns], [164 ns],
    [Fortran], [51 ns], [107 ns], [123 ns], [146 ns], [167 ns],
    [Python], [14,700 ns], [21,900 ns], [21,200 ns], [18,500 ns], [16,300 ns],
  ),
  caption: [Single-point oscillation probability computation times (ns per call). N refers to the number of Newton-Raphson iterations.]
) <tab:results>

== Key Findings

#figure(
  table(
    columns: (auto, auto, auto, auto),
    inset: 8pt,
    align: horizon,
    table.header(
      [*Comparison*], [*Vacuum*], [*Matter (N=0)*], [*Notes*],
    ),
    [Rust vs C++], [+24%], [*−27%*], [Rust faster for matter],
    [Rust vs Fortran], [+20%], [−11%], [Rust faster for matter],
    [Rust vs Python], [×241], [×230], [Compiled vs interpreted],
  ),
  caption: [Relative performance (negative = Rust is faster). Rust shows a significant advantage for matter calculations.]
) <tab:comparison>

=== Rust is Faster for Matter Calculations

The most surprising result is that Rust outperforms C++ by *27%* for matter
oscillations. This is likely due to:

1. *Better loop optimization*: LLVM's optimization of the Newton iteration
2. *Stricter aliasing rules*: Rust's ownership model enables more aggressive optimization
3. *Modern code patterns*: The Rust implementation uses idiomatic patterns that optimize well

=== Vacuum Performance

For vacuum calculations, C++ and Fortran are ~20% faster than Rust. This is a
smaller and simpler computation where traditional numerical languages have
an edge.

=== Python Performance

As expected, Python is approximately 240× slower than compiled languages.
However, for interactive exploration or small-scale calculations, this
overhead is acceptable.

== Throughput

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 8pt,
    align: horizon,
    table.header(
      [*Language*], [*Vacuum*], [*Matter (N=0)*],
    ),
    [Rust], [17.5 M/s], [10.5 M/s],
    [C++], [20.3 M/s], [7.7 M/s],
    [Fortran], [19.7 M/s], [9.4 M/s],
    [Python], [0.07 M/s], [0.05 M/s],
  ),
  caption: [Throughput in millions of probability calculations per second.]
) <tab:throughput>

= Applications

== Interactive Visualization

The `nufast-wasm` module powers "Imagining the Neutrino," an interactive
web-based visualization of neutrino oscillations:
#link("https://planckeon.github.io/itn/")

With WASM, the visualization computes 400-point energy spectra in real-time
as users adjust parameters.

== Integration with Analysis Frameworks

The crate can be integrated with Rust-based physics analysis pipelines, or
called from Python via PyO3 bindings (future work).

= Conclusion

We have demonstrated that Rust provides a viable and performant platform for
computational neutrino physics. The `nufast` crate achieves:

- *27% speedup* over C++ for matter calculations
- *Competitive performance* for vacuum calculations
- *Memory safety* guarantees without runtime overhead
- *WebAssembly support* for browser applications
- *Published on crates.io* for easy integration

The combination of performance, safety, and modern tooling makes Rust an
attractive choice for future neutrino physics software development.

#heading(numbering: none)[Acknowledgments]

The author thanks Dr. Peter B. Denton (Brookhaven National Laboratory) for
recommending the NuFast algorithm and providing guidance on optimal approaches
to neutrino oscillation calculations. This work builds on the NuFast algorithm
developed by Denton and Parke, with original implementations available at
#link("https://github.com/PeterDenton/NuFast").

The author also acknowledges the neutrino physics education received during
undergraduate research at Krea University under Dr. Sushant Raut, which
provided the theoretical foundation for this implementation.

#heading(numbering: none)[Code Availability]

The `nufast` crate is open source under the MIT license:
- Crates.io: #link("https://crates.io/crates/nufast")
- GitHub: #link("https://github.com/planckeon/nufast")
- Documentation: #link("https://docs.rs/nufast")

All benchmark implementations (Rust, C++, Fortran, Python) are included in the
repository under `benchmarks/`.

#heading(numbering: none)[References]

+ P. B. Denton, "Neutrino Oscillation Probabilities: A Compact Multi-algorithm Approach," arXiv:2405.02400 (2024).

+ P. B. Denton and S. J. Parke, "Simple and precise factorization of the Jarlskog invariant for neutrino oscillations in matter," Phys. Rev. D 100, 053004 (2019). arXiv:1902.07185

+ DUNE Collaboration, "Long-Baseline Neutrino Facility (LBNF) and Deep Underground Neutrino Experiment (DUNE) Conceptual Design Report," arXiv:1601.05471

+ Hyper-Kamiokande Proto-Collaboration, "Hyper-Kamiokande Design Report," arXiv:1805.04163

+ JUNO Collaboration, "Neutrino Physics with JUNO," J. Phys. G 43, 030401 (2016). arXiv:1507.05613
