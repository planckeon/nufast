//! NuFast CLI - Quick probability calculations
//!
//! Run with: `cargo run --release`

use nufast::{MatterParameters, VacuumParameters, probability_matter_lbl, probability_vacuum_lbl};
use std::f64::consts::PI;

fn main() {
    println!("NuFast - Fast Neutrino Oscillation Probabilities");
    println!("================================================\n");

    // DUNE-like parameters
    let vacuum = VacuumParameters {
        s12sq: 0.307,
        s13sq: 0.0218,
        s23sq: 0.545,
        delta: 1.36 * PI,
        Dmsq21: 7.42e-5,
        Dmsq31: 2.517e-3,
        L: 1300.0,
        E: 2.5,
    };

    let matter = MatterParameters {
        s12sq: vacuum.s12sq,
        s13sq: vacuum.s13sq,
        s23sq: vacuum.s23sq,
        delta: vacuum.delta,
        Dmsq21: vacuum.Dmsq21,
        Dmsq31: vacuum.Dmsq31,
        L: vacuum.L,
        E: vacuum.E,
        rho: 2.848,
        Ye: 0.5,
        N_Newton: 0,
    };

    println!("Parameters (DUNE-like, Normal Ordering):");
    println!("  Baseline: {} km", vacuum.L);
    println!("  Energy: {} GeV", vacuum.E);
    println!("  Density: {} g/cm³", matter.rho);
    println!();

    let probs_vac = probability_vacuum_lbl(&vacuum);
    let probs_mat = probability_matter_lbl(&matter);

    println!("Vacuum Oscillation Probabilities:");
    print_matrix(&probs_vac);

    println!("\nMatter Oscillation Probabilities (N_Newton={}):", matter.N_Newton);
    print_matrix(&probs_mat);

    println!("\nMatter Effect on P(νμ → νe):");
    println!("  Vacuum: {:.6}", probs_vac[1][0]);
    println!("  Matter: {:.6}", probs_mat[1][0]);
    println!("  Difference: {:.6} ({:+.1}%)", 
        probs_mat[1][0] - probs_vac[1][0],
        (probs_mat[1][0] / probs_vac[1][0] - 1.0) * 100.0
    );
}

fn print_matrix(probs: &[[f64; 3]; 3]) {
    println!("         e          μ          τ");
    let labels = ['e', 'μ', 'τ'];
    for (i, row) in probs.iter().enumerate() {
        println!("  {} → {:>9.6}  {:>9.6}  {:>9.6}", 
            labels[i], row[0], row[1], row[2]);
    }
}
