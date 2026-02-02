//! Benchmarks for NuFast oscillation probability calculations
//!
//! Run with: cargo bench
//!
//! This benchmark compares the Rust implementation against:
//! - Vacuum vs Matter calculations
//! - Different N_Newton iterations
//! - Batch calculations (energy spectrum)
//! - VacuumBatch optimization

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use nufast::{VacuumParameters, MatterParameters, VacuumBatch, probability_vacuum_lbl, probability_matter_lbl};
use std::f64::consts::PI;

/// Standard DUNE-like parameters for benchmarking
fn dune_vacuum_params(e: f64) -> VacuumParameters {
    VacuumParameters {
        s12sq: 0.31,
        s13sq: 0.02,
        s23sq: 0.55,
        delta: -0.7 * PI,
        Dmsq21: 7.5e-5,
        Dmsq31: 2.5e-3,
        L: 1300.0,
        E: e,
    }
}

fn dune_matter_params(e: f64, n_newton: u8) -> MatterParameters {
    MatterParameters {
        s12sq: 0.31,
        s13sq: 0.02,
        s23sq: 0.55,
        delta: -0.7 * PI,
        Dmsq21: 7.5e-5,
        Dmsq31: 2.5e-3,
        L: 1300.0,
        E: e,
        rho: 3.0,
        Ye: 0.5,
        N_Newton: n_newton,
    }
}

/// Single-point vacuum oscillation benchmark
fn bench_vacuum_single(c: &mut Criterion) {
    let params = dune_vacuum_params(2.5);
    
    c.bench_function("vacuum_single", |b| {
        b.iter(|| probability_vacuum_lbl(black_box(&params)))
    });
}

/// Single-point matter oscillation benchmark at various N_Newton levels
fn bench_matter_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("matter_single");
    
    for n_newton in [0u8, 1, 2, 3] {
        let params = dune_matter_params(2.5, n_newton);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("N_Newton={}", n_newton)),
            &params,
            |b, params| b.iter(|| probability_matter_lbl(black_box(params))),
        );
    }
    
    group.finish();
}

/// Batch calculation: 1000-point energy spectrum (DUNE-like)
fn bench_energy_spectrum(c: &mut Criterion) {
    let e_min = 0.5;
    let e_max = 5.0;
    let n_points = 1000;
    
    c.bench_function("vacuum_spectrum_1000", |b| {
        b.iter(|| {
            for i in 0..n_points {
                let e = e_min + (e_max - e_min) * (i as f64 / n_points as f64);
                let params = dune_vacuum_params(e);
                black_box(probability_vacuum_lbl(&params));
            }
        })
    });
    
    // Optimized batch version using pre-computed mixing elements
    c.bench_function("vacuum_batch_spectrum_1000", |b| {
        let batch = VacuumBatch::new(0.31, 0.02, 0.55, -0.7 * PI, 7.5e-5, 2.5e-3);
        b.iter(|| {
            for i in 0..n_points {
                let e = e_min + (e_max - e_min) * (i as f64 / n_points as f64);
                black_box(batch.probability_at(1300.0, e));
            }
        })
    });
    
    c.bench_function("matter_N0_spectrum_1000", |b| {
        b.iter(|| {
            for i in 0..n_points {
                let e = e_min + (e_max - e_min) * (i as f64 / n_points as f64);
                let params = dune_matter_params(e, 0);
                black_box(probability_matter_lbl(&params));
            }
        })
    });
    
    c.bench_function("matter_N1_spectrum_1000", |b| {
        b.iter(|| {
            for i in 0..n_points {
                let e = e_min + (e_max - e_min) * (i as f64 / n_points as f64);
                let params = dune_matter_params(e, 1);
                black_box(probability_matter_lbl(&params));
            }
        })
    });
}

/// Throughput benchmark: how many calculations per second
fn bench_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("throughput");
    group.throughput(criterion::Throughput::Elements(1_000_000));
    
    group.bench_function("vacuum_1M", |b| {
        let params = dune_vacuum_params(2.5);
        b.iter(|| {
            for _ in 0..1_000_000 {
                black_box(probability_vacuum_lbl(&params));
            }
        })
    });
    
    group.bench_function("matter_N0_1M", |b| {
        let params = dune_matter_params(2.5, 0);
        b.iter(|| {
            for _ in 0..1_000_000 {
                black_box(probability_matter_lbl(&params));
            }
        })
    });
    
    group.finish();
}

criterion_group!(
    benches,
    bench_vacuum_single,
    bench_matter_single,
    bench_energy_spectrum,
    bench_throughput,
);

criterion_main!(benches);
