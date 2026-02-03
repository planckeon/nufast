"""
Tests for NuFast Python bindings.

Run with: python test_nufast.py
"""

import math
import sys
from pathlib import Path

# Add parent directory to path to import nufast
sys.path.insert(0, str(Path(__file__).parent))

import nufast


def test_vacuum_probability_conservation():
    """Test that vacuum probability rows and columns sum to 1."""
    print("Testing vacuum probability conservation...")
    
    probs = nufast.vacuum_probability(1300.0, 2.5)
    
    # Check rows sum to 1
    row_e = probs.Pee + probs.Pem + probs.Pet
    row_m = probs.Pme + probs.Pmm + probs.Pmt
    row_t = probs.Pte + probs.Ptm + probs.Ptt
    
    assert abs(row_e - 1.0) < 1e-10, f"Row e sums to {row_e}, expected 1.0"
    assert abs(row_m - 1.0) < 1e-10, f"Row μ sums to {row_m}, expected 1.0"
    assert abs(row_t - 1.0) < 1e-10, f"Row τ sums to {row_t}, expected 1.0"
    
    # Check columns sum to 1
    col_e = probs.Pee + probs.Pme + probs.Pte
    col_m = probs.Pem + probs.Pmm + probs.Ptm
    col_t = probs.Pet + probs.Pmt + probs.Ptt
    
    assert abs(col_e - 1.0) < 1e-10, f"Column e sums to {col_e}, expected 1.0"
    assert abs(col_m - 1.0) < 1e-10, f"Column μ sums to {col_m}, expected 1.0"
    assert abs(col_t - 1.0) < 1e-10, f"Column τ sums to {col_t}, expected 1.0"
    
    print("  ✓ Row sums: e={:.10f}, μ={:.10f}, τ={:.10f}".format(row_e, row_m, row_t))
    print("  ✓ Column sums: e={:.10f}, μ={:.10f}, τ={:.10f}".format(col_e, col_m, col_t))


def test_matter_probability_conservation():
    """Test that matter probability rows and columns sum to 1."""
    print("Testing matter probability conservation...")
    
    probs = nufast.matter_probability(1300.0, 2.5, rho=2.848)
    
    # Check rows sum to 1
    row_e = probs.Pee + probs.Pem + probs.Pet
    row_m = probs.Pme + probs.Pmm + probs.Pmt
    row_t = probs.Pte + probs.Ptm + probs.Ptt
    
    assert abs(row_e - 1.0) < 1e-10, f"Row e sums to {row_e}, expected 1.0"
    assert abs(row_m - 1.0) < 1e-10, f"Row μ sums to {row_m}, expected 1.0"
    assert abs(row_t - 1.0) < 1e-10, f"Row τ sums to {row_t}, expected 1.0"
    
    print("  ✓ Row sums: e={:.10f}, μ={:.10f}, τ={:.10f}".format(row_e, row_m, row_t))


def test_vacuum_known_values():
    """Test vacuum probabilities against known reference values."""
    print("Testing vacuum probabilities against reference values...")
    
    # DUNE-like: L=1300 km, E=2.5 GeV
    probs = nufast.vacuum_probability(1300.0, 2.5)
    
    # Reference values from NuFast Zig tests (which match Python reference)
    expected_Pee = 0.9120570461953496
    expected_Pme = 0.05859330228592972
    expected_Pmm = 0.004726368319321272
    
    assert abs(probs.Pee - expected_Pee) < 1e-10, \
        f"Pee = {probs.Pee}, expected {expected_Pee}"
    assert abs(probs.Pme - expected_Pme) < 1e-10, \
        f"Pme = {probs.Pme}, expected {expected_Pme}"
    assert abs(probs.Pmm - expected_Pmm) < 1e-10, \
        f"Pmm = {probs.Pmm}, expected {expected_Pmm}"
    
    print(f"  ✓ Pee = {probs.Pee:.10f} (expected {expected_Pee})")
    print(f"  ✓ Pme = {probs.Pme:.10f} (expected {expected_Pme})")
    print(f"  ✓ Pmm = {probs.Pmm:.10f} (expected {expected_Pmm})")


def test_matter_known_values():
    """Test matter probabilities against known reference values."""
    print("Testing matter probabilities against reference values...")
    
    # DUNE-like: L=1300 km, E=2.5 GeV, rho=2.848 g/cm³
    probs = nufast.matter_probability(1300.0, 2.5, rho=2.848)
    
    # Reference values from NuFast Zig tests
    expected_Pee = 8.727627013058699e-01
    expected_Pem = 5.731530983046839e-02
    expected_Pme = 8.249293182552475e-02
    expected_Pmm = 5.260049583365356e-02
    
    assert abs(probs.Pee - expected_Pee) < 1e-10, \
        f"Pee = {probs.Pee}, expected {expected_Pee}"
    assert abs(probs.Pem - expected_Pem) < 1e-10, \
        f"Pem = {probs.Pem}, expected {expected_Pem}"
    assert abs(probs.Pme - expected_Pme) < 1e-10, \
        f"Pme = {probs.Pme}, expected {expected_Pme}"
    assert abs(probs.Pmm - expected_Pmm) < 1e-10, \
        f"Pmm = {probs.Pmm}, expected {expected_Pmm}"
    
    print(f"  ✓ Pee = {probs.Pee:.10f} (expected {expected_Pee})")
    print(f"  ✓ Pem = {probs.Pem:.10f} (expected {expected_Pem})")
    print(f"  ✓ Pme = {probs.Pme:.10f} (expected {expected_Pme})")
    print(f"  ✓ Pmm = {probs.Pmm:.10f} (expected {expected_Pmm})")


def test_quick_pme():
    """Test quick Pme functions match full matrix."""
    print("Testing quick Pme functions...")
    
    L, E = 1300.0, 2.5
    
    # Vacuum
    probs_full = nufast.vacuum_probability(L, E)
    pme_quick = nufast.vacuum_Pme(L, E)
    
    assert abs(probs_full.Pme - pme_quick) < 1e-15, \
        f"vacuum_Pme mismatch: {pme_quick} vs {probs_full.Pme}"
    print(f"  ✓ vacuum_Pme = {pme_quick:.10f}")
    
    # Matter
    probs_full = nufast.matter_probability(L, E, rho=2.848)
    pme_quick = nufast.matter_Pme(L, E, rho=2.848)
    
    assert abs(probs_full.Pme - pme_quick) < 1e-15, \
        f"matter_Pme mismatch: {pme_quick} vs {probs_full.Pme}"
    print(f"  ✓ matter_Pme = {pme_quick:.10f}")


def test_batch_processing():
    """Test batch processing functions."""
    print("Testing batch processing...")
    
    L = 1300.0
    energies = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    
    # Vacuum batch
    results = nufast.vacuum_Pme_batch(L, energies)
    
    # Compare with individual calls
    for i, E in enumerate(energies):
        expected = nufast.vacuum_Pme(L, E)
        assert abs(results[i] - expected) < 1e-15, \
            f"Batch mismatch at E={E}: {results[i]} vs {expected}"
    
    print(f"  ✓ Vacuum batch: {len(results)} energies computed correctly")
    
    # Matter batch
    results = nufast.matter_Pme_batch(L, energies, rho=2.848)
    
    for i, E in enumerate(energies):
        expected = nufast.matter_Pme(L, E, rho=2.848)
        assert abs(results[i] - expected) < 1e-15, \
            f"Batch mismatch at E={E}: {results[i]} vs {expected}"
    
    print(f"  ✓ Matter batch: {len(results)} energies computed correctly")


def test_antineutrino():
    """Test antineutrino mode."""
    print("Testing antineutrino mode...")
    
    L, E = 1300.0, 2.5
    
    # Vacuum: Pme(ν) should equal Pem(ν̄) (CPT theorem)
    probs_nu = nufast.vacuum_probability(L, E, antineutrino=False)
    probs_nubar = nufast.vacuum_probability(L, E, antineutrino=True)
    
    assert abs(probs_nu.Pme - probs_nubar.Pem) < 1e-10, \
        f"CPT violation: Pme(ν)={probs_nu.Pme} != Pem(ν̄)={probs_nubar.Pem}"
    assert abs(probs_nu.Pem - probs_nubar.Pme) < 1e-10, \
        f"CPT violation: Pem(ν)={probs_nu.Pem} != Pme(ν̄)={probs_nubar.Pme}"
    
    print(f"  ✓ CPT: Pme(ν) = {probs_nu.Pme:.6f}, Pem(ν̄) = {probs_nubar.Pem:.6f}")
    
    # Matter: neutrinos and antineutrinos should differ
    probs_nu = nufast.matter_probability(L, E, rho=2.848, antineutrino=False)
    probs_nubar = nufast.matter_probability(L, E, rho=2.848, antineutrino=True)
    
    diff = abs(probs_nu.Pme - probs_nubar.Pme)
    assert diff > 0.01, f"Matter effect too small: |Pme(ν) - Pme(ν̄)| = {diff}"
    print(f"  ✓ Matter: Pme(ν) = {probs_nu.Pme:.6f}, Pme(ν̄) = {probs_nubar.Pme:.6f}")


def test_custom_params():
    """Test custom oscillation parameters."""
    print("Testing custom oscillation parameters...")
    
    # Inverted mass ordering
    params = nufast.OscillationParams(
        s12sq=0.307,
        s13sq=0.0220,
        s23sq=0.546,
        delta=-0.7 * math.pi,
        Dmsq21=7.53e-5,
        Dmsq31=-2.498e-3,  # Inverted ordering (negative)
    )
    
    probs = nufast.vacuum_probability(1300.0, 2.5, params=params)
    
    # Reference from Zig tests
    expected_Pee = 0.9126759218475518
    expected_Pme = 0.05722581249769201
    expected_Pmm = 0.01733314218311000
    
    assert abs(probs.Pee - expected_Pee) < 1e-10, \
        f"Pee = {probs.Pee}, expected {expected_Pee}"
    assert abs(probs.Pme - expected_Pme) < 1e-10, \
        f"Pme = {probs.Pme}, expected {expected_Pme}"
    assert abs(probs.Pmm - expected_Pmm) < 1e-10, \
        f"Pmm = {probs.Pmm}, expected {expected_Pmm}"
    
    print(f"  ✓ Inverted ordering: Pme = {probs.Pme:.10f}")


def test_experiments():
    """Test experiment presets."""
    print("Testing experiment presets...")
    
    for name in ["T2K", "NOvA", "DUNE"]:
        probs = nufast.experiment_probability(name)
        row_sum = probs.Pme + probs.Pmm + probs.Pmt
        assert abs(row_sum - 1.0) < 1e-10, f"{name}: row sum = {row_sum}"
        print(f"  ✓ {name}: Pme = {probs.Pme:.6f}")


def test_edge_cases():
    """Test edge cases."""
    print("Testing edge cases...")
    
    # Near-zero baseline (should approach identity)
    probs = nufast.vacuum_probability(0.001, 2.5)
    assert abs(probs.Pee - 1.0) < 1e-10, f"L→0: Pee = {probs.Pee}, expected 1.0"
    assert abs(probs.Pme - 0.0) < 1e-10, f"L→0: Pme = {probs.Pme}, expected 0.0"
    print(f"  ✓ L→0: Pee = {probs.Pee:.10f}, Pme = {probs.Pme:.10f}")
    
    # Very high energy (oscillations suppressed)
    probs = nufast.vacuum_probability(1300.0, 100.0)
    assert probs.Pee > 0.99, f"High E: Pee = {probs.Pee}, expected > 0.99"
    print(f"  ✓ High E: Pee = {probs.Pee:.10f}")
    
    # Zero matter density should equal vacuum
    probs_vac = nufast.vacuum_probability(1300.0, 2.5)
    probs_mat = nufast.matter_probability(1300.0, 2.5, rho=0.0)
    
    assert abs(probs_vac.Pme - probs_mat.Pme) < 1e-10, \
        f"rho=0: matter Pme = {probs_mat.Pme}, vacuum = {probs_vac.Pme}"
    print(f"  ✓ rho=0: vacuum Pme = matter Pme = {probs_vac.Pme:.10f}")


def test_default_params():
    """Test default parameter access."""
    print("Testing default parameters...")
    
    params = nufast.OscillationParams.default()
    
    assert abs(params.s12sq - 0.307) < 1e-10
    assert abs(params.s13sq - 0.0220) < 1e-10
    assert abs(params.s23sq - 0.546) < 1e-10
    assert abs(params.delta - (-0.7 * math.pi)) < 1e-10
    assert abs(params.Dmsq21 - 7.53e-5) < 1e-15
    assert abs(params.Dmsq31 - 2.453e-3) < 1e-15
    
    print(f"  ✓ s12sq = {params.s12sq}")
    print(f"  ✓ s13sq = {params.s13sq}")
    print(f"  ✓ s23sq = {params.s23sq}")
    print(f"  ✓ delta = {params.delta:.6f} ({params.delta/math.pi:.2f}π)")
    print(f"  ✓ Δm²₂₁ = {params.Dmsq21} eV²")
    print(f"  ✓ Δm²₃₁ = {params.Dmsq31} eV²")


def main():
    """Run all tests."""
    print("=" * 60)
    print("NuFast Python Bindings Test Suite")
    print("=" * 60)
    print()
    
    tests = [
        test_default_params,
        test_vacuum_probability_conservation,
        test_matter_probability_conservation,
        test_vacuum_known_values,
        test_matter_known_values,
        test_quick_pme,
        test_batch_processing,
        test_antineutrino,
        test_custom_params,
        test_experiments,
        test_edge_cases,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            print()
            passed += 1
        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            print()
            failed += 1
    
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
