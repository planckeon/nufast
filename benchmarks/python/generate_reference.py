#!/usr/bin/env python3
"""Generate reference values for cross-validation tests - robust version."""

import numpy as np
from benchmark import probability_matter_lbl

PI = np.pi

# NuFit 5.2 parameters
s12sq = 0.307
s13sq = 0.0220
s23sq = 0.546
delta = -0.7 * PI
Dmsq21 = 7.53e-5
Dmsq31 = 2.453e-3

print("NuFit 5.2 Reference Values for Cross-Validation")
print("=" * 60)
print(f"\nParameters:")
print(f"  s12sq = {s12sq}")
print(f"  s13sq = {s13sq}")
print(f"  s23sq = {s23sq}")
print(f"  delta = {delta} rad ({delta/PI:.1f}*pi)")
print(f"  Dmsq21 = {Dmsq21}")
print(f"  Dmsq31 = {Dmsq31}")

# DUNE-like: L=1300 km, E=2.5 GeV, rho=2.848, Ye=0.5, N_Newton=0
print("\n" + "-" * 60)
print("DUNE-like: L=1300 km, E=2.5 GeV, rho=2.848, Ye=0.5, N=0")
print("-" * 60)

L, E, rho, Ye = 1300.0, 2.5, 2.848, 0.5
probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                               L, E, rho, Ye, 0)
print(f"\nN_Newton = 0:")
print(f"  P(e->e)  = {probs[0, 0]:.15e}")
print(f"  P(e->mu) = {probs[0, 1]:.15e}")
print(f"  P(e->tau)= {probs[0, 2]:.15e}")
print(f"  P(mu->e) = {probs[1, 0]:.15e}")
print(f"  P(mu->mu)= {probs[1, 1]:.15e}")
print(f"  P(mu->tau)={probs[1, 2]:.15e}")

# T2K-like: L=295 km, E=0.6 GeV, rho=2.6, Ye=0.5, N_Newton=0
print("\n" + "-" * 60)
print("T2K-like: L=295 km, E=0.6 GeV, rho=2.6, Ye=0.5, N=0")
print("-" * 60)

L, E, rho, Ye = 295.0, 0.6, 2.6, 0.5
probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                               L, E, rho, Ye, 0)
print(f"\nN_Newton = 0:")
print(f"  P(e->e)  = {probs[0, 0]:.15e}")
print(f"  P(e->mu) = {probs[0, 1]:.15e}")
print(f"  P(e->tau)= {probs[0, 2]:.15e}")
print(f"  P(mu->e) = {probs[1, 0]:.15e}")
print(f"  P(mu->mu)= {probs[1, 1]:.15e}")
print(f"  P(mu->tau)={probs[1, 2]:.15e}")

# Also test at different energies where Newton iterations might be stable
print("\n" + "-" * 60)
print("DUNE-like: L=1300 km, E=3.0 GeV, rho=2.848, Ye=0.5")
print("-" * 60)

L, E, rho, Ye = 1300.0, 3.0, 2.848, 0.5
for N_Newton in [0, 1, 2]:
    probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                                   L, E, rho, Ye, N_Newton)
    # Check if probabilities are valid
    valid = all(0 <= probs[i,j] <= 1 for i in range(3) for j in range(3))
    print(f"\nN_Newton = {N_Newton}: {'VALID' if valid else 'INVALID'}")
    print(f"  P(e->e)  = {probs[0, 0]:.15e}")
    print(f"  P(mu->e) = {probs[1, 0]:.15e}")
    print(f"  P(mu->mu)= {probs[1, 1]:.15e}")

# Try another energy
print("\n" + "-" * 60)
print("DUNE-like: L=1300 km, E=1.5 GeV, rho=2.848, Ye=0.5")
print("-" * 60)

L, E, rho, Ye = 1300.0, 1.5, 2.848, 0.5
for N_Newton in [0, 1, 2]:
    probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                                   L, E, rho, Ye, N_Newton)
    valid = all(0 <= probs[i,j] <= 1 for i in range(3) for j in range(3))
    print(f"\nN_Newton = {N_Newton}: {'VALID' if valid else 'INVALID'}")
    print(f"  P(e->e)  = {probs[0, 0]:.15e}")
    print(f"  P(mu->e) = {probs[1, 0]:.15e}")
    print(f"  P(mu->mu)= {probs[1, 1]:.15e}")

print("\n" + "=" * 60)
print("ZIG TEST FORMAT (N_Newton=0 only, always stable):")
print("=" * 60)

# DUNE-like
L, E, rho, Ye = 1300.0, 2.5, 2.848, 0.5
probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                               L, E, rho, Ye, 0)
print(f"""
// DUNE-like: L=1300 km, E=2.5 GeV, rho=2.848, Ye=0.5, N_Newton=0
// Python reference (NuFit 5.2):
//   P(e->e)  = {probs[0, 0]:.15e}
//   P(e->mu) = {probs[0, 1]:.15e}
//   P(mu->e) = {probs[1, 0]:.15e}
//   P(mu->mu)= {probs[1, 1]:.15e}
try std.testing.expectApproxEqAbs(probs[0][0], {probs[0, 0]:.15e}, 1e-10);
try std.testing.expectApproxEqAbs(probs[0][1], {probs[0, 1]:.15e}, 1e-10);
try std.testing.expectApproxEqAbs(probs[1][0], {probs[1, 0]:.15e}, 1e-10);
try std.testing.expectApproxEqAbs(probs[1][1], {probs[1, 1]:.15e}, 1e-10);""")

# T2K-like
L, E, rho, Ye = 295.0, 0.6, 2.6, 0.5
probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                               L, E, rho, Ye, 0)
print(f"""
// T2K-like: L=295 km, E=0.6 GeV, rho=2.6, Ye=0.5, N_Newton=0
// Python reference (NuFit 5.2):
//   P(e->e)  = {probs[0, 0]:.15e}
//   P(e->mu) = {probs[0, 1]:.15e}
//   P(mu->e) = {probs[1, 0]:.15e}
//   P(mu->mu)= {probs[1, 1]:.15e}
try std.testing.expectApproxEqAbs(probs[0][0], {probs[0, 0]:.15e}, 1e-10);
try std.testing.expectApproxEqAbs(probs[0][1], {probs[0, 1]:.15e}, 1e-10);
try std.testing.expectApproxEqAbs(probs[1][0], {probs[1, 0]:.15e}, 1e-10);
try std.testing.expectApproxEqAbs(probs[1][1], {probs[1, 1]:.15e}, 1e-10);""")
