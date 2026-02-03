#!/usr/bin/env python3
"""Debug probability sums."""

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

# DUNE-like at E=2.5
L, E, rho, Ye = 1300.0, 2.5, 2.848, 0.5
probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                               L, E, rho, Ye, 0)

print("DUNE-like E=2.5 GeV, N=0:")
print(f"Full matrix:")
print(probs)
print(f"\nRow sums: {probs.sum(axis=1)}")
print(f"Col sums: {probs.sum(axis=0)}")

# Check T2K-like
L, E, rho, Ye = 295.0, 0.6, 2.6, 0.5
probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                               L, E, rho, Ye, 0)

print("\nT2K-like E=0.6 GeV, N=0:")
print(f"Full matrix:")
print(probs)
print(f"\nRow sums: {probs.sum(axis=1)}")
print(f"Col sums: {probs.sum(axis=0)}")

# All probabilities in [0,1]?
print("\nAll in [0,1]:", all(0 <= p <= 1 for p in probs.flat))
