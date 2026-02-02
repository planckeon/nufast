#!/usr/bin/env python3
"""
NuFast Python Benchmark

Original algorithm by Peter Denton, adapted for benchmarking.

Run: python benchmark.py
"""

import numpy as np
import time

PI = np.pi
YerhoE2a = 1.52e-4
eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4.0


def probability_vacuum_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E):
    """Calculate vacuum oscillation probabilities."""
    c13sq = 1.0 - s13sq
    Ue3sq = s13sq
    Ue2sq = c13sq * s12sq
    Um3sq = c13sq * s23sq
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1.0 - s12sq) * (1.0 - s23sq)
    
    Jrr = np.sqrt(Um2sq * Ut2sq)
    sind = np.sin(delta)
    cosd = np.cos(delta)
    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd
    Jvac = 8.0 * Jrr * c13sq * sind
    
    Ue1sq = 1.0 - Ue3sq - Ue2sq
    Um1sq = 1.0 - Um3sq - Um2sq
    Ut3sq = 1.0 - Um3sq - Ue3sq
    Ut2sq = 1.0 - Um2sq - Ue2sq
    Ut1sq = 1.0 - Um1sq - Ue1sq
    
    Lover4E = eVsqkm_to_GeV_over4 * L / E
    D21 = Dmsq21 * Lover4E
    D31 = Dmsq31 * Lover4E
    
    sinD21 = np.sin(D21)
    sinD31 = np.sin(D31)
    sinD32 = np.sin(D31 - D21)
    
    triple_sin = sinD21 * sinD31 * sinD32
    sinsqD21_2 = 2.0 * sinD21 * sinD21
    sinsqD31_2 = 2.0 * sinD31 * sinD31
    sinsqD32_2 = 2.0 * sinD32 * sinD32
    
    Pme_CPC = ((Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
               (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
               (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2)
    Pme_CPV = -Jvac * triple_sin
    
    Pmm = 1.0 - 2.0 * (Um2sq * Um1sq * sinsqD21_2 +
                       Um3sq * Um1sq * sinsqD31_2 +
                       Um3sq * Um2sq * sinsqD32_2)
    Pee = 1.0 - 2.0 * (Ue2sq * Ue1sq * sinsqD21_2 +
                       Ue3sq * Ue1sq * sinsqD31_2 +
                       Ue3sq * Ue2sq * sinsqD32_2)
    
    probs = np.zeros((3, 3))
    probs[0, 0] = Pee
    probs[0, 1] = Pme_CPC - Pme_CPV
    probs[0, 2] = 1.0 - Pee - probs[0, 1]
    probs[1, 0] = Pme_CPC + Pme_CPV
    probs[1, 1] = Pmm
    probs[1, 2] = 1.0 - probs[1, 0] - Pmm
    probs[2, 0] = 1.0 - Pee - probs[1, 0]
    probs[2, 1] = 1.0 - probs[0, 1] - Pmm
    probs[2, 2] = 1.0 - probs[0, 2] - probs[1, 2]
    
    return probs


def probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                           L, E, rho, Ye, N_Newton):
    """Calculate matter oscillation probabilities."""
    c13sq = 1.0 - s13sq
    Ue2sq = c13sq * s12sq
    Ue3sq = s13sq
    Um3sq = c13sq * s23sq
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1.0 - s12sq) * (1.0 - s23sq)
    
    Jrr = np.sqrt(Um2sq * Ut2sq)
    sind = np.sin(delta)
    cosd = np.cos(delta)
    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd
    Jmatter = 8.0 * Jrr * c13sq * sind
    Amatter = Ye * rho * E * YerhoE2a
    Dmsqee = Dmsq31 - s12sq * Dmsq21
    
    A_sum = Dmsq21 + Dmsq31
    See = A_sum - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq
    Tmm_base = Dmsq21 * Dmsq31
    Tee = Tmm_base * (1.0 - Ue3sq - Ue2sq)
    C = Amatter * Tee
    A = A_sum + Amatter
    
    xmat = Amatter / Dmsqee
    tmp = 1.0 - xmat
    lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1.0 + np.sqrt(tmp * tmp + 4.0 * s13sq * xmat))
    
    B = Tmm_base + Amatter * See
    for _ in range(N_Newton):
        lambda3 = (lambda3 * lambda3 * (lambda3 - A) + C) / (lambda3 * (2.0 * lambda3 - A) + B)
    
    tmp = A - lambda3
    Dlambda21 = np.sqrt(tmp * tmp - 4.0 * C / lambda3)
    lambda2 = 0.5 * (A - lambda3 + Dlambda21)
    Dlambda32 = lambda3 - lambda2
    Dlambda31 = Dlambda32 + Dlambda21
    
    PiDlambdaInv = 1.0 / (Dlambda31 * Dlambda32 * Dlambda21)
    Xp3 = PiDlambdaInv * Dlambda21
    Xp2 = -PiDlambdaInv * Dlambda31
    
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2
    
    Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq
    Tmm = Tmm_base * (1.0 - Um3sq - Um2sq) + Amatter * (See + Smm - A_sum)
    
    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2
    
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv
    
    Ue1sq = 1.0 - Ue3sq - Ue2sq
    Um1sq = 1.0 - Um3sq - Um2sq
    Ut3sq = 1.0 - Um3sq - Ue3sq
    Ut2sq = 1.0 - Um2sq - Ue2sq
    Ut1sq = 1.0 - Um1sq - Ue1sq
    
    Lover4E = eVsqkm_to_GeV_over4 * L / E
    D21 = Dlambda21 * Lover4E
    D32 = Dlambda32 * Lover4E
    
    sinD21 = np.sin(D21)
    sinD31 = np.sin(D32 + D21)
    sinD32 = np.sin(D32)
    
    triple_sin = sinD21 * sinD31 * sinD32
    sinsqD21_2 = 2.0 * sinD21 * sinD21
    sinsqD31_2 = 2.0 * sinD31 * sinD31
    sinsqD32_2 = 2.0 * sinD32 * sinD32
    
    Pme_CPC = ((Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 +
               (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 +
               (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2)
    Pme_CPV = -Jmatter * triple_sin
    
    Pmm = 1.0 - 2.0 * (Um2sq * Um1sq * sinsqD21_2 +
                       Um3sq * Um1sq * sinsqD31_2 +
                       Um3sq * Um2sq * sinsqD32_2)
    Pee = 1.0 - 2.0 * (Ue2sq * Ue1sq * sinsqD21_2 +
                       Ue3sq * Ue1sq * sinsqD31_2 +
                       Ue3sq * Ue2sq * sinsqD32_2)
    
    probs = np.zeros((3, 3))
    probs[0, 0] = Pee
    probs[0, 1] = Pme_CPC - Pme_CPV
    probs[0, 2] = 1.0 - Pee - probs[0, 1]
    probs[1, 0] = Pme_CPC + Pme_CPV
    probs[1, 1] = Pmm
    probs[1, 2] = 1.0 - probs[1, 0] - Pmm
    probs[2, 0] = 1.0 - Pee - probs[1, 0]
    probs[2, 1] = 1.0 - probs[0, 1] - Pmm
    probs[2, 2] = 1.0 - probs[0, 2] - probs[1, 2]
    
    return probs


def main():
    n = 100000  # Fewer iterations for Python (it's slow)
    Emin, Emax = 0.5, 5.0
    
    s12sq, s13sq, s23sq = 0.31, 0.02, 0.55
    delta = -0.7 * PI
    Dmsq21, Dmsq31 = 7.5e-5, 2.5e-3
    L, rho, Ye = 1300.0, 3.0, 0.5
    
    print(f"NuFast Python Benchmark (n={n} iterations)")
    print("=" * 44)
    print()
    print("Single-point calculations:")
    
    # Vacuum
    total = 0.0
    t1 = time.perf_counter()
    for i in range(n):
        E = Emin + (Emax - Emin) * (i % 1000) / 1000.0
        probs = probability_vacuum_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E)
        total += probs[1, 0]
    t2 = time.perf_counter()
    ns_vac = (t2 - t1) * 1e9 / n
    print(f"  Vacuum:            {ns_vac:10.2f} ns")
    
    # Matter N=0,1,2,3
    for N in range(4):
        total = 0.0
        t1 = time.perf_counter()
        for i in range(n):
            E = Emin + (Emax - Emin) * (i % 1000) / 1000.0
            probs = probability_matter_lbl(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, 
                                          L, E, rho, Ye, N)
            total += probs[1, 0]
        t2 = time.perf_counter()
        ns = (t2 - t1) * 1e9 / n
        print(f"  Matter N_Newton={N}: {ns:10.2f} ns")
    
    print()
    print("Throughput:")
    print(f"  Vacuum:            {n / (ns_vac * 1e-9) / 1e6:10.4f} M/s")
    print(f"  Matter N_Newton=0: {n / (ns * 1e-9) / 1e6:10.4f} M/s")


if __name__ == "__main__":
    main()
