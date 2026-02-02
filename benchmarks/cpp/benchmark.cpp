// NuFast C++ Benchmark
// Original algorithm by Peter Denton, adapted for benchmarking
//
// Compile: g++ -O3 -march=native -o benchmark benchmark.cpp
// Run: ./benchmark

#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>

const double PI = 3.14159265358979323846;
const double YerhoE2a = 1.52e-4;
const double eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4.0;

// Vacuum oscillation probabilities
void Probability_Vacuum_LBL(
    double s12sq, double s13sq, double s23sq, double delta,
    double Dmsq21, double Dmsq31, double L, double E,
    double probs[3][3]
) {
    double c13sq = 1.0 - s13sq;
    double Ue3sq = s13sq;
    double Ue2sq = c13sq * s12sq;
    double Um3sq = c13sq * s23sq;
    double Ut2sq = s13sq * s12sq * s23sq;
    double Um2sq = (1.0 - s12sq) * (1.0 - s23sq);
    
    double Jrr = std::sqrt(Um2sq * Ut2sq);
    double sind = std::sin(delta);
    double cosd = std::cos(delta);
    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd;
    double Jvac = 8.0 * Jrr * c13sq * sind;
    
    double Ue1sq = 1.0 - Ue3sq - Ue2sq;
    double Um1sq = 1.0 - Um3sq - Um2sq;
    double Ut3sq = 1.0 - Um3sq - Ue3sq;
    Ut2sq = 1.0 - Um2sq - Ue2sq;
    double Ut1sq = 1.0 - Um1sq - Ue1sq;
    
    double Lover4E = eVsqkm_to_GeV_over4 * L / E;
    double D21 = Dmsq21 * Lover4E;
    double D31 = Dmsq31 * Lover4E;
    
    double sinD21 = std::sin(D21);
    double sinD31 = std::sin(D31);
    double sinD32 = std::sin(D31 - D21);
    
    double triple_sin = sinD21 * sinD31 * sinD32;
    double sinsqD21_2 = 2.0 * sinD21 * sinD21;
    double sinsqD31_2 = 2.0 * sinD31 * sinD31;
    double sinsqD32_2 = 2.0 * sinD32 * sinD32;
    
    double Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
                   + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
                   + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    double Pme_CPV = -Jvac * triple_sin;
    
    double Pmm = 1.0 - 2.0 * (Um2sq * Um1sq * sinsqD21_2 
                            + Um3sq * Um1sq * sinsqD31_2 
                            + Um3sq * Um2sq * sinsqD32_2);
    double Pee = 1.0 - 2.0 * (Ue2sq * Ue1sq * sinsqD21_2 
                            + Ue3sq * Ue1sq * sinsqD31_2 
                            + Ue3sq * Ue2sq * sinsqD32_2);
    
    probs[0][0] = Pee;
    probs[0][1] = Pme_CPC - Pme_CPV;
    probs[0][2] = 1.0 - Pee - probs[0][1];
    probs[1][0] = Pme_CPC + Pme_CPV;
    probs[1][1] = Pmm;
    probs[1][2] = 1.0 - probs[1][0] - Pmm;
    probs[2][0] = 1.0 - Pee - probs[1][0];
    probs[2][1] = 1.0 - probs[0][1] - Pmm;
    probs[2][2] = 1.0 - probs[0][2] - probs[1][2];
}

// Matter oscillation probabilities
void Probability_Matter_LBL(
    double s12sq, double s13sq, double s23sq, double delta,
    double Dmsq21, double Dmsq31, double L, double E,
    double rho, double Ye, int N_Newton,
    double probs[3][3]
) {
    double c13sq = 1.0 - s13sq;
    double Ue2sq = c13sq * s12sq;
    double Ue3sq = s13sq;
    double Um3sq = c13sq * s23sq;
    double Ut2sq = s13sq * s12sq * s23sq;
    double Um2sq = (1.0 - s12sq) * (1.0 - s23sq);
    
    double Jrr = std::sqrt(Um2sq * Ut2sq);
    double sind = std::sin(delta);
    double cosd = std::cos(delta);
    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd;
    double Jmatter = 8.0 * Jrr * c13sq * sind;
    double Amatter = Ye * rho * E * YerhoE2a;
    double Dmsqee = Dmsq31 - s12sq * Dmsq21;
    
    double A_sum = Dmsq21 + Dmsq31;
    double See = A_sum - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq;
    double Tmm_base = Dmsq21 * Dmsq31;
    double Tee = Tmm_base * (1.0 - Ue3sq - Ue2sq);
    double C = Amatter * Tee;
    double A = A_sum + Amatter;
    
    double xmat = Amatter / Dmsqee;
    double tmp = 1.0 - xmat;
    double lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1.0 + std::sqrt(tmp * tmp + 4.0 * s13sq * xmat));
    
    double B = Tmm_base + Amatter * See;
    for (int i = 0; i < N_Newton; i++) {
        lambda3 = (lambda3 * lambda3 * (lambda3 - A) + C) / (lambda3 * (2.0 * lambda3 - A) + B);
    }
    
    tmp = A - lambda3;
    double Dlambda21 = std::sqrt(tmp * tmp - 4.0 * C / lambda3);
    double lambda2 = 0.5 * (A - lambda3 + Dlambda21);
    double Dlambda32 = lambda3 - lambda2;
    double Dlambda31 = Dlambda32 + Dlambda21;
    
    double PiDlambdaInv = 1.0 / (Dlambda31 * Dlambda32 * Dlambda21);
    double Xp3 = PiDlambdaInv * Dlambda21;
    double Xp2 = -PiDlambdaInv * Dlambda31;
    
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;
    
    double Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq;
    double Tmm = Tmm_base * (1.0 - Um3sq - Um2sq) + Amatter * (See + Smm - A_sum);
    
    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;
    
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;
    
    double Ue1sq = 1.0 - Ue3sq - Ue2sq;
    double Um1sq = 1.0 - Um3sq - Um2sq;
    double Ut3sq = 1.0 - Um3sq - Ue3sq;
    Ut2sq = 1.0 - Um2sq - Ue2sq;
    double Ut1sq = 1.0 - Um1sq - Ue1sq;
    
    double Lover4E = eVsqkm_to_GeV_over4 * L / E;
    double D21 = Dlambda21 * Lover4E;
    double D32 = Dlambda32 * Lover4E;
    
    double sinD21 = std::sin(D21);
    double sinD31 = std::sin(D32 + D21);
    double sinD32 = std::sin(D32);
    
    double triple_sin = sinD21 * sinD31 * sinD32;
    double sinsqD21_2 = 2.0 * sinD21 * sinD21;
    double sinsqD31_2 = 2.0 * sinD31 * sinD31;
    double sinsqD32_2 = 2.0 * sinD32 * sinD32;
    
    double Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
                   + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
                   + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    double Pme_CPV = -Jmatter * triple_sin;
    
    double Pmm = 1.0 - 2.0 * (Um2sq * Um1sq * sinsqD21_2 
                            + Um3sq * Um1sq * sinsqD31_2 
                            + Um3sq * Um2sq * sinsqD32_2);
    double Pee = 1.0 - 2.0 * (Ue2sq * Ue1sq * sinsqD21_2 
                            + Ue3sq * Ue1sq * sinsqD31_2 
                            + Ue3sq * Ue2sq * sinsqD32_2);
    
    probs[0][0] = Pee;
    probs[0][1] = Pme_CPC - Pme_CPV;
    probs[0][2] = 1.0 - Pee - probs[0][1];
    probs[1][0] = Pme_CPC + Pme_CPV;
    probs[1][1] = Pmm;
    probs[1][2] = 1.0 - probs[1][0] - Pmm;
    probs[2][0] = 1.0 - Pee - probs[1][0];
    probs[2][1] = 1.0 - probs[0][1] - Pmm;
    probs[2][2] = 1.0 - probs[0][2] - probs[1][2];
}

// Prevent optimization
volatile double sink = 0;

int main() {
    const int n = 10000000;
    const double Emin = 0.5, Emax = 5.0;
    
    double s12sq = 0.31, s13sq = 0.02, s23sq = 0.55;
    double delta = -0.7 * PI;
    double Dmsq21 = 7.5e-5, Dmsq31 = 2.5e-3;
    double L = 1300.0, rho = 3.0, Ye = 0.5;
    
    double probs[3][3];
    double total;
    
    std::cout << "NuFast C++ Benchmark (n=" << n << " iterations)\n";
    std::cout << "============================================\n\n";
    std::cout << std::fixed << std::setprecision(2);
    
    // Warmup
    for (int i = 0; i < 10000; i++) {
        double E = Emin + (Emax - Emin) * (i % 1000) / 1000.0;
        Probability_Vacuum_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, probs);
        sink += probs[1][0];
    }
    
    std::cout << "Single-point calculations:\n";
    
    // Vacuum
    total = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        double E = Emin + (Emax - Emin) * (i % 1000) / 1000.0;
        Probability_Vacuum_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, probs);
        total += probs[1][0];
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    sink += total;
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    std::cout << "  Vacuum:            " << std::setw(10) << (double)ns / n << " ns\n";
    
    // Matter N=0
    for (int N = 0; N <= 3; N++) {
        total = 0;
        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; i++) {
            double E = Emin + (Emax - Emin) * (i % 1000) / 1000.0;
            Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N, probs);
            total += probs[1][0];
        }
        t2 = std::chrono::high_resolution_clock::now();
        sink += total;
        ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        std::cout << "  Matter N_Newton=" << N << ": " << std::setw(10) << (double)ns / n << " ns\n";
    }
    
    std::cout << "\nThroughput:\n";
    std::cout << "  Vacuum:            " << std::setw(10) << std::setprecision(4) << n * 1e-6 / (ns / 1e9) << " M/s\n";
    
    return 0;
}
