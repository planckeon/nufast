module Parameters
    implicit none
    real(8), parameter :: PI = 3.14159265358979323846264338327950288419716939937510d0
    real(8), parameter :: YerhoE2a = 1.52588e-4_8
    real(8), parameter :: eVsqkm_to_GeV_over4 = 1e-9_8 / 1.97327e-7_8 * 1e3_8 / 4
end module Parameters

subroutine Probability_Vacuum_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, probs_returned)
    use Parameters
    implicit none
    real(8), intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E
    real(8), intent(out) :: probs_returned(3, 3)

    real(8) :: c13sq, sind, cosd, Jrr, Jvac
    real(8) :: Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq
    real(8) :: Lover4E, D21, D31
    real(8) :: sinD21, sinD31, sinD32
    real(8) :: sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin
    real(8) :: Pme_CPC, Pme_CPV, Pmm, Pee

    c13sq = 1 - s13sq
    Ue3sq = s13sq
    Ue2sq = c13sq * s12sq
    Um3sq = c13sq * s23sq
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1 - s12sq) * (1 - s23sq)
    Jrr = sqrt(Um2sq * Ut2sq)
    sind = sin(delta)
    cosd = cos(delta)
    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd
    Jvac = 8 * Jrr * c13sq * sind
    Ue1sq = 1 - Ue3sq - Ue2sq
    Um1sq = 1 - Um3sq - Um2sq
    Ut3sq = 1 - Um3sq - Ue3sq
    Ut2sq = 1 - Um2sq - Ue2sq
    Ut1sq = 1 - Um1sq - Ue1sq
    Lover4E = eVsqkm_to_GeV_over4 * L / E
    D21 = dmsq21 * Lover4E
    D31 = dmsq31 * Lover4E
    sinD21 = sin(D21)
    sinD31 = sin(D31)
    sinD32 = sin(D31-D21)
    triple_sin = sinD21 * sinD31 * sinD32
    sinsqD21_2 = 2 * sinD21 * sinD21
    sinsqD31_2 = 2 * sinD31 * sinD31
    sinsqD32_2 = 2 * sinD32 * sinD32
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 &
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 &
            + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2
    Pme_CPV = -Jvac * triple_sin
    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 + Um3sq * Um1sq * sinsqD31_2 + Um3sq * Um2sq * sinsqD32_2)
    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 + Ue3sq * Ue1sq * sinsqD31_2 + Ue3sq * Ue2sq * sinsqD32_2)
    probs_returned(1, 1) = Pee
    probs_returned(1, 2) = Pme_CPC - Pme_CPV
    probs_returned(1, 3) = 1 - Pee - probs_returned(1, 2)
    probs_returned(2, 1) = Pme_CPC + Pme_CPV
    probs_returned(2, 2) = Pmm
    probs_returned(2, 3) = 1 - probs_returned(2, 1) - Pmm
    probs_returned(3, 1) = 1 - Pee - probs_returned(2, 1)
    probs_returned(3, 2) = 1 - probs_returned(1, 2) - Pmm
    probs_returned(3, 3) = 1 - probs_returned(1, 3) - probs_returned(2, 3)
end subroutine Probability_Vacuum_LBL

subroutine Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
    use Parameters
    implicit none
    real(8), intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye
    integer, intent(in) :: N_Newton
    real(8), intent(out) :: probs_returned(3, 3)

    real(8) :: c13sq, sind, cosd, Jrr, Jmatter, Dmsqee, Amatter
    real(8) :: Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq
    real(8) :: A, B, C
    real(8) :: See, Tee, Smm, Tmm
    real(8) :: xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32
    real(8) :: Xp2, Xp3, PiDlambdaInv
    real(8) :: Lover4E, D21, D32
    real(8) :: sinD21, sinD31, sinD32
    real(8) :: sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin
    real(8) :: Pme_CPC, Pme_CPV, Pmm, Pee
    integer :: i

    c13sq = 1 - s13sq
    Ue2sq = c13sq * s12sq
    Ue3sq = s13sq
    Um3sq = c13sq * s23sq
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1 - s12sq) * (1 - s23sq)
    Jrr = sqrt(Um2sq * Ut2sq)
    sind = sin(delta)
    cosd = cos(delta)
    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd
    Jmatter = 8 * Jrr * c13sq * sind
    Amatter = Ye * rho * E * YerhoE2a
    Dmsqee = Dmsq31 - s12sq * Dmsq21
    A = Dmsq21 + Dmsq31
    See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq
    Tmm = Dmsq21 * Dmsq31
    Tee = Tmm * (1 - Ue3sq - Ue2sq)
    C = Amatter * Tee
    A = A + Amatter
    xmat = Amatter / Dmsqee
    lambda3 = Dmsq31 + 0.5d0 * Dmsqee * (xmat - 1 + sqrt((1 - xmat) ** 2 + 4 * s13sq * xmat))
    B = Tmm + Amatter * See
    do i = 1, N_Newton
        lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B)
    enddo
    Dlambda21 = sqrt((A - lambda3) ** 2 - 4 * C / lambda3)
    lambda2 = 0.5d0 * (A - lambda3 + Dlambda21)
    Dlambda32 = lambda3 - lambda2
    Dlambda31 = Dlambda32 + Dlambda21
    PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21)
    Xp3 = PiDlambdaInv * Dlambda21
    Xp2 = -PiDlambdaInv * Dlambda31
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2
    Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq
    Tmm = Tmm * (1 - Um3sq - Um2sq) + Amatter * (See + Smm - A)
    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv
    Ue1sq = 1 - Ue3sq - Ue2sq
    Um1sq = 1 - Um3sq - Um2sq
    Ut3sq = 1 - Um3sq - Ue3sq
    Ut2sq = 1 - Um2sq - Ue2sq
    Ut1sq = 1 - Um1sq - Ue1sq
    Lover4E = eVsqkm_to_GeV_over4 * L / E
    D21 = Dlambda21 * Lover4E
    D32 = Dlambda32 * Lover4E
    sinD21 = sin(D21)
    sinD31 = sin(D32 + D21)
    sinD32 = sin(D32)
    triple_sin = sinD21 * sinD31 * sinD32
    sinsqD21_2 = 2 * sinD21 * sinD21
    sinsqD31_2 = 2 * sinD31 * sinD31
    sinsqD32_2 = 2 * sinD32 * sinD32
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 &
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 &
            + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2
    Pme_CPV = -Jmatter * triple_sin
    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 + Um3sq * Um1sq * sinsqD31_2 + Um3sq * Um2sq * sinsqD32_2)
    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 + Ue3sq * Ue1sq * sinsqD31_2 + Ue3sq * Ue2sq * sinsqD32_2)
    probs_returned(1, 1) = Pee
    probs_returned(1, 2) = Pme_CPC - Pme_CPV
    probs_returned(1, 3) = 1 - Pee - probs_returned(1, 2)
    probs_returned(2, 1) = Pme_CPC + Pme_CPV
    probs_returned(2, 2) = Pmm
    probs_returned(2, 3) = 1 - probs_returned(2, 1) - Pmm
    probs_returned(3, 1) = 1 - Pee - probs_returned(2, 1)
    probs_returned(3, 2) = 1 - probs_returned(1, 2) - Pmm
    probs_returned(3, 3) = 1 - probs_returned(1, 3) - probs_returned(2, 3)
end subroutine Probability_Matter_LBL

! Use result to prevent optimization
real(8) function use_result(probs)
    real(8), intent(in) :: probs(3,3)
    use_result = probs(1,1) + probs(2,1) + probs(3,3)
end function use_result

program benchmark
    use Parameters
    implicit none
    external :: Probability_Matter_LBL, Probability_Vacuum_LBL
    real(8), external :: use_result
    real(8) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, rho, Ye, E
    real(8) :: probs_returned(3, 3)
    real(8) :: Emin, Emax, total
    integer :: N_Newton, n, i
    integer(8) :: t1, t2, count_rate
    real(8) :: elapsed
    
    n = 10000000  ! 10M iterations
    Emin = 0.5d0
    Emax = 5.0d0
    
    s12sq = 0.31d0
    s13sq = 0.02d0
    s23sq = 0.55d0
    delta = -0.7d0 * PI
    Dmsq21 = 7.5e-5_8
    Dmsq31 = 2.5e-3_8
    L = 1300.0d0
    rho = 3.0d0
    Ye = 0.5d0
    
    call system_clock(count_rate=count_rate)
    
    write (*,"(A,I0,A)") "NuFast Fortran Benchmark (n=", n, " iterations)"
    write (*,*) "============================================"
    write (*,*)
    
    ! Warm up
    total = 0.0d0
    do i = 1, 100000
        E = Emin + (Emax - Emin) * mod(i, 1000) / 1000.0d0
        call Probability_Vacuum_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, probs_returned)
        total = total + use_result(probs_returned)
    end do
    if (total < -1d100) write(*,*) "never"
    
    write (*,*) "Single-point calculations:"
    
    ! Vacuum benchmark
    total = 0.0d0
    call system_clock(t1)
    do i = 1, n
        E = Emin + (Emax - Emin) * mod(i, 1000) / 1000.0d0
        call Probability_Vacuum_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, probs_returned)
        total = total + use_result(probs_returned)
    end do
    call system_clock(t2)
    elapsed = real(t2 - t1, 8) / real(count_rate, 8) * 1e9 / n
    if (total < -1d100) write(*,*) "never"
    write (*,"(A,F10.2,A)") "  Vacuum:            ", elapsed, " ns"
    
    ! Matter N_Newton=0
    N_Newton = 0
    total = 0.0d0
    call system_clock(t1)
    do i = 1, n
        E = Emin + (Emax - Emin) * mod(i, 1000) / 1000.0d0
        call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
        total = total + use_result(probs_returned)
    end do
    call system_clock(t2)
    elapsed = real(t2 - t1, 8) / real(count_rate, 8) * 1e9 / n
    if (total < -1d100) write(*,*) "never"
    write (*,"(A,F10.2,A)") "  Matter N_Newton=0: ", elapsed, " ns"
    
    ! Matter N_Newton=1
    N_Newton = 1
    total = 0.0d0
    call system_clock(t1)
    do i = 1, n
        E = Emin + (Emax - Emin) * mod(i, 1000) / 1000.0d0
        call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
        total = total + use_result(probs_returned)
    end do
    call system_clock(t2)
    elapsed = real(t2 - t1, 8) / real(count_rate, 8) * 1e9 / n
    if (total < -1d100) write(*,*) "never"
    write (*,"(A,F10.2,A)") "  Matter N_Newton=1: ", elapsed, " ns"
    
    ! Matter N_Newton=2
    N_Newton = 2
    total = 0.0d0
    call system_clock(t1)
    do i = 1, n
        E = Emin + (Emax - Emin) * mod(i, 1000) / 1000.0d0
        call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
        total = total + use_result(probs_returned)
    end do
    call system_clock(t2)
    elapsed = real(t2 - t1, 8) / real(count_rate, 8) * 1e9 / n
    if (total < -1d100) write(*,*) "never"
    write (*,"(A,F10.2,A)") "  Matter N_Newton=2: ", elapsed, " ns"
    
    ! Matter N_Newton=3
    N_Newton = 3
    total = 0.0d0
    call system_clock(t1)
    do i = 1, n
        E = Emin + (Emax - Emin) * mod(i, 1000) / 1000.0d0
        call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
        total = total + use_result(probs_returned)
    end do
    call system_clock(t2)
    elapsed = real(t2 - t1, 8) / real(count_rate, 8) * 1e9 / n
    if (total < -1d100) write(*,*) "never"
    write (*,"(A,F10.2,A)") "  Matter N_Newton=3: ", elapsed, " ns"
    
end program benchmark
