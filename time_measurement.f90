program main
    use sparse_ir
    use sparse_ir_io
    use sparse_ir_preset
    implicit none
    integer, parameter :: num = 1000

    call time_fermion(.false., .false., num)
    call time_boson  (.false., .false., num)
    call time_fermion(.false., .true., num)
    call time_boson  (.false., .true., num)

    contains

    ! fermion
    subroutine time_fermion(preset, positive_only, num)
        logical, intent(in) :: preset
        logical, intent(in) :: positive_only
        integer, intent(in) :: num
        type(IR) :: ir_obj
        integer, parameter :: ndigit = 10, nlambda = 4
        double precision, parameter :: lambda = 1.d1 ** nlambda
        double precision, parameter :: wmax = 1.d0
        double precision :: PI

        double precision, parameter :: beta = lambda/wmax, omega0 = 1.d0/beta
        double precision, parameter :: eps = 1.d-1**ndigit

        complex(kind(0d0)),allocatable :: giv(:,:), gl_matsu(:, :), gl_tau(:, :), gtau(:, :), &
            gtau_reconst(:, :), giv_reconst(:, :)
        double precision,allocatable :: gl_matsu_d(:, :), gl_tau_d(:, :), gtau_d(:, :), &
            gtau_reconst_d(:, :)
        integer n, t, l, i
        integer:: time_begin_c, time_end_c, CountPerSec, CountMax

        PI =4.D0*DATAN(1.D0)

        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta, positive_only)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta, positive_only)
            close(99)
        end if

        if (abs(ir_obj%beta - beta) > 1d-10) then
            write(*,*) "beta does not match"
            stop 1
        end if
        if (abs(ir_obj%wmax - wmax) > 1d-10) then
            write(*,*) "wmax does not match"
            stop 1
        end if

        ! With ω0 = 1/β,
        !   G(iv) = 1/(iv - ω0),
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        allocate(giv(1, ir_obj%nfreq_f))
        allocate(giv_reconst(1, ir_obj%nfreq_f))
        allocate(gtau(1, ir_obj%ntau))
        allocate(gl_matsu(1, ir_obj%size))
        allocate(gl_tau(1, ir_obj%size))
        allocate(gtau_reconst(1, ir_obj%ntau))
        allocate(gtau_d(1, ir_obj%ntau))
        allocate(gl_matsu_d(1, ir_obj%size))
        allocate(gl_tau_d(1, ir_obj%size))
        allocate(gtau_reconst_d(1, ir_obj%ntau))

        ! From Matsubara
        do n = 1, ir_obj%nfreq_f
            giv(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_f(n)/beta, kind(0d0)) - omega0)
        end do

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 + exp(-beta * omega0))
            gtau_d(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 + exp(-beta * omega0))
        end do

        call system_clock(time_begin_c, CountPerSec, CountMax)
        if (positive_only) then
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_matsubara_f(ir_obj, giv, gl_matsu_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_matsubara_f_zd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_matsubara_f(ir_obj, giv, gl_matsu)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_matsubara_f_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau_d, gl_tau_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_dd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau_d, gl_tau)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_dz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau, gl_tau_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_zd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau, gl_tau)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_zz)"
            call sleep(2)
        
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_matsubara_f(ir_obj, gl_matsu_d, giv_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_matsubara_f_dz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_matsubara_f(ir_obj, gl_matsu, giv_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_matsubara_f_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_dd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_dz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau, gtau_reconst_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_zd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_zz)"
            call sleep(2)
        else
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_matsubara_f(ir_obj, giv, gl_matsu)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_matsubara_f_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau, gl_tau)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_matsubara_f(ir_obj, gl_matsu, giv_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_matsubara_f_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_zz)"
            call sleep(2)
        end if

        deallocate(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        deallocate(gtau_d, gl_matsu_d, gl_tau_d, gtau_reconst_d)

        call finalize_ir(ir_obj)
    end subroutine


    ! boson
    subroutine time_boson(preset, positive_only, num)
        logical, intent(in) :: preset
        logical, intent(in) :: positive_only
        integer, intent(in) :: num
        type(IR) :: ir_obj
        integer, parameter :: ndigit = 10, nlambda = 4
        double precision, parameter :: lambda = 1.d1 ** nlambda
        double precision, parameter :: wmax = 1.d0
        double precision :: PI
    
        double precision, parameter :: beta = lambda/wmax, omega0 = 1.d0/beta
        double precision, parameter :: eps = 1.d-1**ndigit
    
        complex(kind(0d0)),allocatable :: giv(:,:), gl_matsu(:, :), gl_tau(:, :), gtau(:, :), &
            gtau_reconst(:, :), giv_reconst(:, :)
        double precision,allocatable :: gl_matsu_d(:, :), gl_tau_d(:, :), gtau_d(:, :), &
            gtau_reconst_d(:, :)
        integer n, t, l, i
        integer:: time_begin_c, time_end_c, CountPerSec, CountMax
  
        PI =4.D0*DATAN(1.D0)
  
        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta, positive_only)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta, positive_only)
            close(99)
        end if
  
        if (abs(ir_obj%beta - beta) > 1d-10) then
            write(*,*) "beta does not match"
            stop 1
        end if
        if (abs(ir_obj%wmax - wmax) > 1d-10) then
            write(*,*) "wmax does not match"
            stop 1
        end if
  
        ! With ω0 = 1/β,
        !   G(iv) = 1/(iv - ω0),
        !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
        allocate(giv(1, ir_obj%nfreq_b))
        allocate(giv_reconst(1, ir_obj%nfreq_b))
        allocate(gtau(1, ir_obj%ntau))
        allocate(gl_matsu(1, ir_obj%size))
        allocate(gl_tau(1, ir_obj%size))
        allocate(gtau_reconst(1, ir_obj%ntau))
        allocate(gtau_d(1, ir_obj%ntau))
        allocate(gl_matsu_d(1, ir_obj%size))
        allocate(gl_tau_d(1, ir_obj%size))
        allocate(gtau_reconst_d(1, ir_obj%ntau))

        ! From Matsubara
        do n = 1, ir_obj%nfreq_b
            giv(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_b(n)/beta, kind(0d0)) - omega0)
        end do

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 - exp(-beta * omega0))
            gtau_d(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 - exp(-beta * omega0))
        end do

        call system_clock(time_begin_c, CountPerSec, CountMax)
        if (positive_only) then
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_matsubara_b(ir_obj, giv, gl_matsu_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_matsubara_b_zd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_matsubara_b(ir_obj, giv, gl_matsu)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_matsubara_b_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau_d, gl_tau_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_dd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau_d, gl_tau)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_dz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau, gl_tau_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_zd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau, gl_tau)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_zz)"
            call sleep(2)
        
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_matsubara_b(ir_obj, gl_matsu_d, giv_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_matsubara_b_dz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_matsubara_b(ir_obj, gl_matsu, giv_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_matsubara_b_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_dd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_dz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau, gtau_reconst_d)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_zd)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_zz)"
            call sleep(2)
        else
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_matsubara_b(ir_obj, giv, gl_matsu)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_matsubara_b_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call fit_tau(ir_obj, gtau, gl_tau)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (fit_tau_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_matsubara_b(ir_obj, gl_matsu, giv_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_matsubara_b_zz)"
            call sleep(2)
        
            call sleep(2)
            call system_clock(time_begin_c)
            do i = 1, num
                call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
            end do
            call system_clock(time_end_c)
            write(*,*) real(time_end_c - time_begin_c)/CountPerSec," sec (evaluate_tau_zz)"
            call sleep(2)
        end if

        deallocate(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        deallocate(gtau_d, gl_matsu_d, gl_tau_d, gtau_reconst_d)
  
        call finalize_ir(ir_obj)
    end subroutine

end program
