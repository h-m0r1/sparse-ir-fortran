program main
    use sparse_ir
    use sparse_ir_io
    use sparse_ir_preset
    implicit none

    call test_fermion(.false., .false., .false., .false.)
    call test_boson  (.false., .false., .false., .false.)
    call test_fermion(.false., .true., .false., .false.)
    call test_boson  (.false., .true., .false., .false.)
    call test_fermion(.false., .true., .false., .true.)
    call test_boson  (.false., .true., .false., .true.)
    call test_fermion(.false., .true., .true., .false.)
    call test_boson  (.false., .true., .true., .false.)
    call test_fermion(.false., .true., .true., .true.)
    call test_boson  (.false., .true., .true., .true.)

    contains

    ! fermion
    subroutine test_fermion(preset, positive_only, lflag_gtau, lflag_gl)
        logical, intent(in) :: preset
        logical, intent(in) :: positive_only
        logical, intent(in) :: lflag_gtau
        logical, intent(in) :: lflag_gl
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
        integer n, t

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
        write(*,*) "test_fermion0"
        ! From Matsubara
        do n = 1, ir_obj%nfreq_f
            giv(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_f(n)/beta, kind(0d0)) - omega0)
        end do
        if (lflag_gl) then
            call fit_matsubara_f(ir_obj, giv, gl_matsu_d)
        else
            call fit_matsubara_f(ir_obj, giv, gl_matsu)
        end if
        write(*,*) "test_fermion1"

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 + exp(-beta * omega0))
            gtau_d(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 + exp(-beta * omega0))
        end do
        if (lflag_gtau .and. lflag_gl) then
            call fit_tau(ir_obj, gtau_d, gl_tau_d)
        elseif ((.not. lflag_gtau) .and. lflag_gl) then
            call fit_tau(ir_obj, gtau, gl_tau_d)
        elseif (lflag_gtau .and. (.not. lflag_gl)) then
            call fit_tau(ir_obj, gtau_d, gl_tau)
        else
            call fit_tau(ir_obj, gtau, gl_tau)
        end if
        write(*,*) "test_fermion2"
        !do l = 1, ir_obj%size
            !write(*,*) real(gl_matsu(1,l)), real(gl_tau(1,l))
        !end do
        if (lflag_gl) then
            do l = 1, ir_obj%size
                write(*,*) gl_matsu_d(:,l)
                write(*,*) gl_tau_d(:,l)
            end do
            if (maxval(abs(gl_matsu_d - gl_tau_d)) > 1d2*eps) then
                write(*,*) "gl_matsu and gl_tau do not match!"
                stop 1
            end if
        else
            do l = 1, ir_obj%size
                write(*,*) real(gl_matsu(:,l))
                write(*,*) real(gl_tau(:,l))
            end do
            if (maxval(abs(gl_matsu - gl_tau)) > 1d2*eps) then
                write(*,*) "gl_matsu and gl_tau do not match!"
                stop 1
            end if
        end if
        if (lflag_gl) then
            call evaluate_matsubara_f(ir_obj, gl_matsu_d, giv_reconst)
        else
            call evaluate_matsubara_f(ir_obj, gl_matsu, giv_reconst)
        end if
        if (maxval(abs(giv - giv_reconst)) > 1d2*eps) then
            write(*,*) "giv do not match!"
            stop 1
        end if
        write(*,*) "test_fermion3"
        if (lflag_gtau .and. lflag_gl) then
            call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst_d)
        elseif ((.not. lflag_gtau) .and. lflag_gl) then
            call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst)
            gtau_reconst_d = real(gtau_reconst, kind(0d0))
        elseif (lflag_gtau .and. (.not. lflag_gl)) then
            call evaluate_tau(ir_obj, gl_tau, gtau_reconst_d)
        else
            call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
            gtau_reconst_d = real(gtau_reconst, kind(0d0))
        end if
        if (maxval(abs(gtau_d - gtau_reconst_d)) > 1d2*eps) then
            write(*,*) "gtau do not match!"
            stop 1
        end if

        write(*,*) "test_fermion"
        write(*,*) "positive_only = ", positive_only
        write(*,*) "lflag_gtau = ", lflag_gtau
        write(*,*) "lflag_gl = ", lflag_gl

        deallocate(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        deallocate(gtau_d, gl_matsu_d, gl_tau_d, gtau_reconst_d)

        call finalize_ir(ir_obj)
    end subroutine


    ! boson
    subroutine test_boson(preset, positive_only, lflag_gtau, lflag_gl)
        logical, intent(in) :: preset
        logical, intent(in) :: positive_only
        logical, intent(in) :: lflag_gtau
        logical, intent(in) :: lflag_gl
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
        integer n, t
    
        PI =4.D0*DATAN(1.D0)
    
        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta)
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
        if (lflag_gl) then
            call fit_matsubara_b(ir_obj, giv, gl_matsu_d)
        else
            call fit_matsubara_b(ir_obj, giv, gl_matsu)
        end if
    
        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 - exp(-beta * omega0))
            gtau_d(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 - exp(-beta * omega0))
        end do
        if (lflag_gtau .and. lflag_gl) then
            call fit_tau(ir_obj, gtau_d, gl_tau_d)
        elseif ((.not. lflag_gtau) .and. lflag_gl) then
            call fit_tau(ir_obj, gtau, gl_tau_d)
        elseif (lflag_gtau .and. (.not. lflag_gl)) then
            call fit_tau(ir_obj, gtau_d, gl_tau)
        else
            call fit_tau(ir_obj, gtau, gl_tau)
        end if
    
        !do l = 1, ir_obj%size
            !write(*,*) real(gl_matsu(1,l)), real(gl_tau(1,l))
        !end do
        if (lflag_gl) then
            if (maxval(abs(gl_matsu_d - gl_tau_d)) > 1d2*eps) then
                write(*,*) "gl_matsu and gl_tau do not match!"
                stop 1
            end if
        else
            if (maxval(abs(gl_matsu - gl_tau)) > 1d2*eps) then
                write(*,*) "gl_matsu and gl_tau do not match!"
                stop 1
            end if
        end if
        if (lflag_gl) then
            call evaluate_matsubara_b(ir_obj, gl_matsu_d, giv_reconst)
        else
            call evaluate_matsubara_b(ir_obj, gl_matsu, giv_reconst)
        end if
        if (maxval(abs(giv - giv_reconst)) > 1d2*eps) then
            write(*,*) "giv do not match!"
            stop 1
        end if
    
        if (lflag_gtau .and. lflag_gl) then
            call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst_d)
        elseif ((.not. lflag_gtau) .and. lflag_gl) then
            call evaluate_tau(ir_obj, gl_tau_d, gtau_reconst)
            gtau_reconst_d = real(gtau_reconst, kind(0d0))
        elseif (lflag_gtau .and. (.not. lflag_gl)) then
            call evaluate_tau(ir_obj, gl_tau, gtau_reconst_d)
        else
            call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
            gtau_reconst_d = real(gtau_reconst, kind(0d0))
        end if
        if (maxval(abs(gtau_d - gtau_reconst_d)) > 1d2*eps) then
            write(*,*) "gtau do not match!"
            stop 1
        end if
    
        write(*,*) "test_boson"
        write(*,*) "positive_only = ", positive_only
        write(*,*) "lflag_gtau = ", lflag_gtau
        write(*,*) "lflag_gl = ", lflag_gl
    
        deallocate(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        deallocate(gtau_d, gl_matsu_d, gl_tau_d, gtau_reconst_d)
    
        call finalize_ir(ir_obj)
    end subroutine

end program
