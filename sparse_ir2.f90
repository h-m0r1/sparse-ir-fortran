module sparse_ir
    implicit none
    complex(kind(0d0)), PARAMETER :: cone  = (1.0d0, 0.0d0)
    complex(kind(0d0)), PARAMETER :: ci  = (0.0d0, 1.0d0)
    complex(kind(0d0)), PARAMETER :: czero  = (0.0d0, 0.0d0)
    double precision, PARAMETER :: one = 1.0d0
    double precision, PARAMETER :: zero = 0.0d0

    interface decompose
        module procedure decompose_z, decompose_d
    end interface decompose

    interface evaluate_tau
        module procedure evaluate_tau_zz, evaluate_tau_dd, evaluate_tau_dz, evaluate_tau_zd
    end interface evaluate_tau

    interface evaluate_matsubara_f
        module procedure evaluate_matsubara_f_zz, evaluate_matsubara_f_dz
    end interface evaluate_matsubara_f

    interface evaluate_matsubara_b
        module procedure evaluate_matsubara_b_zz, evaluate_matsubara_b_dz
    end interface evaluate_matsubara_b

    interface fit_tau
        module procedure fit_tau_zz, fit_tau_dd, fit_tau_dz, fit_tau_zd
    end interface fit_tau

    interface fit_matsubara_f
        module procedure fit_matsubara_f_zz, fit_matsubara_f_zd
    end interface fit_matsubara_f

    interface fit_matsubara_b
        module procedure fit_matsubara_b_zz, fit_matsubara_b_zd
    end interface fit_matsubara_b

    ! Matrix decomposed in SVD for fitting
    type DecomposedMatrix_z
        complex(kind(0d0)), allocatable :: a(:, :) ! Original matrix
        double precision, allocatable :: a_real(:, :) ! Real part of original matrix
        double precision, allocatable :: a_imag(:, :) ! Imaginary part of original matrix
        double precision, allocatable :: a_odd(:, :) ! having the elements from odd rows of original matrix
        double precision, allocatable :: a_even(:, :) ! having the elements from even rows of original matrix
        double precision, allocatable :: inv_s_dl(:) ! Inverse of dimensionless singular values
        double precision, allocatable :: inv_s(:) ! Inverse of singular values
        complex(kind(0d0)), allocatable :: ut(:, :), v(:, :)
        double precision, allocatable :: ut_real(:, :), v_real(:, :)
        double precision, allocatable :: ut_imag(:, :), v_imag(:, :)
        integer :: m, n, ns
    end type

    type DecomposedMatrix_d
        complex(kind(0d0)), allocatable :: a(:, :) ! Original matrix
        double precision, allocatable :: a_real(:, :) ! Original matrix
        double precision, allocatable :: inv_s_dl(:) ! Inverse of dimensionless singular values
        double precision, allocatable :: inv_s(:) ! Inverse of singular values
        complex(kind(0d0)), allocatable :: ut(:, :), v(:, :)
        double precision, allocatable :: ut_real(:, :), v_real(:, :)
        integer :: m, n, ns
    end type

    ! Sampling points, basis functions
    type IR
        integer :: size, ntau, nfreq_f, nfreq_b, nomega
        double precision :: beta, lambda, wmax, eps, eps_svd
        double precision, allocatable :: s(:), tau(:), x(:)
        double precision, allocatable :: omega(:), y(:)
        integer, allocatable :: freq_f(:), freq_b(:)
        double precision, allocatable :: u_data(:,:)
        complex(kind(0d0)), allocatable :: uhat_f_data(:,:), uhat_b_data(:,:)
        double precision, allocatable :: v_data(:,:), dlr_data(:,:)
        type(DecomposedMatrix_d) :: u
        type(DecomposedMatrix_z) :: uhat_f, uhat_b
        !type(DecomposedMatrix_d) :: dlr
        logical :: positive_only
    end type

    contains

    subroutine init_ir(obj, beta, lambda, eps, s, x, freq_f, freq_b, u, uhat_f, uhat_b, y, v, dlr, eps_svd, positive_only)
        type(IR), intent(inout) :: obj
        double precision, intent(in) :: beta, lambda, eps, s(:), x(:), y(:), eps_svd
        double precision, intent(in) :: u(:,:), v(:, :), dlr(:, :)
        complex(kind(0d0)), intent(in) :: uhat_f(:, :), uhat_b(:, :)
        integer, intent(in) :: freq_f(:), freq_b(:)
        logical, intent(in), optional :: positive_only

        if ((.not. present(positive_only)) .or. (.not. positive_only)) then
            obj%positive_only = .false.
        else
            obj%positive_only = .true.
        end if

        if (allocated(obj%x)) then
            stop 'IR%x is already allocated. You should call finalize_ir before recalling init_ir.'
        end if

        obj%size = size(s)
        obj%ntau = size(x)
        if (.not. obj%positive_only) then
            obj%nfreq_f = size(freq_f)
            obj%nfreq_b = size(freq_b)
        else
            obj%nfreq_f = size(freq_f) / 2
            obj%nfreq_b = (size(freq_b)-1) / 2 + 1
        end if
        obj%lambda = lambda
        obj%eps = eps
        obj%eps_svd = eps_svd
        obj%nomega = size(y)

        allocate(obj%x(obj%ntau))
        obj%x = x

        allocate(obj%tau(obj%ntau))

        allocate(obj%s(obj%size))
        obj%s = sqrt(5.0d-1*obj%lambda) * s

        allocate(obj%freq_f(obj%nfreq_f))

        allocate(obj%freq_b(obj%nfreq_b))

        if (.not. obj%positive_only) then
            obj%freq_f = freq_f
            obj%freq_b = freq_b
        else
            obj%freq_f(1:obj%nfreq_f) = freq_f((size(freq_f) / 2 + 1):size(freq_f))
            obj%freq_b(1:obj%nfreq_b) = freq_b(((size(freq_b) + 1) / 2):size(freq_b))
        end if

        allocate(obj%u_data(obj%ntau, obj%size))
        obj%u_data = u

        allocate(obj%uhat_f_data(obj%nfreq_f, obj%size))

        allocate(obj%uhat_b_data(obj%nfreq_b, obj%size))

        if (.not. obj%positive_only) then
            obj%uhat_f_data = uhat_f
            obj%uhat_b_data = uhat_b
        else
            obj%uhat_f_data(1:obj%nfreq_f, :) = uhat_f((size(freq_f) / 2 + 1):size(freq_f), :)
            obj%uhat_b_data(1:obj%nfreq_b, :) = uhat_b(((size(freq_b) + 1) / 2):size(freq_b), :)
        end if

        allocate(obj%y(obj%nomega))
        obj%y = y

        allocate(obj%omega(obj%nomega))

        allocate(obj%v_data(obj%nomega, obj%size))
        obj%v_data = v

        allocate(obj%dlr_data(obj%size, obj%nomega))
        obj%dlr_data = transpose(dlr)

        obj%u = decompose(obj%u_data, obj%eps_svd, .false.)
        if (.not. obj%positive_only) then
            obj%uhat_f = decompose(obj%uhat_f_data, obj%eps_svd, .false.)
            obj%uhat_b = decompose(obj%uhat_b_data, obj%eps_svd, .false.)
        else
            obj%uhat_f = split_decompose(obj%uhat_f_data, .false., obj%eps_svd, .false.)
            obj%uhat_b = split_decompose(obj%uhat_b_data, .true., obj%eps_svd, .false.)
        end if

        !obj%dlr = decompose(obj%dlr_data, obj%eps_svd, .true.)

        ! Here we define basis sets for the input value of beta. 
        call set_beta(obj, beta)
    end subroutine

    subroutine set_beta(obj, beta)
        type(IR), intent(inout) :: obj
        double precision, intent(in) :: beta

        obj%beta = beta
        obj%wmax = obj%lambda / beta

        obj%tau = 5.0d-1 * beta * (obj%x + 1.d0)
        obj%omega = obj%y * obj%wmax

        obj%u%a_real(:, :) = sqrt(2.0d0/beta)*obj%u_data(:, :)
        obj%uhat_f%a(:, :) = sqrt(beta) * obj%uhat_f_data(:, :)
        obj%uhat_b%a(:, :) = sqrt(beta) * obj%uhat_b_data(:, :)
        !obj%dlr%a(:, :) = sqrt(5.0d-1*beta)*obj%dlr_data(:, :)

        obj%u%a = cmplx(obj%u%a_real, zero, kind(0d0))
        obj%uhat_f%a_real = REAL(obj%uhat_f%a, KIND(0d0))
        obj%uhat_f%a_imag = AIMAG(obj%uhat_f%a)
        obj%uhat_b%a_real = REAL(obj%uhat_b%a, KIND(0d0))
        obj%uhat_b%a_imag = AIMAG(obj%uhat_b%a)
        obj%uhat_f%a_odd(:, :) = obj%uhat_f%a_imag(:, 1:(obj%uhat_f%n - 1):2)
        obj%uhat_f%a_even(:, :) = obj%uhat_f%a_real(:, 2:obj%uhat_f%n:2)
        obj%uhat_b%a_odd(:, :) = obj%uhat_b%a_real(:, 1:(obj%uhat_b%n - 1):2)
        obj%uhat_b%a_even(:, :) = obj%uhat_b%a_imag(:, 2:obj%uhat_b%n:2)

        obj%u%inv_s(:) = sqrt(5.0d-1*beta) * obj%u%inv_s_dl(:)
        obj%uhat_f%inv_s(:) = (1.0d0 / sqrt(beta)) * obj%uhat_f%inv_s_dl(:)
        obj%uhat_b%inv_s(:) = (1.0d0 / sqrt(beta)) * obj%uhat_b%inv_s_dl(:)
        !obj%dlr%inv_s(:) = sqrt(2.0d0 / beta) * obj%dlr%inv_s_dl(:)

    end subroutine

    subroutine finalize_ir(obj)
        type(IR) :: obj
    
        if (allocated(obj%x)) deallocate(obj%x)
        if (allocated(obj%tau)) deallocate(obj%tau)
        if (allocated(obj%s)) deallocate(obj%s)
        if (allocated(obj%freq_f)) deallocate(obj%freq_f)
        if (allocated(obj%freq_b)) deallocate(obj%freq_b)
        if (allocated(obj%u_data)) deallocate(obj%u_data)
        if (allocated(obj%uhat_f_data)) deallocate(obj%uhat_f_data)
        if (allocated(obj%uhat_b_data)) deallocate(obj%uhat_b_data)
        if (allocated(obj%y)) deallocate(obj%y)
        if (allocated(obj%omega)) deallocate(obj%omega)
        if (allocated(obj%v_data)) deallocate(obj%v_data)
        if (allocated(obj%dlr_data)) deallocate(obj%dlr_data)
    
        call finalize_dmat_d(obj%u)
        call finalize_dmat_z(obj%uhat_f)
        call finalize_dmat_z(obj%uhat_b)
        !call finalize_dmat_d(obj%dlr)
    end subroutine

    subroutine finalize_dmat_z(dmat)
        type(DecomposedMatrix_z) :: dmat
    
        if (allocated(dmat%a)) deallocate(dmat%a)
        if (allocated(dmat%a_real)) deallocate(dmat%a_real)
        if (allocated(dmat%a_imag)) deallocate(dmat%a_imag)
        if (allocated(dmat%a_odd)) deallocate(dmat%a_odd)
        if (allocated(dmat%a_even)) deallocate(dmat%a_even)
        if (allocated(dmat%inv_s)) deallocate(dmat%inv_s)
        if (allocated(dmat%inv_s)) deallocate(dmat%inv_s_dl)
        if (allocated(dmat%ut)) deallocate(dmat%ut)
        if (allocated(dmat%v)) deallocate(dmat%v)
        if (allocated(dmat%ut_real)) deallocate(dmat%ut_real)
        if (allocated(dmat%v_real)) deallocate(dmat%v_real)
        if (allocated(dmat%ut_imag)) deallocate(dmat%ut_imag)
        if (allocated(dmat%v_imag)) deallocate(dmat%v_imag)
    end subroutine

    subroutine finalize_dmat_d(dmat)
        type(DecomposedMatrix_d) :: dmat
    
        if (allocated(dmat%a)) deallocate(dmat%a)
        if (allocated(dmat%a_real)) deallocate(dmat%a_real)
        if (allocated(dmat%inv_s)) deallocate(dmat%inv_s)
        if (allocated(dmat%inv_s)) deallocate(dmat%inv_s_dl)
        if (allocated(dmat%ut)) deallocate(dmat%ut)
        if (allocated(dmat%v)) deallocate(dmat%v)
        if (allocated(dmat%ut_real)) deallocate(dmat%ut_real)
        if (allocated(dmat%v_real)) deallocate(dmat%v_real)
    end subroutine

    function decompose_z(a, eps, ill_conditioned) result(dmat)
        complex(kind(0d0)), intent(in) :: a(:, :)
        double precision, intent(in) :: eps
        logical, intent(in) :: ill_conditioned
      
        integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns
        complex(kind(0d0)), allocatable :: a_copy(:, :), u(:, :), &
            vt(:, :), work(:)
        double precision, allocatable :: rwork(:), s(:)
        integer, allocatable :: iwork(:)
        type(DecomposedMatrix_z)::dmat
      
        if (allocated(dmat%a)) then
            stop 'DMAT%a is already allocated. You should call finalize_dmat before recalling decompose.'
        end if
      
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        lda = m
        ldu = m
        ldvt = n
        if (.not. ill_conditioned) then
            lwork = mn*mn + 3*mn
            allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork((5*mn+7)*mn), iwork(8*mn))
        else
            lwork = 2*mn + m + n
            allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork(5*n))
        endif
      
        a_copy(1:m, 1:n) = a(1:m, 1:n)
        if (.not. ill_conditioned) then
            lwork = mn*mn + 3*mn
            call zgesdd('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
            if (info /= 0) then
                stop 'Failure in ZGESDD.'
            end if
        else
            lwork = 2*mn + m + n
            call zgesvd('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
            if (info /= 0) then
                stop 'Failure in ZGESVD.'
            end if
        endif
      
        ! Number of relevant singular values s(i)/s(1) >= eps
        ns = 0
        do i = 1, mn
            if (s(i)/s(1) < eps) then
                exit
            end if
            ns = ns + 1
        end do
      
        allocate(dmat%a(m, n))
        allocate(dmat%a_real(m, n))
        allocate(dmat%a_imag(m, n))
        allocate(dmat%inv_s_dl(ns))
        allocate(dmat%inv_s(ns))
        allocate(dmat%ut(ns, m))
        allocate(dmat%ut_real(ns, m))
        allocate(dmat%ut_imag(ns, m))
        allocate(dmat%v(n, ns))
        allocate(dmat%v_real(n, ns))
        allocate(dmat%v_imag(n, ns))
      
        ! dmat%a temporarily stores the same data of input a
        dmat%a = a
        dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
        ! inv_s temporarily stores the same data of inv_s_dl
        dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
        dmat%ut(1:ns, 1:m) = conjg(transpose(u(1:m, 1:ns)))
        dmat%v(1:n, 1:ns) = conjg(transpose(vt(1:ns, 1:n)))
        dmat%m = size(a, 1)
        dmat%n = size(a, 2)
        dmat%ns = ns
    
        dmat%a_real = REAL(dmat%a, KIND(0d0))
        dmat%a_imag = AIMAG(dmat%a)
        dmat%ut_real = REAL(dmat%ut, KIND(0d0))
        dmat%ut_imag = AIMAG(dmat%ut)
        dmat%v_real = REAL(dmat%v, KIND(0d0))
        dmat%v_imag = AIMAG(dmat%v)
      
        if (.not. ill_conditioned) then
            deallocate(work, a_copy, s, u, vt, rwork, iwork)
        else
            deallocate(work, a_copy, s, u, vt, rwork)
        endif
    end function
    
    function decompose_d(a, eps, ill_conditioned) result(dmat)
        double precision, intent(in) :: a(:, :)
        double precision, intent(in) :: eps
        logical, intent(in) :: ill_conditioned
      
        integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns
        double precision, allocatable :: a_copy(:, :), u(:, :), &
            vt(:, :), work(:), s(:)
        integer, allocatable :: iwork(:)
        type(DecomposedMatrix_d)::dmat
      
        if (allocated(dmat%a)) then
            stop 'DMAT%a is already allocated. You should call finalize_dmat before recalling decompose.'
        end if
      
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        lda = m
        ldu = m
        ldvt = n
        if (.not. ill_conditioned) then
            lwork = 4*mn*mn + 7*mn
            allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), iwork(8*mn))
        else
            lwork = 5*mn + m + n
            allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n))
        endif
      
        a_copy(1:m, 1:n) = a(1:m, 1:n)
        if (.not. ill_conditioned) then
            lwork = 4*mn*mn + 7*mn
            call dgesdd('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            if (info /= 0) then
                stop 'Failure in DGESDD.'
            end if
        else
            lwork = 5*mn + m + n
            call dgesvd('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            if (info /= 0) then
                stop 'Failure in DGESVD.'
            end if
        endif
      
        ! Number of relevant singular values s(i)/s(1) >= eps
        ns = 0
        do i = 1, mn
            if (s(i)/s(1) < eps) then
                exit
            end if
            ns = ns + 1
        end do
      
        allocate(dmat%a(m, n))
        allocate(dmat%a_real(m, n))
        allocate(dmat%inv_s_dl(ns))
        allocate(dmat%inv_s(ns))
        allocate(dmat%ut(ns, m))
        allocate(dmat%ut_real(ns, m))
        allocate(dmat%v(n, ns))
        allocate(dmat%v_real(n, ns))
      
        ! dmat%a temporarily stores the same data of input a
        dmat%a_real = a
        dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
        ! inv_s temporarily stores the same data of inv_s_dl
        dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
        dmat%ut_real(1:ns, 1:m) = transpose(u(1:m, 1:ns))
        dmat%v_real(1:n, 1:ns) = transpose(vt(1:ns, 1:n))
        dmat%m = size(a, 1)
        dmat%n = size(a, 2)
        dmat%ns = ns
    
        dmat%a = dmat%a_real
        dmat%ut = dmat%ut_real
        dmat%v = dmat%v_real
      
        if (.not. ill_conditioned) then
            deallocate(work, a_copy, s, u, vt, iwork)
        else
            deallocate(work, a_copy, s, u, vt)
        endif
    end function
    
    function split_decompose(a, has_zero, eps,ill_conditioned) result(dmat)
        complex(kind(0d0)), intent(in) :: a(:, :)
        logical, intent(in) :: has_zero
        double precision, intent(in) :: eps
        logical, intent(in) :: ill_conditioned
      
        integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns, m_half
        double precision, allocatable :: a_copy(:, :), u(:, :), &
            vt(:, :), work(:), s(:)
        complex(kind(0d0)), allocatable :: u_copy(:, :)
        integer, allocatable :: iwork(:)
        type(DecomposedMatrix_z)::dmat
      
        if (allocated(dmat%a)) then
            stop 'DMAT%a is already allocated. You should call finalize_dmat before recalling decompose.'
        end if
        
        m_half = size(a, 1)
        if (has_zero) then
            m = 2 * m_half - 1
        else
            m = 2 * m_half
        end if
    
        n = size(a, 2)
        mn = min(m, n)
        lda = m
        ldu = m
        ldvt = n
        if (.not. ill_conditioned) then
            lwork = 4*mn*mn + 7*mn
            allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), iwork(8*mn))
        else
            lwork = 5*mn + m + n
            allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n))
        endif
    
        a_copy(1:m_half, 1:n) = real(a(1:m_half, 1:n))
    
        if (has_zero) then
            a_copy(m_half+1:m, 1:n) = aimag(a(2:m_half, 1:n))
        else
            a_copy(m_half+1:m, 1:n) = aimag(a(1:m_half, 1:n))
        end if
    
        if (.not. ill_conditioned) then
            lwork = 4*mn*mn + 7*mn
            call dgesdd('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            if (info /= 0) then
                stop 'Failure in DGESDD.'
            end if
        else
            lwork = 5*mn + m + n
            call dgesvd('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            if (info /= 0) then
                stop 'Failure in DGESVD.'
            end if
        end if
      
        ! Number of relevant singular values s(i)/s(1) >= eps
        ns = 0
        do i = 1, mn
            if (s(i)/s(1) < eps) then
                exit
            end if
            ns = ns + 1
        end do
        
        allocate(u_copy(m_half, ns))
    
        if (has_zero) then
            u_copy(1, 1:ns) = cmplx(u(1, 1:ns), 0.0d0, kind(0d0))
            u_copy(2:m_half, 1:ns) = cmplx(u(2:m_half, 1:ns), u(m_half+1:m, 1:ns), kind(0d0))
        else
            u_copy(1:m_half, 1:ns) = cmplx(u(1:m_half, 1:ns), u(m_half+1:m, 1:ns), kind(0d0))
        end if
    
        allocate(dmat%a(m_half, n))
        allocate(dmat%a_real(m_half, n))
        allocate(dmat%a_imag(m_half, n))
        allocate(dmat%inv_s_dl(ns))
        allocate(dmat%inv_s(ns))
        allocate(dmat%ut(ns, m_half))
        allocate(dmat%ut_real(ns, m_half))
        allocate(dmat%ut_imag(ns, m_half))
        allocate(dmat%v(n, ns))
        allocate(dmat%v_real(n, ns))
        allocate(dmat%v_imag(n, ns))
        allocate(dmat%a_odd(m_half, n/2))
        allocate(dmat%a_even(m_half, n/2))
      
        ! dmat%a temporarily stores the same data of input a
        dmat%a = a
        dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
        ! inv_s temporarily stores the same data of inv_s_dl
        dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
        dmat%ut(1:ns, 1:m_half) = conjg(transpose(u_copy(1:m_half, 1:ns)))
        dmat%v_real(1:n, 1:ns) = transpose(vt(1:ns, 1:n))
        dmat%m = size(a, 1)
        dmat%n = size(a, 2)
        dmat%ns = ns
        dmat%a_odd = zero
        dmat%a_even = zero
    
        dmat%a_real = real(dmat%a, kind(0d0))
        dmat%a_imag = aimag(dmat%a)
        dmat%ut_real = real(dmat%ut, kind(0d0))
        dmat%ut_imag = aimag(dmat%ut)
        dmat%v_imag = zero
        dmat%v = cmplx(dmat%v_real, zero, kind(0d0))
      
        if (.not. ill_conditioned) then
            deallocate(work, a_copy, s, u, vt, iwork, u_copy)
        else
            deallocate(work, a_copy, s, u, vt, u_copy)
        endif
    end function
    
    subroutine evaluate_tau_zz(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: arr_tmp(:, :)
        double precision, allocatable :: res_r(:, :)
        double precision, allocatable :: res_i(:, :)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of output array.'
        res(:, :) = czero
        !
        IF (.NOT. obj%positive_only) THEN
            !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
            !           l2, obj%u%a, m, czero, res(:, :), l2)
            !
            allocate(arr_tmp(l1, n))
            allocate(res_r(l2, m))
            allocate(res_i(l2, m))
            ! calculate the real part
            arr_tmp = REAL(arr, KIND(0d0))
            CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
                       l2, obj%u%a_real, m, zero, res_r(:, :), l2)
            !
            ! calculate the imaginary part
            arr_tmp = AIMAG(arr)
            CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
                       l2, obj%u%a_real, m, zero, res_i(:, :), l2)
            res = cmplx(res_r, res_i, kind(0d0))
            deallocate(arr_tmp)
            deallocate(res_r)
            deallocate(res_i)
        ELSE
            allocate(arr_tmp(l1, n))
            allocate(res_r(l2, m))
            ! only calculate the real part
            arr_tmp = REAL(arr, KIND(0d0))
            CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
                       l2, obj%u%a_real, m, zero, res_r(:, :), l2)
            res = cmplx(res_r, zero, kind(0d0))
            deallocate(arr_tmp)
            deallocate(res_r)
        ENDIF
    end subroutine
    
    subroutine evaluate_tau_dd(obj, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent (in) :: arr(:, :)
        double precision, intent(out) :: res(:, :)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'input and output arrays should be complex arrays.'
        res(:, :) = zero
        !
        CALL DGEMM('n', 't', l1, m, n, one, arr(:,:), &
                   l2, obj%u%a_real, m, zero, res(:, :), l2)
    end subroutine
    
    subroutine evaluate_tau_dz(obj, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: res_tmp(:, :)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'input array should be a complex array.'
        res(:, :) = czero
        !
        allocate(res_tmp(l2, m))
        CALL DGEMM('n', 't', l1, m, n, one, arr(:,:), &
                   l2, obj%u%a_real, m, zero, res_tmp(:, :), l2)
        res = cmplx(res_tmp, zero, kind(0d0))
        deallocate(res_tmp)
    end subroutine
    
    subroutine evaluate_tau_zd(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        double precision, intent(out) :: res(:, :)
        double precision, allocatable :: arr_tmp(:, :)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'output array should be a complex array.'
        res(:, :) = zero
        !
        allocate(arr_tmp(l1, n))
        arr_tmp = REAL(arr, KIND(0d0))
        CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
                   l2, obj%u%a_real, m, zero, res(:, :), l2)
        deallocate(arr_tmp)
    end subroutine
    
    subroutine evaluate_matsubara_f_zz(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: res_r(:, :)
        double precision, allocatable :: res_i(:, :)
        double precision, allocatable :: arr_half(:, :)
        integer :: m, n, l1, l2, n_half
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of output array.'
        res(:, :) = czero
        !
        !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
        !           l2, obj%uhat_f%a, m, czero, res(:, :), l2)
        !
        n_half = n / 2
        !
        allocate(res_r(l2, m))
        allocate(res_i(l2, m))
        allocate(arr_half(l1, n_half))
        !
        IF (.NOT. obj%positive_only) THEN
            ! take partial summations over even l
            arr_half(:,:) = real(arr(:,2:n:2))
            ! calculate the real part
            ! Re(arr) * Re(uhat_f)**T
            res_r(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_f%a_even, m, zero, res_r(:, :), l2)
            !
            arr_half(:,:) = aimag(arr(:,2:n:2))
            ! calculate the imaginary part
            ! Im(arr) * Re(uhat_f)**T
            res_i(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_f%a_even, m, zero, res_i(:, :), l2)
            !
            res = cmplx(res_r, res_i, kind(0d0))
            !
            ! take partial summations over odd l
            arr_half(:,:) = aimag(arr(:,1:(n-1):2))
            ! calculate the real part
            ! Im(arr) * Im(uhat_f)**T
            res_r(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_f%a_odd, m, zero, res_r(:, :), l2)
            !
            arr_half(:,:) = real(arr(:,1:(n-1):2))
            ! calculate the imaginary part
            ! Re(arr) * Im(uhat_f)**T
            res_i(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_f%a_odd, m, zero, res_i(:, :), l2)
            !
            res = res + cmplx(-res_r, res_i, kind(0d0))
        ELSE
            ! take partial summations over even l
            arr_half(:,:) = real(arr(:,2:n:2))
            ! calculate the real part
            ! Re(arr) * Re(uhat_f)**T
            res_r(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_f%a_even, m, zero, res_r(:, :), l2)
            !
            ! take partial summations over odd l
            arr_half(:,:) = real(arr(:,1:(n-1):2))
            ! calculate the imaginary part
            ! Re(arr) * Im(uhat_f)**T
            res_i(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_f%a_odd, m, zero, res_i(:, :), l2)
            !
            res = cmplx(res_r, res_i, kind(0d0))
        ENDIF
        deallocate(res_r)
        deallocate(res_i)
    end subroutine
    
    subroutine evaluate_matsubara_f_dz(obj, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: res_r(:, :)
        double precision, allocatable :: res_i(:, :)
        double precision, allocatable :: arr_half(:, :)
        integer :: m, n, l1, l2, n_half
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'input array should be a complex array.'
        res(:, :) = czero
        !
        !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
        !           l2, obj%uhat_f%a, m, czero, res(:, :), l2)
        !
        n_half = n / 2
        !
        allocate(res_r(l2, m))
        allocate(res_i(l2, m))
        allocate(arr_half(l1, n_half))
        !
        ! take partial summations over even l
        arr_half(:,:) = arr(:,2:n:2)
        ! calculate the real part
        ! Re(arr) * Re(uhat_f)**T
        res_r(:, :) = zero
        CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                   l2, obj%uhat_f%a_even, m, zero, res_r(:, :), l2)
        !
        ! take partial summations over odd l
        arr_half(:,:) = arr(:,1:(n-1):2)
        ! calculate the imaginary part
        ! Re(arr) * Im(uhat_f)**T
        res_i(:, :) = zero
        CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                   l2, obj%uhat_f%a_odd, m, zero, res_i(:, :), l2)
        !
        res = cmplx(res_r, res_i, kind(0d0))
        deallocate(res_r)
        deallocate(res_i)
    end subroutine
    
    subroutine evaluate_matsubara_b_zz(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: res_r(:, :)
        double precision, allocatable :: res_i(:, :)
        double precision, allocatable :: arr_half(:, :)
        integer :: m, n, l1, l2, n_half
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of output array.'
        res(:, :) = czero
        !
        !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
        !           l2, obj%uhat_b%a, m, czero, res(:, :), l2)
        !
        n_half = n / 2
        !
        allocate(res_r(l2, m))
        allocate(res_i(l2, m))
        allocate(arr_half(l1, n_half))
        !
        IF (.NOT. obj%positive_only) THEN
            ! take partial summations over odd l
            arr_half(:,:) = real(arr(:,1:(n-1):2))
            ! calculate the real part
            ! Re(arr) * Re(uhat_b)**T
            res_r(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_b%a_odd, m, zero, res_r(:, :), l2)
            !
            arr_half(:,:) = aimag(arr(:,1:(n-1):2))
            ! calculate the imaginary part
            ! Im(arr) * Re(uhat_b)**T
            res_i(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_b%a_odd, m, zero, res_i(:, :), l2)
            !
            res = cmplx(res_r, res_i, kind(0d0))
            !
            ! take partial summations over even l
            arr_half(:,:) = aimag(arr(:,2:n:2))
            ! calculate the real part
            ! Im(arr) * Im(uhat_b)**T
            res_r(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_b%a_even, m, zero, res_r(:, :), l2)
            !
            arr_half(:,:) = real(arr(:,2:n:2))
            ! calculate the imaginary part
            ! Re(arr) * Im(uhat_b)**T
            res_i(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_b%a_even, m, zero, res_i(:, :), l2)
            !
            res = res + cmplx(-res_r, res_i, kind(0d0))
        ELSE
            ! take partial summations over odd l
            arr_half(:,:) = real(arr(:,1:(n-1):2))
            ! calculate the real part
            ! Re(arr) * Re(uhat_b)**T
            res_r(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_b%a_odd, m, zero, res_r(:, :), l2)
            !
            ! take partial summations over even l
            arr_half(:,:) = real(arr(:,2:n:2))
            ! calculate the imaginary part
            ! Re(arr) * Im(uhat_b)**T
            res_i(:, :) = zero
            CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                       l2, obj%uhat_b%a_even, m, zero, res_i(:, :), l2)
            !
            res = cmplx(res_r, res_i, kind(0d0))
        ENDIF
        deallocate(res_r)
        deallocate(res_i)
    end subroutine
    
    subroutine evaluate_matsubara_b_dz(obj, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: res_r(:, :)
        double precision, allocatable :: res_i(:, :)
        double precision, allocatable :: arr_half(:, :)
        integer :: m, n, l1, l2, n_half
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'input array should be a complex array.'
        res(:, :) = czero
        !
        !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
        !           l2, obj%uhat_b%a, m, czero, res(:, :), l2)
        !
        n_half = n / 2
        !
        allocate(res_r(l2, m))
        allocate(res_i(l2, m))
        allocate(arr_half(l1, n_half))
        !
        ! take partial summations over odd l
        arr_half(:,:) = arr(:,1:(n-1):2)
        ! calculate the real part
        ! Re(arr) * Re(uhat_b)**T
        res_r(:, :) = zero
        CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                   l2, obj%uhat_b%a_odd, m, zero, res_r(:, :), l2)
        !
        ! take partial summations over even l
        arr_half(:,:) = arr(:,2:n:2)
        ! calculate the imaginary part
        ! Re(arr) * Im(uhat_b)**T
        res_i(:, :) = zero
        CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
                   l2, obj%uhat_b%a_even, m, zero, res_i(:, :), l2)
        !
        res = cmplx(res_r, res_i, kind(0d0))
        deallocate(res_r)
        deallocate(res_i)
    end subroutine
    
    subroutine fit_tau_zz(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: ut_arr(:, :)
        double precision, allocatable :: arr_tmp(:, :)
        double precision, allocatable :: res_r(:, :)
        double precision, allocatable :: res_i(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%u%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of output array.'
        allocate(ut_arr(ns, l1))
    
        IF (.NOT. obj%positive_only) THEN
            allocate(arr_tmp(l1, m))
            allocate(res_r(l2, n))
            allocate(res_i(l2, n))
            ! calculate the real part
            arr_tmp = REAL(arr, KIND(0d0))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
            do j = 1, ns
                do i = 1, l1
                    ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
                end do
            end do
    
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res_r(:, :) = zero
            call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_r, l2)
            !
            ! calculate the imaginary part
            arr_tmp = AIMAG(arr)
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
            do j = 1, ns
                do i = 1, l1
                    ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
                end do
            end do
    
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res_r(:, :) = zero
            call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_i, l2)
            res = cmplx(res_r, res_i, kind(0d0))
            deallocate(arr_tmp)
            deallocate(res_r)
            deallocate(res_i)
        ELSE
            allocate(arr_tmp(l1, m))
            allocate(res_r(l2, n))
            ! calculate the real part
            arr_tmp = REAL(arr, KIND(0d0))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
            do j = 1, ns
                do i = 1, l1
                    ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
                end do
            end do
    
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res_r(:, :) = zero
            call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_r, l2)
            res = cmplx(res_r, zero, kind(0d0))
            deallocate(arr_tmp)
            deallocate(res_r)
        ENDIF
    
        deallocate(ut_arr)
    end subroutine
    
    subroutine fit_tau_dz(obj, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision, allocatable :: ut_arr(:, :)
        double precision, allocatable :: res_tmp(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%u%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'input array should be a complex array.'
        allocate(res_tmp(l2, n))
        allocate(ut_arr(ns, l1))
    
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr, l1, zero, ut_arr, ns)
        do j = 1, ns
            do i = 1, l1
                ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
            end do
        end do
    
        ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
        res_tmp(:, :) = zero
        call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v, n, zero, res_tmp, l2)
    
        res(:, :) = cmplx(res_tmp(:, :), zero, kind(0d0)) 
        deallocate(ut_arr, res_tmp)
    end subroutine
    
    subroutine fit_tau_zd(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        double precision, intent(out) :: res(:, :)
        double precision, allocatable :: ut_arr(:, :)
        double precision, allocatable :: arr_tmp(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%u%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'output array should be a complex array.'
        allocate(arr_tmp(l1, m))
        arr_tmp(:, :) = real(arr(:, :), kind(0d0))
        allocate(ut_arr(ns, l1))
    
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
        do j = 1, ns
            do i = 1, l1
                ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
            end do
        end do
    
        ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
        call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v, n, zero, res, l2)
    
        deallocate(ut_arr, arr_tmp)
    end subroutine
    
    subroutine fit_tau_dd(obj, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent (in) :: arr(:, :)
        double precision, intent(out) :: res(:, :)
        double precision, allocatable :: ut_arr(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%u%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'input and output arrays should be complex arrays.'
        allocate(ut_arr(ns, l1))
    
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr, l1, zero, ut_arr, ns)
        do j = 1, ns
            do i = 1, l1
                ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
            end do
        end do
    
        ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
        call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res, l2)
    
        deallocate(ut_arr)
    end subroutine
    
    subroutine fit_matsubara_f_zz(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), allocatable :: ut_arr(:, :)
        double precision, allocatable :: arr_tmp(:, :)
        double precision, allocatable :: res_tmp(:, :)
        double precision, allocatable :: ut_arr_r(:, :)
        double precision, allocatable :: ut_arr_tmp(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%uhat_f%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of output array.'
        
        IF (.not. obj%positive_only) then
            allocate(ut_arr(ns, l1))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr(:, :) = czero
            call zgemm("n", "t", ns, l1, m, cone, obj%uhat_f%ut, ns, arr, l1, czero, ut_arr, ns)
            do j = 1, ns
                do i = 1, l1
                    ut_arr(j, i) = ut_arr(j, i) * obj%uhat_f%inv_s(j)
                end do
            end do
    
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res(:, :) = czero
            call zgemm("t", "t", l1, n, ns, cone, ut_arr, ns, obj%uhat_f%v, n, czero, res, l2)
            deallocate(ut_arr)
        ELSE
            allocate(arr_tmp(l1, m))
            allocate(res_tmp(l2, n))
            allocate(ut_arr_r(ns, l1))
            allocate(ut_arr_tmp(ns, l1))
            arr_tmp(:, :) = real(arr(:, :), kind(0d0))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr_r(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%uhat_f%ut_real, ns, arr_tmp, l1, zero, ut_arr_r, ns)
        
            arr_tmp(:, :) = aimag(arr(:, :))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr_tmp(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%uhat_f%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
        
            ut_arr_r(:, :) = ut_arr_r(:, :) + ut_arr_tmp(:, :)
            do j = 1, ns
                do i = 1, l1
                    ut_arr_r(j, i) = ut_arr_r(j, i) * obj%uhat_f%inv_s(j)
                end do
            end do
        
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res_tmp(:, :) = zero
            call dgemm("t", "t", l1, n, ns, one, ut_arr_r, ns, obj%uhat_f%v_real, n, zero, res_tmp, l2)

            res = cmplx(res_tmp, zero, kind(0d0))
        
            deallocate(ut_arr_r, ut_arr_tmp, arr_tmp, res_tmp)
        ENDIF
    
    end subroutine

    subroutine fit_matsubara_f_zd(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        double precision, intent(out) :: res(:, :)
        double precision, allocatable :: ut_arr(:, :)
        double precision, allocatable :: ut_arr_tmp(:, :)
        double precision, allocatable :: arr_tmp(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%uhat_f%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'output array should be a complex array.'
        allocate(arr_tmp(l1, m))
        allocate(ut_arr(ns, l1))
        allocate(ut_arr_tmp(ns, l1))
    
        arr_tmp(:, :) = real(arr(:, :), kind(0d0))
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%uhat_f%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    
        arr_tmp(:, :) = aimag(arr(:, :))
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr_tmp(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%uhat_f%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
    
        ut_arr(:, :) = ut_arr(:, :) + ut_arr_tmp(:, :)
        do j = 1, ns
            do i = 1, l1
                ut_arr(j, i) = ut_arr(j, i) * obj%uhat_f%inv_s(j)
            end do
        end do
    
        ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
        res(:, :) = zero
        call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%uhat_f%v_real, n, zero, res, l2)
    
        deallocate(ut_arr, ut_arr_tmp, arr_tmp)
    end subroutine

    subroutine fit_matsubara_b_zz(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), allocatable :: ut_arr(:, :)
        double precision, allocatable :: arr_tmp(:, :)
        double precision, allocatable :: res_tmp(:, :)
        double precision, allocatable :: ut_arr_r(:, :)
        double precision, allocatable :: ut_arr_tmp(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%uhat_b%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of output array.'
        
        IF (.not. obj%positive_only) then
            allocate(ut_arr(ns, l1))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr(:, :) = czero
            call zgemm("n", "t", ns, l1, m, cone, obj%uhat_b%ut, ns, arr, l1, czero, ut_arr, ns)
            do j = 1, ns
                do i = 1, l1
                    ut_arr(j, i) = ut_arr(j, i) * obj%uhat_b%inv_s(j)
                end do
            end do
    
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res(:, :) = czero
            call zgemm("t", "t", l1, n, ns, cone, ut_arr, ns, obj%uhat_b%v, n, czero, res, l2)
            deallocate(ut_arr)
        ELSE
            allocate(arr_tmp(l1, m))
            allocate(res_tmp(l2, n))
            allocate(ut_arr_r(ns, l1))
            allocate(ut_arr_tmp(ns, l1))
            arr_tmp(:, :) = real(arr(:, :), kind(0d0))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr_r(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%uhat_b%ut_real, ns, arr_tmp, l1, zero, ut_arr_r, ns)
        
            arr_tmp(:, :) = aimag(arr(:, :))
            !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
            ut_arr_tmp(:, :) = zero
            call dgemm("n", "t", ns, l1, m, one, obj%uhat_b%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
        
            ut_arr_r(:, :) = ut_arr_r(:, :) + ut_arr_tmp(:, :)
            do j = 1, ns
                do i = 1, l1
                    ut_arr_r(j, i) = ut_arr_r(j, i) * obj%uhat_b%inv_s(j)
                end do
            end do
        
            ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
            res_tmp(:, :) = zero
            call dgemm("t", "t", l1, n, ns, one, ut_arr_r, ns, obj%uhat_b%v_real, n, zero, res_tmp, l2)

            res = cmplx(res_tmp, zero, kind(0d0))
        
            deallocate(ut_arr_r, ut_arr_tmp, arr_tmp, res_tmp)
        ENDIF
    
    end subroutine

    subroutine fit_matsubara_b_zd(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        double precision, intent(out) :: res(:, :)
        double precision, allocatable :: ut_arr(:, :)
        double precision, allocatable :: ut_arr_tmp(:, :)
        double precision, allocatable :: arr_tmp(:, :)
    
        integer :: m, n, l1, l2, ns, i, j
    
        ! ut(ns, m)
        ! v(n, ns)
        ! arr(l1, m)
        ! mat(m, n)
        ! ut_arr(ns, l1)
        ! res(l2, n)
        l1 = size(arr, 1)
        m = size(arr, 2)
        l2 = size(res, 1)
        n = size(res, 2)
        ns = obj%uhat_b%ns
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of input array.'
        IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of output array.'
        IF (.not. obj%positive_only) stop 'output array should be a complex array.'
        allocate(arr_tmp(l1, m))
        allocate(ut_arr(ns, l1))
        allocate(ut_arr_tmp(ns, l1))
    
        arr_tmp(:, :) = real(arr(:, :), kind(0d0))
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%uhat_b%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    
        arr_tmp(:, :) = aimag(arr(:, :))
        !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
        ut_arr_tmp(:, :) = zero
        call dgemm("n", "t", ns, l1, m, one, obj%uhat_b%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
    
        ut_arr(:, :) = ut_arr(:, :) + ut_arr_tmp(:, :)
        do j = 1, ns
            do i = 1, l1
                ut_arr(j, i) = ut_arr(j, i) * obj%uhat_b%inv_s(j)
            end do
        end do
    
        ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
        res(:, :) = zero
        call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%uhat_b%v_real, n, zero, res, l2)
    
        deallocate(ut_arr, ut_arr_tmp, arr_tmp)
    end subroutine

end module