module sparse_ir_modules2
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

    contains

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
        vt(:, :), work(:)
    double precision, allocatable :: rwork(:), s(:)
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
        lwork = mn*mn + 3*mn
        allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork((5*mn+7)*mn), iwork(8*mn))
    else
        lwork = 2*mn + m + n
        allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork(5*n))
    endif
  
    a_copy(1:m, 1:n) = a(1:m, 1:n)
    if (.not. ill_conditioned) then
        lwork = mn*mn + 3*mn
        call dgesdd('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
        if (info /= 0) then
            stop 'Failure in DGESDD.'
        end if
    else
        lwork = 2*mn + m + n
        call dgesvd('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
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
    dmat%ut_real(1:ns, 1:m) = (transpose(u(1:m, 1:ns)))
    dmat%v_real(1:n, 1:ns) = (transpose(vt(1:ns, 1:n)))
    dmat%m = size(a, 1)
    dmat%n = size(a, 2)
    dmat%ns = ns

    dmat%a = dmat%a_real
    dmat%ut = dmat%ut_real
    dmat%v = dmat%v_real
  
    if (.not. ill_conditioned) then
        deallocate(work, a_copy, s, u, vt, rwork, iwork)
    else
        deallocate(work, a_copy, s, u, vt, rwork)
    endif
end function

function split_decompose(a, has_zero, eps,ill_conditioned) result(dmat)
    complex(kind(0d0)), intent(in) :: a(:, :)
    logical, intent(in) :: has_zero
    double precision, intent(in) :: eps
    logical, intent(in) :: ill_conditioned
  
    integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns, m_half
    double precision, allocatable :: a_copy(:, :), u(:, :), &
        vt(:, :), work(:)
    complex(kind(0d0)) :: u_copy(:, :)
    double precision, allocatable :: rwork(:), s(:)
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
        lwork = mn*mn + 3*mn
        allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork((5*mn+7)*mn), iwork(8*mn))
    else
        lwork = 2*mn + m + n
        allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork(5*n))
    end if

    a_copy(1:m_half, 1:n) = real(a(1:m_half, 1:n))

    if (has_zero) then
        a_copy(m_half+1:m, 1:n) = aimag(a(2:m, 1:n))
    else
        a_copy(m_half+1:m, 1:n) = aimag(a(1:m, 1:n))
    end if

    if (.not. ill_conditioned) then
        lwork = mn*mn + 3*mn
        call dgesdd('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
        if (info /= 0) then
            stop 'Failure in DGESDD.'
        end if
    else
        lwork = 2*mn + m + n
        call dgesvd('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
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
    allocate(dmat%v_imag(ns, m))
    allocate(dmat%a_odd(m, n/2))
    allocate(dmat%a_even(m, n/2))
  
    ! dmat%a temporarily stores the same data of input a
    dmat%a = a
    dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
    ! inv_s temporarily stores the same data of inv_s_dl
    dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
    dmat%ut(1:ns, 1:m) = (transpose(u_copy(1:m, 1:ns)))
    dmat%v_real(1:n, 1:ns) = (transpose(vt(1:ns, 1:n)))
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
        deallocate(work, a_copy, s, u, vt, rwork, iwork)
    else
        deallocate(work, a_copy, s, u, vt, rwork)
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
    ! res(l1, n)
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
        allocate(arr_tmp(l1, n))
        allocate(res_r(l2, m))
        allocate(res_i(l2, m))
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
        allocate(arr_tmp(l1, n))
        allocate(res_r(l2, m))
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
    ! res(l1, n)
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
    call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v, n, zero, res, l2)

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
    ! res(l1, n)
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
    ! res(l1, n)
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
    call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_tmp, l2)

    deallocate(ut_arr)
end subroutine

subroutine fit_matsubara_f_zz(obj, arr, res)
    type(IR), intent(in) :: obj
    complex(kind(0d0)), intent (in) :: arr(:, :)
    complex(kind(0d0)), intent(out) :: res(:, :)
    complex(kind(0d0)), allocatable :: ut_arr(:, :)

    integer :: m, n, l1, l2, ns, i, j

    ! ut(ns, m)
    ! v(n, ns)
    ! arr(l1, m)
    ! mat(m, n)
    ! ut_arr(ns, l1)
    ! res(l1, n)
    l1 = size(arr, 1)
    m = size(arr, 2)
    l2 = size(res, 1)
    n = size(res, 2)
    ns = obj%uhat_f%ns
    IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
    IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of input array.'
    IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of output array.'
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
    ! res(l1, n)
    l1 = size(arr, 1)
    m = size(arr, 2)
    l2 = size(res, 1)
    n = size(res, 2)
    ns = obj%uhat_f%ns
    IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
    IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of input array.'
    IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of output array.'
    IF (.not. obj%positive_only) stop 'input array should be a complex array.'
    allocate(arr_tmp(l1, m))
    allocate(ut_arr(ns, l1))
    allocate(ut_arr_tmp(ns, l1))

    arr_tmp(:, :) = real(arr(:, :), kind(0d0))
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    call dgemm("n", "t", ns, l1, m, one, obj%uhat_f%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    do j = 1, ns
        do i = 1, l1
            ut_arr(j, i) = ut_arr(j, i) * obj%uhat_f%inv_s(j)
        end do
    end do

    arr_tmp(:, :) = aimag(arr(:, :))
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr_tmp(:, :) = zero
    call dgemm("n", "t", ns, l1, m, one, obj%uhat_f%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
    do j = 1, ns
        do i = 1, l1
            ut_arr_tmp(j, i) = ut_arr_tmp(j, i) * obj%uhat_f%inv_s(j)
        end do
    end do

    ut_arr(:, :) = ut_arr(:, :) + ut_arr_tmp(:, :)

    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_tmp(:, :) = czero
    call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%uhat_f%v_real, n, zero, res, l2)

    deallocate(ut_arr, ut_arr_tmp, arr_tmp)
end subroutine

subroutine fit_matsubara_b_zz(obj, arr, res)
    type(IR), intent(in) :: obj
    complex(kind(0d0)), intent (in) :: arr(:, :)
    complex(kind(0d0)), intent(out) :: res(:, :)
    complex(kind(0d0)), allocatable :: ut_arr(:, :)

    integer :: m, n, l1, l2, ns, i, j

    ! ut(ns, m)
    ! v(n, ns)
    ! arr(l1, m)
    ! mat(m, n)
    ! ut_arr(ns, l1)
    ! res(l1, n)
    l1 = size(arr, 1)
    m = size(arr, 2)
    l2 = size(res, 1)
    n = size(res, 2)
    ns = obj%uhat_b%ns
    IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
    IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of input array.'
    IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of output array.'
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
    ! res(l1, n)
    l1 = size(arr, 1)
    m = size(arr, 2)
    l2 = size(res, 1)
    n = size(res, 2)
    ns = obj%uhat_b%ns
    IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
    IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of input array.'
    IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of output array.'
    IF (.not. obj%positive_only) stop 'input array should be a complex array.'
    allocate(arr_tmp(l1, m))
    allocate(ut_arr(ns, l1))
    allocate(ut_arr_tmp(ns, l1))

    arr_tmp(:, :) = real(arr(:, :), kind(0d0))
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    call dgemm("n", "t", ns, l1, m, one, obj%uhat_b%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    do j = 1, ns
        do i = 1, l1
            ut_arr(j, i) = ut_arr(j, i) * obj%uhat_b%inv_s(j)
        end do
    end do

    arr_tmp(:, :) = aimag(arr(:, :))
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr_tmp(:, :) = zero
    call dgemm("n", "t", ns, l1, m, one, obj%uhat_b%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
    do j = 1, ns
        do i = 1, l1
            ut_arr_tmp(j, i) = ut_arr_tmp(j, i) * obj%uhat_b%inv_s(j)
        end do
    end do

    ut_arr(:, :) = ut_arr(:, :) + ut_arr_tmp(:, :)

    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_tmp(:, :) = czero
    call dgemm("t", "t", l1, n, ns, one, ut_arr, ns, obj%uhat_b%v_real, n, zero, res, l2)

    deallocate(ut_arr, ut_arr_tmp, arr_tmp)
end subroutine

end module