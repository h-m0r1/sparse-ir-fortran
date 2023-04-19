module sparse_ir
    use sparse_ir_modules_orig
    implicit none

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
        complex(kind(0d0)), allocatable :: u_data(:,:)
        complex(kind(0d0)), allocatable :: uhat_f_data(:,:), uhat_b_data(:,:)
        complex(kind(0d0)), allocatable :: v_data(:,:), dlr_data(:,:)
        type(DecomposedMatrix_d) :: u
        type(DecomposedMatrix_z) :: uhat_f, uhat_b
        !type(DecomposedMatrix_d) :: dlr
        logical :: positive_only
    end type

    contains

    subroutine init_ir(obj, beta, lambda, eps, s, x, freq_f, freq_b, u, uhat_f, uhat_b, y, v, dlr, eps_svd, positive_only)
        type(IR), intent(inout) :: obj
        double precision, intent(in) :: beta, lambda, eps, s(:), x(:), y(:), eps_svd
        complex(kind(0d0)), intent(in) :: u(:,:), uhat_f(:, :), uhat_b(:, :), v(:, :), dlr(:, :)
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
            obj%nfreq_b = size(freq_b - 1) / 2 + 1
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
            obj%freq_f = freq_f((size(freq_f) / 2):size(freq_f))
            obj%freq_b = freq_b((size(freq_b - 1) / 2):size(freq_b))
        end if

        allocate(obj%u_data(obj%ntau, obj%size))
        obj%u_data = u

        allocate(obj%uhat_f_data(obj%nfreq_f, obj%size))

        allocate(obj%uhat_b_data(obj%nfreq_b, obj%size))

        if (.not. obj%positive_only) then
            obj%uhat_f_data = uhat_f
            obj%uhat_b_data = uhat_b
        else
            obj%uhat_f_data = uhat_f((size(freq_f) / 2):size(freq_f), :)
            obj%uhat_b_data = uhat_b((size(freq_b - 1) / 2):size(freq_b), :)
        end if

        allocate(obj%y(obj%nomega))
        obj%y = y

        allocate(obj%omega(obj%nomega))

        allocate(obj%v_data(obj%nomega, obj%size))
        obj%v_data = v

        allocate(obj%dlr_data(obj%size, obj%nomega))
        obj%dlr_data = transpose(dlr)

        obj%u = decompose(obj%u_data, obj%eps_svd)
        if (.not. obj%positive_only) then
            obj%uhat_f = decompose(obj%uhat_f_data, obj%eps_svd, .false.)
            obj%uhat_b = decompose(obj%uhat_b_data, obj%eps_svd, .false.)
        else
            obj%uhat_f = split_decompose(obj%uhat_f_data, .false., obj%eps_svd, .false.)
            obj%uhat_b = split_decompose(obj%uhat_b_data, .true., obj%eps_svd, .false.)
        end if
        obj%uhat_f%a_odd(:, :) = obj%uhat_f%a_imag(:, 1:(obj%uhat_f%n - 1):2)
        obj%uhat_f%a_even(:, :) = obj%uhat_f%a_real(:, 2:obj%uhat_f%n:2)
        obj%uhat_b%a_odd(:, :) = obj%uhat_b%a_real(:, 1:(obj%uhat_b%n - 1):2)
        obj%uhat_b%a_even(:, :) = obj%uhat_b%a_imag(:, 2:obj%uhat_b%n:2)

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

        obj%u%a(:, :) = sqrt(2.0d0/beta)*obj%u_data(:, :)
        obj%uhat_f%a(:, :) = sqrt(beta) * obj%uhat_f_data(:, :)
        obj%uhat_b%a(:, :) = sqrt(beta) * obj%uhat_b_data(:, :)
        !obj%dlr%a(:, :) = sqrt(5.0d-1*beta)*obj%dlr_data(:, :)

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
    
        call finalize_dmat(obj%u)
        call finalize_dmat(obj%uhat_f)
        call finalize_dmat(obj%uhat_b)
        !call finalize_dmat(obj%dlr)
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

end module