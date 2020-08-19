subroutine cholesky(m, mat_a, sigma)
    use argumt
    implicit none

    !Global variables
    integer(KIND = wInt), intent(in) :: m
    real(KIND = wReal), intent(in), dimension(m, m) :: mat_a
    real(KIND = wReal), intent(out), dimension(m, m):: sigma

    !Local variable
    integer(KIND = wInt) :: i, j

    sigma = 0.0_wReal

    do j = 1, m
        sigma(j, j) = sqrt( mat_a(j, j) - dot_product(sigma(j,1:j-1),sigma(j,1:j-1)) )
        if ( sigma(j, j) .le. 0.000001_wReal ) then
            print*, 'Numerical instability'
            print*, sigma(j, j)
            stop
        end if
        do i = j+1, m
            sigma(i, j)  = (mat_a(i, j) - dot_product(sigma(i, 1: j-1), &
                        sigma(j, 1: j-1)) ) / sigma(j, j)
       end do
    end do

end subroutine cholesky

!! gasdev : Returns a normally distributed deviate with zero mean and unit variance from Numerical recipes
double precision function gasdev()
    use argumt

    ! Global variables
    use box_pot_parameters

! Global variables
    implicit none

    real(KIND = wReal):: v1, v2, fac, rsq
    real(KIND = wReal), save :: gset

    logical, save :: available = .false.

    if (available) then
        gasdev = gset
        available = .false.
    else
        do
            call random_number(v1)
            call random_number(v2)
            v1 = 2.0_wReal*v1-1.0_wReal
            v2 = 2.0_wReal*v2-1.0_wReal
            rsq = v1**2+v2**2
            if (rsq > 0.0_wReal .and. rsq < 1.0_wReal) exit
        end do
        fac = sqrt(-2.0_wReal*log(rsq)/rsq)
        gasdev = v1 * fac
        gset = v2 * fac
        available = .true.
    end if
return
end function gasdev

subroutine show_m(A, n, m)
    !! Subrutina simple para imprimir un arreglo de nxm
    !! como matriz en la terminal
    use iso_fortran_env

    implicit none

    integer, intent(in) :: n, m
    real(real64), intent(in), dimension(n, m) :: A
    integer :: ix

    do ix = 1, n
        write(*,*) A(ix, :)
    end do
end subroutine show_m

subroutine check_unity(A, n)
    use iso_fortran_env
    implicit none

    real(real64), dimension(n,n), intent(inout) :: A
    integer, intent(in) :: n
    integer :: i

    do i = 1, n
        if ( (A(i,i)-1.0d0) .le. 0.000001 ) then
            print*, 'Not unity!'
            stop
        end if
    end do

    print*, 'This matrix is unity!'
end subroutine check_unity

subroutine matrix_file(A, n)
    use iso_fortran_env
    implicit none

    real(real64), dimension(n,n), intent(inout) :: A
    integer, intent(in) :: n
    integer :: i,j

    open(1, file = 'matrix.dat', status='unknown')  
    do i = 1, n
        do j = 1, n
            write(1,'(20G12.4)',advance='no') A(i, j)
        end do
        write(1,*) ''
    end do  
    close(1)
end subroutine matrix_file

subroutine check_nan(x, n)
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    integer, intent(in) :: n
    real(real64), dimension(n), intent(in) :: x
    integer :: i

    do i = 1,n
        if ( ieee_is_nan(x(i)) ) then
            print*, 'Found NaN. Stopping now.'
            stop
        end if
    end do
end subroutine check_nan