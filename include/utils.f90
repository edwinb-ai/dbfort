module utils
implicit none
private
public unit_matrix, cholesky, check_nan,save_timeseries

contains

subroutine unit_matrix(mat) ! matrix dimension

    ! Sets matrix to be a unit matrix.
    use iso_fortran_env, only: real64, int32

    real(real64), DIMENSION(:, :), INTENT(OUT) :: mat
    ! Local:
    integer(int32) :: i, n

    n = min(size(mat, 1), size(mat, 2))

    mat(:, :) = 0.0d0

    do i=1, n
        mat(i, i) = 1.0d0
    end do

end subroutine unit_matrix

subroutine cholesky(m, mat_a, sigma)
    use iso_fortran_env, only: real64, int32
    implicit none

    !Global variables
    integer(int32), intent(in) :: m
    real(real64), intent(in), dimension(m, m) :: mat_a
    real(real64), intent(out), dimension(m, m):: sigma

    !Local variable
    integer(int32) :: i, j

    sigma = 0.0d0

    do j = 1, m
        sigma(j, j) = sqrt( mat_a(j, j) - dot_product(sigma(j,1:j-1),sigma(j,1:j-1)) )
        if ( sigma(j, j) <= 0.000001d0 ) then
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
subroutine show_m(A, n, m)
    !! Subrutina simple para imprimir un arreglo de nxm
    !! como matriz en la terminal
    use iso_fortran_env, only: real64

    implicit none

    integer, intent(in) :: n, m
    real(real64), intent(in), dimension(n, m) :: A
    integer :: ix

    do ix = 1, n
        write(*,*) A(ix, :)
    end do
end subroutine show_m

subroutine check_unity(A)
    use iso_fortran_env, only: real64
    implicit none

    real(real64), intent(inout) :: A(:,:)
    integer :: i,n

    n = size(A,1)

    do i = 1, n
        if ( (A(i,i)-1.0d0) <= 0.000001 ) then
            print*, 'Not unity!'
            stop
        end if
    end do

    print*, 'This matrix is unity!'
end subroutine check_unity

subroutine matrix_file(A)
    use iso_fortran_env, only: real64
    implicit none

    real(real64), intent(inout) :: A(:,:)
    integer :: i,j,n

    n = size(A,1)

    open(1, file = 'matrix.dat', status='unknown')  
    do i = 1, n
        do j = 1, n
            write(1,'(20G12.4)',advance='no') A(i, j)
        end do
        write(1,*) ''
    end do  
    close(1)
end subroutine matrix_file

subroutine check_nan(x)
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    real(kind=real64), intent(in) :: x(:)
    integer :: i,n

    n = size(x)

    do i = 1,n
        if ( ieee_is_nan(x(i)) ) then
            print*, 'Found NaN. Stopping now.'
            stop
        end if
    end do
end subroutine check_nan

subroutine save_timeseries(filename,x,y,z)
    use iso_fortran_env, only: real64, int32

    real(kind=real64), intent(in) :: x(:,:),y(:,:),z(:,:)
    character(len=*), intent(in) :: filename
    character(len=1024) :: newname
    character(len=8) :: fmt
    character(len=8) :: x1
    integer(kind=int32) :: i, j, n, m

    n = size(x,1)
    m = size(x,2)
    fmt = '(I5.1)'

    ! Ciclar sobre todas las partículas
    do i = 1, m
        write(x1,fmt) i
        newname = filename//trim(adjustl(x1))//'.dat'
        open(60, file=newname, status="new")
        do j = 1, n
            write(60,'(3f16.8)') x(j, i),y(j, i),z(j, i)
        end do
        close(60)
    end do

end subroutine save_timeseries

end module utils