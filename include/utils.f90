module utils
    use types
    use box_pot_parameters
    implicit none
    private
    public unit_matrix, cholesky, check_nan, save_timeseries, &
            iniconfig, show_m, save_msd

contains

subroutine iniconfig(xc, yc, zc, d)
! defining three vector of mp dimension, it indicate that only are out variables
    real(dp), intent(out) :: xc(:), yc(:), zc(:)
    real(dp), intent(in) :: d
    ! Local variables
    integer :: i

    xc(1) = -(boxl-d)/2.0_dp
    yc(1) = -(boxl-d)/2.0_dp
    zc(1) = -(boxl-d)/2.0_dp

    do i = 2,np
        xc(i) = xc(i-1) + d
        yc(i) = yc(i-1)
        zc(i) = zc(i-1)
        if (xc(i) > rc) then
            xc(i) = xc(1)
            yc(i) = yc(i-1) + d
            if (yc(i) > rc) then
                xc(i) = xc(1)
                yc(i) = yc(1)
                zc(i) = zc(i-1) + d
            end if
        end if
    end do
end subroutine iniconfig

subroutine unit_matrix(mat) ! matrix dimension
    real(dp), intent(inout) :: mat(:, :)
    ! Local:
    integer :: i, n

    n = min(size(mat, 1), size(mat, 2))

    mat = 0.0_dp

    forall ( i = 1:n ) mat(i, i) = 1.0_dp

end subroutine unit_matrix

subroutine cholesky(mat_a, sigma)
    real(dp), intent(in) :: mat_a(:,:)
    real(dp), intent(out):: sigma(:,:)

    !Local variable
    integer :: i, j, m

    m = size(mat_a, 1)

    sigma = 0.0_dp

    do j = 1, m
        sigma(j, j) = dsqrt( mat_a(j, j) - dot_product(sigma(j,1:j-1),sigma(j,1:j-1)) )
        if ( sigma(j, j) <= 0.000001_dp ) then
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
subroutine show_m(A)
    real(dp), intent(in) :: A(:,:)
    integer :: ix, n

    n = size(A, 1)

    do ix = 1, n
        write(*,*) A(ix, :)
    end do
end subroutine show_m

subroutine check_unity(A)
    real(dp), intent(inout) :: A(:,:)
    integer :: i,n

    n = size(A,1)

    do i = 1, n
        if ( (A(i,i)-1.0_dp) <= 0.000001 ) then
            print*, 'Not unity!'
            stop
        end if
    end do

    print*, 'This matrix is unity!'
end subroutine check_unity

subroutine matrix_file(A)
    real(dp), intent(inout) :: A(:,:)
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
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    real(dp), intent(in) :: x(:)
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
    real(dp), intent(in) :: x(:,:), y(:,:), z(:,:)
    character(len=*), intent(in) :: filename
    character(len=1024) :: newname
    character(len=8) :: fmt
    character(len=8) :: x1
    integer :: i, j, n, m, u

    n = size(x,1)
    m = size(x,2)
    fmt = '(I5.1)'

    ! Ciclar sobre todas las partículas
    do i = 1, m
        write(x1,fmt) i
        newname = filename//trim(adjustl(x1))//'.dat'
        open(newunit=u, file=newname, status="new")
        do j = 1, n
            write(u,'(f11.8,A,f12.8,A,f12.8)') x(j, i),',',y(j, i),',',z(j, i)
        end do
        close(u)
    end do

end subroutine save_timeseries

subroutine save_msd(t, wt, ft, nprom, filename)
    ! Variables de entrada/salida
    real(dp), intent(in) :: t(:), wt(:), ft(:)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nprom

    ! Variables locales
    integer :: u, i

    open(newunit=u, file=filename, status='unknown')
    do i=1,nprom-1
        write(u,'(f15.11,A,f15.11,A,f14.11)') t(i+1),',',wt(i),',',ft(i)
    end do
    close(u)
end subroutine save_msd

end module utils