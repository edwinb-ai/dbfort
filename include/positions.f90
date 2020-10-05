module positions
    use iso_fortran_env, only: real64, int32
    use box_pot_parameters, only: np, sqtwodt, boxl, deltat

    implicit none
    private

    public :: position, position_ih

    contains

    subroutine position_ih(x, y, z, fx, fy, fz, dij, Rz, pbc)
        real(real64), intent(in), dimension(np) :: fx, fy, fz
        real(real64), intent(inout), dimension(np) :: x, y, z
        integer(int32), intent(in):: pbc
        real(real64), intent(in), dimension(3*np, 3*np) :: dij
        real(real64), intent(in), dimension(3*np) :: Rz
        ! Local variables
        integer(int32) :: i, j, ij
        integer(int32), parameter :: k = 3 ! Dimensión espacial del sistema
        real(real64), dimension(k, np) :: fuerzas
        real(real64), dimension(np*k) :: temp
        real(real64), dimension(k) :: mulout

        fuerzas(1, :) = fx
        fuerzas(2, :) = fy
        fuerzas(3, :) = fz

        temp = 0.0d0

        do i = 1, k*np, k
            do j = 1, np-1
                ij = (k*j) - k + 1
                call dgemv( 'n',k,k,1.0d0,dij(i:i+2,ij:ij+2),k,fuerzas(:, j),1,0.d0,mulout,1 )
                ! temp(i:i+2) = temp(i:i+2) + matmul( dij(i:i+2,ij:ij+2), fuerzas(:, j) )
                temp(i:i+2) = temp(i:i+2) + mulout
            end do
        end do

        do i = 1, np
            ij = (k * i) - k + 1
            x(i) = x(i) + temp(ij)*deltat + sqtwodt*Rz(ij)
            y(i) = y(i) + temp(ij+1)*deltat + sqtwodt*Rz(ij+1)
            z(i) = z(i) + temp(ij+2)*deltat + sqtwodt*Rz(ij+2)
            if (pbc == 1) then
                x(i) = x(i) - boxl*dnint(x(i)/boxl)
                y(i) = y(i) - boxl*dnint(y(i)/boxl)
                z(i) = z(i) - boxl*dnint(z(i)/boxl)
            end if
        end do
    end subroutine position_ih

    subroutine position(x, y, z, fx, fy, fz, pbc)
        use randomm, only: gasdev

        real(real64), intent(in), dimension(np) :: fx, fy, fz
        real(real64), intent(inout), dimension(np) :: x, y, z
        integer(int32), intent(in):: pbc
        ! Local variables
        integer(int32) :: i

        do i = 1, np
            x(i) = x(i) + (fx(i)*deltat) + (sqtwodt*gasdev())
            y(i) = y(i) + (fy(i)*deltat) + (sqtwodt*gasdev())
            z(i) = z(i) + (fz(i)*deltat) + (sqtwodt*gasdev())
            if (pbc == 1) then
                x(i) = x(i) - boxl*dnint(x(i)/boxl)
                y(i) = y(i) - boxl*dnint(y(i)/boxl)
                z(i) = z(i) - boxl*dnint(z(i)/boxl)
            end if
        end do
    end subroutine position

end module positions