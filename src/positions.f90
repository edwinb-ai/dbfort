subroutine position_ih(x, y, z, fx, fy, fz, dij, Rz, pbc)
    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl, dT

    implicit none

    real(real64), intent(in), dimension(np) :: fx, fy, fz
    real(real64), intent(inout), dimension(np) :: x, y, z
    real(real64), intent(in):: pbc
    real(real64), intent(in), dimension(3*np, 3*np) :: dij
    real(real64), intent(in), dimension(3*np) :: Rz
    ! Local variables
    integer(int32) :: i, j, ij
    integer(int32), parameter :: k = 3 ! Dimensión espacial del sistema
    real(real64), dimension(k, np) :: fuerzas
    real(real64), dimension(np*k) :: temp

    fuerzas(1, :) = fx
    fuerzas(2, :) = fy
    fuerzas(3, :) = fz

    temp = 0.0d0

    do i = 1, k*np, k
        do j = 1, np-1
            ij = (k*j) - k + 1
            temp(i:i+2) = temp(i:i+2) + matmul( dij(i:i+2,ij:ij+2), fuerzas(:, j) )
        enddo
    enddo

    do i = 1, np
        ij = (k * i) - k + 1
        x(i) = x(i) + temp(ij)*deltat + sqrt(2.0d0*deltat)*Rz(ij)
        y(i) = y(i) + temp(ij+1)*deltat + sqrt(2.0d0*deltat)*Rz(ij+1)
        z(i) = z(i) + temp(ij+2)*deltat + sqrt(2.0d0*deltat)*Rz(ij+2)
        if (pbc > 0.0d0) then
            x(i) = x(i) - boxl*anint(x(i)/boxl)
            y(i) = y(i) - boxl*anint(y(i)/boxl)
            z(i) = z(i) - boxl*anint(z(i)/boxl)
        endif
    enddo
return
end subroutine position_ih

subroutine position(x, y, z, fx, fy, fz, pbc)
    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl, dT
    use randomm, only: gasdev
    implicit none

    real(real64), intent(in), dimension(np) :: fx, fy, fz
    real(real64), intent(inout), dimension(np) :: x, y, z
    real(real64), intent(in):: pbc
    ! Local variables
    integer(int32) :: i
    real(real64) :: sigma, dx, dy, dz
    ! real(real64), external :: gasdev ! function

    sigma = dsqrt(2.d0*deltat)
    do i = 1, np
        dx = sigma*gasdev()
        dy = sigma*gasdev()
        dz = sigma*gasdev()
        x(i) = x(i) + fx(i)*deltat + dx
        y(i) = y(i) + fy(i)*deltat + dy
        z(i) = z(i) + fz(i)*deltat + dz
        if (pbc > 0.d0) then
            x(i) = x(i) - boxl*dnint(x(i)/boxl)
            y(i) = y(i) - boxl*dnint(y(i)/boxl)
            z(i) = z(i) - boxl*dnint(z(i)/boxl)
        endif
    enddo
return
end subroutine position