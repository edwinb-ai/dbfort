module energies
    use types
    use box_pot_parameters

    implicit none
    
    save
    public force

contains

    subroutine  force(x, y, z, fx, fy, fz, ener)
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(out) :: fx(:), fy(:), fz(:)
        real(dp), intent(out) :: ener

        ! Local variables
        integer :: i, j
        real(dp) :: rij, xij, yij, zij, uij
        real(dp) :: fij, fxij, fyij, fzij

        ! Inicializar arreglos y variables
        ener = 0.0_dp
        fx(:) = 0.0_dp
        fy(:) = 0.0_dp
        fz(:) = 0.0_dp

        do i = 1, np-1
            do j = i+1, np
                uij = 0._dp
                fij = 0._dp
                fxij = 0._dp
                fyij = 0._dp
                fzij = 0._dp
                
                xij = x(i)-x(j)
                yij = y(i)-y(j)
                zij = z(i)-z(j)
                
                xij = xij-boxl*idnint(xij/boxl)
                yij = yij-boxl*idnint(yij/boxl)
                zij = zij-boxl*idnint(zij/boxl)
                
                rij = norm2( [xij, yij, zij] )
                
                if (rij < rc) then
                    call hardsphere( rij,uij,xij,yij,zij,fxij,fyij,fzij )
                    
                    ener = ener + uij
                    fx(i) = fx(i) + fxij
                    fy(i) = fy(i) + fyij
                    fz(i) = fz(i) + fzij
                    
                    fx(j) = fx(j) - fxij
                    fy(j) = fy(j) - fyij
                    fz(j) = fz(j) - fzij
                end if
            end do
        end do
    end subroutine force

    subroutine hardsphere(rij, uij, xij, yij, zij, fxij, fyij, fzij)
        real(dp), intent(in) :: rij, xij, yij, zij
        real(dp), intent(inout) :: fxij, fyij, fzij, uij

        ! Local variables
        real(dp) ::  fij

        if (rij < bpot) then
            uij = (a2/dT)*((1.0_dp/rij)**dlr-(1.0_dp/rij)**dla)
            uij = uij + 1.0_dp/dT
            fij = dlr*(1.0_dp/rij)**(dlr+1.0_dp)-dla*(1.0_dp/rij)**(dla+1.0_dp)
            fij = a2*fij/dT
        else
            uij = 0.0_dp
            fij = 0.0_dp
        end if

        fxij = fij*xij/rij
        fyij = fij*yij/rij
        fzij = fij*zij/rij
    end subroutine hardsphere
end module energies