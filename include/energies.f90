module energies
    use types
    use box_pot_parameters

    implicit none
    
    private
    public force

contains

    subroutine  force(r,f,ener)
        real(dp), intent(in) :: r(:,:)
        real(dp), intent(inout) :: f(:,:)
        real(dp), intent(out) :: ener

        ! Local variables
        integer, parameter :: n = 3 ! Dimensión espacial
        integer :: i, j
        real(dp) :: rij, rposij(n), uij
        real(dp) :: fij(n), fxij, fyij, fzij
        real(dp) :: dnrm2 ! Norma euclidiana de BLAS

        ! Inicializar arreglos y variables
        ener = 0.0_dp
        f = 0.0_dp
        fij = 0.0_dp

        do i = 1, np-1
            do j = i+1, np
                uij = 0._dp
                fxij = 0._dp
                fyij = 0._dp
                fzij = 0._dp
                
                ! xij = x(i)-x(j)
                ! yij = y(i)-y(j)
                ! zij = z(i)-z(j)
                rposij = r(:,i) - r(:,j)
                
                ! xij = xij-boxl*idnint(xij/boxl)
                ! yij = yij-boxl*idnint(yij/boxl)
                ! zij = zij-boxl*idnint(zij/boxl)
                rposij = rposij - boxl*idnint(rposij/boxl)
                
                ! rij = norm2( [xij, yij, zij] )
                ! Mejor rendimiento con BLAS
                ! rij = dnrm2( 3,[xij, yij, zij],1 )
                rij = dnrm2( n,rposij,1 )
                
                if (rij < rc) then
                    call hardsphere( rij,uij,rposij,fxij,fyij,fzij )
                    fij = [fxij,fyij,fzij]
                    
                    ener = ener + uij
                    ! fx(i) = fx(i) + fxij
                    ! fy(i) = fy(i) + fyij
                    ! fz(i) = fz(i) + fzij
                    f(:,i) = f(:,i) + fij
                    
                    ! fx(j) = fx(j) - fxij
                    ! fy(j) = fy(j) - fyij
                    ! fz(j) = fz(j) - fzij
                    f(:,j) = f(:,j) - fij
                end if
            end do
        end do
    end subroutine force

    subroutine hardsphere(rij, uij, rposij, fxij, fyij, fzij)
        real(dp), intent(in) :: rij, rposij(:)
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

        fxij = fij*rposij(1)/rij
        fyij = fij*rposij(2)/rij
        fzij = fij*rposij(3)/rij
    end subroutine hardsphere
end module energies