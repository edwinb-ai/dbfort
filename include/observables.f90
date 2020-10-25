module observables
    use types
    use box_pot_parameters
    implicit none
    save
    public gr, difusion
contains
    subroutine gr(x,y,z,g,dr,pbc)
        real(dp), intent(in) ::  x(:), y(:), z(:)
        real(dp), intent(inout) :: g(:)
        real(dp), intent(in) :: dr
        logical, intent(in) :: pbc
        real(dp) :: dnrm2
        
        ! Local variables
        integer :: i, j, nbin
        real(dp) :: xij, yij, zij, rij

        do i = 1,np-1
            do j = i+1,np

                xij = x(i)-x(j)
                yij = y(i)-y(j)
                zij = z(i)-z(j)

                if (pbc) then
                    xij = xij-boxl*idnint(xij/boxl)
                    yij = yij-boxl*idnint(yij/boxl)
                    zij = zij-boxl*idnint(zij/boxl)
                end if

                ! rij = norm2( [xij,yij,zij] )
                ! Mejor rendimiento con LAPACK
                rij = dnrm2( 3,[xij, yij, zij],1 )
                if (rij < rc) then
                    nbin = idnint(rij/dr) + 1
                    if (nbin <= mr) then
                        g(nbin) = g(nbin) + 2.0_dp
                    end if
                end if
            end do
        end do
    end subroutine gr

    ! This sobroutine computes the mean-square displacement
    subroutine difusion(nprom, cfx, cfy, cfz, wt, ft)
        real(dp), intent(in) :: cfx(:,:), cfy(:,:), cfz(:,:)
        real(dp), intent(inout) :: wt(:), ft(:)
        integer, intent(in):: nprom
        ! Local variables
        integer :: i, j, k
        real(dp) :: dif2, dself, dx, dy, dz, dk, aux, aux2

        !Mean-square displacement and intermediate scattering function
        dk = 6.6_dp
        do i = 1, nprom-1
            dif2 = 0.0_dp
            dself = 0.0_dp
            do j = 1, nprom-i
                do k = 1, np
                    dx = cfx(j+i,k)-cfx(j,k)
                    dy = cfy(j+i,k)-cfy(j,k)
                    dz = cfz(j+i,k)-cfz(j,k)
                    aux = dx*dx + dy*dy + dz*dz
                    dif2 = dif2 + aux
                    aux2 = dk*dsqrt(aux)
                    aux = dsin(aux2)/aux2
                    dself = dself+aux
                end do
            end do
            aux2 = real(np*(nprom-i))
            dif2 = dif2/aux2
            dself = dself/aux2
            wt(i) = wt(i)+dif2
            ft(i) = ft(i)+dself
        end do
    end subroutine difusion
end module observables