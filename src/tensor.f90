subroutine IH( x, y, z, k, mat_a, R )
    use argumt
    use outerprod_m
    use unit_matrix_m
    use box_pot_parameters

    implicit none

    !declaracion de variables
    integer(KIND = wInt), parameter :: n = 3 !la dimension de la matriz (x,y,z)
    integer(KIND = wInt), intent(in) :: k   ! numero de submatrices (particulas)

    !real(KIND = wReal), dimension(n, n) :: ident0
    real(KIND = wReal), intent(in), dimension(k) :: x !
    real(KIND = wReal), intent(in), dimension(k) :: y !
    real(KIND = wReal), intent(in), dimension(k) :: z !
    real(KIND = wReal), dimension(n, k) :: datos ! (3, particulas)

    ! i dif j
    real(KIND = wReal), DIMENSION(n, n) :: Dij
    real(KIND = wReal), dimension(n) :: part1 !
    real(KIND = wReal), dimension(n) :: part2 !
    real(KIND = wReal), dimension(n, n) :: ident

    real(KIND = wReal), dimension(n*k) :: Xr !
    real(KIND = wReal), dimension(n*k), intent(out) :: R !

    integer(KIND = wInt) :: ix, ij, pos, il   ! indices
    integer(kind=wInt) :: info

    real(KIND = wReal), dimension(n*k, n*k), intent(out) :: mat_a
    real(KIND = wReal), dimension(n*k, n*k) :: sigma

    real(KIND = wReal), external :: gasdev ! function

    ! Concatenarlos para construir la matriz
    datos(1, :) = x
    datos(2, :) = y
    datos(3, :) = z

    call unit_matrix(ident)
    R = 0.0_wReal
    Xr = 0.0_wReal
    mat_a = 0.0_wReal
    Dij = 0.0_wReal

    ! Distribucion normal estándar
    do ix = 1, k*n
        Xr(ix) = gasdev()
    end do

    do ix = 1, k*n, n
        mat_a(ix:ix+n-1,ix:ix+n-1) = ident
    end do

    do ix = 1, k
        il = (ix * n) - n + 1
        do ij = (ix+1), k
            pos = (ij*n) - n + 1

            part1 = datos(:, ix)
            part2 = datos(:, ij)

            call matrizDij(part1, part2, Dij, n)
            ! call dpotrf( 'L',n,Dij,n,info )
            ! if ( info .ne. 0 ) then
            !     print*, 'Decomposition not possible'
            !     print*,''
            !     call show_m( Dij, n, n )
            !     ! call matrix_file( mat_a, n*k )
            !     print*, part1
            !     print*,''
            !     print*, part2
            !     print*,''
            !     print*, ix
            !     print*,''
            !     print*, ij
            !     ! print*, info
            !     stop
            ! end if

            mat_a(il:il+n-1,pos:pos+n-1) = Dij
            mat_a(pos:pos+n-1,il:il+n-1) = Dij
            
        end do
    end do

    mat_a = mat_a + 0.00005_wReal
    call cholesky(n*k, mat_a, sigma)
    ! sigma = mat_a + 0.00005_wReal
    ! call dpotrf( 'L',n*k,sigma,n*k,info )
    ! if ( info .ne. 0 ) then
    !     print*, 'Decomposition not possible'
    !     print*,''
    !     ! call show_m( sigma, n*k, n*k )
    !     call matrix_file( mat_a, n*k )
    !     print*, info
    !     stop
    ! end if

    
    
    R = matmul(sigma, Xr)
end subroutine IH


subroutine matrizDij(part1, part2, Dij, n)
    use argumt
    use outerprod_m
    use box_pot_parameters
    use unit_matrix_m

    implicit none

    ! Variables de entrada y salida
    integer(KIND = wInt), intent(in) :: n
    real(KIND = wReal), dimension(n), intent(in) :: part1 !
    real(KIND = wReal), dimension(n), intent(in) :: part2 !
    real(KIND = wReal), dimension(n, n), intent(out) :: Dij

    ! Variables locales
    real(KIND = wReal) :: rij,sqrddist,xij,yij,zij
    ! integer(KIND = wInt) :: ix
    real(KIND = wReal), dimension(n, n) :: ident, prodout
    real(KIND = wReal), dimension(n) :: temp

    !! Calcular las distancias
    xij = part1(1) - part2(1)
    yij = part1(2) - part2(2)
    zij = part1(3) - part2(3)
    ! condiciones periódicas a la frontera
    ! xij = xij - boxl*dnint(xij/boxl)
    ! yij = yij - boxl*dnint(yij/boxl)
    ! zij = zij - boxl*dnint(zij/boxl)
    rij = xij**2 + yij**2 + zij**2
    temp = (/ xij, yij, zij /)
    ! print*,''
    ! print*,temp

    ! calculando Dij
    Dij = 0.0_wReal
    prodout = 0.0_wReal

    sqrddist = sqrt(rij)
    call unit_matrix(ident)

    ! prodout = outerprod_d( temp, temp ) / rij
    ! print*, prodout
    ! Dij = ident + prodout
    ! Dij = Dij + ( (ident/3.0_wReal) - prodout )/(2.0_wReal*rij)
    ! Dij = 3.0_wReal*Dij/(8.0_wReal*sqrddist)

    if (sqrddist .ge. 1.0_wReal) then
        prodout = outerprod_d( temp, temp ) / rij
        Dij = ident + prodout
        
        Dij = Dij + ( (ident/3.0_wReal) - prodout )/(2.0_wReal*rij)
        
        Dij = 3.0_wReal*Dij/(8.0_wReal*sqrddist)
    else
        Dij = (1.0_wReal-(9.0_wReal*sqrddist/16.0_wReal))*ident
        Dij = Dij + 3.0_wReal*outerprod_d( temp, temp )/(16.0_wReal*sqrddist)
    end if
end subroutine matrizDij