module tensor
use types
use outerprod
use utils, only: unit_matrix
use randomm, only: gasdev

implicit none
private
public ih

contains

    subroutine ih(x, y, z, k, mat_a, R)
        !declaracion de variables
        integer, parameter :: n = 3 !la dimension de la matriz (x,y,z)
        integer, intent(in) :: k   ! numero de submatrices (particulas)
        real(dp), intent(in) :: x(:), y(:), z(:)

        ! Variables locales
        real(dp) :: datos(n, k)
        real(dp) :: Dij(n, n)
        real(dp) :: part1(n)
        real(dp) :: part2(n)
        real(dp) :: ident(n, n)

        real(dp) :: Xr(n*k)
        real(dp), intent(out) :: R(n*k)

        integer :: ix, ij, pos, il   ! indices
        integer :: info

        real(dp), intent(out) :: mat_a(n*k, n*k)
        real(dp) :: sigma(n*k, n*k)

        ! Concatenarlos para construir la matriz
        datos(1, :) = x
        datos(2, :) = y
        datos(3, :) = z

        call unit_matrix( ident )
        R = 0.0_dp
        Xr = 0.0_dp
        mat_a = 0.0_dp
        Dij = 0.0_dp

        ! Distribucion normal estándar
        do ix = 1, k*n
            Xr(ix) = gasdev()
        end do

        do ix = 1, k*n, n
            mat_a(ix:ix+n-1,ix:ix+n-1) = ident
        end do

        do ix = 1, k
            il = (ix * n) - n + 1
            do ij = ix+1, k
                pos = (ij*n) - n + 1

                part1 = datos(:, ix)
                part2 = datos(:, ij)

                call matrizDij( part1,part2,Dij,n )

                mat_a(il:il+n-1,pos:pos+n-1) = Dij
                mat_a(pos:pos+n-1,il:il+n-1) = Dij
            end do
        end do

        sigma = mat_a + 0.00001_dp
        ! Descomposición de Cholesky
        call dpotrf( 'L',n*k,sigma,n*k,info )
        ! Multiplicar L*Xr para obtener el vector de números aleatorios
        ! R = matmul( sigma, Xr )
        call dgemv( 'n',n*k,n*k,1.0_dp,sigma,n*k,Xr,1,0.0_dp,R,1 )
    end subroutine ih


    subroutine matrizDij(part1, part2, Dij, n)
        ! Variables de entrada y salida
        integer, intent(in) :: n
        real(dp), intent(in) :: part1(:)
        real(dp), intent(in) :: part2(:)
        real(dp), intent(out) :: Dij(:,:)

        ! Variables locales
        real(dp) :: rij,sqrddist,xij,yij,zij
        real(dp) :: ident(n, n), prodout(n, n)
        real(dp) :: temp(n)

        ! Calcular las distancias
        xij = part1(1) - part2(1)
        yij = part1(2) - part2(2)
        zij = part1(3) - part2(3)
        rij = xij**2 + yij**2 + zij**2
        temp = [xij, yij, zij]

        ! calculando Dij
        Dij = 0.0_dp
        prodout = 0.0_dp

        sqrddist = sqrt( rij )
        call unit_matrix( ident )

        if (sqrddist >= 1.0_dp) then
            prodout = outerprod_d( temp, temp ) / rij
            Dij = ident + prodout
            
            Dij = Dij + ( (ident/3.0_dp) - prodout )/(2.0_dp*rij)
            
            Dij = 3.0_dp*Dij/(8.0_dp*sqrddist)
        else
            Dij = (1.0_dp-(9.0_dp*sqrddist/16.0_dp))*ident
            Dij = Dij + 3.0_dp*outerprod_d( temp, temp )/(16.0_dp*sqrddist)
        end if
    end subroutine matrizDij
end module tensor