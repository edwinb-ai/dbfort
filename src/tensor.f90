module tensor
use iso_fortran_env, only: real64, int32
use outerprod
use utils, only: unit_matrix
use randomm, only: gasdev

implicit none
private
public ih

contains

    subroutine IH( x, y, z, k, mat_a, R )
        !declaracion de variables
        integer, parameter :: n = 3 !la dimension de la matriz (x,y,z)
        integer, intent(in) :: k   ! numero de submatrices (particulas)

        !real(real64)(KIND = wreal(real64)), dimension(n, n) :: ident0
        real(real64), intent(in), dimension(k) :: x !
        real(real64), intent(in), dimension(k) :: y !
        real(real64), intent(in), dimension(k) :: z !
        real(real64), dimension(n, k) :: datos ! (3, particulas)

        ! i dif j
        real(real64), DIMENSION(n, n) :: Dij
        real(real64), dimension(n) :: part1 !
        real(real64), dimension(n) :: part2 !
        real(real64), dimension(n, n) :: ident

        real(real64), dimension(n*k) :: Xr !
        real(real64), dimension(n*k), intent(out) :: R !

        integer :: ix, ij, pos, il   ! indices
        integer :: info

        real(real64), dimension(n*k, n*k), intent(out) :: mat_a
        real(real64), dimension(n*k, n*k) :: sigma

        ! real(real64), external :: gasdev ! function

        ! Concatenarlos para construir la matriz
        datos(1, :) = x
        datos(2, :) = y
        datos(3, :) = z

        call unit_matrix( ident )
        R = 0.0d0
        Xr = 0.0d0
        mat_a = 0.0d0
        Dij = 0.0d0

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

                call matrizDij( part1, part2, Dij, n )

                mat_a(il:il+n-1,pos:pos+n-1) = Dij
                mat_a(pos:pos+n-1,il:il+n-1) = Dij
                
            end do
        end do

        sigma = mat_a + 0.00005d0
        ! Descomposición de Cholesky
        call dpotrf( 'L',n*k,sigma,n*k,info )
        ! Multiplicar L*Xr para obtener el vector de números aleatorios
        ! R = matmul( sigma, Xr )
        call dgemv( 'n',n*k,n*k,1.0d0,sigma,n*k,Xr,1,0.0d0,R,1 )
    end subroutine IH


    subroutine matrizDij(part1, part2, Dij, n)
        ! Variables de entrada y salida
        integer(int32), intent(in) :: n
        real(real64), dimension(n), intent(in) :: part1 !
        real(real64), dimension(n), intent(in) :: part2 !
        real(real64), dimension(n, n), intent(out) :: Dij

        ! Variables locales
        real(real64) :: rij,sqrddist,xij,yij,zij
        real(real64), dimension(n, n) :: ident, prodout
        real(real64), dimension(n) :: temp

        !! Calcular las distancias
        xij = part1(1) - part2(1)
        yij = part1(2) - part2(2)
        zij = part1(3) - part2(3)
        rij = xij**2 + yij**2 + zij**2
        temp = [xij, yij, zij]

        ! calculando Dij
        Dij = 0.0d0
        prodout = 0.0d0

        sqrddist = sqrt( rij )
        call unit_matrix( ident )

        if (sqrddist .ge. 1.0d0) then
            prodout = outerprod_d( temp, temp ) / rij
            Dij = ident + prodout
            
            Dij = Dij + ( (ident/3.0d0) - prodout )/(2.0d0*rij)
            
            Dij = 3.0d0*Dij/(8.0d0*sqrddist)
        else
            Dij = (1.0d0-(9.0d0*sqrddist/16.0d0))*ident
            Dij = Dij + 3.0d0*outerprod_d( temp, temp )/(16.0d0*sqrddist)
        end if
    end subroutine matrizDij
end module tensor