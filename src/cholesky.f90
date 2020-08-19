subroutine cholesky(n, A, G)
    ! Argumentos:
    ! n (integer): la dimension de la matriz
    ! A (flotante): arreglo de nxn
    ! G (flotante): arreglo de nxn, la matriz diagonal
    !               inferior
    implicit none
    integer,    intent(in)    :: n
    real*8,     intent(in)    :: A(n,n)
    real*8,     intent(out)   :: G(n,n)
    integer    :: i,j
 
    G = 0.0d0

   do j = 1, n

      G(j,j) = dsqrt( A(j,j) - dot_product(G(j,1:j-1),G(j,1:j-1)) )
      
      do i = j+1, n
         G(i,j)  = ( A(i,j) - dot_product(G(i,1:j-1),G(j,1:j-1)) ) / G(j,j)
      end do

   end do
end subroutine cholesky