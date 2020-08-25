module outerprod
    use iso_fortran_env, only: real64

    implicit none
    private
    public outerprod_d

contains

function outerprod_d(a,b)
    real(real64), DIMENSION(:), INTENT(IN) :: a,b
    real(real64), DIMENSION(size(a), size(b)) :: outerprod_d

    outerprod_d = spread(a, dim=2, ncopies=size(b)) * &
                    spread(b, dim=1, ncopies=size(a))

end function outerprod_d

end module outerprod