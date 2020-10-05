module randomm
use iso_fortran_env, only: real64

implicit none
private
public :: gasdev

contains

    real(real64) function gasdev()
        real(real64) :: v1, v2, fac, rsq
        real(real64), save :: gset

        logical, save :: available = .false.

        if (available) then
            gasdev = gset
            available = .false.
        else
            do
                call random_number(v1)
                call random_number(v2)
                v1 = 2.0d0*v1-1.0d0
                v2 = 2.0d0*v2-1.0d0
                rsq = v1**2+v2**2
                if (rsq > 0.0d0 .and. rsq < 1.0d0) exit
            end do
            fac = sqrt(-2.0d0*log(rsq)/rsq)
            gasdev = v1 * fac
            gset = v2 * fac
            available = .true.
        end if
    end function gasdev

end module randomm