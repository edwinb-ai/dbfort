module randomm
use types

implicit none
private
public gasdev

contains

    real(dp) function gasdev()
        real(dp) :: v1, v2, fac, rsq
        real(dp), save :: gset

        logical, save :: available = .false.

        if (available) then
            gasdev = gset
            available = .false.
        else
            do
                call random_number(v1)
                call random_number(v2)
                v1 = 2.0_dp*v1-1.0_dp
                v2 = 2.0_dp*v2-1.0_dp
                rsq = v1**2+v2**2
                if (rsq > 0.0_dp .and. rsq < 1.0_dp) exit
            end do
            fac = dsqrt(-2.0_dp*log(rsq)/rsq)
            gasdev = v1 * fac
            gset = v2 * fac
            available = .true.
        end if
    end function gasdev

end module randomm