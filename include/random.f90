module randomm
use types
use box_pot_parameters, only: pi

implicit none
private
public gasdev, normal

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
        fac = dsqrt(-2.0_dp*dlog(rsq)/rsq)
        gasdev = v1 * fac
        gset = v2 * fac
        available = .true.
    end if
end function gasdev

subroutine normal(a)
real(dp), intent(inout) :: a(:)
integer :: i, n
real(dp) :: temp, mean = 0.0_dp, sd = 1.0_dp

n = size(a)

call random_number(a)

! Now convert to normal distribution
do i = 1, n-1, 2
    temp = sd * dsqrt(-2.0_dp*dlog(a(i))) * dcos(2.0_dp*pi*a(i+1)) + mean
    a(i+1) = sd * dsqrt(-2.0_dp*dlog(a(i))) * dsin(2.0_dp*pi*a(i+1)) + mean
    a(i) = temp
end do
end subroutine normal

end module randomm