module box_pot_parameters
    use iso_fortran_env, only: real64, int32
    implicit none
    save
! CONSTANT PARAMETERS
    ! parameters arguments
                                      !packing fraction
    real(real64), parameter :: phi = 0.05d0, &
                                     pi = 4.d0*datan(1.0d0), &
                                     rho = 6.d0*phi/pi, &
                                     diam = 1.0d0
! phi = filling fraction , diam = diameter of the particles
! pi = value with double presicion , rho = reduced density

    ! box arguments
    integer(int32), parameter :: np = 6**3! number of particles

    real(real64), parameter :: boxl = (np/rho)**(1.d0/3.d0), &
                                     rc = boxl/2.d0

  ! np = number of particule , boxl = box length, rc = radio corte

    ! pot arguments
    real(real64), parameter :: dlr = 50.d0, dT = 1.4737d0
    real(real64), parameter :: dla = 49.d0, deltat = 0.00001d0
    real(real64), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))

    ! mp and mr arguments
    integer(int32), parameter :: mp = 1024, mr = 2**9, &
                                        mt = 1000
end module box_pot_parameters