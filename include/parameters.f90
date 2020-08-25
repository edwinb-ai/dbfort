module box_pot_parameters
    use iso_fortran_env, only: real64, int32
    implicit none
    save
! CONSTANT PARAMETERS
    ! pot arguments
    real(real64), parameter :: dlr = 50.d0, dT = 1.4737d0
    real(real64), parameter :: dla = 49.d0, deltat = 0.000001d0
    real(real64), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))

    ! parameters arguments
    real(real64), parameter :: phi = 0.50d0
    real(real64), parameter :: pi = 4.d0*datan(1.0d0)
    real(real64), parameter :: rho = 6.d0*phi/pi
    real(real64), parameter :: diam = 1.0d0
    real(real64), parameter :: sqtwodt = sqrt(2.0d0*deltat)

    ! box arguments
    integer(int32), parameter :: np = 6**3! number of particles

    real(real64), parameter :: boxl = (np/rho)**(1.d0/3.d0)
    real(real64), parameter :: rc = boxl/2.d0

    ! mp and mr arguments
    integer(int32), parameter :: mp = 1024, mr = 2**9
    integer(int32), parameter :: mt = 1000
end module box_pot_parameters