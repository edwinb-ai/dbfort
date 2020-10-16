module box_pot_parameters
    use types
    implicit none
    public
! CONSTANT PARAMETERS
    ! pot arguments
    real(dp), parameter :: dlr = 50._dp, dT = 1.4737_dp
    real(dp), parameter :: dla = 49._dp, deltat = 0.00001_dp
    real(dp), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
    real(dp), parameter :: bpot = (dlr/dla)**(1.0_dp/(dlr-dla))

    ! parameters arguments
    real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)
    
    real(dp), parameter :: diam = 1.0_dp
    real(dp), parameter :: sqtwodt = dsqrt(2.0_dp*deltat)

    ! Dependen del usuario
    real(dp) :: phi, rho, boxl, rc
    integer :: np
    logical :: with_ih ! Para prender o apagar las IH

    ! mp and mr arguments
    integer, parameter :: mp = 1024, mr = 2**9
    integer, parameter :: mt = 300000
end module box_pot_parameters