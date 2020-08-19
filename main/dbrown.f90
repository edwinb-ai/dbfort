! codigo
! DEFINING THE MODULES
module argumt
   use, intrinsic :: iso_fortran_env, only: wInt => INT32, wReal => real64, &
                                                            xReal => real32
! las variables q se declaren con wp trabajar en doble precision las q no con simple precision
end module argumt

module outerprod_m
    use argumt

    implicit NONE

    interface outerprod
        ! Returns the outer product of two vectors.
        module procedure outerprod_r
        module procedure outerprod_d
    end interface

    real(KIND = wReal), parameter :: h = 3

    !private outerprod_r, outerprod_d
    public outerprod_r, outerprod_d

contains

    function outerprod_r(a,b)
        real(KIND = xReal), DIMENSION(:), INTENT(IN) :: a,b
        real(KIND = xReal), DIMENSION(size(a), size(b)) :: outerprod_r

        outerprod_r = spread(a, dim=2, ncopies=size(b)) * &
                      spread(b, dim=1, ncopies=size(a))

    end function outerprod_r

    function outerprod_d(a,b)
        real(KIND = wReal), DIMENSION(:), INTENT(IN) :: a,b
        real(KIND = wReal), DIMENSION(size(a), size(b)) :: outerprod_d

        outerprod_d = spread(a, dim=2, ncopies=size(b)) * &
                      spread(b, dim=1, ncopies=size(a))

   end function outerprod_d

end module outerprod_m

module unit_matrix_m
    implicit none

contains

    subroutine unit_matrix(mat) ! matrix dimension

        ! Sets matrix to be a unit matrix.
        use argumt

        real(KIND = wReal), DIMENSION(:, :), INTENT(OUT) :: mat
        ! Local:
        integer(KIND = wInt) :: i, n

        n = min(size(mat, 1), size(mat, 2))

        mat(:, :) = 0.0_wReal

        do i=1, n
            mat(i, i) = 1.0_wReal
        end do

    end subroutine unit_matrix

end module unit_matrix_m


! definiendo common/box/boxl,rc,np   common/pot/dl,dT   common/parameters/diam,rho
! en este modulo incluimos mp y mr pq es mas facil cambiar su valor en este modulo y no en todas las subrutinas
module box_pot_parameters
    use argumt ! calling the module argumt becuse we use wReal, wInt
    implicit none
    save
! CONSTANT PARAMETERS
    ! parameters arguments
                                      !packing fraction
    real(KIND = wReal), parameter :: phi = 0.35_wReal, &
                                     pi = 4.d0*datan(1.0_wReal), &
                                     rho = 6._wReal*phi/pi, &
                                     diam = 1.0_wReal
! phi = filling fraction , diam = diameter of the particles
! pi = value with double presicion , rho = reduced density

    ! box arguments
    integer(KIND = wInt), parameter :: np = 5**3_wInt! number of particles

    real(KIND = wReal), parameter :: boxl = (dfloat(np)/rho)**(1._wReal/3._wReal), &
                                     rc = boxl/2._wReal

  ! np = number of particule , boxl = box length, rc = radio corte

    ! pot arguments
    real(KIND = wReal), parameter :: dlr = 50._wReal, dT = 1.4737_wReal
    real(KIND = wReal), parameter :: dla = 49._wReal, deltat = 0.00001_wReal
    real(KIND = wReal), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))

    ! mp and mr arguments
    integer(KIND = wInt), parameter :: mp = 1024_wInt, mr = 2**9_wInt, &
                                        mt = 1000_wInt  !!! si lo aumento  se le acaba la memoria

end module box_pot_parameters

! END OF MODULES


PROGRAM PRINCIPAL
    use argumt
    use box_pot_parameters

    implicit none

    ! Local variables, note that somes variables was initialized
    real(KIND = wReal), dimension(np) :: x, y, z, fx, fy, fz
    real(KIND = wReal), dimension(mr) :: r, g, q, sq, h ! five vector of mr dimension
    real(KIND = wReal), dimension(mt) :: t, wt, ft
    real(KIND = wReal), dimension(mt, np) :: cfx, cfy, cfz

    real(KIND = wReal), parameter :: dr = rc/mr, dq = pi/rc !ancho de un punto antro dentrodelradio de corte
    real(KIND = wReal), parameter :: d = (1.0_wReal/rho)**(1._wReal/3._wReal)
    integer(KIND = wInt), parameter :: k = 3

    real(KIND = wReal), dimension(np*k, np*k) :: dij
    real(KIND = wReal), dimension(np*k) :: Rz


    integer(KIND = wInt) :: i, istep, nprom, j, ncep, nconf, ncp
    real(KIND = wReal) :: ener, enerpot, epotn, dv, fnorm, &
                        graux, hraux, pbc = 1.0_wReal

    integer(KIND = wInt), parameter :: limT = 50000

    print*, 'The length of the box is: ', boxl
    print*, 'The mean interparticle distance is: ', d
    print*, 'Cut radius: ', rc

    call iniconfig(x, y, z, d) ! this subroutine create the data

    open(20, file='iniconfBD.dat', status='unknown')

    do i=1, np
       write(20, '(3f16.8)') x(i), y(i), z(i)
    enddo

    close(20)

    !initial force on the particles
    call force( x, y, z, fx, fy, fz, ener )
    call IH( x, y, z, np, dij, Rz )

    ! Cálculo inicial de interacciones hidrodinámicas, inicializar arreglos
    dij = 0.0_wReal
    Rz = 0.0_wReal

    !Energy of the initial configuration
    print*, 'E/N = ', ener/np

     !The system will thermalize

     !Periodic boundary conditions; pbc > 0
    open(51, file = 'energy_BD.dat', status = 'unknown')
    do istep = 1, limT
        call position_ih( x, y, z, fx, fy, fz, dij, Rz, pbc )
        ! call position( x, y, z, fx, fy, fz, pbc )
        call force( x, y, z, fx, fy, fz, enerpot )
        ! Revisar posibles indeterminaciones
        call check_nan(x, np)
        call check_nan(y, np)
        call check_nan(z, np)
        call check_nan(fx, np)
        call check_nan(fy, np)
        call check_nan(fz, np)
        call IH( x, y, z, np, dij, Rz )
        epotn = enerpot/dfloat(np)
        if (mod(istep, 10000) .eq. 0) then
            print*, istep, epotn, 'Thermal'
        end if
        if (mod(istep, 1000) .eq. 0) then 
            write(51,'(3f16.8)') dfloat(istep), epotn
        end if
    enddo
    close(51)

    print*, 'The system has thermalized'

    open(60, file = 'finalconBD.dat', status = 'unknown')
    do i = 1,np
        write(60,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
    enddo
    close(60)
    !Thermalization ends

    !cycle to calculate the g(r), h, sq, t, wt, ft
    g = 0.0_wReal
    h = 0.0_wReal
    sq = 0.0_wReal
    t = 0.0_wReal
    wt = 0.0_wReal
    ft = 0.0_wReal

    ncep = 1 ! es el paso a la hora de guardar los datos
    ncp = 100000 !cantidad de configuraciones
    nprom = 0
    nconf = ncp
    pbc = 0.0_wReal

    ! call IH( x, y, z, np, dij, Rz )

    do i = 1, nconf
        call position_ih( x, y, z, fx, fy, fz, dij, Rz, pbc )
        ! call position(x, y, z, fx, fy, fz, pbc)
        call force(x, y, z, fx, fy, fz, enerpot)
        ! Revisar posibles indeterminaciones
        call check_nan(x, np)
        call check_nan(y, np)
        call check_nan(z, np)
        call check_nan(fx, np)
        call check_nan(fy, np)
        call check_nan(fz, np)
        call IH( x, y, z, np, dij, Rz )
        if ( mod(i,10000) .eq. 0 ) then
            print*, i, enerpot/np, 'Average'
        end if
        if ( mod(i, ncep) .eq. 0 ) then
            nprom = nprom + 1
            ! t(nprom) = deltat*ncep*(nprom-1)
            ! do j = 1, np
            !     cfx(nprom, j) = x(j)
            !     cfy(nprom, j) = y(j)
            !     cfz(nprom, j) = z(j)
            ! enddo
            call gr(x,y,z,g,dr)
        endif
    enddo

    open(65,file='gr_BD.dat',status='unknown')
      write(65,'(3f16.8)') r(1), g(1)

!      print*,dr,nprom

      do i=2,mr
         r(i)=(i-1)*dr
         q(i)=(i-1)*dq
         dv=4.d0*pi*r(i)**2.*dr
         fnorm=boxl**3./(np**2*nprom*dv)
         graux=g(i)*fnorm
         hraux=graux-1.d0
         g(i)=graux
         h(i)=hraux
         write(65,'(3f16.8)')r(i),graux,hraux
      enddo
      close(65)

END PROGRAM PRINCIPAL


!--------------------------------
! START DEFINITIONS OF SUBROUTINES
!---------------------------------

! This subroutine calculates the initial configuration in 3D.
subroutine iniconfig(xc, yc, zc, d)
    use argumt

! Global variables
    use box_pot_parameters  ! gives diam, rho, rc, np, mp
    implicit none

! defining three vector of mp dimension, it indicate that only are out variables
    real(KIND = wReal), intent(out), dimension(np) :: xc, yc, zc
    real(KIND = wReal), intent(in) :: d
    ! Local variables
    integer(KIND = wInt) :: i ! it is neccesary for restart i

    xc(1) = -(boxl-d)/2.0_wReal
    yc(1) = -(boxl-d)/2.0_wReal
    zc(1) = -(boxl-d)/2.0_wReal

    do i = 2,np ! note that 'i' was defined as integer in argumt
        xc(i) = xc(i-1) + d
        yc(i) = yc(i-1)
        zc(i) = zc(i-1)
        if (xc(i) .gt. rc) then  ! .gt. it is >
            xc(i) = xc(1)
            yc(i) = yc(i-1) + d
            if (yc(i) .gt. rc) then
                xc(i) = xc(1)
                yc(i) = yc(1)
                zc(i) = zc(i-1) + d
            endif
        endif
    enddo

return
end subroutine iniconfig

subroutine  force(x, y, z, fx, fy, fz, ener)
    use argumt

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl

    implicit none

    real(KIND = wReal), intent(in), dimension(np) :: x, y, z
    real(KIND = wReal), intent(out), dimension(np) :: fx, fy, fz
    real(KIND = wReal), intent(out) :: ener

! Local variables
    integer(KIND = wInt) :: i, j
    real(KIND = wReal) :: rij, xij, yij, zij,uij,fij,fxij,fyij,fzij, rij2

    ener = 0._wReal ! initializing
    do i = 1, np
       fx(i) = 0._wReal
       fy(i) = 0._wReal
       fz(i) = 0._wReal
    enddo

    !pair contribution
    do i = 1, np-1
        do j = i + 1, np
            uij = 0._wReal
            fij = 0._wReal
            fxij = 0._wReal
            fyij = 0._wReal
            fzij = 0._wReal
            xij = x(i)-x(j)
            yij = y(i)-y(j)
            zij = z(i)-z(j)
            xij = xij-boxl*anint(xij/boxl)
            yij = yij-boxl*anint(yij/boxl)
            zij = zij-boxl*anint(zij/boxl)
            rij2 = xij*xij+yij*yij+zij*zij
            rij = sqrt(rij2)  !sqrt -> square root
            if (rij .lt. rc) then ! .lt. it is <
                call hardsphere(rij,uij,xij,yij,zij,fxij,fyij,fzij)
                ener = ener+uij
                fx(i) = fx(i)+fxij
                fy(i) = fy(i)+fyij
                fz(i) = fz(i)+fzij
                fx(j) = fx(j)-fxij
                fy(j) = fy(j)-fyij
                fz(j) = fz(j)-fzij
            endif
        enddo
    enddo
return
end subroutine force

subroutine hardsphere(rij, uij, xij, yij, zij, fxij, fyij, fzij)
    use argumt

    ! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl

    implicit none

    real(KIND = wReal), intent(in) :: rij, xij, yij, zij
    real(KIND = wReal), intent(out) :: fxij, fyij, fzij, uij

! Local variables
    real(KIND = wReal) ::  fij

    if (rij .lt. (dlr/dla)**(1._wReal/(dlr-dla))) then
       uij = (a2/dT)*((1./rij)**dlr-(1./rij)**dla) + 1./dT
       fij = dlr*(1._wReal/rij)**(dlr+1._wReal)-dla*(1._wReal/rij)**(dla+1._wReal)
       fij = (a2/dT)*fij
    else
       uij = 0._wReal
       fij = 0._wReal
    endif

    fxij = fij*xij/rij
    fyij = fij*yij/rij
    fzij = fij*zij/rij

return
end subroutine hardsphere

! This sobroutine computes the mean-square displacement

subroutine difusion(nprom, cfx, cfy, cfz, t, wt, ft)
    use argumt

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl, dT

    implicit none

    real(KIND = wReal), intent(in), dimension(mt, np) :: cfx, cfy, cfz
    real(KIND = wReal), intent(inout), dimension(mt) :: wt, ft
    real(KIND = wReal), intent(in), dimension(mt) :: t
    integer(KIND = wInt), intent(in):: nprom
! Local variables
    integer(KIND = wInt) :: i, j, k
    real(KIND = wReal) :: dif2, dself, dx, dy, dz, dk, aux, aux2

!Mean-square displacement and intermediate scattering function
    dk = 6.6_wReal
    do i = 1, nprom-1
        dif2 = 0._wReal
        dself = 0._wReal
        do j = 1, nprom-i
            do k = 1,np
                dx = cfx(j+i,k)-cfx(j,k)
                dy = cfy(j+i,k)-cfy(j,k)
                dz = cfz(j+i,k)-cfz(j,k)
                dif2 = dif2+dx*dx+dy*dy+dz*dz
                aux = sqrt(dx*dx+dy*dy+dz*dz)
                aux = sin(dk*aux)/(dk*aux)
                dself = dself+aux
            enddo
        enddo
        aux2 = (np*(nprom-i))
        dif2 = dif2/aux2
        dself = dself/aux2
        wt(i) = wt(i)+dif2
        ft(i) = ft(i)+dself
    enddo
return
end subroutine difusion

! This subroutine calculates the g(r)

subroutine gr(x,y,z,g,dr)

    use argumt

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl, dT

    implicit none

    real(KIND = wReal), intent(in), dimension(np) ::  x, y, z
    real(KIND = wReal), intent(inout), dimension(mr) :: g
    real(KIND = wReal), intent(in) :: dr
! Local variables
    integer(KIND = wInt) :: i, j, nbin
    real(KIND = wReal) :: xij, yij, zij, rij2, rij

    do i = 1,np-1
        do j = i+1,np
            xij = x(j)-x(i)
            yij = y(j)-y(i)
            zij = z(j)-z(i)
            xij = xij-boxl*anint(xij/boxl)
            yij = yij-boxl*anint(yij/boxl)
            zij = zij-boxl*anint(zij/boxl)
            rij2 = xij*xij+yij*yij+zij*zij
            rij = sqrt(rij2)
            if (rij .lt. rc) then
                nbin = nint(rij/dr)+1
                if (nbin .le. mr) then
                    g(nbin) = g(nbin)+2.0_wReal
                endif
            endif
        enddo
    enddo
return
end subroutine gr
