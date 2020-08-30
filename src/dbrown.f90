! definiendo common/box/boxl,rc,np   common/pot/dl,dT   common/parameters/diam,rho
! en este modulo incluimos mp y mr pq es mas facil cambiar su valor en este modulo y no en todas las subrutinas


! END OF MODULES


PROGRAM PRINCIPAL
    use iso_fortran_env, only: real64, int32
    use box_pot_parameters
    use tensor, only: ih
    use utils, only: save_timeseries

    implicit none

    ! Variables locales
    real(real64), dimension(np) :: x, y, z, fx, fy, fz
    real(real64), dimension(mr) :: r, g, q, sq, h
    real(real64), dimension(mt) :: t, wt, ft
    ! Deben ser `allocatable` dado que son arreglos grandes
    real(real64), allocatable :: cfx(:,:), cfy(:,:), cfz(:,:)

    real(real64), parameter :: dr = rc/mr, dq = pi/rc
    real(real64), parameter :: d = (1.0d0/rho)**(1.d0/3.d0)
    integer(int32), parameter :: k = 3

    real(real64), dimension(np*k, np*k) :: dij
    real(real64), dimension(np*k) :: Rz

    integer(int32) :: i, istep, nprom, j, ncep, ncp
    real(real64) :: ener, enerpot, epotn, dv, fnorm
    real(real64) :: graux, hraux, pbc = 1.0d0

    integer(int32), parameter :: limT = 20000

    print*, 'The length of the box is: ', boxl
    print*, 'The mean interparticle distance is: ', d
    print*, 'Cut radius: ', rc

    ! call iniconfig(x, y, z, d)

    open(60, file = 'finalconBD_040.dat', status = 'unknown')
    do i = 1,np
        read(60,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
    enddo
    close(60)

    !initial force on the particles
    call force( x, y, z, fx, fy, fz, ener )
    ! call IH( x, y, z, np, dij, Rz )

    ! Cálculo inicial de interacciones hidrodinámicas, inicializar arreglos
    dij = 0.0d0
    Rz = 0.0d0

    !Energy of the initial configuration
    print*, 'E/N = ', ener/np

     !The system will thermalize

     !Periodic boundary conditions; pbc > 0
    open(51, file = 'energy_BD.dat', status = 'unknown')
    do istep = 1, limT
        ! call position_ih( x, y, z, fx, fy, fz, dij, Rz, pbc )
        call position( x, y, z, fx, fy, fz, pbc )
        call force( x, y, z, fx, fy, fz, enerpot )
        ! call check_nan(fz, np)
        ! call IH( x, y, z, np, dij, Rz )
        epotn = enerpot/real(np)
        if (mod(istep, 10000) == 0) then
            print*, istep, epotn, 'Thermal'
        end if
        if (mod(istep, 1000) == 0) then 
            write(51,'(3f16.8)') real(istep), epotn
        end if
    enddo
    close(51)

    print*, 'The system has thermalized'

    open(60, file = 'finalconBD_040.dat', status = 'unknown')
    do i = 1,np
        write(60,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
    enddo
    close(60)
    !Thermalization ends

    !cycle to calculate the g(r), h, sq, t, wt, ft
    g = 0.0d0
    h = 0.0d0
    sq = 0.0d0
    t = 0.0d0
    wt = 0.0d0
    ft = 0.0d0

    ncep = 100 ! Modificar este parámetro
    ncp = 3000000
    nprom = 0
    pbc = 0.0d0

    allocate( cfx(mt,np),cfy(mt,np),cfz(mt,np) )

    do i = 1, ncp
        ! call position_ih( x, y, z, fx, fy, fz, dij, Rz, pbc )
        call position( x, y, z, fx, fy, fz, pbc )
        call force( x, y, z, fx, fy, fz, enerpot )
        ! call IH( x, y, z, np, dij, Rz )
        if ( mod(i,100000) == 0 ) then
            print*, i, enerpot/np, 'Average'
        end if
        if ( mod(i, ncep) == 0 ) then
            nprom = nprom + 1
            t(nprom) = deltat*ncep*(nprom-1)
            do j = 1, np
                cfx(nprom, j) = x(j)
                cfy(nprom, j) = y(j)
                cfz(nprom, j) = z(j)
            enddo
            ! call gr( x,y,z,g,dr )
        end if
    end do

!     open(65,file='gr_BD_040_ih.dat',status='old')
!       write(65,'(3f16.8)') r(1), g(1)

! !      print*,dr,nprom

!     do i=2,mr
!         r(i)=(i-1)*dr
!         q(i)=(i-1)*dq
!         dv=4.d0*pi*r(i)**2.0d0*dr
!         fnorm=boxl**3.d0/( np**2.0d0 * nprom*dv )
!         graux=g(i)*fnorm
!         hraux=graux-1.d0
!         g(i)=graux
!         h(i)=hraux
!         write(65,'(3f16.8)')r(i),graux,hraux
!     enddo
!     close(65)

    
    call difusion( nprom,cfx,cfy,cfz,t,wt,ft )
    open(80,file='wt_p512_phi040.dat',status='unknown')

    do i=1,(ncp/ncep)-1
        write(80,"(3f16.8)") t(i+1),wt(i),ft(i)
    enddo

    close(80)

    print*, "Saving MSD to files..."
    call save_timeseries( 'msd_data/msd_',cfx,cfy,cfz )
    deallocate( cfx,cfy,cfz )
    print*, "Done!"

END PROGRAM PRINCIPAL


!--------------------------------
! START DEFINITIONS OF SUBROUTINES
!---------------------------------

! This subroutine calculates the initial configuration in 3D.
subroutine iniconfig(xc, yc, zc, d)
    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters  ! gives diam, rho, rc, np, mp
    implicit none

! defining three vector of mp dimension, it indicate that only are out variables
    real(real64), intent(out), dimension(np) :: xc, yc, zc
    real(real64), intent(in) :: d
    ! Local variables
    integer(int32) :: i

    xc(1) = -(boxl-d)/2.0d0
    yc(1) = -(boxl-d)/2.0d0
    zc(1) = -(boxl-d)/2.0d0

    do i = 2,np
        xc(i) = xc(i-1) + d
        yc(i) = yc(i-1)
        zc(i) = zc(i-1)
        if (xc(i) > rc) then
            xc(i) = xc(1)
            yc(i) = yc(i-1) + d
            if (yc(i) > rc) then
                xc(i) = xc(1)
                yc(i) = yc(1)
                zc(i) = zc(i-1) + d
            endif
        endif
    enddo

return
end subroutine iniconfig

subroutine  force(x, y, z, fx, fy, fz, ener)
    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl

    implicit none

    real(real64), intent(in), dimension(np) :: x, y, z
    real(real64), intent(out), dimension(np) :: fx, fy, fz
    real(real64), intent(out) :: ener

! Local variables
    integer(int32) :: i, j
    real(real64) :: rij, xij, yij, zij,uij,fij,fxij,fyij,fzij, rij2

    ener = 0.d0 ! initializing
    do i = 1, np
       fx(i) = 0.d0
       fy(i) = 0.d0
       fz(i) = 0.d0
    enddo

    !pair contribution
    do i = 1, np-1
        do j = i + 1, np
            uij = 0.d0
            fij = 0.d0
            fxij = 0.d0
            fyij = 0.d0
            fzij = 0.d0
            xij = x(i)-x(j)
            yij = y(i)-y(j)
            zij = z(i)-z(j)
            xij = xij-boxl*dnint(xij/boxl)
            yij = yij-boxl*dnint(yij/boxl)
            zij = zij-boxl*dnint(zij/boxl)
            rij2 = xij*xij+yij*yij+zij*zij
            rij = sqrt(rij2)
            if (rij < rc) then
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
    use iso_fortran_env, only: real64, int32

    ! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl

    implicit none

    real(real64), intent(in) :: rij, xij, yij, zij
    real(real64), intent(out) :: fxij, fyij, fzij, uij

! Local variables
    real(real64) ::  fij

    if (rij < (dlr/dla)**(1.d0/(dlr-dla))) then
       uij = (a2/dT)*((1./rij)**dlr-(1./rij)**dla) + 1./dT
       fij = dlr*(1.d0/rij)**(dlr+1.d0)-dla*(1.d0/rij)**(dla+1.d0)
       fij = (a2/dT)*fij
    else
       uij = 0.d0
       fij = 0.d0
    endif

    fxij = fij*xij/rij
    fyij = fij*yij/rij
    fzij = fij*zij/rij

return
end subroutine hardsphere

! This sobroutine computes the mean-square displacement

subroutine difusion(nprom, cfx, cfy, cfz, t, wt, ft)
    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl, dT

    implicit none

    real(real64), intent(in), dimension(mt, np) :: cfx, cfy, cfz
    real(real64), intent(inout), dimension(mt) :: wt, ft
    real(real64), intent(in), dimension(mt) :: t
    integer(int32), intent(in):: nprom
! Local variables
    integer(int32) :: i, j, k
    real(real64) :: dif2, dself, dx, dy, dz, dk, aux, aux2

!Mean-square displacement and intermediate scattering function
    dk = 6.6d0
    do i = 1, nprom-1
        dif2 = 0.d0
        dself = 0.d0
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

    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters ! gives diam, rho, rc, np, mp, pi, boxl, dl, dT

    implicit none

    real(real64), intent(in), dimension(np) ::  x, y, z
    real(real64), intent(inout), dimension(mr) :: g
    real(real64), intent(in) :: dr
! Local variables
    integer(int32) :: i, j, nbin
    real(real64) :: xij, yij, zij, rij2, rij

    do i = 1,np-1
        do j = i+1,np
            xij = x(j)-x(i)
            yij = y(j)-y(i)
            zij = z(j)-z(i)
            xij = xij-boxl*dnint(xij/boxl)
            yij = yij-boxl*dnint(yij/boxl)
            zij = zij-boxl*dnint(zij/boxl)
            rij2 = xij*xij+yij*yij+zij*zij
            rij = sqrt(rij2)
            if (rij < rc) then
                nbin = dnint(rij/dr)+1
                if (nbin <= mr) then
                    g(nbin) = g(nbin)+2.0d0
                endif
            endif
        enddo
    enddo
return
end subroutine gr
