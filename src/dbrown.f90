! definiendo common/box/boxl,rc,np   common/pot/dl,dT   common/parameters/diam,rho
! en este modulo incluimos mp y mr pq es mas facil cambiar su valor en este modulo y no en todas las subrutinas


! END OF MODULES


PROGRAM PRINCIPAL
    use iso_fortran_env, only: real64, int32
    use box_pot_parameters
    ! use tensor, only: ih
    use utils, only: save_timeseries
    use positions, only: position

    implicit none

    ! Variables locales
    real(real64), dimension(np) :: x, y, z, fx, fy, fz
    real(real64), dimension(mr) :: r, g, h
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
    real(real64) :: graux, hraux
    integer(int32) :: pbc = 1

    integer(int32), parameter :: limT = 1000000

    print*, 'The length of the box is: ', boxl
    print*, 'The mean interparticle distance is: ', d
    print*, 'Cut radius: ', rc

    call iniconfig(x, y, z, d)

    ! open(60, file = 'finalconBD_035.dat', status = 'unknown')
    ! do i = 1,np
    !     read(60,'(3f16.8)') x(i), y(i), z(i)
    ! end do
    ! close(60)

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
        if (mod(istep, 100000) == 0) then
            print*, istep, epotn, 'Thermal'
        end if
        if (mod(istep, 10000) == 0) then 
            write(51,'(3f16.8)') real(istep), epotn
        end if
    enddo
    close(51)

    print*, 'The system has thermalized'

    open(60, file = 'finalconBD_045.dat', status = 'unknown')
    do i = 1,np
        write(60,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
    enddo
    close(60)
    !Thermalization ends

    !cycle to calculate the g(r), h, sq, t, wt, ft
    g = 0.0d0
    h = 0.0d0
    t = 0.0d0
    wt = 0.0d0
    ft = 0.0d0

    ncep = 100
    ncp = 3000000
    nprom = 0
    pbc = 0

    allocate( cfx(mt,np),cfy(mt,np),cfz(mt,np) )
    cfx = 0.0d0
    cfy = 0.0d0
    cfz = 0.0d0

    do i = 1, ncp
        ! call position_ih( x, y, z, fx, fy, fz, dij, Rz, pbc )
        call position( x, y, z, fx, fy, fz, pbc )
        call force( x, y, z, fx, fy, fz, enerpot )
        ! call IH( x, y, z, np, dij, Rz )
        if ( mod(i, 100000) == 0 ) then
            print*, i, enerpot/np, 'Average'
        end if
        if ( mod(i, ncep) == 0 ) then
            nprom = nprom + 1
            t(nprom) = deltat*ncep*(nprom-1)
            do j = 1, np
                cfx(nprom, j) = x(j)
                cfy(nprom, j) = y(j)
                cfz(nprom, j) = z(j)
            end do
            call gr( x,y,z,g,dr,pbc ) ! Con PBC
        end if
    end do

    open(65,file='gr_BD_045_ih.dat',status='unknown')
      write(65,'(3f16.8)') r(1), g(1)

    do i=2,mr
        r(i)=(i-1)*dr
        dv=4.d0*pi*r(i)**2.0d0*dr
        fnorm=boxl**3.d0/( np**2.0d0 * nprom*dv )
        graux=g(i)*fnorm
        hraux=graux-1.d0
        g(i)=graux
        h(i)=hraux
        write(65,'(3f16.8)')r(i),graux,hraux
    enddo
    close(65)

    
    call difusion( nprom,cfx,cfy,cfz,t,wt,ft )
    open(80,file='wt_p125_phi045.dat',status='unknown')

    do i=1,(ncp/ncep)-1
        write(80,"(3f16.8)") t(i+1),wt(i),ft(i)
    enddo

    close(80)

    ! print*, "Saving MSD to files..."
    ! call save_timeseries( 'msd_data/msd_',cfx,cfy,cfz )
    deallocate( cfx,cfy,cfz )
    ! print*, "Done!"

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
    real(real64) :: rij, xij, yij, zij, uij
    real(real64) :: fij, fxij, fyij, fzij

    ! Inicializar arreglos y variables
    ener = 0.0d0
    fx(:) = 0.0d0
    fy(:) = 0.0d0
    fz(:) = 0.0d0

    do i = 1, np-1
        do j = i+1, np
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
            
            rij = norm2( [xij, yij, zij] )
            ! rij = dsqrt(rij2)
            
            if (rij < rc) then
                call hardsphere( rij,uij,xij,yij,zij,fxij,fyij,fzij )
                
                ener = ener+uij
                fx(i) = fx(i)+fxij
                fy(i) = fy(i)+fyij
                fz(i) = fz(i)+fzij
                
                fx(j) = fx(j)-fxij
                fy(j) = fy(j)-fyij
                fz(j) = fz(j)-fzij
            end if
        end do
    end do
end subroutine force

subroutine hardsphere(rij, uij, xij, yij, zij, fxij, fyij, fzij)
    use iso_fortran_env, only: real64, int32

    ! Global variables
    use box_pot_parameters, only: dlr, dla, dT, a2, bpot

    implicit none

    real(real64), intent(in) :: rij, xij, yij, zij
    real(real64), intent(inout) :: fxij, fyij, fzij, uij

! Local variables
    real(real64) ::  fij

    if (rij < bpot) then
       uij = (a2/dT)*((1.0d0/rij)**dlr-(1.0d0/rij)**dla)
       uij = uij + 1.0d0/dT
       fij = dlr*(1.0d0/rij)**(dlr+1.0d0)-dla*(1.0d0/rij)**(dla+1.0d0)
       fij = a2*fij/dT
    else
       uij = 0.0d0
       fij = 0.0d0
    end if

    fxij = fij*xij/rij
    fyij = fij*yij/rij
    fzij = fij*zij/rij
end subroutine hardsphere

! This sobroutine computes the mean-square displacement

subroutine difusion(nprom, cfx, cfy, cfz, t, wt, ft)
    use iso_fortran_env, only: real64, int32

! Global variables
    use box_pot_parameters, only: mt, np

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
        dif2 = 0.0d0
        dself = 0.0d0
        do j = 1, nprom-i
            do k = 1,np
                dx = cfx(j+i,k)-cfx(j,k)
                dy = cfy(j+i,k)-cfy(j,k)
                dz = cfz(j+i,k)-cfz(j,k)
                aux = dx*dx + dy*dy + dz*dz
                dif2 = dif2 + aux
                aux2 = dk*dsqrt(aux)
                aux = dsin(aux2)/aux2
                dself = dself+aux
            end do
        end do
        aux2 = real(np*(nprom-i))
        dif2 = dif2/aux2
        dself = dself/aux2
        wt(i) = wt(i)+dif2
        ft(i) = ft(i)+dself
    end do
end subroutine difusion

! This subroutine calculates the g(r)

subroutine gr(x,y,z,g,dr,pbc)
    use iso_fortran_env, only: real64, int32
    use box_pot_parameters, only: np, mr, boxl, rc
    implicit none

    real(real64), intent(in), dimension(np) ::  x, y, z
    real(real64), intent(inout), dimension(mr) :: g
    real(real64), intent(in) :: dr
    integer, intent(in) :: pbc
    
    ! Local variables
    integer(int32) :: i, j, nbin
    real(real64) :: xij, yij, zij, rij

    do i = 1,np-1
        do j = i+1,np

            xij = x(i)-x(j)
            yij = y(i)-y(j)
            zij = z(i)-z(j)

            if ( pbc == 1 ) then
                xij = xij-boxl*nint(xij/boxl)
                yij = yij-boxl*nint(yij/boxl)
                zij = zij-boxl*nint(zij/boxl)
            end if

            rij = norm2( [xij,yij,zij] )
            if (rij < rc) then
                nbin = nint(rij/dr) + 1
                if (nbin <= mr) then
                    g(nbin) = g(nbin) + 2.0d0
                end if
            end if
        end do
    end do
end subroutine gr
