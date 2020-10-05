program principal
    use types
    use box_pot_parameters
    ! use tensor, only: ih
    use utils, only: save_timeseries, iniconfig
    use positions, only: position
    use energies
    use observables

    implicit none

    ! Variables locales
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable :: fx(:), fy(:), fz(:)
    real(dp) :: r(mr), g(mr), h(mr)
    real(dp) :: t(mt), wt(mt), ft(mt)
    ! Deben ser `allocatable` dado que son arreglos grandes
    real(dp), allocatable :: cfx(:,:), cfy(:,:), cfz(:,:)

    real(dp) :: dr, dq
    real(dp) :: d
    integer, parameter :: k = 3

    real(dp), allocatable :: dij(:, :)
    real(dp), allocatable :: Rz(:)

    integer :: i, istep, nprom, j, ncep, ncp
    real(dp) :: ener, enerpot, epotn, dv, fnorm
    real(dp) :: graux, hraux
    integer :: pbc = 1

    integer, parameter :: limT = 1000000

    ! Leer de un archivo de entrada los valores del usuario
    integer :: u
    open(newunit=u, file = 'entrada.in', status = 'old')
    read(u, *) phi, np
    close(u)
    ! Actualizar los parámetros de simulación
    rho = 6.0_dp * phi / pi
    boxl = (np/rho)**(1._dp/3._dp)
    rc = boxl / 2.0_dp
    d = (1.0_dp/rho)**(1._dp/3._dp)
    dr = rc/mr
    dq = pi/rc

    print*, 'The length of the box is: ', boxl
    print*, 'The mean interparticle distance is: ', d
    print*, 'Cut radius: ', rc

    ! Inicializar la memoria de los arreglos
    allocate( x(np), y(np), z(np), fx(np), fy(np), fz(np) )
    allocate( dij(np*k, np*k) )
    allocate( Rz(np*k) )
    allocate( cfx(mt,np),cfy(mt,np),cfz(mt,np) )

    ! Crear la configuración inicial de malla
    call iniconfig(x, y, z, d)

    ! Calcular la energía y fuerza iniciales
    call force( x, y, z, fx, fy, fz, ener )
    ! call IH( x, y, z, np, dij, Rz )

    ! Cálculo inicial de interacciones hidrodinámicas, inicializar arreglos
    dij = 0.0_dp
    Rz = 0.0_dp

    !Energy of the initial configuration
    print*, 'E/N = ', ener/np

    ! Comienza la termalización
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
    end do
    close(51)

    print*, 'The system has thermalized'

    open(60, file = 'finalconBD.dat', status = 'unknown')
    do i = 1,np
        write(60,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
    end do
    close(60)
    !Thermalization ends

    ! Aquí se inicializan los arreglos para el cálculo de las observables
    g = 0.0_dp
    h = 0.0_dp
    t = 0.0_dp
    wt = 0.0_dp
    ft = 0.0_dp

    ncep = 100
    ncp = 3000000
    nprom = 0
    pbc = 0
    
    cfx = 0.0_dp
    cfy = 0.0_dp
    cfz = 0.0_dp

    ! Ciclos de promediación para observables
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
            call gr( x,y,z,g,dr,1 ) ! Siempre con PBC
        end if
    end do

    open(newunit=u, file='gr_BD.dat', status='unknown')
    write(u,'(3f16.8)') r(1), g(1)

    do i=2,mr
        r(i)=(i-1)*dr
        dv=4._dp*pi*r(i)**2.0_dp*dr
        fnorm=boxl**3._dp/( np**2.0_dp * nprom*dv )
        graux=g(i)*fnorm
        hraux=graux-1._dp
        g(i)=graux
        h(i)=hraux
        write(u,'(3f16.8)')r(i),graux,hraux
    end do
    close(u)

    
    call difusion( nprom,cfx,cfy,cfz,t,wt,ft )
    open(newunit=u,file='wt_fself.dat',status='unknown')

    do i=1,(ncp/ncep)-1
        write(u,"(3f16.8)") t(i+1),wt(i),ft(i)
    end do

    close(u)

    ! print*, "Saving MSD to files..."
    ! call save_timeseries( 'msd_data/msd_',cfx,cfy,cfz )
    deallocate( cfx,cfy,cfz,x,y,z,fx,fy,fz )
    deallocate( dij,Rz )
    ! print*, "Done!"

end program principal
