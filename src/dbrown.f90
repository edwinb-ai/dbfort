program principal
use types
use box_pot_parameters
use tensor, only: ih
use utils, only: save_timeseries, iniconfig, show_m
use positions, only: position, position_ih
use energies
use observables

implicit none

! Variables locales
real(dp), allocatable :: x(:), y(:), z(:)
real(dp), allocatable :: fx(:), fy(:), fz(:)
real(dp) :: r(mr), g(mr), h(mr)
real(dp), allocatable :: t(:), wt(:), ft(:)
! Deben ser `allocatable` dado que son arreglos grandes
real(dp), allocatable :: cfx(:,:), cfy(:,:), cfz(:,:)

real(dp) :: dr, dq
real(dp) :: d
integer, parameter :: k = 3

real(dp), allocatable :: dij(:, :)
real(dp), allocatable :: Rz(:)

integer :: i, istep, nprom, j, ncep, ncp
real(dp) :: ener, enerpot, epotn, dv, fnorm
! real(dp) :: graux, hraux
logical :: pbc = .true.
integer :: u

integer, parameter :: limT = 300000

! Leer de un archivo de entrada los valores del usuario
open(newunit=u, file = 'entrada.in', status = 'old')
read(u, *) phi, np, with_ih, ncp, ncep
close(u)
! Actualizar los parámetros de simulación
rho = 6.0_dp * real(phi) / pi
boxl = (np/rho)**(1._dp/3._dp)
rc = boxl / 2.0_dp
d = (1.0_dp/rho)**(1._dp/3._dp)
dr = rc/mr
dq = pi/rc

print*, 'The length of the box is: ', boxl
print*, 'The mean interparticle distance is: ', d
print*, 'Cut radius: ', rc

if (with_ih) then
    print*, 'With hydrodynamic interactions'
else
    print*, 'Without hydrodynamic interactions'
end if

! Inicializar la memoria de los arreglos
allocate( x(np), y(np), z(np), fx(np), fy(np), fz(np) )
allocate( dij(np*k, np*k) )
allocate( Rz(np*k) )

! Crear la configuración inicial de malla
call iniconfig(x, y, z, d)

! open(newunit=u, file = 'finalconBD.dat', status = 'unknown')
! do i = 1,np
!     read(u,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
! end do
! close(u)

! Cálculo inicial de interacciones hidrodinámicas, inicializar arreglos
dij = 0.0_dp
Rz = 0.0_dp

! Calcular la energía y fuerza iniciales
call force( x, y, z, fx, fy, fz, ener )
if (with_ih) call ih( x, y, z, np, dij, Rz )

!Energy of the initial configuration
print*, 'E/N = ', ener/np

! Comienza la termalización
!Periodic boundary conditions; pbc > 0
open(newunit=u, file = 'energy_BD.dat', status = 'unknown')
do istep = 1, limT
    call position( x, y, z, fx, fy, fz, pbc )
    call force( x, y, z, fx, fy, fz, enerpot )
    epotn = enerpot/real(np)
    if (mod(istep, 100000) == 0) then
        print*, istep, epotn, 'Thermal'
    end if
    if (mod(istep, 10000) == 0) then 
        write(u,'(2f16.8)') real(istep), epotn
    end if
end do
close(u)

print*, 'The system has thermalized'

open(newunit=u, file = 'finalconBD.dat', status = 'unknown')
do i = 1,np
    write(u,'(3f16.8)') x(i), y(i), z(i) !guardo mi foto final
end do
close(u)
!Thermalization ends

! Aquí se inicializan los arreglos para el cálculo de las observables
g = 0.0_dp
h = 0.0_dp

! ncep = 10
! ncp = 200000
mt = ncp / ncep
nprom = 0
pbc = .false.

allocate( cfx(mt,np),cfy(mt,np),cfz(mt,np) )
allocate( t(mt), wt(mt), ft(mt) )

cfx = 0.0_dp
cfy = 0.0_dp
cfz = 0.0_dp
t = 0.0_dp
wt = 0.0_dp
ft = 0.0_dp

! Se inicializa el tensor de difusión y los números aleatorios
if (with_ih) call ih( x, y, z, np, dij, Rz )

! Ciclos de promediación para observables
do i = 1, ncp
    if (with_ih) then
        call position_ih( x, y, z, fx, fy, fz, dij, Rz, pbc )
    else
        call position( x, y, z, fx, fy, fz, pbc )
    end if
    ! Siempre se calcula la energía de la misma forma
    call force( x, y, z, fx, fy, fz, enerpot )
    if (with_ih) call ih( x, y, z, np, dij, Rz )

    if ( mod(i, 100000) == 0 ) then
        print*, i, enerpot/np, 'Average'
    end if
    if ( mod(i, ncep) == 0 ) then
        nprom = nprom + 1
        t(nprom) = deltat*ncep*(nprom-1)
        cfx(nprom, :) = x
        cfy(nprom, :) = y
        cfz(nprom, :) = z
        ! call gr( x,y,z,g,dr,1 ) ! Siempre con PBC
    end if
end do

! open(newunit=u, file='gr_BD.dat', status='unknown')
! write(u,'(3f16.8)') r(1), g(1)

! do i=2,mr
!     r(i)=(i-1)*dr
!     dv=4._dp*pi*r(i)**2.0_dp*dr
!     fnorm=boxl**3._dp/( np**2.0_dp * nprom*dv )
!     graux=g(i)*fnorm
!     hraux=graux-1._dp
!     g(i)=graux
!     h(i)=hraux
!     write(u,'(3f16.8)')r(i),graux,hraux
! end do
! close(u)


call difusion( nprom,cfx,cfy,cfz,wt,ft )

open(newunit=u,file='wt_fself.dat',status='unknown')
do i=1,nprom-1
    write(u,'(f10.8,A,f10.8,A,f11.8)') t(i+1),',',wt(i),',',ft(i)
end do
close(u)

print*, "Saving MSD to files..."
call save_timeseries( 'msd_data/msd_',cfx,cfy,cfz )
print*, "Done!"

! Desalojar toda la memoria utilizada
deallocate( cfx,cfy,cfz,x,y,z )
deallocate( fx,fy,fz,dij,Rz )

end program principal
