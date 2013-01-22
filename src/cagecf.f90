! This program computes the cage correlation function as defined in
! XX
! Original program by Mathieu Salanne
! Rewritten in fortran 2008 with improved performance by Maximilien Levesque

MODULE MYMOD
double precision, allocatable, dimension(:) :: x, y, z
double precision :: cutoff, width
double precision :: boxlen, halfbox, halfboxrec
double precision :: dtime, ndump
integer, dimension(:), allocatable :: ntype
integer, dimension(:), allocatable :: nfilecf
integer, dimension(:), allocatable :: nion
integer :: ncorrtime, ncorr, ncorrcall
integer :: nrun,npereng,nstep,npervel,nperfri,nskip
integer :: num, ncut
integer :: nspecies, ncoordav,ntot
integer, allocatable, dimension(:,:) :: ncfin, ncfout, nthetainout, nthetaout, norm
integer(kind=1), allocatable, dimension(:,:,:) :: istorederfunc ! kind=1 really improves performance cause derfunc only contains 0 or 1.
END MODULE MYMOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM CAGECF
use mymod
implicit none
integer :: nummax, nspmax, nfiles, ncnt, nfilecnt
 character(len=80), allocatable, dimension(:) :: posfile
integer :: i, j, k
double precision :: du

ncoordav = 0

open(10,file='cagecf.inpt',status='old')

read(10,*) nspecies
if( nspecies < 1 ) stop 'nspecies should be higher than 0'
allocate( nion(nspecies) )
nspmax = nspecies

! read number of atoms of each species
do i=1,nspecies
  read(10,*) nion(i)
  if( nion(i)<1 ) stop 'nb of ions should be >0 for all ions'
end do

! allocate accordingly
nummax = sum(nion)
allocate( x(nummax), y(nummax), z(nummax) ) ! coordinates of all atoms
allocate( ntype(nummax) ) ! nummax is the total number of ions in the simulation box

! fullfill the type of each atom in the box
num=0
do i=1,nspecies
  do j=num+1,num+nion(i)
    ntype(j)=i
  end do
  num=num+nion(i)
end do

read(10,*)nskip
if( nskip < 1) stop 'nskip should be > 0'
read(10,*)boxlen
if( boxlen < 0. ) stop 'box length should be positive'
read(10,*)nfiles ! number of position files
if( nfiles < 1 ) stop 'nb of position files should be >=1'

allocate( posfile(nfiles), nfilecf(nfiles) )

do i=1,nfiles
  read(10,'(a)')posfile(i)
  read(10,*)nfilecf(i)
end do

read(10,*)ncorr ! Time-span of the correlation function (in steps)
allocate( istorederfunc(nummax,nummax,ncorr) )
allocate( ncfin(nummax,0:ncorr) )
allocate( ncfout(nummax,0:ncorr) )
allocate( nthetainout(nummax,0:ncorr) )
allocate( nthetaout(nummax,0:ncorr) )
allocate( norm(nummax,0:ncorr) )


read(10,*)cutoff ! cutoff for bond
read(10,*)ncut ! number of ions for change in cage
read(10,*)dtime ! MD timestep
read(10,*)ndump ! number of configs between MD dumps

 close(10)

! Auxillary variables.
ncorrtime=1
ncorrcall=0

halfbox = boxlen/2.0d0
halfboxrec = 1.0d0/halfbox

ntot = sum(nfilecf)
ncnt = 0
nfilecnt = 1

open(10,file=posfile(1),status='old')

loopj: do j=1,ntot
  print*, j
  if ((mod(j,nskip))==0) then
    do k=1,nummax
      read(10,*) x(k), y(k), z(k)
    end do
    call ircfcalc
    ncorrtime = mod(ncorrtime,ncorr)
    ncorrcall = ncorrcall+1
    ncorrtime = ncorrtime+1
    ncnt=ncnt+1
  else
    do k=1,nummax
      read(10,*)du,du,du
    end do
    ncnt=ncnt+1
  endif

  if (ncnt == nfilecf(nfilecnt)) then
    if (nfilecnt == nfiles) exit loopj
    nfilecnt=nfilecnt+1
    close(10)
    open(10,file=posfile(nfilecnt),status='old')
    ncnt = 1
  endif
end do loopj

deallocate(istorederfunc, ncfin, ncfout, x, y, z)


call ircfdump


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! where cage correlation function is calculated
SUBROUTINE IRCFCALC
use mymod
implicit none
integer(kind=1), dimension(nummax,nummax) :: derfunc
integer :: i, j, m, nerfc, nt, imin, imax, kmin,kmax,mmin,mmax
double precision :: xi, yi, zi, dx, dy, dz, dr
integer(kind=1) :: t_istorederfunc_ncorrtime, t_istorederfunc_m, t_mult
integer(kind=1), allocatable, dimension(:,:) :: i_istorederfunc

derfunc = 0 ! init

imin=nion(1)+1
imax=nion(1)+nion(2)
kmin=1
kmax=nion(1)

do i=imin,imax
  xi=x(i) ; yi=y(i) ; zi=z(i)

  do k=kmin,kmax
    dx=xi-x(k) ; dy=yi-y(k) ; dz=zi-z(k)
    ! apply minimum image convention
    dx=dx-boxlen*int(dx*halfboxrec)
    dy=dy-boxlen*int(dy*halfboxrec)
    dz=dz-boxlen*int(dz*halfboxrec)
    dr = norm2((/dx,dy,dz/))
    if( dr <= cutoff ) then
      derfunc(i,k) = 1
      ncoordav = ncoordav+1
    else
      derfunc(i,k) = 0
    endif
  end do

end do

mmin=1
mmax=ncorrtime

istorederfunc(imin:imax,kmin:kmax,ncorrtime) = derfunc(imin:imax,kmin:kmax)
allocate(i_istorederfunc(kmin:kmax,mmin:mmax))

! Store anion and cation correlation functions.
do i=imin,imax

  i_istorederfunc(kmin:kmax,mmin:mmax) = istorederfunc(i,kmin:kmax,mmin:mmax)

  do m=mmin,mmax
    nt=ncorrtime-m
    norm(i,nt)=norm(i,nt)+1
    ncfout(i,nt)=0
    ncfin(i,nt)=0
    
    do k=kmin,kmax
      t_istorederfunc_ncorrtime = derfunc(i,k) ! seems wierd but is correct
      t_istorederfunc_m = i_istorederfunc(k,m)
      t_mult = t_istorederfunc_ncorrtime * t_istorederfunc_m
      ! Auto-correlation terms.
      ncfout(i,nt) = ncfout(i,nt) + t_istorederfunc_m**2 - t_mult
      ncfin(i,nt)  = ncfin(i,nt)  + t_istorederfunc_ncorrtime**2 - t_mult
     end do

     if ((ncfout(i,nt))<ncut) nthetaout(i,nt)=nthetaout(i,nt)+1
     if (((ncfin(i,nt))<ncut).and.((ncfout(i,nt))<ncut)) nthetainout(i,nt)=nthetainout(i,nt) +1

   end do
end do


! Secondary loop.
if (ncorrcall > ncorr) then

  do i=imin,imax

    do m=ncorrtime+1,ncorr
      nt=ncorrtime-m+ncorr
      norm(i,nt)=norm(i,nt)+1
      ncfout(i,nt)=0
      ncfin(i,nt)=0

      do k=kmin,kmax
        t_istorederfunc_ncorrtime = istorederfunc(i,k,ncorrtime)
        t_istorederfunc_m = istorederfunc(i,k,m)
        t_mult = t_istorederfunc_ncorrtime * t_istorederfunc_m
        ! Auto-correlation terms.
        ncfout(i,nt)=ncfout(i,nt) + t_istorederfunc_m**2 - t_mult
        ncfin(i,nt)=ncfin(i,nt) + t_istorederfunc_ncorrtime**2 - t_mult
      end do

      if ((ncfout(i,nt))<ncut) nthetaout(i,nt)=nthetaout(i,nt)+1
      if (((ncfin(i,nt))<ncut) .and. ((ncfout(i,nt))<ncut)) nthetainout(i,nt)=nthetainout(i,nt)+1
    end do
  end do
endif

END SUBROUTINE IRCFCALC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! where everything gets printed
SUBROUTINE IRCFDUMP
use mymod
implicit none
double precision, allocatable, dimension(:) :: cfinouttot, cfouttot
integer :: i, j, ns
double precision :: coordav, time, t_norm

allocate( cfinouttot(nspecies) )
allocate( cfouttot(nspecies) )
coordav=dble(ncoordav*nskip)/dble(ntot*nion(2))

open (60,file='cagecfinout.dat')
open (61,file='cagecfout.dat')

do j=0,ncorr-1

  cfinouttot = 0.0d0
  cfouttot = 0.0d0

  do i=nion(1)+1,nion(1)+nion(2)
    ns=ntype(i)
    t_norm = dble(norm(i,j))
    if( t_norm/=0 ) then
      cfouttot(ns)=cfouttot(ns)+dble(nthetaout(i,j)) / t_norm
      cfinouttot(ns)=cfinouttot(ns)+dble(nthetainout(i,j)) / t_norm
    end if
  end do

  time=dble(j)*dble(nskip)*dble(ndump)*dtime*2.418d-5
  write(60,*)time,cfinouttot(ns)/dble(nion(ns))
  write(61,*)time,cfouttot(ns)/dble(nion(ns))

end do

 close (60)
 close (61)

print*,
print*, 'Average coordination no. is',coordav
print*,
print*, '**** Finished. Correlation functions written in ****'
print*, 'cagecfout.dat   and   cagecfinout.dat'
print*,
print*,

END SUBROUTINE IRCFDUMP

END PROGRAM CAGECF
