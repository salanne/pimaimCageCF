! This program computes the cage correlation function as defined in
! Original program by Mathieu Salanne
! Rewritten quickly in fortran 2008 with improved performance by Maximilien Levesque

MODULE MYMOD
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), allocatable, dimension(:) :: x, y, z,xelec,yelec,zelec
    real(dp) :: cutoff, width
    real(dp) :: boxlenx,boxleny,boxlenz
    real(dp) :: dtime, ndump
    integer(1), dimension(:,:), allocatable :: derfunc
    integer, dimension(:), allocatable :: ntype
    integer, dimension(:), allocatable :: nion
    integer :: ncorrtime, ncorr, ncorrcall
    integer :: nrun,npereng,nstep,npervel,nperfri,nskip
    integer :: num, ncut
    integer :: nspecies, ncoordav,nanalyze,nelectrode
    integer, allocatable, dimension(:,:) :: ncfin, ncfout, nthetainout, nthetaout, norm
    integer(1), allocatable, dimension(:,:,:) :: istorederfunc ! kind=1 really improves performance cause derfunc only contains 0 or 1.
    logical :: externalelectrode,dumpinfo,fullelectrode
    character(len=80) :: poselecfile
    type FileWithPositions
        character(len=80) :: name ! file name
        integer :: nbSets ! number of velocity sets in file
    end type
    type (FileWithPositions), allocatable, dimension(:) :: PosFile
END MODULE MYMOD


PROGRAM CAGECF
    use mymod
    use iso_fortran_env, only: dp => real64
    implicit none

    integer :: nummax, ncnt, nfilecnt
    integer :: i, j, k
    integer, parameter :: inputUnit = 10, positionsUnit = 11, inputUnit2 = 20
    
    call checkInputFileExists( "cagecf.inpt")
    call openInputFile
    call readInputFile
    call closeInputFile
!~     call testInputFile

    ncoordav = 0
    ncorrtime=1
    ncorrcall=0

    ncnt = 0
    nfilecnt = 1

    open(inputUnit,file=posFile(1)%name,status='old')
    if(externalelectrode)then
       open(inputUnit2,file=poselecfile,status='old')
       do i=1,nelectrode
          read(inputUnit2,*)xelec(i),yelec(i),zelec(i)
       enddo
    endif
    
    dumpinfo=.true.
    loopj: do j=1, sum(posFile(:)%nbSets) ! sum over all configurations in all position files
        call tellUserHowItAdvances (j)

        if ((mod(j,nskip))==0) then
            do k=1,nummax
                read(inputUnit,*) x(k), y(k), z(k)
            end do
            call ircfcalc
            ncorrtime = mod(ncorrtime,ncorr)
            ncorrcall = ncorrcall+1
            ncorrtime = ncorrtime+1
            ncnt = ncnt+1
        else
            do k=1,nummax
                read(inputUnit,*)
            end do
            ncnt = ncnt+1
        endif

        if (ncnt == posFile(nfilecnt)%nbSets) then
            if (nfilecnt == size(posFile)) exit loopj
            nfilecnt=nfilecnt+1
            close(inputUnit)
            open(positionsUnit,file=posFile(nfilecnt)%name,status='old')
            ncnt = 1
        endif
    end do loopj
    
    deallocate(istorederfunc, ncfin, ncfout, x, y, z)
    if(externalelectrode)then
       deallocate(xelec,yelec,zelec)
    endif
    call ircfdump

    CONTAINS

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! where cage correlation function is calculated
        SUBROUTINE IRCFCALC

            integer :: i, m, nt, imin, imax, kmin,kmax,mmin,mmax,kmin2,kmax2
            double precision :: xi, yi, zi, dx, dy, dz, dr
            integer(1) :: t_istorederfunc_ncorrtime, t_istorederfunc_m, t_mult
            integer(1), allocatable, dimension(:,:) :: i_istorederfunc
            integer :: tmp_ncfout_i_nt, tmp_ncfin_i_nt
            logical :: attrib

            !derfunc = 0 ! init
            if(externalelectrode)then
               if(fullelectrode)then
                  if(nanalyze.eq.1)then
                     if(dumpinfo)write(*,*)'We look at species number 1 - A'
                     imin=1
                     imax=nion(1)
                  elseif(nanalyze.eq.2)then
                     if(dumpinfo)write(*,*)'We look at species number 2 - Im1'
                     imin=nion(1)+1
                     imax=nion(1)+nion(2)
                  elseif(nanalyze.eq.3)then
                     if(dumpinfo)write(*,*)'We look at species number 3 - ACN'
                     imin=nion(1)+nion(2)+1
                     imax=nion(1)+nion(2)+nion(3)
                  endif
               else
                  imin=1
                  imax=nelectrode
                  if(dumpinfo)write(*,*)'The reference species is the electrode atoms'
               endif
            else
               if(nanalyze.eq.1)then
                  if(dumpinfo)write(*,*)'The program cannot calculate the anion-anion cagecf'
                  stop
               elseif(nanalyze.eq.2)then
                  if(dumpinfo)write(*,*)'The reference species is cation 1'
                  imin=nion(1)+1
                  imax=nion(1)+nion(2)
               elseif(nanalyze.eq.3)then
                  if(dumpinfo)write(*,*)'The reference species is cation 2'
                  imin=nion(1)+nion(2)+1
                  imax=nion(1)+nion(2)+nion(3)
               else
                  if(dumpinfo)write(*,*)'The program cannot calculate the cagecf for cation 3 and above'
                  stop
               endif
            endif
            if(externalelectrode)then
               if(fullelectrode)then
                  kmin=1
                  kmax=1
                  kmin2=1
                  kmax2=nelectrode
                  if(dumpinfo)write(*,*)'We would like to know how fast it escapes from the surface'
               else
                  if(nanalyze.eq.1)then
                     if(dumpinfo)write(*,*)'The species which escapes the cage is number 1 - A'
                     kmin=1
                     kmax=nion(1)
                  elseif(nanalyze.eq.2)then
                     if(dumpinfo)write(*,*)'The species which escapes the cage is number 2 - Im'
                     kmin=nion(1)+1
                     kmax=nion(1)+nion(2)
                  elseif(nanalyze.eq.3)then
                     if(dumpinfo)write(*,*)'The species which escapes the cage is number 3 - ACN'
                     kmin=nion(1)+nion(2)+1
                     kmax=nion(1)+nion(2)+nion(3)
                  else
                     if(dumpinfo)write(*,*)'Can work only with species 1, 2 or 3'
                     stop
                  endif
               endif
            else
               if(dumpinfo)write(*,*)'The species which escapes the cage is the anion (1)'
               kmin=1
               kmax=nion(1)
            endif
            dumpinfo=.false.

            !allocate(derfunc(imin:imax,kmin:kmax), )
            if(.not.allocated(derfunc)) allocate(derfunc(imin:imax,kmin:kmax))
            derfunc=0
            if(fullelectrode)then
               do i=imin,imax
                  attrib=.false.
                  xi=x(i) ; yi=y(i) ; zi=z(i)
                  do k=kmin2,kmax2
                     if(attrib)cycle
                     dx=xi-xelec(k) ; dy=yi-yelec(k) ; dz=zi-zelec(k)
                     dx=dx-boxlenx*int(dx*2.0d0/boxlenx)
                     dy=dy-boxleny*int(dy*2.0d0/boxleny)
                     dr = norm2((/dx,dy,dz/))
                     if( dr <= cutoff ) then
                        derfunc(i,1) = 1
                        ncoordav = ncoordav+1
                        attrib=.true.
                     endif
                  enddo
               enddo
            else
               do i=imin,imax
                  if(externalelectrode)then
                     xi=xelec(i) ; yi=yelec(i) ; zi=zelec(i)
                  else
                     xi=x(i) ; yi=y(i) ; zi=z(i)
                  endif
                  do k=kmin,kmax
                     dx=xi-x(k) ; dy=yi-y(k) ; dz=zi-z(k)
                     dx=dx-boxlenx*int(dx*2.0d0/boxlenx)
                     dy=dy-boxleny*int(dy*2.0d0/boxleny)
                     if(.not.externalelectrode)then
                       ! apply minimum image convention
                       dz=dz-boxlenz*int(dz*2.0d0/boxlenz)
                     endif
                     dr = norm2((/dx,dy,dz/))
                     if( dr <= cutoff ) then
                        derfunc(i,k) = 1
                        ncoordav = ncoordav+1
                     else
                        derfunc(i,k) = 0
                     endif
                  end do
               end do
            endif

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
                    tmp_ncfout_i_nt = 0
                    tmp_ncfin_i_nt = 0
                    do k=kmin,kmax
                        t_istorederfunc_ncorrtime = derfunc(i,k) ! seems wierd but is correct
                        t_istorederfunc_m = i_istorederfunc(k,m)
                        t_mult = t_istorederfunc_ncorrtime * t_istorederfunc_m
                        ! Auto-correlation terms.
                        tmp_ncfout_i_nt = tmp_ncfout_i_nt + t_istorederfunc_m**2 - t_mult
                        tmp_ncfin_i_nt = tmp_ncfin_i_nt  + t_istorederfunc_ncorrtime**2 - t_mult
                        !      ncfout(i,nt) = ncfout(i,nt) + t_istorederfunc_m**2 - t_mult
                        !      ncfin(i,nt)  = ncfin(i,nt)  + t_istorederfunc_ncorrtime**2 - t_mult
                    end do
                    ncfout(i,nt) = tmp_ncfout_i_nt
                    ncfin(i,nt) = tmp_ncfin_i_nt
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

        END SUBROUTINE

        ! where everything gets printed
        SUBROUTINE IRCFDUMP
            real(dp) :: cfinouttot, cfouttot
            integer :: i, j, ns,imin,imax
            real(dp) :: averageCoordinationNb ! average coordination number
            real(dp) :: time, t_norm

!           allocate( cfinouttot(nspecies) )
!           allocate( cfouttot(nspecies) )
            averageCoordinationNb = dble(ncoordav*nskip)/dble(sum(posFile(:)%nbSets)*nion(2))

            open (60,file='cagecfinout.dat')
            open (61,file='cagecfout.dat')

            if(externalelectrode)then
               if(fullelectrode)then
                  if(nanalyze.eq.1)then
                     imin=1
                     imax=nion(1)
                  elseif(nanalyze.eq.2)then
                     imin=nion(1)+1
                     imax=nion(1)+nion(2)
                  elseif(nanalyze.eq.3)then
                     imin=nion(1)+nion(2)+1
                     imax=nion(1)+nion(2)+nion(3)
                  endif
               else
                  imin=1
                  imax=nelectrode
               endif
            else
               if(nanalyze.eq.2)then
                  imin=nion(1)+1
                  imax=nion(1)+nion(2)
               elseif(nanalyze.eq.3)then
                  imin=nion(1)+nion(2)+1
                  imax=nion(1)+nion(2)+nion(3)
               endif
            endif
            do j=0,ncorr-1
                cfinouttot = 0.0d0
                cfouttot = 0.0d0
                do i=imin,imax
                    t_norm = dble(norm(i,j))
                    if( t_norm/=0 ) then
                        cfouttot=cfouttot+dble(nthetaout(i,j)) / t_norm
                        cfinouttot=cfinouttot+dble(nthetainout(i,j)) / t_norm
                    end if
                end do
                time=dble(j)*dble(nskip)*dble(ndump)*dtime*2.418d-5

                if(externalelectrode)then
                   if(fullelectrode)then
                      ns=ntype(imin)
                      write(60,*)time,cfinouttot/dble(nion(ns))
                      write(61,*)time,cfouttot/dble(nion(ns))
                   else
                      write(60,*)time,cfinouttot/dble(nelectrode)
                      write(61,*)time,cfouttot/dble(nelectrode)
                   endif
                else
                   ns=ntype(imin)
                   write(60,*)time,cfinouttot/dble(nion(ns))
                   write(61,*)time,cfouttot/dble(nion(ns))
                endif
            end do

            close (60)
            close (61)

            print*, 'Average coordination no. is', averageCoordinationNb
            print*, 'cagecfout.dat cagecfinout.dat      written'
            print*, '**** Finished correctly ****'

        END SUBROUTINE

        SUBROUTINE checkInputFileExists( filename)
            character(*), intent(in) :: filename
            logical :: isOK
            inquire (file=filename, exist=isOK)
            if (.not. isOK) stop "STOP. Input file called cagecf.inpt not found."
        END SUBROUTINE
        
        SUBROUTINE openInputFile
            character(*), parameter :: inputFilename = 'cagecf.inpt' ! historical name
            open(inputUnit, file=inputFilename)
        END SUBROUTINE
        
        SUBROUTINE readInputFile
            read(inputUnit,*) nspecies ! number of species
            if( nspecies < 1 ) stop 'STOP. nspecies should be higher than 0'
            allocate( nion(nspecies) )
        
            ! read number of atoms of each species
            do i= 1, size(nion)
                read(inputUnit,*) nion(i)
                if( nion(i)<1 ) stop 'STOP. nb of ions should be >0 for all ions'
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
        
            read(inputUnit,*)nskip
            if( nskip < 1) stop 'STOP. nskip should be > 0'
            read(inputUnit,*)nanalyze
            fullelectrode=.false.
            read(inputUnit,*)externalelectrode
            if(externalelectrode)then
               read(inputUnit,'(a)') poselecfile
               read(inputUnit,*)nelectrode
               allocate( xelec(nelectrode), yelec(nelectrode), zelec(nelectrode) ) ! coordinates of electrode atoms
               read(inputUnit,'(a)') fullelectrode
            endif
            read(inputUnit,*)boxlenx
            read(inputUnit,*)boxleny
            read(inputUnit,*)boxlenz
            if( boxlenx <= 0. ) stop 'STOP. box length should be strictly positive'
            call readNumberOfPositionFiles (i)
            allocate( PosFile(i) )

            do i= 1, size(PosFile)
                read(inputUnit,'(a)') posFile(i)%name
                read(inputUnit,*) posFile(i)%nbSets
            end do

            read(inputUnit,*)ncorr ! Time-span of the correlation function (in steps)
            allocate( istorederfunc(nummax,nummax,ncorr) )
            allocate( ncfin(nummax,0:ncorr) )
            allocate( ncfout(nummax,0:ncorr) )
            allocate( nthetainout(nummax,0:ncorr) )
            allocate( nthetaout(nummax,0:ncorr) )
            allocate( norm(nummax,0:ncorr) )

            read(inputUnit,*)cutoff ! cutoff for bond
            read(inputUnit,*)ncut ! number of ions for change in cage
            read(inputUnit,*)dtime ! MD timestep
            read(inputUnit,*)ndump ! number of configs between MD dumps
        END SUBROUTINE
        
        SUBROUTINE closeInputFile
            close(inputUnit)
        END SUBROUTINE

        SUBROUTINE tellUserHowItAdvances (j)
            integer, intent(in) :: j
            integer :: maxi
            maxi = sum(posFile(:)%nbSets)
            call eraserStandardOutputLastCharacters(100)
            write(*,'(a,i10,a,i10)',advance='no')'step',j,' /',maxi
            if( j == maxi ) call eraserStandardOutputLastCharacters(100)
        END SUBROUTINE
        
        SUBROUTINE eraserStandardOutputLastCharacters (imax) ! erase the last imax characters in the standard output (terminal) in the current line
            integer, intent(in) :: imax
            integer :: i
            do i = 1, imax
                write(*,'(a)',advance='no')8
            end do
        END SUBROUTINE
        
        SUBROUTINE readNumberOfPositionFiles (nbOfPositionFiles)
            integer, intent(out) :: nbOfPositionFiles ! number of files in which positions may be found. Given in input file.
            read(inputUnit,*) nbOfPositionFiles ! number of position files
            if( nbOfPositionFiles < 1 ) stop 'STOP. nb of position files should be >=1'
        END SUBROUTINE

END PROGRAM
