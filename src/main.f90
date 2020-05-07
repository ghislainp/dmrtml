!----------------------------------------------------------------------------!
!									     !
! DMRT-ML (Dense Media Radiative Transfer - Multi-Layer)		     !
! Copyright (c), all rights reserved, 2007-2011,  Ghislain Picard	     !
! email: Ghislain.Picard@ujf-grenoble.fr				     !
!									     !
! Main contributors: Ludovic Brucker, Alexandre Roy, Florent Dupont	     !
!									     !
!									     !
! Licensed under the terms of the GNU General Public License version 3:	     !
! http://www.opensource.org/licenses/gpl-3.0.html 			     !
!									     !
! dmrtml'url: http://lgge.osug.fr/~picard/dmrtml/                            !
!                                                                            !
! Recommended citation:							     !
!    G. Picard, L. Brucker, A. Roy, F. Dupont, M. Fily, and A. Royer,        !
!    Simulation of the microwave emission of multi-layered snowpacks using   !
!    the Dense Media Radiative Transfer theory: the DMRT-ML model,           !
!    Geoscientific Model Development, 6, 1061-1078, 2013                     !
!    doi:10.5194/gmd-6-1061-2013                                             !
!    http://www.geosci-model-dev.net/6/1061/2013/gmd-6-1061-2013.html        !
!									     !
!----------------------------------------------------------------------------!

program dmrt_prog
  use mod_dmrtml
  use mod_options
  use mod_soil

  implicit none

  !------------------------------------------------------------
  integer,parameter        :: lmax=200         ! maximum number of layer

  real*8,parameter         :: pi=3.14159265

  integer                  :: n                 ! number of stream
  integer                  :: l                 ! number of layer
  integer                  :: i,ith0,ifreq
  integer,dimension(1)     :: ith0_arr

  ! Input paramaters
  real*8,dimension(lmax)   :: temp,height,radius,density

  ! Command line
  Type(options)            :: opt
  !type(soilparams)         :: soilp             ! soil parameters

  real*8,dimension(:),allocatable              :: tau
  ! outputs
  real*8,dimension(:),allocatable              :: TbV,TbH,outmhu

  !------------------------------------------------------------
  real*8                 :: freq
  integer                :: ierr

  character(len=255)     :: str
  !-------------------------------------------------------------
  ! read from the command line

  call readoptions(opt)

  !-------------------------------------------------------------
  ! read the input file

  open(unit=9,file=opt%infilename)
  l=1
  do while (l<=lmax)
     read(9,'(A)',end=999) str
     if (str(1:1)=='#') cycle
     read(str,*,IOSTAT=ierr) height(l),density(l),radius(l),temp(l)
     if (ierr/=0) then
        print *,'Error in '// trim(opt%infilename)// ' on line '//trim(str)
        stop
     endif
     l=l+1
  enddo
  print *,'to many layer, increase lmax in the source code'
  stop


999  close(9)

  l=l-1
  if (opt%verbose) then
     write (*,'(A,I10)') '#number of layers: ',l
     write (*,'(A,F10.2)') '#total height (m): ',sum(height(1:l))
     write (*,'(A,F10.2)') '#mean temperature (K)',sum(temp(1:l))/dble(l)
     write (*,'(A,F10.2)') '#mean density (kg/m3)',sum(density(1:l))/dble(l)
     write (*,'(A,F10.2)') '#mean grain radius (microns)',sum(radius(1:l))/dble(l)
  endif

  ! convert radius to meter
  radius=radius*1e-6

  n=opt%nstream
  allocate(outmhu(n))
  allocate(TbV(n))
  allocate(TbH(n))

  allocate(tau(l))
  tau(:)=opt%tau

  ! loop over frequencies

  if (opt%verbose) write (*,'(A)') '# freq trueinc TbV TbH'

  do ifreq=1,opt%ifreq
     freq=opt%freqlist(ifreq)*1e9

     call dmrtml(freq,l,n,height,density,radius,temp,TbV,TbH,outmhu,tau=tau)

     ! outputs
     if (opt%iinc>0) then
        do i=1,opt%iinc
           ith0_arr=minloc(abs(outmhu(:)-cos(opt%inclist(i)*pi/180))) ! search for the closest angle
           ith0=ith0_arr(1)
           write (*,'(F6.1,F5.1,2F6.1)') freq/1e9,180/pi*acos(outmhu(ith0)),TbV(ith0),TbH(ith0)
        enddo
     else
        i=1
        do while (outmhu(i)>0)
           write (*,'(F6.1,1X,F5.1,1X,2F7.1)') freq/1e9,180/pi*acos(outmhu(i)),TbV(i),TbH(i)
           i=i+1
        enddo
     endif
  enddo ! frequencies

end program dmrt_prog
