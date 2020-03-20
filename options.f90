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

Module mod_options

  integer,parameter     :: fmax=100 ! maximum number of frequency

  Type options
     character(len=1024)    :: infilename
     real*8,dimension(fmax) :: freqlist
     integer                :: ifreq
     real*8,dimension(fmax) :: inclist
     integer                :: iinc
     integer                :: nstream
     logical                :: verbose
     real*8                 :: tau ! stickiness
  end Type options


contains

subroutine readoptions(opt)

  Type(options),intent(out)  :: opt


  character(len=255)         :: str
  integer                    :: iparam,iargc
  logical                    :: e
  !------------------------------

  opt%nstream=64
  opt%infilename=''
  opt%iinc=0
  opt%ifreq=0
  opt%verbose=.TRUE.
  opt%tau=1e6 ! Nonsticky

  iparam=1
  do while (iargc()>=iparam) 
     call getarg(iparam,str)
     select case(str)
     case ('-n')
        iparam=iparam + 1
        call getarg(iparam,str)
        !read (str,'(I)') opt%nstream
        read (str,*) opt%nstream

     case ('-i','--input')
        iparam=iparam + 1
        call getarg(iparam,opt%infilename)
        inquire(file=opt%infilename,exist=e)
        if (.not. e) then
           print *,'input filename ' // trim(opt%infilename) // ' doesn''t exist'
           stop
        endif
     case ('-f','--freq')

        call readlist(iparam,opt%freqlist,opt%ifreq)
        if (opt%ifreq<0) then
           print *,'error parsing option --freq'
           stop
        endif
     case ('-t','--stickiness')
        iparam=iparam + 1
        call getarg(iparam,str)
        read (str,*) opt%tau
        if (opt%tau<0) then
           print *,'error parsing option --tau'
           stop
        endif

     case ('-inc','--inc')

        call readlist(iparam,opt%inclist,opt%iinc)
        if (opt%ifreq<0) then
           print *,'error parsing option --inc'
           stop
        endif
        
     case ('-v0','--silent')
        opt%verbose=.FALSE.

     case('-s','--sensor')
        iparam=iparam+1
        call getarg(iparam,str)
        select case(str)
        case ('ssmi')
           opt%freqlist(1)=19.3
           opt%freqlist(2)=37.0
           opt%freqlist(3)=85.5
           opt%ifreq=3
           opt%inclist(1)=53.1
           opt%iinc=1
        case('amsre','amsr')
           opt%freqlist(1)=6.925
           opt%freqlist(2)=10.65
           opt%freqlist(3)=18.7
           opt%freqlist(4)=36.5
           opt%freqlist(5)=89.0
           opt%ifreq=5
           opt%inclist(1)=55
           opt%iinc=1
        case default
           print *,'Unknown sensor: ',str
           stop
        end select
     case default
        print *,'Unknown argument: ',str
        stop
     end select
     iparam=iparam+1
  enddo

  if (opt%infilename==''.or.opt%ifreq==0) then
     print *,'Usage:'
     print *,'dmrt -i mediumfile -f frequencies [-n] [-inc incidence angles] [-s sensor] [-v0]'
     print *,'-i,--input    input file describing the medium'
     print *,'-f,--freq     list or range of frequencies in (GHz)'
     print *,'-n            number of stream (default 64)'
     print *,'-s,--sensor   sensor name to automatically fill frequencies and incidence angle'
     print *,'-v0,--silent  no comment output'
     stop
  endif

end subroutine readoptions



subroutine readlist(iparam,list,size)
  implicit none

  integer,intent(inout)                 :: iparam
  real*8,intent(out),dimension(fmax)      :: list
  integer,intent(out)                   :: size

  !---------------------------------
  character(len=255)         :: str
  integer                    :: i,iargc
  logical                    :: range
  real*8                       :: x,x1,step

  range=.false.
  iparam=iparam+1
  call getarg(iparam,str)
  i=0
  do while (isnumeric(str).or.str==':')
     if (str==':') then
        range=.true.
     else
        i=i+1
        read(str,*) list(i)
        if (i>fmax) then
           print *,'to many value'
           size=-1
           return
        endif
     endif
     iparam=iparam+1
     if (iparam>iargc()) exit
     call getarg(iparam,str)
  enddo


  iparam=iparam-1

  if (range) then
     if (i<2.or.i>3) then
        print *,'invalid range'
        size=-1
        return
     endif
     x=min(list(1),list(2))
     x1=max(list(1),list(2))
     if (i==3) then
        step=abs(list(3))
     else
        step=1
     endif
     i=0
     do while (x<=x1)
        i=i+1
        list(i)=x
        x=x+step
     enddo
  end if
  size=i
end subroutine readlist



FUNCTION isnumeric(string)
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: string
  LOGICAL :: isnumeric
  real*8 :: x
  INTEGER :: e
  READ(string,*,IOSTAT=e) x
  isnumeric = e == 0
END FUNCTION isnumeric


  !function isnumeric(str)
  !  implicit none
  !  character(len=*),intent(in)          :: str
  !  logical                              :: isnumeric
  !  !--------------------
  !  real :: x
  !  integer :: stat
  !
  !  read(str,'(E)',iostat=stat) x
  !
  !  isnumeric=(stat==0)
  !
  !end function isnumeric

end Module mod_options
