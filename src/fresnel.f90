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

module mod_fresnel


contains

!!-----------------------------------------------------------------------------!
!!-----------------------------------------------------------------------------!
!subroutine fresnel_coefficients(E,theta0,rv,rh)
!  ! calculate the Fresnel coefficient for a single incidence angle
!  implicit none
!  real*8,intent(in)       :: theta0
!  complex*16,intent(in)    :: E
!  real*8,intent(out)      :: rv,rh
!  !----------------------------------
!  real*8                  :: n,theta1
!  stop "incorrecte... a reecrire"
!
!  n=real(sqrt(E))
!  theta1=asin(sin(theta0)/n)
!
!  rv=abs( (n*cos(theta0)-cos(theta1))/(n*cos(theta0)+cos(theta1)) )**2
!  rh=abs( (cos(theta0)-n*cos(theta1))/(cos(theta0)+n*cos(theta1)) )**2
!
!end subroutine fresnel_coefficients




!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine fresnel_coefficients_multiangle(E0,E1,costheta0,rv,rh)
  ! calculate the Fresnel coefficient for an array of incidence angle
  ! in medium 0 with transmission into the medium 1
  implicit none
  real*8,intent(in),dimension(:)       :: costheta0
  complex*16,intent(in)                 :: E0,E1
  real*8,intent(out),dimension(:)      :: rv,rh
  !----------------------------------
  complex*16                            :: costheta1
  complex*16                            :: n,b
  integer                              :: i
  
!

  n=sqrt(E1/E0)

  do i=1,size(costheta0)
     if (real(costheta0(i)**2 - (1.0-n**2)) <=0) then

        rv(i)=1
        rh(i)=1
     else
        b = 1.0-(1.0-costheta0(i)**2)/(n**2)
        costheta1=sqrt(b)
        rv(i)=abs( (n*costheta0(i)-costheta1)/(n*costheta0(i)+costheta1) )**2
        rh(i)=abs( (costheta0(i)-n*costheta1)/(costheta0(i)+n*costheta1) )**2

    end if
  end do

end subroutine fresnel_coefficients_multiangle

end module mod_fresnel
