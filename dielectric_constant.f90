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

module mod_dielectric_constant

contains
        
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine icedielectric(frequency,temperature,Ei)
  ! calculate the complex ice dielectric constante depending on frequency and 
  ! temperature.
  ! code from MEMLS
  !
  ! Inputs:
  ! frequency (Hz)
  ! temperature (K)
  ! Outputs:
  ! Ew dielectric constant of ice
  !
  implicit none
  real*8,intent(in)       :: frequency,temperature 
  complex*16,intent(out)   :: Ei

  !----------------------------
  real*8                  :: freqGHz,theta, alpha, beta
  real*8                  :: B1,B2,b,deltabeta,betam
  real*8                  :: Ereal, Eimag


  freqGHz=frequency/1e9

  Ereal=3.1884 + 9.1e-4 * (temperature-273.0)

  theta=300.0/temperature - 1.0
  alpha=(0.00504+0.0062*theta) * exp(-22.1*theta)

  B1 = 0.0207
  B2 = 1.16e-11
  b = 335
  deltabeta = exp(-9.963 + 0.0372 * (temperature-273.16))
  betam = (B1/temperature) * ( exp(b/temperature)/ ((exp(b/temperature)-1)**2) ) + B2*freqGHz**2
  beta = betam + deltabeta

  Eimag=alpha / freqGHz + beta * freqGHz

  Ei=cmplx(Ereal,Eimag)

end subroutine icedielectric


!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine waterdielectric (frequency,temperature,Ew)
  ! calculate the complex water dielectric constante depending on frequency 
  ! and temperature
  ! code from MEMLS

  ! Inputs:
  ! frequency (Hz)
  ! temperature (K)
  ! Outputs:
  ! Ew dielectric constant of water
  !
  implicit none
  real*8,intent(in)       :: frequency,temperature
  complex*16,intent(out)   :: Ew

  !----------------------------
  real*8                  :: freqGHz,theta
  real*8                  :: e0,e1,e2,f1,f2


  freqGHz=frequency/1e9

  theta=1 - 300.0/temperature

  e0=77.66-103.3*theta
  e1=0.0671*e0

  f1=20.2+146.4*theta+316*theta**2
  e2=3.52+7.52*theta
  !!% version of Liebe MPM 1993 uses: e2=3.52 
  f2=39.8*f1

  !print *,'weps',e0,e1,e2,f1,f2
  Ew=e2 + (e1-e2)/cmplx(1,-freqGHz/f2) + (e0-e1)/cmplx(1,-freqGHz/f1)

end subroutine waterdielectric


!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine weticedielectric (frequency,temperature,fwetness,Ewi)
  ! calculate the dielectric constant of wet particule of ice
  ! using Bohren and Huffman 1983 according to Ya Qi Jin, eq 8-69, 1996 p282
  !
  ! Inputs:
  ! frequency (Hz)
  ! temperature (K)
  ! wetness (fractional volume of water with respect to ice)
  ! Outputs:
  ! Ewi dielectric constant of the ice sphere coated with water
  !
  implicit none
  real*8,intent(in)       :: frequency,temperature,fwetness
  complex*16               :: Ewi

  !----------------------------
  complex*16               :: epsice,epswater,Cminus,Cplus
  real*8                  :: S
  ! from http://books.google.com/books?id=Y_-113zIvgkC&pg=PA142&lpg=PA142&dq=effective+dielectric+constant+sphere+coated+small&source=bl&ots=ZVfwvkA0K1&sig=P7fHb0Jff8C-7-GrlEnWRZkkxY8&hl=en&ei=RHfDTrmjJYXj8AO3v7ScCw&sa=X&oi=book_result&ct=result&resnum=3&ved=0CDYQ6AEwAg#v=onepage&q=effective%20dielectric%20constant%20sphere%20coated%20small&f=false

  ! see also: K L CHOPRA and G B REDDY, Praman.a- Optically selective coatings, J. Phys., Vol. 27, Nos 1 & 2, July & August 1986, pp. 193-217.

  call icedielectric(frequency,temperature,epsice)
  if (fwetness<=0.0) then
     Ewi=epsice
     return
  endif

  call waterdielectric(frequency,temperature,epswater)

  S = 1-fwetness

  Cplus=epsice+2*epswater
  Cminus=(epsice-epswater)*S

  !!alpha = ( (epswater-eps0)*Cplus + Cminus*(eps0+2*epswater) ) /  &
  !!     ( (epswater+2*eps0)*Cplus + 2*(epswater-eps0)*Cminus )
  !! and
  !! alpha= (epseff-eps0) / (epseff + 2*eps0)
  !! after some calculation:

  Ewi = (Cplus + 2*Cminus) / (Cplus - Cminus) *epswater

  !print *,'Ewi=',Ewi

end subroutine weticedielectric


end module mod_dielectric_constant
