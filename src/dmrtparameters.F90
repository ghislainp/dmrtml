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

!> Compute dense media radiative transfer parameters (extinction and single scattering albedo) 
!! under various assumptions (monodispese or polydisperse, Grody's limitation, ...)
!!

module mod_dmrtparameters

#ifdef MKL
  use mkl95_lapack
#else
#ifdef SUNPERF
  use sunperf
#endif
#endif
  implicit none
  real*8,parameter :: NONESTICKY=1000

contains


!> Compute radiative transfer parameters (extinction and single scattering albedo) using 
!! DMRT theory under the approximations: QCA-CP, small particules, 
!! sticky short range, single species. In addition, an empirical correction for large particules 
!! is applied based on  Grody, 2008. This correction was further improved by G. Picard, A. Roy and F. Dupont but it is still far from perfect. 
!! Its main purpose is to prevent DMRT divergence for large particules when ONLY a few layers contain coarse grains.
!! However, it should not be considered sufficiently accurate to investigate the scattering behavior of large particules.
!! To disable this ad-hoc correction use the subroutine \ref dmrtparameters instead of this one.
!!
!! REFERENCES:
!!
!! Gordy, N. Relationship between snow parameters and microwave satellite measurements: Theory compared with Advanced Microwave Sounding Unit observations from 23 to 150 GHz. J. Geophys. Res. 113, 22108–+ (2008).
!!
!! @param[in] frequency frequency in Hertz
!! @param[in] temperature temperature in Kelvin
!! @param[in] f fractional volume
!! @param[in] a radius in meter
!! @param[in] tau stickiness
!! @param[in] fwetness liquid water content in m3/m3
!! @param[in] medium type of medium (valid value is S=Snow or I=ice with air bubbles)
!! @param[out] Eeff0 effective dielectric constant, order 0 (=no scattering)
!! @param[out] Eeff effective dielectric constant
!! @param[out] albedo single scattering albedo
!! @param[out] beta extinction coefficient (also known as Ke)
!! @param[in] Eice optional dielectric constant of ice (don't use it)
!! @todo add W=Water in the list of supported medium
!!

subroutine dmrtparameters_grodyapproach(frequency,temperature,f,a,tau,fwetness,medium,&
     Eeff0,Eeff,albedo,beta,Eice)


  implicit none

  real*8,intent(in)       :: frequency,temperature
  real*8,intent(in)       :: f,a, fwetness,tau ! fractional volume, radius, fwetness and stickiness
  character,intent(in)    :: medium ! ice sphere (.FALSE.) or air bubble (.TRUE.)
  complex*16,intent(out)   :: Eeff0,Eeff
  real*8,intent(out)      :: albedo,beta
  complex*16,intent(in),optional :: Eice
  !-----------------------------
  real*8 :: aol,diff,lastdiff,fb_dyn,Qe,Qs!,kr !aol = a for optic limit
  real*8,parameter :: pi = 3.14159265

  real*8,parameter :: optic_limit = 2.0 ! parameter for ad hoc correction based on Grody 2008
  real*8,parameter :: res = 1e-3 ! res = resolution in the internal optimization loops 

  call dmrtparameters(frequency,temperature,f,a,tau,fwetness,medium,&
       Eeff0,Eeff,albedo,beta,Eice)


  ! compute extinction and scattering efficiency
  Qe=(beta*a*4) / (3*f)
  Qs=albedo*Qe

  if (Qs>optic_limit) then ! if True, the particule size is too large.
    
     ! decrease particle size until Qs==2
    aol = a
    diff = Qs-optic_limit
    lastdiff=diff

    fb_dyn = aol / 2

    do while (abs(diff)>res)

      if (diff<0) then
	aol = aol + fb_dyn
      else
	aol = aol - fb_dyn
      endif
      ! decrease the step only if around the optic_limit
      if (diff*lastdiff<0) then
         fb_dyn = fb_dyn/2
      endif
      lastdiff=diff

      call dmrtparameters(frequency,temperature,f,aol,tau,fwetness,medium,&
           Eeff0,Eeff,albedo,beta,Eice)

      beta=beta*aol/a
      Qe=(beta*a*4) / (3*f)
      Qs=albedo*Qe

      diff = Qs-optic_limit
    end do
    !the last  dmrtparameters call provides the value to be return.
  endif
end subroutine dmrtparameters_grodyapproach


function sticky_t_parameter(f,tau)
  ! calculate the t parameter: equation 8.4.22  Tsang vol II
  implicit none
  real*8,intent(in) :: f,tau
  real*8            :: sticky_t_parameter
  !--------------------------------------
  real*8            :: a,b,c,discr,mhu,mhulim

  ! calculate t
  if (tau<(2-sqrt(2.0))/6) then
     print *,"tau must be higher than ",(2-sqrt(2.0))/6
     stop
  endif
  ! solve equation 8.4.22  vol II
  a=f/12
  b=-(tau+f/(1-f))
  c=(1+f/2)/(1-f)**2

  discr=sqrt(b**2-4*a*c)

  sticky_t_parameter=(-b-discr)/(2*a)

  ! check mhu<1+2f
  mhu=sticky_t_parameter*f*(1-f)
  mhulim=1+2*f
    
  if (mhu>mhulim) then
     print *,"largest"
     sticky_t_parameter=(-b+discr)/(2*a)
     mhu=sticky_t_parameter*f*(1-f)

     if (mhu>mhulim) then
        print *,"no tau solution"
        stop
     endif
  else
     print *,'smallest'
  endif
end function sticky_t_parameter


!> Compute radiative transfer parameters (extinction and single scattering albedo) using DMRT theory
!! under the assumption: QCA-CP, small particules, sticky short range, single species
!!
!! Warning: this only works for short range, 
!! i.e. moderate stickiness / small grains (high tau value)!!
!!
!! @param[in] frequency frequency in Hertz
!! @param[in] temperature temperature in Kelvin
!! @param[in] f fractional volume
!! @param[in] a radius in meter
!! @param[in] tau stickiness
!! @param[in] fwetness liquid water content in m3/m3
!! @param[in] medium type of medium (valid value is S=Snow or I=ice with air bubbles)
!! @param[out] Eeff0 effective dielectric constant, order 0 (=no scattering)
!! @param[out] Eeff effective dielectric constant
!! @param[out] albedo single scattering albedo
!! @param[out] beta extinction coefficient (also known as Ke)
!! @param[in] Eice optional dielectric constant of ice (don't use it)
!! @todo add W=Water in the list of supported medium
subroutine dmrtparameters(frequency,temperature,f,a,tau,fwetness,medium, &
     Eeff0,Eeff,albedo,beta,Eice)
!
  use mod_dielectric_constant
!
!
  implicit none
  real*8,intent(in)       :: frequency,temperature ! Frequency in Hz, Temperature in K
  real*8,intent(in)       :: f,a,tau,fwetness  ! fractional volume, radius, stickiness and fwetness
  character,intent(in)    :: medium ! ice sphere 'S' or air bubble 'I'
  complex*16,intent(out)  :: Eeff0,Eeff
  real*8,intent(out)      :: albedo,beta
  complex*16,intent(in),optional   :: Eice

  !-----------------------------
  ! local variables
  complex*16               :: Ei,es,e0
  complex*16               :: b,c,discriminant
  real*8                   :: lambda,t, snow_moisture
  real*8,parameter         :: pi=3.14159265
  complex*16,parameter     :: cone=cmplx(1.0,0.0)

  !-----------------------------

  if (PRESENT(Eice).and.real(Eice)>=1.0) then
     Ei=Eice
  else
     call weticedielectric(frequency,temperature, fwetness, Ei)
  endif

  if (medium=='S') then ! SNOW
     ! ice spheres in air
     es=Ei
     e0=cone
  else if (medium=='I') then ! ICE
     ! air spheres in ice
     es=cone
     e0=Ei
  else if (medium=='B') then ! Bridging
     call dmrtparameters_bridge(frequency,temperature,f,a,tau,fwetness,Eeff0,Eeff,albedo,beta,Ei)
     return
  end if
  
  !Solve the 0th-order solution: Eeff0
  ! Eeff0^2 + Eeff0 *[ (Ei-eo)/3*(1-4f)-eo] - eo (Ei-1)/3*(1-f) = 0

  b=(es-e0)*(1.0-4.0*f)/3.0 - e0
  c=-e0*(es-e0)*(1.0-f)/3.0

  discriminant= b**2 - 4*c

  !solution

  Eeff0 = 0.5*(-b+sqrt(discriminant))

  if (real(Eeff0)<1) then
     Eeff0 = 0.5*(-b-sqrt(discriminant))
     print *,'attention.... par la!!'
  endif

  if (tau>=NONESTICKY) then
     t=0
  else
     t=sticky_t_parameter(f,tau)
  endif


  !Solve 1st-order solution: E
  lambda=3e8/frequency

  Eeff= e0 + (Eeff0 - e0) * ( cone + cmplx(0.0,2.0/9.0) * &
   (2*pi*a/lambda)**3 * sqrt(Eeff0)*(es-e0) / (cone + (es-e0)/(3*Eeff0)*(1.0-f) ) * &
   (1.0-f)**4 /(1.0+2*f-t*f*(1.0-f))**2 )
    
  albedo=2.0/9.0 * (2*pi*a/lambda)**3 * f / (2*aimag(sqrt(Eeff))) *  &
       abs( (es-e0)/(cone+(es-e0)/(3*Eeff0)*(1.0-f) ) )**2 * &
       (1.0-f)**4 /(1.0+2*f-t*f*(1.0-f))**2   
    
  beta=2*pi/lambda * 2*aimag(sqrt(Eeff))

end subroutine dmrtparameters




!> Compute radiative transfer parameters (extinction and single scattering albedo) using DMRT theory
!! under the assumption: QCA-CP, small particules, sticky short range, single species
!! apply the bridging correction proposed by Dierking et al., 2012
!!
!! Warning: this only works for short range, 
!! i.e. moderate stickiness / small grains (high tau value)!!
!!
!! @param[in] frequency frequency in Hertz
!! @param[in] temperature temperature in Kelvin
!! @param[in] f fractional volume
!! @param[in] a radius in meter
!! @param[in] tau stickiness
!! @param[in] fwetness liquid water content in m3/m3
!! @param[out] Eeff0 effective dielectric constant, order 0 (=no scattering)
!! @param[out] Eeff effective dielectric constant
!! @param[out] albedo single scattering albedo
!! @param[out] beta extinction coefficient (also known as Ke)
!! @param[in] Eice optional dielectric constant of ice (don't use it)
!! @todo add W=Water in the list of supported medium

subroutine dmrtparameters_bridge(frequency,temperature,f,a,tau,fwetness, &
     Eeff0,Eeff,albedo,beta,Eice)

  use mod_dielectric_constant
  implicit none
  real*8,intent(in)       :: frequency,temperature ! Frequency in Hz, Temperature in K
  real*8,intent(in)       :: f,a,tau,fwetness  ! fractional volume, radius, stickiness and fwetness
  complex*16,intent(out)   :: Eeff0,Eeff
  real*8,intent(out)      :: albedo,beta
  complex*16,intent(in),optional   :: Eice

  !--------------------------------
  complex*16              :: Ei
  real*8,dimension(4)     :: Ks,Ka
  complex*16,dimension(4) :: Eeff0_pt,Eeff_pt
  real*8,dimension(3)     :: b,d
  real*8,dimension(4)     :: c,l,mu,z,fpt
  real*8                  :: Ks_bridge

  integer :: i


  if (PRESENT(Eice).and.real(Eice)>=1.0) then
     Ei=Eice
  else
     call weticedielectric(frequency,temperature,fwetness,Ei)
  endif

  ! tie-point for the spline
  !fpt(1)=0.2
  !fpt(2)=0.3
  !fpt(3)=0.7
  !fpt(4)=0.8

  fpt(1)=0.3
  fpt(2)=0.4
  fpt(3)=0.6
  fpt(4)=0.7
  
  if (f<=fpt(2)) then
     call dmrtparameters(frequency,temperature,f,a,tau,fwetness,'S',Eeff0,Eeff,albedo,beta,Ei)

  else if (f>=fpt(3)) then
     call dmrtparameters(frequency,temperature,1-f,a,tau,fwetness,'I',Eeff0,Eeff,albedo,beta,Ei)

  else

     do i=1,4
        if (fpt(i)<=0.5) then
           call dmrtparameters(frequency,temperature,fpt(i),a,tau,fwetness,'S',Eeff0_pt(i),Eeff_pt(i),albedo,beta,Ei)
        else
           call dmrtparameters(frequency,temperature,1-fpt(i),a,tau,fwetness,'I',Eeff0_pt(i),Eeff_pt(i),albedo,beta,Ei)
        endif
        Ks(i)=albedo*beta
        Ka(i)=(1.0-albedo)*beta
     enddo

     Ks_bridge = cubicspline4(Ks,fpt,f)
     beta = cubicspline4(Ka,fpt,f) + Ks_bridge
  
     albedo = Ks_bridge/beta

     Eeff0 = cmplx(cubicspline4(real(Eeff0_pt),fpt,f),cubicspline4(imag(Eeff0_pt),fpt,f))
     Eeff = cmplx(cubicspline4(real(Eeff_pt),fpt,f),cubicspline4(imag(Eeff_pt),fpt,f))
  end if

  return

end subroutine dmrtparameters_bridge


function cubicspline4(ys,xs,x)

  real*8 :: cubicspline4
  real*8,intent(in)       :: x
  real*8,intent(in),dimension(4) :: ys,xs

  real*8,dimension(3) :: h
  real*8,dimension(4) :: c,mu,z
  real*8,dimension(2) :: alpha
  real*8              :: l,d,b
  integer i,j

  do i=1,3
     h(i)=xs(i+1)-xs(i)
  end do

  do i=1,2
     alpha(i)=3*( (ys(i+2)-ys(i+1))/h(i+1) - (ys(i+1)-ys(i))/h(i) )
  end do

  mu(1)=0
  z(1)=0
  c(4)=0

  do i=2,3
     l = 2*(xs(i+1)-xs(i-1))-h(i-1)*mu(i-1)
     mu(i) = h(i)/l
     z(i) = (alpha(i-1)-h(i-1)*z(i-1))/l
  end do

  do j=3,1,-1
     c(j) = z(j) - mu(j)*c(j+1)
     b = (ys(j+1) -ys(j))/h(j) - h(j)*(c(j+1)+2*c(j))/3
     d = (c(j+1)-c(j))/(3*h(j))

     if (x>=xs(j)) then
        cubicspline4 = ys(j)+b*(x-xs(j))+c(j)*(x-xs(j))**2+d*(x-xs(j))**3
        return
     end if
  end do
end function cubicspline4


!------------------------------
function conjg(c)
  complex*16,intent(in) :: c
  complex*16 :: conjg

  conjg=cmplx(real(c),-aimag(c))
end function conjg


subroutine computeHmultispecies(n,a,voldens,H)
  integer,intent(in)             :: n
  real*8,dimension(n),intent(in)   :: a,voldens
  real*8,dimension(n,n),intent(out):: H
  !----------
  real*8,dimension(n,n)            :: C,I_C
  real*8                           :: Ri,Rj,Mi,Mj,Ni,Nj,one_zeta3
  real*8,dimension(0:3)            :: zeta
  integer,dimension(n)           :: ipiv
  integer :: i,j
  real*8,parameter        :: pi=3.14159265

#ifndef MKL
  integer :: info
  real*8,dimension(5*n) :: work
#endif
  !----------

  zeta=0.0
  I_C=0.0

  do i=1,n
     !print *,'voldens=',voldens(i)
     zeta(0)=zeta(0)+voldens(i)
     zeta(1)=zeta(1)+voldens(i)*(2*a(i))
     zeta(2)=zeta(2)+voldens(i)*(2*a(i))**2
     zeta(3)=zeta(3)+voldens(i)*(2*a(i))**3

     I_C(i,i)=1.0
  enddo
  zeta=zeta*(pi/6.0)
  one_zeta3=1.0-zeta(3)

  !print *,'zeta0=',zeta(0)

  do i=1,n
     Ri=2*a(i)
     Mi=Ri**3
     Ni=Ri**2

     do j=1,n
        Rj=2*a(j)
        Mj=Rj**3
        Nj=Rj**2

        C(i,j)=-pi/6.0*sqrt(voldens(i)*voldens(j))/one_zeta3 * (                   &
              Mj * (1+ 3*(zeta(2)*Ri + zeta(1)*Ni)/one_zeta3 + 9*zeta(2)**2*Ni/one_zeta3**2) &
             +Mi * (1+ 3*(zeta(2)*Rj + zeta(1)*Nj)/one_zeta3 + 9*zeta(2)**2*Nj/one_zeta3**2) &
             +Mi*Mj * (zeta(0)/one_zeta3 + 6*zeta(1)*zeta(2)/one_zeta3**2 + 9*zeta(2)**3/one_zeta3**3) &
             +3*Ni*Rj + 3*Nj*Ri + 9*zeta(2)*Ni*Nj/one_zeta3   ) 
     enddo
  enddo

  I_C=I_C-C

#ifdef MKL
  call getrf(I_C,ipiv)
  call getri(I_C,ipiv)
#else
#ifdef SUNPERF
#error "a implementer, inversion pour SUNPERF"
#else
  call dgetrf(n,n,I_C,n,ipiv,info)
  call dgetri(n,I_C,n,ipiv,work,size(work),info)
#endif
#endif


  H=matmul(I_C,C)

  do i=1,n
     do j=1,n
        H(i,j)=H(i,j)/sqrt(voldens(i)*voldens(j))/(2*pi)**3
     enddo
  enddo
end subroutine computeHmultispecies

!
! don't converge when used without Grody-like correction.
!
!subroutine computeLogNormalDistribution(a,sigma,n,aa,number)
!  real,intent(in)              :: a,sigma
!  integer,intent(in)           :: n
!  real,intent(out),dimension(n) :: aa,number
!  !----------------------
!  real                         :: rm,amin,amax
!  integer                      :: i
!
!  ! a est le rayon optique et rm est le rayon moyen
!  rm=a*exp(-5.0/2.0*log(sigma)**2)   !(attentio, erreur dans le papier de Charlie)
!
!  amin=0.1*rm
!  amax=30.0*rm
!  do i=1,n
!     aa(i)=amin+i*(amax-amin)/n
!     number(i)=exp(-0.5*(log(aa(i)/rm)/log(sigma))**2)/aa(i)
!  enddo
!
!  print *,'#meanaa=',sum(number*aa)/sum(number),rm
!  print *,'#amax=',amax
!
!  if (.False.) then
!     print *,'#pdf'
!     do i=1,n
!        print *,aa(i)*1e6,number(i)
!     enddo
!     stop
!  endif
!  if (.False.) then
!     print *,'#pdf aa3'
!     do i=1,n
!        print *,aa(i)*1e6,aa(i)**3*number(i)
!     enddo
!     stop
!  endif
!end subroutine computeLogNormalDistribution


subroutine computeRayleighDistribution(a,n,aa,number)
  real*8,intent(in)               :: a
  integer,intent(in)              :: n
  real*8,intent(out),dimension(n) :: aa,number
  !-------------------------------------------------
  real*8                          :: rm,amin,amax,sigma
  integer                         :: i
  real*8,parameter                :: pi=3.14159265
  !-------------------------------------------------

  ! a is the optical radius and rm is the mean radius
  ! rm=a*exp(-5.0/2.0*log(sigma)**2)   !(warning error in Charlie's paper)
  sigma=a/(sqrt(2.0)*1.329340) !! calculation done by Ghislain
  rm=sigma*sqrt(pi/2)

  amin=0.1*rm
  amax=10.0*rm
  do i=1,n
     aa(i)=amin+i*(amax-amin)/n
     number(i)=exp(-0.5*(aa(i)/sigma)**2)*aa(i)/sigma**2
  enddo

  if (.False.) then
     print *,'#pdf'
     do i=1,n
        print *,aa(i)*1e6,number(i)
     enddo
     stop
  endif
  if (.False.) then
     print *,'#pdf aa6'
     do i=1,n
        print *,aa(i)*1e6,aa(i)**6*number(i)
     enddo
     stop
  endif
end subroutine computeRayleighDistribution



!> Compute radiative transfer parameters (extinction and single scattering albedo) 
!! using DMRT theory under the assumptions: QCA-CP, small particules, *none-sticky*, multispecies.
!! Grain size repartition is given by a Rayleigh distribution
!!
!! @param[in] frequency frequency in Herz.
!! @param[in] temperature temperature in Kelvin.
!! @param[in] f fractional volume.
!! @param[in] a radius in meter.
!! @param[in] tau stickiness (must be zero).
!! @param[in] fwetness liquid water content in m3/m3.
!! @param[in] medium type of medium (valid value is S=Snow or I=ice with air bubbles).
!! @param[out] Eeff0 effective dielectric constant, order 0 (=no scattering).
!! @param[out] Eeff effective dielectric constant.
!! @param[out] albedo single scattering albedo.
!! @param[out] beta extinction coefficient (also known as Ke).
!! @param[in] Eice optional dielectric constant of ice (don't use it).
!! @todo implement sticky multispecies... not easy.

subroutine dmrtparameters_dist(frequency,temperature,f,a,tau,fwetness,medium,&
     Eeff0,Eeff,albedo,beta,Eice)

  use mod_dielectric_constant

  implicit none
  real*8,intent(in)          :: frequency,temperature
  real*8,intent(in)          :: f,a,tau,fwetness   ! fractional volume, radius, stickiness and fwetness
  character,intent(in)       :: medium ! ice sphere (.FALSE.) or air bubble (.TRUE.)
  complex*16,intent(out)     :: Eeff0,Eeff
  real*8,intent(out)         :: albedo,beta
  complex*16,intent(in),optional   :: Eice

  !-----------------------------
  complex*16              :: Ei,es,e0
  complex*16              :: b,c,discriminant
  real*8                  :: lambda,k0
  real*8,parameter        :: pi=3.14159265
  complex*16,parameter    :: cone=cmplx(1.0,0.0)

  integer,parameter       :: ll=100
  real*8,dimension(ll,ll) :: H
  real*8,dimension(ll)    :: aa,ff,voldens,number

  integer                 :: l,j
  real*8                  :: suml2,norm
  complex*16              :: sumj,suml1,D,y
  !-----------------------------

  if (tau<NONESTICKY) then
     print *,'Mutlispecies none sticky version is not implemented yet'
     print *,'tau=',tau
     stop
  endif

  if (PRESENT(Eice).and.real(Eice)>=1.0) then
     Ei=Eice
  else
     ! get the ice dielectric constant
     call weticedielectric(frequency,temperature,fwetness,Ei)
     !!print *,'#Ei=',Ei
     !!print *,'#sqrt(Ei)',sqrt(EI)
     !Ei=cmplx(3.2,0.016)
  endif

  if (medium=='S') then ! SNOW
     ! ice spheres in air
     es=Ei
     e0=cone
  else ! ICE
     ! air spheres in ice
     es=cone
     e0=Ei
  end if
  
  !Solve the 0th-order solution: Eeff0
  !Eeff0^2 + Eeff0 *[ (Ei-eo)/3*(1-4f)-eo] - eo (Ei-1)/3*(1-f) = 0

  b=(es-e0)*((1.0-4.0*f)/3.0)-e0
  c=-e0*(es-e0)*(1.0-f)/3.0

  discriminant= b**2 - 4*c

  !solution

  Eeff0 = 0.5*(-b+sqrt(discriminant))
 
  if (real(Eeff0)<1) then
     Eeff0 = 0.5*(-b-sqrt(discriminant))
     print *,'attention.... par la!!'
  endif

  call computeRayleighDistribution(a,ll,aa,number)

  if (.False.) then
     aa(:)=a    !! Test monodisperse
     number(:)=1
  endif

  ff=4.0/3.0*pi*aa**3 * number
  ! normalize
  norm=f/sum(ff)
  voldens=number*norm
  ff=ff*norm

  call computeHmultispecies(ll,aa,voldens,H)

  !Solve 1st-order solution: E   (equation 3-39 Ya Qiu Jin et 5.3.145 Tsang III)
  lambda=3e8/frequency
  k0=2*pi/lambda

  Eeff=Eeff0
  y=(es-e0)/(3*Eeff0 + es-e0)
  D=cone-f*y ! f (=sum(ff)) because all our particules have the same dielectric constant here

  suml1=0.0
  suml2=0.0
  do l=1,ll
     sumj=aa(l)**3*y
     do j=1,ll
        sumj= sumj + (2*pi*aa(j))**3*voldens(j)*y*H(j,l)  ! (5.3.145, Tsang III)
     enddo
     !!print *,aa(l),ff(l)*aa(l)**3,ff(l)*abs(sumj/y)
     suml1= suml1 + ff(l)*y*(cone + cmplx(0.0,2.0/3.0)*Eeff0*sqrt(Eeff0)*k0**3/ D * sumj) !! to compute E
     suml2= suml2 + real(ff(l)*y*conjg(sumj)) !! to compute albedo (7.4.6 Tsang III)
  enddo
     
  Eeff=e0 + 3*Eeff0/D*suml1

  ! albedo equation: 7.4.6 Tsang III
  beta=k0 * 2*aimag(sqrt(Eeff))

  albedo=2*k0**4*abs(Eeff0)**2/beta/abs(D)**2 * real(suml2)

end subroutine dmrtparameters_dist


end module mod_dmrtparameters
