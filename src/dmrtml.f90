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

!> Main entry point for the DMRT-ML model. Most users should use the subroutine \ref dmrtml.
!! This module includes subroutines to be called from Fortran as well as those exposed from Python (see dmrtml.py for the main python module)
!!
module mod_dmrtml

contains

!> Compute the brightness temperatures (at H and V polarisations and for several angles) for a given 
!! frequency and layered medium.
!!
!! The medium is characterized by the number of layers, the depth of each layer, the type of medium in each layer
!! (Snow or bubbly Ice) and the lower interface (called soil but can also be pure ice, etc...).
!! In the case of snow layers, the grain size, density and temperature are required parameters in every layer. Stickiness and wetness inputs 
!! are optional. In the case of bubbly ice, the bubble size, density and temperature are required and the stickiness 
!! is optional. In the future, medium of type Water will be added.
!! 
!! @param[in] frequency frequency (Hz).
!! @param[in] l number of layers.
!! @param[in] n number of stream in the most refringent layer (see \ref mod_disort for details).
!! @param[in] depth depth/thickness (in meter) of each layer.
!! @param[in] density (kg/m3) of Snow or Ice in each layer.
!! @param[in] radius (m) of snow grains or air bubbles in each layer.
!! @param[in] temperature temperature (K) in each layer.
!!
!! @param[out] TbV outgoing V-polarized brightness temperature (K) for every emerging streams.
!! @param[out] TbH outgoing H-polarized brightness temperature (K) for every emerging streams.
!! @param[out] mhu cosine angles in the air of the emerging streams.
!!
!! @param[in] tau stickiness factor in each layer. If not present, none-sticky calculation (tau=infinity).
!! @param[in] fwetness fractional volume of water with respect to volume of ice. Default is 0.
!! @param[in] medium type of medium in each layer (valid value is S=Snow or I=ice with air bubbles or A=automatic selection between Snow and Ice or B=bridging). Default is S=Snow.
!!
!! @param[in] dist if True a poly-disperse distribution of particule size is considered. Radius is taken from a Rayleigh distribution. The radius parameter in input corresponds to the optical radius (equivalent S/V grain size). If False a mono-disperse distribution is used, this is the default.
!! @param[in] soilp description of the soil. See \ref mod_soil for details.
!! @param[in] Tbatmodown atmosphere downwelling brightness temperature (assumed isotropic and unpolarized). Can be used to compute snow emissivity (see \ref mod_disort for details).
!! @param[in] Eice ice dielectric constant (optional). To presribe a dielectric constant instead of using the Matlzer&Wegmuller formula.
!! @param[out] profileV if present and allocated, the computed profile of upwelling V-polarized Tb at the top of each layer (before being transmitted).
!! @param[out] profileH if present and allocated, the computed profile of upwelling H-polarized Tb at the top of each layer (before being transmitted).
!!

!!
!!
!! This subroutine calculates for each layer the DMRT parameters (exctinction and single scattering albedo) using 
!! one of the subroutines in the \ref mod_dmrtparameters module and then solves the 
!! radiative transfer equation with the DMRT parameters of each layer, the soil boundary conditions (if present)
!! and the downwelling atmospheric contribution (if present) using the disort subroutine in the \ref mod_disort module.
!! It is recommanded to read the documentation of both modules to understand the details.
!! See also \ref mod_soil for describing the soil (or for modeling any kind of interfaces lying below the layers).
!!
!! The subroutine returns the brightness temperature at H and V polarization and at several angles given by acos(mhu).
!! It is not possible to choose these angles in input because they are automatically determined in disortml
!! to optimize the solution (i.e. they result from a gaussian quadrature in the most refringent layer and refraction in other layers (see \ref mod_disort for details). To compute the brightness temperature at a particular angle, it is recommand to use interpolation and/or 
!! increase the number of stream to refine the quadrature. 
!! Fortran users have to implement their own interpolation method.
!! Python users can use the dmrtml.py module that implements linear interpolation.
!!

subroutine dmrtml(frequency,l,n,depth,density,radius,temperature,&
     TbV,TbH,mhu,tau,fwetness,medium,dist,soilp,Tbatmodown,Eice,profileV,profileH)
  use mod_dmrtparameters
  use mod_disort
  use mod_soil

  implicit none

  !--------------------------------------
  integer,intent(in)                         :: l,n   ! l number of layer, n number of stream in the most refringent layer
  real*8,intent(in)                          :: frequency ! microwave frequency (in Hz)
  real*8,intent(in),dimension(l)             :: temperature,depth,radius,density ! temperature (K), depth (m), radius (m) and density (kg/m3) in each layer
  ! outputs                                  
  real*8,intent(out),dimension(n)            :: TbV,TbH,mhu ! brightness temperature (K)
  ! optional arguments                       
  real*8,intent(in),dimension(l),optional    :: tau ! stickiness factor
  real*8,intent(in),dimension(l),optional    :: fwetness ! fwetness (kg/kg)
  character,intent(in),dimension(l),optional :: medium ! ice sphere (.FALSE.) or air bubble (.TRUE.)

  logical,intent(in),optional                :: dist        ! use rayleigh distribution instead of monodispere particules
  real*8,intent(in),optional                 :: Tbatmodown ! atmosphere brightness temperature 
  complex*16,intent(in),optional             :: Eice
  type(soilparams),intent(in),optional       :: soilp ! soil parameters (see mod_soil)
  real*8,intent(out),dimension(n,l),optional :: profileV,profileH ! ! brightness temperature profile

  !--------------------------------------
  ! DMRT paramaters
  real*8,dimension(l)             :: ke,albedo
  complex*16,dimension(l)         :: eps
  character                       :: medium_
  integer                         :: k
  complex*16                      :: Eo
  real*8                          :: f,tau_,fwetness_
 
  if (frequency<1e9) then
     print *,'frequency must be given in Hz.',frequency
     stop
  endif

  fwetness_=0.0
  medium_='S'
  !print *,'#dmrtfunc',l,freq,depth,temperature,radius,density
  ! calculate DMRT parameters
  do k=1,l
     !print *,'f=',f,tau(k)
     if (present(tau)) then
        tau_=tau(k)
     else
        tau_=1000
     endif

    if (present(fwetness)) then
        fwetness_=fwetness(k)
     endif

     if (present(medium)) then
        if (medium(k)=='A') then
           if (density(k)>0.5*917) then
              medium_='I'
           else
              medium_='S'
           endif
        else
           medium_=medium(k)
        endif
     endif

     if (medium_=='I') then
        f = 1.0- density(k)/917.0
     else
        f = density(k)/917.0
     endif

     if (present(dist) .and. dist) then
        call dmrtparameters_dist(frequency,temperature(k),f,radius(k),&
             tau_,fwetness_,medium_,Eo,eps(k),albedo(k),ke(k),Eice)
     else
        call dmrtparameters_grodyapproach(frequency,temperature(k),f,radius(k),&
             tau_,fwetness_,medium_,Eo,eps(k),albedo(k),ke(k),Eice)
     endif
  enddo

  call mldisort(l,n,depth,temperature,Ke,albedo,eps,&
       TbV,TbH,mhu,soilp,frequency,Tbatmodown,profileV,profileH)

end subroutine dmrtml



!> Compute the extinction and the single scattering albedo.
!!
!! This subroutine is a convenience routine that calls one of the subroutine in the mod_dmrtparameters module 
!! and has a list of arguments similar to those in dmrtml. It is recommended to read the documentation of 
!! the \ref mod_dmrtparameters module
!!
!! @param[in] frequency frequency in Hertz.
!! @param[in] density density of Snow or bubbly Ice (in kg/m3).
!! @param[in] radius radius (in meter) of grains (for Snow) or bubbles (for Ice).
!! @param[in] temperature temperature (in Kelvin).
!!
!! @param[out] albedo single scattering albedo.
!! @param[out] beta exctinction coefficient, also written Ke in some books.
!!
!! @param[in] tau stickiness. If not present, none-sticky calculation (tau=infinity).
!! @param[in] fwetness fractional volume of water with respect to ice. Default is 0.
!! @param[in] medium type of medium in each layer (valid value is S=Snow or I=ice with air bubbles). Default is Snow.
!!
!! @param[in] distribution.  If set to False a mono-disperse distribution is used, this is the default. If set to True a poly-disperse distribution of particules is used: Radii are taken from a Rayleigh distribution. The radius parameter in input corresponds to the optical radius (equivalent S/V grain size).
!! @param[in] grodyapproach.  If set to False, the Grody's approach is not used.

subroutine albedobeta(frequency,density,radius,temperature,&
     albedo,beta,tau,fwetness,medium,dist,grodyapproach)
  use mod_dmrtparameters

  implicit none

  !--------------------------------------
  real*8,intent(in)               :: frequency
  real*8,intent(in)               :: temperature,radius,density
  real*8,intent(out)              :: albedo,beta

  real*8,intent(in),optional      :: tau,fwetness
  character,intent(in),optional   :: medium
  logical,intent(in),optional     :: dist        ! use rayleigh distribution instead of monodispere particules
  logical,intent(in),optional     :: grodyapproach        ! don't use grody approximation if False
  !--------------------------------------
  complex*16                  :: Eo,E
  character                   :: medium_
  real*8                      :: tau_,f,fwetness_
  logical                     :: grodyapproach_

  if (present(tau)) then
     tau_=tau
  else
     tau_=10000
  endif
  f=density/917.0

  if (present(grodyapproach)) then
     grodyapproach_=grodyapproach
  else
     grodyapproach_=.true.
  endif

  if (present(medium)) then
     if (medium=='A') then
        if (f>0.5) then
           medium_='I'
        else
           medium_='S'
        endif
     else
        medium_=medium
     endif
     if (medium_=='I') then
        f=1-f
     endif
  else
     medium_='S'
  endif

  if (present(fwetness)) then
     fwetness_=fwetness
  else
     fwetness_=0.0
  endif

  if (present(dist) .and. dist) then
     call dmrtparameters_dist(frequency,temperature,f,radius,tau_,fwetness_,medium_,&
          Eo,E,albedo,beta)
  else if (grodyapproach_) then
     call dmrtparameters_grodyapproach(frequency,temperature,f,radius,&
          tau_,fwetness_,medium_,&
          Eo,E,albedo,beta)
  else
     call dmrtparameters(frequency,temperature,f,radius,&
          tau_,fwetness_,medium_,&
          Eo,E,albedo,beta)
  endif
end subroutine albedobeta

end module mod_dmrtml



!
! F2PY does not support optional argument and user-defined type...
! F2PY does not support modules
!
subroutine dmrtml_pywrapper(frequency,l,n,depth,density,radius,temperature,&
     TbV,TbH,mhu,tau,fwetness,medium,dist, &
     soilp_imodel,soilp_temperature,soilp_eps_r,soilp_eps_i,soilp_sigma,soilp_SM,&
     soilp_sand,soilp_clay,soilp_dm_rho,soilp_Q,soilp_N,Tbatmodown,eps_ice_r,eps_ice_i)


  use mod_dmrtml
  use mod_soil

  implicit none
  integer,intent(in)             :: l,n 
  real*8,intent(in)              :: frequency 
  real*8,intent(in),dimension(l) :: temperature,depth,radius,density,tau,fwetness
  character,intent(in),dimension(l):: medium
  real*8,intent(out),dimension(n):: TbV,TbH,mhu 
  logical,intent(in)             :: dist     
  integer,intent(in)             :: soilp_imodel
  real*8,intent(in)              :: soilp_eps_r,soilp_eps_i
  real*8,intent(in)              :: soilp_temperature,soilp_sigma,soilp_SM
  real*8,intent(in)              :: soilp_sand,soilp_clay,soilp_dm_rho
  real*8,intent(in)              :: soilp_Q,soilp_N
  real*8,intent(in)              :: Tbatmodown
  real*8,intent(in)              :: eps_ice_r,eps_ice_i
  !------------------------------------------------------------------------
  type(soilparams) :: soilp
  complex*16       :: Eice

  soilp%imodel = soilp_imodel
  soilp%temp   = soilp_temperature
  soilp%eps    = cmplx(soilp_eps_r,soilp_eps_i)
  soilp%sigma  = soilp_sigma
  soilp%Q      = soilp_Q
  soilp%N     = soilp_N
  soilp%SM     = soilp_SM
  soilp%sand   = soilp_sand
  soilp%clay   = soilp_clay
  soilp%dm_rho = soilp_dm_rho
  soilp%N = soilp_N
  soilp%Q = soilp_Q

  Eice = cmplx(eps_ice_r,eps_ice_i)
  
  call dmrtml(frequency,l,n,depth,density,radius,temperature,&
       TbV,TbH,mhu,tau=tau,fwetness=fwetness,&
       medium=medium,dist=dist,soilp=soilp,Tbatmodown=Tbatmodown,&
       Eice=Eice)

end subroutine dmrtml_pywrapper



subroutine albedobeta_pywrapper(frequency,density,radius,temperature,&
     albedo,beta,tau,fwetness,medium,dist,grodyapproach)
  use mod_dmrtml

  implicit none

  !--------------------------------------
  real*8,intent(in)     :: frequency
  real*8,intent(in)     :: temperature,radius,density
  real*8,intent(out)    :: albedo,beta
  real*8,intent(in)     :: tau,fwetness
  character,intent(in)  :: medium
  logical,intent(in)    :: dist        ! use rayleigh distribution instead of monodispere particules
  logical,intent(in)    :: grodyapproach

  call albedobeta(frequency,density,radius,temperature,&
       albedo,beta,tau=tau,fwetness=fwetness,medium=medium,dist=dist,grodyapproach=grodyapproach)
end subroutine albedobeta_pywrapper





subroutine compute_streams_pywrapper(frequency,l,n,density,radius,temperature,&
     tau,fwetness,medium,dist,mhu,weight,ns,outmhu,ns0)

  use mod_dmrtparameters
  use mod_disort
  real*8,intent(in)                    :: frequency
  integer,intent(in)                   :: l
  integer,intent(in)                   :: n
  real*8,intent(in),dimension(l)       :: temperature,radius,density 
  real*8,intent(in),dimension(l)       :: tau
  real*8,intent(in),dimension(l)       :: fwetness
  character,intent(in),dimension(l)    :: medium
  logical,intent(in)                   :: dist 

  real*8,intent(out),dimension(n,l)    :: mhu,weight ! cosin of all layers
  integer,intent(out),dimension(l)     :: ns ! number of stream in each layer 
  real*8,intent(out),dimension(n)      :: outmhu      ! cosine of the angles
  integer,intent(out)                  :: ns0
  !!-------------------------------------------------------------------------
  integer :: k
  real*8  :: f,albedo,ke
  complex*16 :: Eo
  complex*16,dimension(l) :: eps

  do k=1,l
     f=density(k)/917.0

     if (medium(k)=='I') then
        f=1-f
     endif
     if (dist) then
        call dmrtparameters_dist(frequency,temperature(k),f,radius(k),&
             tau(k),fwetness(k),medium(k),Eo,eps(k),albedo,ke)
     else
        call dmrtparameters_grodyapproach(frequency,temperature(k),f,radius(k),&
             tau(k),fwetness(k),medium(k),Eo,eps(k),albedo,ke)
     endif
  enddo

  call compute_streams(l,n,eps,mhu,weight,ns,outmhu,ns0)

end subroutine compute_streams_pywrapper


subroutine icedielectric_pywrapper(frequency,temperature,eps_ice_r,eps_ice_i)
  use mod_dielectric_constant

  implicit none
  real*8,intent(in)       :: frequency,temperature 
  real*8,intent(out)      :: eps_ice_r
  real*8,intent(out)      :: eps_ice_i

  complex*16 :: Ei

  call icedielectric(frequency,temperature,Ei)
  eps_ice_r=real(Ei)
  eps_ice_i=aimag(Ei)

end subroutine icedielectric_pywrapper



subroutine mldisort_pywrapper(l,n,depth,temp,Ke,albedo,reps,Tbatmodown,TbV,TbH,outmhu)   !!,soilp,frequency,Tbatmodown)

    use mod_disort


    implicit none
    integer,intent(in)                   :: l           ! number of layers
    integer,intent(in)                   :: n           ! number of stream in the most refringent layer
    real*8,intent(in),dimension(l)       :: depth       ! depth of the layer (first is at the top)
    real*8,intent(in),dimension(l)       :: temp        ! temperature of the layer (first is at the top)
    real*8,intent(in),dimension(l)       :: Ke,albedo   ! effective extinction and single scatt albedo
    real*8,intent(in),dimension(l)       :: reps         ! effective dielectric constant of the snow
    real*8,intent(out),dimension(n)      :: TbV,TbH     ! brightness temperature
    real*8,intent(out),dimension(n)      :: outmhu      ! cosine of the angles
    !!!!!!!!!!!!type(soilparams),intent(in)          :: soilp  ! soil parameters
    !!!!!!!!!!!!real*8,intent(in)                    :: frequency ! necessary for some soil...
    real*8,intent(in)                    :: Tbatmodown  ! atmosphere brightness temperature
    !---------------------
    integer    :: i
    complex*16, dimension(l) :: ceps


    ceps(:) = reps(:)

    call mldisort(l,n,depth,temp,Ke,albedo,ceps,TbV,TbH,outmhu,Tbatmodown=Tbatmodown) !!!,soilp,frequency,Tbatmodown)

  end subroutine mldisort_pywrapper
