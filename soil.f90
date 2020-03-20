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

!> Description of soil properties and calculation of the reflection coefficients

module mod_soil


  !> This type provides information on the soil (type of reflection model, soil properties, ...) used as input of the subroutine \ref soil_coefficients_multiangle.
  !! The imodel and temp variables are required. Other variables are used or not depending on imodel 
  !! (see the code of \ref soil_coefficients_multiangle for details)
  !!
  type soilparams
     integer     :: imodel  !< Identification number of the soil reflection model
     real*8      :: temp    !< Soil temperature (K)

     !-------------- the following parameters are used or not depending on imodel ------------
     complex*16  :: eps               !< Effective dielectric constant of the soil 
     real*8      :: sigma             !< Surface roughness rms height (m)
     real*8      :: SM                !< Soil moisture [0..1]
     real*8      :: sand              !< Sand fraction [0..1]
     real*8      :: clay              !< Clay fraction [0..1]
     real*8      :: dm_rho            !< Dry matter density (kg/m3)
     real*8      :: Q                 !< Polarization mixing effects
     real*8      :: N                 !< Modulation as a function of angle


  end type soilparams



contains

  !> Compute the reflection coefficients of the soil at several angles.
  !!
  !! This subroutine dispatches the calculation depending on the type of interface and/or soil models.
  !! Available models include:
  !!
  !! - imodel=0 no soil (rh=rv=0)
  !! - imodel=1 flat surface, fresnel coefficient with prescribed dieletric constant
  !! - imodel=2 flat surface, fresnel coefficient with dieletric constant calculated with HUT EPSS.m code (Pulliainen et al., 1999)
  !! - imodel=3 flat surface, fresnel coefficient with dieletric constant calculated with Dobson et al. (1985)
  !! - imodel=4 flat surface, fresnel coefficient with dieletric constant calculated with Mironov
  !!
  !! - imodel=101 Wegmüller & Mätzler (1999) reflectivity for rough surface with prescribed dieletric constant
  !! - imodel=102 Wegmüller & Mätzler (1999) reflectivity for rough surface with dieletric constant calculated with HUT EPSS.m code (Pulliainen et al., 1999)
  !! - imodel=103 Wegmüller & Mätzler (1999) reflectivity for rough surface with dieletric constant calculated with Dobson et al. (1985)
  !! - imodel=104 Wegmüller & Mätzler (1999) reflectivity for rough surface with dieletric constant calculated with Mironov
  !!
  !! - imodel=201 ice flat surface, fresnel coefficient with prescribed dielectric constant (equivalent to imodel=1)
  !! - imodel=202 ice flat surface, fresnel coefficient with ice dielectric constant model.
  !!
  !! - imodel=301 QHmodel (Wang et al. 1983) for rough surface with prescribed dieletric constant
  !! - imodel=302 QHmodel (Wang et al. 1983) for rough surface with dieletric constant calculated with HUT EPSS.m code (Pulliainen et al., 1999)
  !! - imodel=303 QHmodel (Wang et al. 1983) for rough surface with dieletric constant calculated with Dobson et al. (1985)
  !! - imodel=304 QHmodel (Wang et al. 1983) for rough surface with dieletric constant calculated with Mironov
  !! - imodel=401 water flat surface, fresnel coefficient with prescribed dielectric constant (equivalent to imodel=1 and imodel=201)
  !! - imodel=402 water flat surface, fresnel coefficient with water dielectric constant model.
  !!
  !! This list can be extended by users by choosing a unique number as model ID (recommended: over 10000 and model familly as multiple of 100).
  !! Please consider to share your piece of code for inclusion in the main distribution of DMRTML (ghislain dot picard at ujf-grenoble fr)
  !!
  !! The soil properties are described by filling the parameter soilp of type(\ref soilparams).
  !!
  !! This subroutine is automatically called by mldisort (see \ref mod_disort).
  !!
  !! REFERENCES:
  !! 
  !! Pulliainen, J.T., J. Grandell and M.T. Hallikainen. 1999. HUT snow emission model and its applicability to snow water equivalent retrieval. IEEE Trans. Geosci. Remote Sens., 37(3), 1378–1390.
  !!
  !! Dobson, M.C., Ulaby, F.T., Hallikainen, M.T., El-Rayes, M.A. 1985. Microwave dielectric behavior of wet soil-part II : dielectric miximg models. IEEE transanctions on geoscience and remote sensing. vol.GE-23 (1). 35-46.
  !!
  !! Wang, J.R., O'Neill, P.E. Jackson, T.J., Engman, E.T. (1983) Multifrequency measurements of the effects of soil moisture, soil texture and surface roughness, IEEE transactions on geoscience and remote sensing, vol.GE-21(1), p.44-51. 
  !! 
  !! Wegmuller, U. & Matzler, C. Rough bare soil reflectivity model. Geoscience and Remote Sensing, IEEE Transactions on 37, 1391–1395 (1999).



  subroutine soil_coefficients_multiangle(eps_upper,soilp,frequency,mhu,rv,rh)
    use mod_fresnel
    use mod_dielectric_constant

    implicit none
    complex*16,intent(in)             :: eps_upper !< Effective dielectric constant of the layer above the soil
    type(soilparams),intent(in)       :: soilp     !< Soil parameters, see \ref soilparams for the documentation
    real*8,intent(in)                 :: frequency !< Frequency (Hz)
    real*8,dimension(:),intent(in)    :: mhu       !< Cosines in the layer above the soil
    real*8,dimension(:),intent(out)   :: rh        !< H-pol reflection coefficients for the given cosines
    real*8,dimension(:),intent(out)   :: rv        !< V-pol reflection coefficients for the given cosines

    !-----------------------------------------------------
    real*8           :: ksigma
    real*8,parameter :: pi=3.14159265
    complex*16       :: soileps ! effective dielectric of the soil
    real*8,dimension(size(rh))        :: rsh,rsv
    complex*16                        :: Ei,Ew
    !-----------------------------------------------------

    soileps=soilp%eps !

    if (soilp%imodel>0.and.soilp%imodel<=4) then !----- flat interface 
       ! basic fresnel coefficient
       if (soilp%imodel==2) then  ! eps calculated from SM,sand,clay and dm_rho (after Pulliainen et al.1999) 
          call soil_dielectric_constant_pulliainen(frequency,soilp%temp,soilp%SM,soilp%sand,soilp%clay,soilp%dm_rho,soileps)
       end if

       if (soilp%imodel==3) then  ! eps calculated from SM,sand,clay(extracted from HUtnlayer) 
          call soil_dielectric_constant_dobson(frequency,soilp%temp,soilp%SM,soilp%sand,soilp%clay,soileps)
       end if

       if (soilp%imodel==4) then  ! eps calculated from SM,clay (after Mironov) 
          call soil_dielectric_constant_mironov(frequency,soilp%clay,soilp%temp,soilp%SM,soileps)
       end if

       call fresnel_coefficients_multiangle(eps_upper,soileps,mhu,rv,rh)

    else if (soilp%imodel>100.and.soilp%imodel<=104) then !----- Wegmüller & Mätzler (1999) family  
       ! inspired from HUT (1 layer) J.Pulliainen/8.12.1997 (Matlab epss.m code) coding soil reflectivity by Wegmüller & Mätzler (1999) 

       if (soilp%imodel==102) then  ! eps calculated form SM,sand,clay and dm_rho (after Pulliainen et al.1999) 
          call soil_dielectric_constant_pulliainen(frequency,soilp%temp,soilp%SM,soilp%sand,soilp%clay,soilp%dm_rho,soileps)
       end if

       if (soilp%imodel==103) then  ! eps calculated form SM,sand,clay(extractedinspired from HUtnlayer) 
          call soil_dielectric_constant_dobson(frequency,soilp%temp,soilp%SM,soilp%sand,soilp%clay,soileps)
       end if

       if (soilp%imodel==104) then  ! eps calculated from SM,clay (after Mironov) 
          call soil_dielectric_constant_mironov(frequency,soilp%clay,soilp%temp,soilp%SM,soileps)
       end if

       call fresnel_coefficients_multiangle(eps_upper,soileps,mhu,rv,rh) ! calculate Fresnel coefficient for a flat surface

       !rough soil reflectivity model of Wegmüller & Mätzler (1999)

       !Calculate ksigma = wavenumber*soilp%sigma(standard deviation of surface height)
       ksigma=real(2*pi*frequency*sqrt((1/2.9979e8)**2*eps_upper))*soilp%sigma 

       !Calculation of rh with ksigma
       rh = rh *exp(-ksigma**(sqrt(0.1*mhu)))

       !calculation of rv with rh (the model is valid for angle between 0-70°
       rv = rh * mhu**0.655 ! <-- * ou ** ??
       WHERE ( mhu< cos(60*pi/180) ) rv = rh*(0.635-0.0014*(acos(mhu)*180/pi-60)) 


    else if (soilp%imodel>200.and.soilp%imodel<=202) then    !Soil = ice     

       if (soilp%imodel==201) then  ! 
          Ei = soilp%eps
       end if
       if (soilp%imodel==202) then  ! 
          call icedielectric(frequency,soilp%temp,Ei)  !Call module ice_dielec for ice dielectric constant (Ei) calculation
       end if

       call fresnel_coefficients_multiangle(eps_upper,Ei,mhu,rv,rh)

    else if (soilp%imodel>300.and.soilp%imodel<=304) then !----- QHmodel (Wang et al. 1983)

       if (soilp%imodel==302) then  ! eps calculated form SM,sand,clay and dm_rho (after Pulliainen et al.1999) 
          call soil_dielectric_constant_pulliainen(frequency,soilp%temp,soilp%SM,soilp%sand,soilp%clay,soilp%dm_rho,soileps)
       end if

       if (soilp%imodel==303) then  ! eps calculated form SM,sand,clay(extracted from HUtnlayer) 
          call soil_dielectric_constant_dobson(frequency,soilp%temp,soilp%SM,soilp%sand,soilp%clay,soileps)
       end if

       if (soilp%imodel==304) then  ! eps calculated from SM,clay (after Mironov) 
          call soil_dielectric_constant_mironov(frequency,soilp%clay,soilp%temp,soilp%SM,soileps)
       end if

       call fresnel_coefficients_multiangle(eps_upper,soileps,mhu,rsv,rsh) !Calculate fresnell coefficient for a plat surface

       rv = ((1-soilp%Q)* rsv + soilp%Q* rsh)*exp(-soilp%sigma*(mhu**soilp%N))
       rh = ((1-soilp%Q)* rsh + soilp%Q *rsv)*exp(-soilp%sigma*(mhu**soilp%N))

    else if (soilp%imodel>400.and.soilp%imodel<=402) then    !Soil = water

       if (soilp%imodel==401) then  ! 
          Ei = soilp%eps
       end if
       if (soilp%imodel==402) then  ! 
          call waterdielectric(frequency,soilp%temp,Ew)  
       end if

       call fresnel_coefficients_multiangle(eps_upper,Ew,mhu,rv,rh)
    else
       print *,'unknown soil model:',soilp%imodel 
       stop
    end if

  end subroutine soil_coefficients_multiangle


  subroutine soil_dielectric_constant_pulliainen(frequency,tempK,SM,sand,clay,dm_rho,soil_dielectric) !(after Pulliainen et al.1999) 
    use mod_dielectric_constant

    implicit none

    real*8,intent(in)      :: frequency,tempK            ! microwave frequency (in Hz)
    real*8,intent(in)      :: SM,sand,clay,dm_rho        ! soil moisture, proportion of sand, clay and dry matter density
    complex*16,intent(out) :: soil_dielectric
    real*8,parameter       :: pi=3.14159265  
    !------------------------------------
    real*8                 :: beta,temp
    real*8                 :: ew_r, ew_i                 !water or ice dielectricity     
    complex*16             :: epsalf,Ei
    !Parameters for soil dielectric constant calculation with water
    real*8,parameter       :: e0=8.854e12,ew_inf=4.9
    real*8                 :: ew0,d,alfa,tw

    temp = tempK-273.15   

    if (temp>0) then
       !calculates real and imag. part of water dielectricity (code HUT 20.12.95 [epsw.m]; K.Tigerstedt)
       ew0 = 87.74 - 0.40008*temp + 9.398e-4 * temp**2 + 1.410e-6 * temp**3
       d = 25 - temp
       alfa = 2.033e-2+1.266e-4*d+2.464e-6*d**2
       tw = 1/(2*pi)*(1.1109e-10 - 3.824e-12 * temp + 6.938e-14 * temp**2 - 5.096e-16*temp**3)

       ew_r = ew_inf + (ew0 - ew_inf)/(1+(2*pi*frequency*tw)**2)
       ew_i = (ew0-ew_inf)*2*pi*frequency*tw/(1+(2*pi*frequency*tw)**2)
       !ok
    else   

       call icedielectric(frequency,tempK,Ei)
       ew_r=real(Ei)
       ew_i=aimag(Ei)

       !option for salt consideration (Mätzler 1987) 
       !iei_S =A/M+B*M**C                 !impure ice 
       !iei_P=Ap/M+Bp*M**Cp                 !pure ice
       !delta_iei = iei_S - iei_P
       !ew_i=ew_i+delta_iei*SS/13

    end if

    beta = 1.09-0.11*sand+0.18*clay
    epsalf = 1 + 0.65*dm_rho/1000.0 + SM**beta*(cmplx(ew_r,ew_i)**0.65 - 1) ! dm_rho is now in SI // Ulaby et al. (1986, p. 2099)
    soil_dielectric = (epsalf)**(1/0.65)

  end subroutine soil_dielectric_constant_pulliainen

  subroutine soil_dielectric_constant_dobson(frequency,tempK,SM,S,C,soil_dielectric) !Dobson et al.,(1985) dielectric constant calculation
    !(extracted from HUTnlayer code [Lemmetyinen et al.,(2010)]) 
    use mod_dielectric_constant

    implicit none

    real*8,intent(in)      :: frequency,tempK            ! microwave frequency (in Hz)
    real*8,intent(in)      :: SM,S,C        ! soil moisture, proportion of sand, clay and dry organic matter
    complex*16,intent(out) :: soil_dielectric
    real*8,parameter       :: pi=3.14159265  
    !------------------------------------
    real*8                 :: beta1,beta2,temp,sigma_eff     
    real*8,parameter       :: e_0=8.854e12,e_w_inf=4.9,e_s=4.7,rho_b=1.3,rho_s=2.664 !e_s=solid soil dielectric constant
    real*8                 :: e_w0,rt_w,e_fw1,e_fw2

    temp = tempK-273.15  

    beta1=1.2748-0.519*S-0.152*C;
    beta2=1.33797-0.603*S-0.166*C;

    sigma_eff=0.0467+0.2204*rho_b-0.4111*S+0.6614*C; 

    e_w0=87.134-1.949e-1*temp-1.276e-2*temp**2+2.491e-4*temp**3;
    rt_w=(1.1109e-10-3.824e-12*temp+6.938e-14*temp**2-5.096e-16*temp**3)/(2*pi);

    e_fw1=e_w_inf+(e_w0-e_w_inf)/(1+(2*pi*frequency*rt_w)**2);
    e_fw2=2*pi*frequency*rt_w*(e_w0-e_w_inf)/(1+(2*pi*frequency*rt_w)**2)+sigma_eff*(rho_s-rho_b)/(2*pi*frequency*e_0*rho_s*SM);

    soil_dielectric=cmplx((1+(rho_b/rho_s)*(e_s**0.65-1)+SM**beta1*e_fw1**0.65-SM)**(1/0.65),(SM**beta2*e_fw2**0.65)**(1/0.65));

  end subroutine soil_dielectric_constant_dobson

  subroutine soil_dielectric_constant_mironov(theF,c,tK,sm,diel)
    ! Translate from a matlab code provides by A. Roy
    ! (initial release J-P Wigneron, H.Lawrence and C.Duffour, and modified by P. Richaume and B. Montpetit)
    !
    ! Input:
    !   theF: double scalar of vector ; frequency in Hz
    !   c:    double scalar or vector ; Clay fraction [0 1]
    !   tK:   double scalar or vector ; temperature in Kelvin
    !   sm:   double scalar or vector ; volumetric soil moisture m3/m3
    !
    ! Output:
    !   diel: double vector ; dielectric constant
    !
    ! Comments:
    ! Three options are not implementes as input in this version and are putted by default:
    ! * Standard (=0) or symetrized (=1) Mironov version
    ! sym=1 by default   
    ! * Substituon polynom function and constraints around SM
    ! sm1=0.02 by default
    ! * M=0 neigborhood where the polynom approximation is done i.e. for |sm| < sm1
    !   tpoly = 0 polynom form: P(X) = aX2 + bX + c       such P(sm1)=eps(sm1) & dP/dsm(0)=0 & P(0)=eps(0)  
    !   tpoly = 1 polynom form: P(X) = aX2 + bX + c       such P(sm1)=eps(sm1) & dP/dsm(0)=0 & dP/dsm(sm1)=deps/dsm(sm1)
    !   tpoly = 3 polynom form: P(X) = aX3 + bX2 + cX + d such P(0)=eps(0) & dP/dsm(0)=0 & P(sm1)=eps(sm1) & dP/dsm(sm1)=deps/dsm(sm1)  
    !   tpoly = 4 polynom form: P(X) = aX4 + bX2 + c      such P(sm1)=eps(sm1) & dP/dsm(0)=0 & P(0)=eps(0)
    ! tpoly=1 by default
    ! 
    !   If used in a retrieval process, sym should be set to 1 to avoid very negative SM values and using the paraboloic appoximation type 1 (default) 
    !
    ! References:
    ! Mironov et al PIERS proc. 2009 

    implicit none

    real*8,intent(in)       :: theF ! frequency in Hz
    real*8,intent(in)       :: c    ! Clay fraction [0 1]
    real*8,intent(in)       :: tK   ! temperature in Kelvin
    real*8,intent(in)       :: sm   ! volumetric soil moisture m¬≥/m¬≥
    complex*16,intent(out)  :: diel ! soil dielectric constant
    !------------------------------------
    real*8,parameter :: pi=3.14159265  
    real*8  :: sym, sm1, tpoly
    real*8  :: thePERMIT0, theEPWI0, theND0, theND1, theND2, theKD0, theKD1, theXMVT0, theXMVT1
    real*8  :: theTF0, theE0PB0, theE0PB1, theE0PB2, theBVB0, theBVB1, theBVB2, theBVB3, theBVB4
    real*8  :: theBSGB0, theBSGB1, theBSGB2, theBSGB3, theBSGB4, theDHBR0, theDHBR1, theDHBR2
    real*8  :: theDSRB0, theDSRB1, theDSRB2, theTAUB0, theSBT0, theSBT1, theE0PU, theBVU0, theBVU1
    real*8  :: theBSGU0, theBSGU1, theDHUR0, theDHUR1, theDSUR0, theDSUR1, theTAUU0, theSUT0, theSUT1
    real*8  :: delta_tk_theTF0
    real*8  :: nd, kd, xmvt, e0pb, Bvb, Bsgb, Fpb, ep0b, dHbR, dSbR, taub, sigmab
    real*8  :: Bvu, Bsgu, Fpu, ep0u, dHuR, dSuR, tauu, sigmau
    real*8  :: twopiFtaub, twopiP0F, cxb, epwbx, epwby, twopiFtauu, cxu, epwux, epwuy
    real*8  :: smabs
    real*8  :: epwbnorm, nb, kb, epwunorm, nu, ku
    real*8  :: epmx, epmy, nm, km
    real*8  :: ex0, ey0, ex1, ey1, dex1, dey1, dx01, dy01
    real*8  :: ax, ay, bx, by, cx, cy, dx, dy
    real*8  :: tsqrt, SMdiel10, SMdiel11, SMdiel20, SMdiel21

    ! Optional inputs
    ! => keep default values
    sym=1    !! sym=0 Mironov standard, sym=1 Symetrized Mironov 
    sm1=0.02 !! substituon polynom function and constraints around SM
    tpoly=1  !! M=0 neigborhood where the polynom approximation is done i.e. for |sm| < sm1 

    ! GPP10
    thePERMIT0=8.854e-12
    ! GPP 11
    theEPWI0=4.9
    ! GPP12
    theND0=1.62847
    ! GPP13
    theND1=-0.70803
    ! GPP14
    theND2=0.4659
    ! GPP15
    theKD0=0.03945
    ! GPP16
    theKD1=-0.03721
    ! GPP17
    theXMVT0=0.00625
    ! GPP18
    theXMVT1=0.33918
    ! GPP19
    theTF0=20.+273.15
    ! GPP20
    theE0PB0=79.82918
    ! GPP21
    theE0PB1=-85.36581
    ! GPP22
    theE0PB2=32.70444
    ! GPP23
    theBVB0=-5.96311e-19
    ! GPP24
    theBVB1=-1.25999e-3
    ! GPP25
    theBVB2=1.83991e-3
    ! GPP26
    theBVB3=-9.77347e-4
    ! GPP27
    theBVB4=-1.39013e-7
    ! GPP28
    theBSGB0=0.0028
    ! GPP29
    theBSGB1=2.37388e-2
    ! GPP30
    theBSGB2=-2.93876e-2
    ! GPP31
    theBSGB3=3.28954e-2
    ! GPP32
    theBSGB4=-2.00582e-2
    ! GPP33
    theDHBR0=1466.80741
    ! GPP34
    theDHBR1=26.97032e2
    ! GPP35
    theDHBR2=-0.09803e4
    ! GPP36
    theDSRB0=0.88775
    ! GPP37
    theDSRB1=0.09697e2
    ! GPP38
    theDSRB2=-4.2622
    ! GPP39
    theTAUB0=48.0e-12
    ! GPP40
    theSBT0=0.29721
    ! GPP41
    theSBT1=0.49
    ! GPP42
    theE0PU=100.
    ! GPP43
    theBVU0=1.10511e-4
    ! GPP44
    theBVU1=-5.16834e-6
    ! GPP45
    theBSGU0=0.00277
    ! GPP46
    theBSGU1=3.58315e-2
    ! GPP47
    theDHUR0=2230.20237
    ! GPP48
    theDHUR1=-0.39234e2
    ! GPP49
    theDSUR0=3.6439
    ! GPP50
    theDSUR1=-0.00134e2
    ! GPP51
    theTAUU0=48.0e-12
    ! GPP52
    theSUT0=0.12799
    ! GPP53
    theSUT1=1.65164

    delta_tk_theTF0=tk-theTF0

    ! MIRONOV formulation

    ! -----------------------------------------------------------
    ! Initializing the GRMDM spectroscopic parameters with clay (fraction)

    ! RI & NAC of dry soils
    !PAR ND0=1.62847 ND1=-0.70803 ND2=0.4659
    nd = theND0 + theND1*c +  theND2*c**2

    !PAR KD0=0.03945 KD1=-0.03721
    kd=theKD0 + theKD1*c

    ! maximum bound water fraction
    !PAR XMVT0=0.00625  XMVT1=0.33918
    xmvt = theXMVT0 + theXMVT1*c

    ! bound water parameters
    !PAR TF0=20
    ! t=20 starting temperature for parameters' fit

    ! ep0b computation
    !PAR E0PB0=79.82918 E0PB1=-85.36581 E0PB2=32.70444
    e0pb = theE0PB0 + theE0PB1*c + theE0PB2*c**2

    !PAR BVB0=-5.96311e-19 BVB1=-1.25999e-3 BVB2=1.83991e-3 BVB3=-9.77347e-4 BVB4=-1.39013e-7
    Bvb = theBVB0 + theBVB1*c + theBVB2*c**2 + theBVB3*c**3 + theBVB4*c**4

    !PAR BSGB0=0.0028 BSGB1=2.37388e-2 BSGB2=-2.93876e-2 BSGB3=3.28954e-2 BSGB4=-2.00582e-2
    Bsgb = theBSGB0 + theBSGB1*c + theBSGB2*c**2 + theBSGB3*c**3 + theBSGB4*c**4

    Fpb = log( (e0pb - 1.) / (e0pb+2.) ) 
    ep0b = (1. + 2.* exp(Fpb - Bvb* delta_tk_theTF0) )/( 1 - exp(Fpb - Bvb* delta_tk_theTF0) )

    !! taub computation
    !PAR DHBR0=1466.80741 DHBR1=26.97032e2 DHBR2=-0.09803e4
    dHbR = theDHBR0 + theDHBR1*c + theDHBR2*c**2

    !PAR DSRB0=0.88775 DSRB1=0.09697e2 DSRB2=4.2622
    dSbR = theDSRB0 + theDSRB1*c + theDSRB2*c**2

    !PAR TAUB0=48e-12
    taub = theTAUB0* exp( dHbR/tK - dSbR )/tK

    !! sigmab computation
    !PAR SBT0=0.29721 SBT1=0.49
    sigmab = theSBT0 + theSBT1*c + Bsgb* delta_tk_theTF0

    ! unbound (free) water parameters
    !! ep0u computation
    !PAR E0PU=100
    !e0pu = 100

    !PAR BVU0=1.10511e-4 BVU1=-5.16834e-6
    Bvu = theBVU0 + theBVU1*c

    !PAR BSGU0=0.00277 BSGU1=3.58315e-2
    Bsgu = theBSGU0 + theBSGU1*c
    Fpu = log( (theE0PU-1.)/(theE0PU+2.) ) - Bvu*delta_tk_theTF0
    ep0u = ( 1. + 2.*exp(Fpu) )/( 1 - exp(Fpu) )

    !!tauu computation
    !PAR DHUR0=2230.20237 DHUR1=-0.39234e2
    dHuR = theDHUR0 + theDHUR1*c
    !PAR DSUR0=3.6439 DSUR1=-0.00134e2
    dSuR = theDSUR0 + theDSUR1*c
    !PAR TAUU0=48e-12
    tauu = theTAUU0* exp( dHuR/tK - dSuR )/tK

    !!sigmau computation
    !PAR SUT0=0.12799 SUT1=1.65164
    sigmau = theSUT0 + theSUT1*c + Bsgu*delta_tk_theTF0

    ! -----------------------------------------------------------
    ! computation of epsilon water (bound & unbound)
    twopiFtaub=2.*pi*theF*taub
    twopiP0F=2.*pi*thePERMIT0*theF

    cxb = (ep0b - theEPWI0) / (1. + twopiFtaub*twopiFtaub)
    epwbx = theEPWI0 + cxb
    epwby =  cxb*twopiFtaub  + sigmab/twopiP0F

    twopiFtauu=2.*pi*theF*tauu
    cxu = (ep0u - theEPWI0) / (1. + twopiFtauu*twopiFtauu)
    epwux = theEPWI0 + cxu
    epwuy =  cxu*twopiFtauu + sigmau/twopiP0F

    ! -----------------------------------------------------------
    ! computation of refractive index of water (bound & unbound):
    epwbnorm=sqrt(epwbx*epwbx+epwby*epwby)

    nb= sqrt( ( epwbnorm + epwbx ) / 2. )
    kb= sqrt( ( epwbnorm - epwbx ) / 2. )

    epwunorm=sqrt(epwux*epwux+epwuy*epwuy)
    nu= sqrt( ( epwunorm + epwux ) / 2. )
    ku= sqrt( ( epwunorm - epwux ) / 2. )

    ! -----------------------------------------------------------
    ! computation of soil refractive index (nm & km):

    ! sym=0
    if (sym==0) then
       !! 'std'
       if (sm >= xmvt) then
          nm = nd + (nb-1.)*xmvt + (nu-1.)*(sm-xmvt)
          km = kd + kb*xmvt + ku*(sm-xmvt)
       else
          nm = nd + (nb-1.)*sm
          km = kd + kb*sm
       end if
       ! -----------------------------------------------------------
       ! computation of soil dielectric constant:

       epmx= nm**2 - km**2    ! real part
       epmy= nm*km*2.      ! imaginary part
    else 
       !! 'symetrized' if (sm <= sm1)
       smabs=abs(sm)
       if (smabs <= sm1) then
          ex0=nd**2. - kd**2.
          ey0=2.*kd*nd

          ex1=(nd + sm1*(nb - 1.))**2. - (kd + kb*sm1)**2.
          ey1=(kd + kb*sm1)*(2.*nd + 2.*sm1*(nb - 1.))

          dex1=2.*(nd + sm1*(nb - 1.))*(nb - 1.) - 2.*kb*(kd + kb*sm1)
          dey1=2.*(kd + kb*sm1)*(nb - 1.) + 2.*kb*(nd + sm1*(nb - 1.))

          dx01=ex0-ex1
          dy01=ey0-ey1

          if (tpoly==0) then
             !! 'V0: aX¬≤ + bx + c'           
             ax=-dx01 / sm1**2.
             ay=-dy01 / sm1**2.

             bx=0.
             by=bx

             cx=ex0
             cy=ey0

             epmx=ax*smabs**2+bx*smabs+cx
             epmy=ay*smabs**2+by*smabs+cy

          else if (tpoly==1) then
             !! 'V1: aX¬≤+ bX + c'
             ax=0.5*dex1/sm1
             ay=0.5*dey1/sm1

             bx=0.
             by=bx

             cx=ex1-0.5*dex1*sm1
             cy=ey1-0.5*dey1*sm1

             epmx=ax*smabs**2+bx*smabs+cx
             epmy=ay*smabs**2+by*smabs+cy

          else if (tpoly==2) then
             !! 'aX¬≥+bX¬≤+cX+d'
             ax=(2*dx01+dex1*sm1)/sm1**3
             ay=(2*dy01+dey1*sm1)/sm1**3

             bx=-(3*dx01+dex1*sm1)/sm1**2
             by=-(3*dy01+dey1*sm1)/sm1**2

             cx=0.
             cy=cx

             dx=ex0
             dy=ey0

             epmx=ax*smabs**3+bx*smabs**2+cx*smabs+dx
             epmy=ay*smabs**3+by*smabs**2+cy*smabs+dy

          else if (tpoly==3) then
             !! 'aX‚?¥+bX¬≤+c'           
             ax=(dx01+0.5*dex1*sm1)/sm1**4
             ay=(dy01+0.5*dey1*sm1)/sm1**4

             bx=-(2*dx01+0.5*dex1*sm1)/sm1**2
             by=-(2*dy01+0.5*dey1*sm1)/sm1**2

             cx=ex0
             cy=ey0

             epmx=ax*smabs**4+bx*smabs**2+cx
             epmy=ay*smabs**4+by*smabs**2+cy

          else 
             print*, 'tpoly (%d) is unknow ',tpoly
          end if
       else
          if (smabs >= xmvt) then
             nm = nd + (nb-1.)*xmvt + (nu-1.)*(smabs-xmvt)
             km = kd + kb*xmvt + ku*(smabs-xmvt)
          else
             nm = nd + (nb-1.)*smabs
             km = kd + kb*smabs
          end if
          ! -----------------------------------------------------------
          ! computation of soil dielectric constant:

          epmx= nm**2 - km**2    ! real part
          epmy= nm*km*2.      ! imaginary part

       end if
       diel=cmplx(epmx, -epmy)
    end if
    ! Value of SM leading to epsr=0, epsi=0
    ! Two roots exists for each
    ! epsr(SM)=0 => SMdiel10, SMdiel11
    tsqrt =(kb**2*nd**2 - kb**2 - 2.*kb*kd*nb*nd + 2.*kb*kd*nd + kd**2*nb**2 - 2.*kd**2*nb + kd**2 + nb**2 - 2.*nb + 1.)**(1./2.)
    SMdiel10= -(nd + kb*kd - nb*nd - tsqrt)/(kb**2 - nb**2 + 2.*nb - 1.)
    SMdiel11= -(nd + kb*kd - nb*nd + tsqrt)/(kb**2 - nb**2 + 2.*nb - 1.)
    ! epsi(SM)=0 => SMdiel20, SMdiel21
    SMdiel20=-nd/(nb-1.)
    SMdiel21=-kd/kb

  end subroutine soil_dielectric_constant_mironov


end module mod_soil
