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

!> Solve the radiative transfer equation in the passive microwave case for 
!! a layered plane-parallel medium
!!
module mod_disort



contains



!> Solve the radiative transfer equation in the passive microwave case for 
!! a layered plane-parallel medium.
!!
!! The solution is based on the Discrete ordinate
!! approach described by Jin, 1994.
!! While the discrete ordinate method is general, the implementation here is specific 
!! to passive microwave configuration with an isotopic medium and considering a Rayleigh phase matrix (i.e. small particules with respect to the wavelength).
!! Under these assumptions the azimuth integration present in the radiative transfer equation can 
!! be calculated  analytically. It means that other phase matrix form would require significant 
!! change in the code.
!! It also implies that using this code for active microwaves is not straightforward, because
!! an explicit treatment of the azimuth is necessary in the active case.
!! (see for instance Picard et al. 2004)
!!
!! The number of streams n specified as input is the number of streams in the most refringent layer 
!! (the one with the largest refractive index). The streams in the other layers are deduced 
!! from the streams in the most refractive layer to ensure continuity of the stream.
!! Since total reflexion can occur 
!! the number of streams in the other layers is lower or equal to n. 
!! Hence, the number of stream emerging from the surface is generally lower than n.
!!
!! The disortml subroutine computes brightness temperatures. To deduce the emissivity e of the medium, the simplest approach is 
!! to use e=Tb/T, but this is only valid when all the layers have the same temperature T (isothermal medium).
!! For non-isothermal media, it is recommended to call disortml twice, once with Tbatmodown=1K
!! and once with Tbatmodown=0K, all the other input parameters being similar. The difference 
!! between both Tb's is equal to the reflectivity r, i.e. r=1-e if the snowpack is semi-infinite 
!! or the soil is opaque.
!!
!! REFERENCES:
!!
!! Jin, Y-Q., 1994, "Electromagnetic Scattering Modeling for Quantitative Remote Sensing", Singapore: World Scientific.
!!
!! Picard, G., T. Le Toan, S. Quegan, Y. Caraglio, T. Castel, 2004, "Radiative Transfer modeling of cross-polarised backscatter from 
!! a pine forest using the discrete ordinate and eigenvalue method". Transaction on Geoscience 
!! and Remote Sensing, Vol 42, No 8, pp.1720-1730
!!
!!
!! @param[in] l number of layer.
!! @param[in] n number of stream in the most refringent layer.
!! @param[in] depth depth (in meter) of each layer.
!! @param[in] temp temperature (in Kelvin) in each layer.
!! @param[in] Ke  exctinction coefficient (in 1/meter) in each layer. Also called beta in other part of the code.
!! @param[in] albedo single scattering albedo in each layer.
!! @param[in] eps effective dielectric constant in each layer.
!! @param[out] TbV outgoing V-polarized brightness temperature (in Kelvin) for each emerging stream.
!! @param[out] TbH outgoing H-polarized brightness temperature (in Kelvin) for each emerging stream.
!! @param[out] outmhu cosine angle in the air for each emerging stream.
!! @param[in] soilp parameters for the surface under the deepest layer (soil, ice, ... see \ref mod_soil).
!! @param[in] frequency frequency (Hz). This parameters is required for some models of soil roughness and for dielectric constant calculation.
!! @param[in] Tbatmodown atmosphere downwelling brightness temperature (assumed isotropic radiation pattern). Can be used to compute snow emissivity.
!! @param[out] profileV if present and allocated, the computed profile of upwelling V-polarized Tb at the top of each layer (before being transmitted).
!! @param[out] profileH if present and allocated, the computed profile of upwelling H-polarized Tb at the top of each layer (before being transmitted).
!!
!!
!!
subroutine mldisort(l,n,depth,temp,Ke,albedo,eps,TbV,TbH,outmhu,soilp,frequency,Tbatmodown,&
     profileV,profileH)
!
#ifdef MKL
  use mkl95_lapack
#else
#ifdef SUNPERF
  use sunperf
#endif
#endif
!
  use mod_fresnel
  use mod_soil
!
  implicit none
  integer,intent(in)                   :: l           ! number of layers
  integer,intent(in)                   :: n           ! number of stream in the most refringent layer
  real*8,intent(in),dimension(l)       :: depth       ! depth of the layer (first is at the top)
  real*8,intent(in),dimension(l)       :: temp        ! temperature of the layer (first is at the top)
  real*8,intent(in),dimension(l)       :: Ke,albedo   ! effective extinction and single scatt albedo
  complex*16,intent(in),dimension(l)    :: eps         ! effective dielectric constant of the snow
  real*8,intent(in),optional           :: Tbatmodown  ! atmosphere brightness temperature
  real*8,intent(out),dimension(n)      :: TbV,TbH     ! brightness temperature
  real*8,intent(out),dimension(n)      :: outmhu      ! cosine of the angles
  type(soilparams),intent(in),optional :: soilp  ! soil parameters
  real*8,intent(in),optional           :: frequency ! necessary for some soil...
  real*8,intent(out),optional,dimension(:,:) :: profileV,profileH

  !-------------------------------------------------
  ! Internal variables
  integer,dimension(l)          :: ns ! number of stream in each layer 
  integer                       :: nsk,nskm1,nskp1,nsmin,ns0
  real*8,dimension(n,l)         :: mhu,weight ! cosin of all layers

  real*8                        :: tempsoil
  integer                       :: j,k
  integer                       :: nboundary
  integer                       :: roff,coff
  complex*16                    :: E1
  complex*16,parameter          :: cone=cmplx(1.0,0)
  complex*16,dimension(l)       :: eps_

  real*8,dimension(2*n)         :: x
  real*8,dimension(:,:),allocatable :: bvector
#ifdef FULLMAT
  real*8,dimension(:,:),allocatable :: BC
#endif
  real*8,dimension(n)           :: rh,rv

  real*8,dimension(:,:),allocatable :: E,Q         ! eigenvector matrix
  real*8                            :: att

  real*8,dimension(:,:),allocatable :: Tbmat     ! matrix to recover birghtness temperature from bvector

  real*8,dimension(:),allocatable   :: alpha
  integer,dimension(:),allocatable  :: ipiv
  integer                           :: info

  ! for banded matrix
  integer                        :: ku,kl,abbase,abbase1,j1
#ifdef FULLMAT
  integer                        :: i0,i1
  real*8,dimension(:,:),allocatable :: AB
#endif
  real*8,dimension(:,:),allocatable :: AB2


  do k=1,l
     if (depth(k)==0) then 
        if (k==1) then
           eps_(k)=cone
        else
           eps_(k)=eps_(k-1)
        endif
     else
        eps_(k)=eps(k)
     endif
  enddo

  call compute_streams(l,n,eps_,mhu,weight,ns,outmhu,ns0)


  !----
  ! solve the eigenvalue problem for each layer

  nboundary=4*sum(ns)
  !print *,'nboundary=',nboundary
#ifdef FULLMAT
  allocate(BC(nboundary,nboundary))
  BC=0.0
#endif

     ! BANDED MATRIX
  ku=6*n ! max(ns) should be n
  kl=ku
  
  allocate(AB2(2*kl+ku+1,nboundary))
  AB2=0.0
  abbase=kl+ku+1

  allocate(bvector(nboundary,1))
  bvector=0.0

  if (present(profileV).or.present(profileH)) then
     allocate(Tbmat(2*n,nboundary))
  else
     allocate(Tbmat(2*n,4*ns(1))) ! save only the top of the matrix 
  endif

  roff=nboundary ! row offset
  coff=nboundary ! col offset

  do k=l,1,-1
     nsk=ns(k)
     coff=coff - 4*nsk

     allocate(E(2*nsk,2*nsk))
     allocate(Q(2*nsk,2*nsk))
     allocate(alpha(2*nsk))

     !-----------------------------
     ! solve the eigenvalue problem
     call solveonelayer(nsk,Ke(k),albedo(k),mhu(1:nsk,k),weight(1:nsk,k),E,Q,alpha)
   
     !-----------------------------
     ! save the matrix to recover Tb
     if (present(profileV).or.present(profileH).or.k==1) then
        Tbmat(1:2*nsk,coff+1:coff+2*nsk)=E+Q
        do j=1,2*nsk
           Tbmat(1:2*nsk,coff+j+2*nsk)=(E(:,j)-Q(:,j))*exp(-alpha(j)*depth(k))
        enddo
     endif
   
     !----------------------------------
     ! calculate boundary condition. Calculate the system to solve: BCx=b and 

     !-------------------------
     !-- bottom boundary condition

     !--
     ! Condition Idown(k+1)=R Iup(k+1) +  T Idown(k)

       if (k<l) then
 
        nskp1=ns(k+1)
        nsmin=min(nsk,nskp1)

        ! calulate fresnel coefficient
        call fresnel_coefficients_multiangle(eps_(k+1),eps_(k),mhu(1:nsmin,k+1),rv,rh)

        roff=roff - 2*nskp1
        !-- -T(E-Q)D x  Eq 5-10a
        do j=1,2*nsk
           att=exp(-alpha(j)*depth(k))
#ifdef FULLMAT
           BC(1+roff:nsmin+roff,j+coff)=-(1-rv(1:nsmin))*(E(1:nsmin,j)-Q(1:nsmin,j))*att
           BC(nskp1+1+roff:nskp1+nsmin+roff,j+coff)=-(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)-Q(nsk+1:nsk+nsmin,j))*att
#endif
           j1=j+coff
           abbase1=abbase-j1
           AB2(abbase1+1+roff:abbase1+nsmin+roff,j1) = -(1-rv(1:nsmin))*(E(1:nsmin,j)-Q(1:nsmin,j))*att
           AB2(abbase1+nskp1+1+roff:abbase1+nskp1+nsmin+roff,j1)=-(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)-Q(nsk+1:nsk+nsmin,j))*att
        enddo

        !--  -T(E+Q)y  Eq 5-10a
        do j=1,2*nsk
#ifdef FULLMAT
           BC(1+roff:nsmin+roff,j+2*nsk+coff)=-(1-rv(1:nsmin))*(E(1:nsmin,j)+Q(1:nsmin,j))
           BC(nskp1+1+roff:nskp1+nsmin+roff,j+2*nsk+coff)=-(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)+Q(nsk+1:nsk+nsmin,j))
#endif
           j1=j+2*nsk+coff
           abbase1=abbase-j1
           AB2(abbase1+1+roff:abbase1+nsmin+roff,j1)=-(1-rv(1:nsmin))*(E(1:nsmin,j)+Q(1:nsmin,j))  ! V 
           AB2(abbase1+nskp1+1+roff:abbase1+nskp1+nsmin+roff,j1)=-(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)+Q(nsk+1:nsk+nsmin,j)) ! H 
        enddo

        !--
        bvector(1+roff:nsmin+roff,1)=bvector(1+roff:nsmin+roff,1)+(1-rv(1:nsmin))*temp(k) ! V    Eq 5-10a
        bvector(nskp1+1+roff:nskp1+nsmin+roff,1)=bvector(nskp1+1+roff:nskp1+nsmin+roff,1)+(1-rh(1:nsmin))*temp(k) ! H
        !! below the bottom layer.
    endif

     !--
     ! Condition Iup(k)=R Idown(k) +  T Iup(k+1)

     roff=roff - 2*nsk
  
     ! calulate fresnel coefficient
     if (k==l) then
        !! below the bottom layer.
        if ( (.not.present(soilp)) .or. (soilp%imodel==0)) then ! not soil, emissivity=1
           rv=0.0
           rh=0.0
           tempsoil=temp(k) !! take the same temperature as the bottom layer
        else
           call soil_coefficients_multiangle(eps_(k),soilp,frequency,mhu(1:nsk,k),rv,rh)
           tempsoil=soilp%temp
        endif

        bvector(1+roff:nsk+roff,1)=bvector(1+roff:nsk+roff,1)-(1-rv(1:nsk))*tempsoil
        bvector(nsk+1+roff:nsk+nsk+roff,1)= bvector(nsk+1+roff:nsk+nsk+roff,1)-(1-rh(1:nsk))*tempsoil
     else
        call fresnel_coefficients_multiangle(eps_(k),eps_(k+1),mhu(1:nsk,k),rv,rh)
     endif

     !-- -(TE+PQ) D x   Eq 5-10c
     do j=1,2*nsk
        att=exp(-alpha(j)*depth(k)) ! U in eq 5-10c
#ifdef FULLMAT
        BC(1+roff:nsk+roff,j+coff)=-( (1-rv(1:nsk))*E(1:nsk,j) + (1+rv(1:nsk))*Q(1:nsk,j) )*att
        BC(nsk+1+roff:nsk+nsk+roff,j+coff)=-( (1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) + (1+rh(1:nsk)) * Q(nsk+1:nsk+nsk,j) )*att
#endif
        j1=j+coff
        abbase1=abbase-j1
        AB2(abbase1+1+roff:abbase1+nsk+roff,j1)=-( (1-rv(1:nsk))*E(1:nsk,j) + (1+rv(1:nsk))*Q(1:nsk,j) )*att
        AB2(abbase1+nsk+1+roff:abbase1+nsk+nsk+roff,j1)=-( (1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) + (1+rh(1:nsk)) &
             * Q(nsk+1:nsk+nsk,j) )*att
     enddo ! could be optimized if rv and rh are merged

     !-- -(TE-PQ) y   Eq 5-10c        
     do j=1,2*nsk
#ifdef FULLMAT
        BC(1+roff:nsk+roff,j+2*nsk+coff)=-( (1-rv(1:nsk))*E(1:nsk,j) - (1+rv(1:nsk))* Q(1:nsk,j) )
        BC(nsk+1+roff:nsk+nsk+roff,j+2*nsk+coff)=-( (1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) - (1+rh(1:nsk)) * Q(nsk+1:nsk+nsk,j) )
#endif
        j1=j+2*nsk+coff
        abbase1=abbase-j1
        AB2(abbase1+1+roff:abbase1+nsk+roff,j1)=-( (1-rv(1:nsk))*E(1:nsk,j) - (1+rv(1:nsk))* Q(1:nsk,j) )
        AB2(abbase1+nsk+1+roff:abbase1+nsk+nsk+roff,j1)=-( (1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) - (1+rh(1:nsk)) * Q(nsk+1:nsk+nsk,j) )
     enddo

     bvector(1+roff:nsk+roff,1)=bvector(1+roff:nsk+roff,1)+(1-rv(1:nsk))*temp(k)   ! Eq 5-10c        
     bvector(nsk+1+roff:nsk+nsk+roff,1)= bvector(nsk+1+roff:nsk+nsk+roff,1)+(1-rh(1:nsk))*temp(k)

     !-------------------------
     ! upper boundary condition

     !--
     ! Condition Idown(k)=R Iup(k) +  T Idown(k-1)  # bottom of layer k
     
     ! calculate fresnel coefficient
     if (k==1) then
        E1=cone
     else
        E1=eps_(k-1)
     endif
     call fresnel_coefficients_multiangle(eps_(k),E1,mhu(1:nsk,k),rv,rh)

     roff=roff - 2*nsk

     !-- (TE-PQ)x   Eq 5-10a 
     do j=1,2*nsk
#ifdef FULLMAT
        BC(1+roff:nsk+roff,j+coff)=(1-rv(1:nsk))*E(1:nsk,j) - (1+rv(1:nsk))* Q(1:nsk,j)
        BC(nsk+1+roff:nsk+nsk+roff,j+coff)=(1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) - (1+rh(1:nsk)) * Q(nsk+1:nsk+nsk,j)
#endif
        j1=j+coff
        abbase1=abbase-j1
        AB2(abbase1+1+roff:abbase1+nsk+roff,j1)=(1-rv(1:nsk))*E(1:nsk,j) - (1+rv(1:nsk))* Q(1:nsk,j)
        AB2(abbase1+nsk+1+roff:abbase1+nsk+nsk+roff,j1)=(1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) - (1+rh(1:nsk)) * Q(nsk+1:nsk+nsk,j)
     enddo ! could be optimized if rv and rh are merged

     !-- (TE+PQ) U y   Eq 5-10a 
     do j=1,2*nsk
        att=exp(-alpha(j)*depth(k)) ! U in eq 5-10c
#ifdef FULLMAT
        BC(1+roff:nsk+roff,j+2*nsk+coff)=( (1-rv(1:nsk))*E(1:nsk,j) + (1+rv(1:nsk))* Q(1:nsk,j) )*att
        BC(nsk+1+roff:nsk+nsk+roff,j+2*nsk+coff)=( (1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) + (1+rh(1:nsk)) * Q(nsk+1:nsk+nsk,j) )*att
#endif
        j1=j+2*nsk+coff
        abbase1=abbase-j1
        AB2(abbase1+1+roff:abbase1+nsk+roff,j1)=( (1-rv(1:nsk))*E(1:nsk,j) + (1+rv(1:nsk))* Q(1:nsk,j) )*att
        AB2(abbase1+nsk+1+roff:abbase1+nsk+nsk+roff,j1)=( (1-rh(1:nsk))*E(nsk+1:nsk+nsk,j) + (1+rh(1:nsk)) &
             * Q(nsk+1:nsk+nsk,j) )*att
     enddo ! could be optimized if rv and rh are merged

     bvector(1+roff:nsk+roff,1)=bvector(1+roff:nsk+roff,1)-(1-rv(1:nsk))*temp(k) !-- -T*Temp(k)   Eq 5-10a 
     bvector(nsk+1+roff:nsk+nsk+roff,1)=bvector(nsk+1+roff:nsk+nsk+roff,1)-(1-rh(1:nsk))*temp(k)

     !--
     ! Condition Iup(k-1)=R Idown(k-1) +  T Iup(k)  ! Top of layer k
     if (k==1) then
        !! above the top layer.
        if (PRESENT(Tbatmodown)) then
           E1=cone
           ! calulate fresnel coefficient between the air and the first snow interface
           call fresnel_coefficients_multiangle(E1,eps_(k),outmhu,rv,rh)
           bvector(1+roff:ns0+roff,1)=bvector(1+roff:ns0+roff,1)+(1-rv(1:ns0))*Tbatmodown
           bvector(nsk+1+roff:nsk+ns0+roff,1)= bvector(nsk+1+roff:nsk+ns0+roff,1)+(1-rh(1:ns0))*Tbatmodown
        endif
     else ! k>1
        nskm1=ns(k-1)
        nsmin=min(nsk,nskm1)

        ! calulate fresnel coefficient
        call fresnel_coefficients_multiangle(eps_(k-1),eps_(k),mhu(1:nsmin,k-1),rv,rh)

        roff=roff-2*nskm1

        !    T*(E+Q)*x   Eq 5-10c
        do j=1,2*nsk
#ifdef FULLMAT
           BC(1+roff:nsmin+roff,j+coff)=(1-rv(1:nsmin))*(E(1:nsmin,j)+Q(1:nsmin,j))
           BC(nskm1+1+roff:nskm1+nsmin+roff,j+coff)=(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)+Q(nsk+1:nsk+nsmin,j))
#endif
           j1=j+coff
           abbase1=abbase-j1
           AB2(abbase1+1+roff:abbase1+nsmin+roff,j1)=(1-rv(1:nsmin))*(E(1:nsmin,j)+Q(1:nsmin,j))
           AB2(abbase1+nskm1+1+roff:abbase1+nskm1+nsmin+roff,j1)=(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)+Q(nsk+1:nsk+nsmin,j))
        enddo
        !    T*(E-Q)*U*y  Eq 5-10c

        do j=1,2*nsk
           att=exp(-alpha(j)*depth(k)) ! U in eq 5-10c
#ifdef FULLMAT
           BC(1+roff:nsmin+roff,j+2*nsk+coff)=(1-rv(1:nsmin))*(E(1:nsmin,j)-Q(1:nsmin,j))*att
           BC(nskm1+1+roff:nskm1+nsmin+roff,j+2*nsk+coff)=(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)-Q(nsk+1:nsk+nsmin,j))*att
#endif
           j1=j+2*nsk+coff
           abbase1=abbase-j1
           AB2(abbase1+1+roff:abbase1+nsmin+roff,j1)=(1-rv(1:nsmin))*(E(1:nsmin,j)-Q(1:nsmin,j))*att
           AB2(abbase1+nskm1+1+roff:abbase1+nskm1+nsmin+roff,j1)=(1-rh(1:nsmin))*(E(nsk+1:nsk+nsmin,j)-Q(nsk+1:nsk+nsmin,j))*att
        enddo

        !
        bvector(1+roff:nsmin+roff,1)=bvector(1+roff:nsmin+roff,1)-(1-rv(1:nsmin))*temp(k)  ! Eq 5-10c
        bvector(nskm1+1+roff:nskm1+nsmin+roff,1)=bvector(nskm1+1+roff:nskm1+nsmin+roff,1)-(1-rh(1:nsmin))*temp(k)

        roff=roff+2*nskm1
     endif
    roff=roff + 2*nsk

    deallocate(E)
    deallocate(Q)
    deallocate(alpha)

 enddo

 roff=roff - 2*nsk

 if (roff/=0) then
    print *,'probleme interne avec roff=',roff
    stop
 endif
 if (coff/=0) then
    print *,'probleme interne avec coff=',coff
    stop
 endif

 !!-----------------------------------------
 !! solve the boundary system BCx=bvector   (x is indeed x and y as in Y-Q Jin)
 
 allocate(ipiv(nboundary))
 if (.FALSE.) then
#ifdef FULLMAT
    print *,'start solve'
    call getrf(BC,ipiv,info=info)
    if (info/=0) then
       print *,'dgetrf info=',info
       stop
    endif
    call getrs(BC,ipiv,bvector,info=info)
    if (info/=0) then
       print *,'dgetrs info=',info
       stop
    endif
    print *,'end solve'
#endif
  else
     ! BANDED MATRIX
#ifdef FULLMAT
     ku=6*n ! max(ns) should be n
     kl=ku        !
     allocate(AB(2*kl+ku+1, nboundary))
     do j=1,nboundary
        i0=max(1,j-ku)
        i1=min(nboundary,j+kl)
        AB(kl+ku+1+i0-j:kl+ku+1+i1-j,j) = BC(i0:i1,j)
     enddo
#ifdef MKL
     call gbsv( AB, bvector, ipiv=ipiv, info=info )
#else
#ifdef SUNPERF
     call gbsv( nboundary,kl,ku,1,AB, nboundary, ipiv, bvector, nboundary,info=info )
#else
     #error "not implemented with LAPACK. FULLMAT is only used for debugging"
#endif
#endif
     deallocate(AB)
#else
#ifdef MKL
     call gbsv( AB2, bvector, ipiv=ipiv, info=info )
#else
#ifdef SUNPERF
     call gbsv( kl=kl,ku=ku,A=AB2, ipivot= ipiv, B=bvector, info=info )
#else
     call dgbsv(nboundary,kl,ku,1, AB2, 2*kl+ku+1,ipiv,bvector,nboundary,info)
#endif
#endif
     deallocate(AB2)
#endif
     if (info/=0) then
        print *,'#gbsv info=',info
        TbV=0
        TbH=0
        outmhu=0
        return
     endif
  endif
  
  deallocate(ipiv)   

 ! calculate the emerging brightness temperatures

 x=matmul(Tbmat(1:2*nsk,1:4*nsk),bvector(1:4*nsk,1))


 call fresnel_coefficients_multiangle(eps_(1),cone,mhu(1:nsk,1),rv,rh)

 if (present(Tbatmodown)) then
    TbV(1:nsk) = (1-rv(1:nsk))*(temp(1)+x(1:nsk)) + rv(1:nsk)*Tbatmodown
    TbH(1:nsk) = (1-rh(1:nsk))*(temp(1)+x(nsk+1:2*nsk)) + rh(1:nsk)*Tbatmodown 
 else
    TbV(1:nsk) = (1-rv(1:nsk))*(temp(1)+x(1:nsk))
    TbH(1:nsk) = (1-rh(1:nsk))*(temp(1)+x(nsk+1:2*nsk))
 endif

 if (nsk<n) then
    TbV(nsk+1:n)=0
    TbH(nsk+1:n)=0
 endif

 if (present(profileV).or.present(profileH)) then

    ! Allocation is done by the calling subroutine.
    !if (present(profileV).and.(.not.(allocated(profileV)))) allocate(profileV(n,l))
    !if (present(profileH).and.(.not.(allocated(profileH)))) allocate(profileH(n,l))

    coff=0
    do k=1,l
       nsk=ns(k)
       x=matmul(Tbmat(1:2*nsk,1+coff:4*nsk+coff),bvector(1+coff:4*nsk+coff,1))
       coff=coff+4*nsk

       if (present(profileV)) profileV(1:nsk,k)=temp(k)+x(1:nsk)
       if (present(profileH)) profileH(1:nsk,k)=temp(k)+x(nsk+1:2*nsk)
    enddo

 endif
#ifdef FULLMAT
 deallocate(BC)
#endif
 deallocate(bvector)
 deallocate(Tbmat)

end subroutine mldisort


!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine solveonelayer(n,Ke,albedo,mhu,weight,E,Q,alpha)

#ifdef MKL
  use mkl95_lapack
#else
#ifdef SUNPERF
  use sunperf
#endif
#endif
  implicit none

#ifndef MKL
  integer :: ihi,ilo
  real*8 :: abnrm
  real*8,dimension(2*n) :: scale, rcone, rconv
  real*8,dimension(2*n,2*n) :: vl
  real*8,dimension(5*2*n) :: work
  real*8,dimension(1) :: iwork
#endif

  !! solve the homogeneous equation for a single layer.
  !! return E,Q and lambdas.

  integer,intent(in)                    :: n ! number of stream
  real*8,intent(in)                     :: Ke,albedo
  real*8,intent(in),dimension(n)        :: mhu,weight
  real*8,intent(out),dimension(2*n,2*n) :: E,Q
  real*8,intent(out),dimension(2*n)     :: alpha

  !! ------------------------------------------------
  real*8                                :: coef
  real*8,dimension(2*n,2*n)             :: A,TD         ! matrix
  real*8,dimension(n)                   :: mhu2
  real*8,dimension(2*n)                 :: lambda_imag
  integer                               :: info,i,j


  !----
  ! calculate the A matrix
  coef=2* 3.0/8.0*Ke*albedo  ! 3/8 comes from Y-Q Jin page 20
                             ! 2 comes from the addition F+B Y-Q Jin page 100

  mhu2=mhu**2

  do j=1,n
     A(1:n,j)=       coef * ( 2*(1-mhu2(:))*(1-mhu2(j)) + mhu2(:)*mhu2(j) )*weight(j)
     A(1:n,j+n)=     coef *  mhu2(:) * weight(j)  
     A(n+1:2*n,j)=   coef *  mhu2(j) * weight(j)
     A(n+1:2*n,j+n)= coef *            weight(j) 
  enddo

  do i=1,2*n
     A(i,i)=A(i,i) - Ke
  enddo

  ! calculate the matrix to diagonalize (TD=mhu^-1 W mhu^-1 A) 
  ! here W is constant

  do i=1,n
     TD(i,:)    =A(i,:)    * (-Ke)/mhu2(i)   ! -Ke c'est W
     TD(i+n,:)  =A(i+n,:)  * (-Ke)/mhu2(i)   ! -Ke c'est W
  enddo

  ! diagonalise the matrix
#ifdef MKL
  call geevx(TD,alpha,lambda_imag,vr=E,info=info,balanc='B') ! here alpha is indeed lambda_real
#else
#ifdef SUNPERF
  call geevx('N','N','V','N',a=TD,wr=alpha,wi=lambda_imag,vl=vl,vr=E,   &
       ilo=ilo,ihi=ihi,scale=scale, abnrm=abnrm,rcone=rcone, rconv=rconv,work=work,ldwork=size(work),info=info) ! here alpha is indeed lambda_real
#else
  ! basic lapack (work with MKL also)
  call dgeevx('N','N','V','N',2*n,TD,2*n, alpha,lambda_imag, vl,2*n, E, 2*n, &
       ilo,ihi,scale,abnrm,rcone,rconv,work,size(work),iwork,info)
#endif
#endif
    if (info/=0) then
     print *,'diagonalisation error. info=',info
     stop
  end if
  alpha=sqrt(alpha) ! now, alpha is really the alpha found in Y-Q Jin

  do j=1,2*n !! "pour l'esthetique"
     if (E(1,j)<0) E(:,j)=-E(:,j)
  enddo

  ! calculate Q

  Q=matmul(A,E)

  ! right multiplication by 1/sqrt(lambda)
  do j=1,2*n
     Q(:,j)=Q(:,j)/alpha(j)
  enddo
  ! left multiplication by 1/mhu
  do j=1,2*n
     Q(1:n,j)=Q(1:n,j)/mhu(:)
     Q(n+1:2*n,j)=Q(n+1:2*n,j)/mhu(:)
  enddo

end subroutine solveonelayer
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!


!> Compute the streams angles in the different layers given the number of streams n
!! in the most refringent layer.
!! The streams in the other layers are deduced from the streams in the most refractive 
!! layer to ensure continuity of the stream. Since total reflexion can occur the number of
!! streams in the other layers is lower or equal to n. Hence, the number of stream 
!! emerging from the surface is generally lower than n.
!!
!! @param[in] l number of layer.
!! @param[in] n number of stream in the most refringent layer.
!! @param[in] eps effective dielectric constant in each layer.
!! @param[out] mhu cosine angles in every layer.
!! @param[out] weight of the quadrature associated with every stream in every layer.
!! @param[out] ns number of streams in every layer.
!! @param[out] outmhu cosine angles in the air for each emerging stream.
!! @param[out] ns0 number of streams in the air.


subroutine compute_streams(l,n,eps,mhu,weight,ns,outmhu,ns0)

  integer,intent(in)                   :: l           ! number of layers
  integer,intent(in)                   :: n           ! number of stream in the most refring
  complex*16,intent(in),dimension(l)    :: eps         ! effective dielectric constant of the

  real*8,intent(out),dimension(n,l)    :: mhu,weight ! cosin of all layers
  integer,intent(out),dimension(l)     :: ns ! number of stream in each layer 
  real*8,intent(out),dimension(n)      :: outmhu      ! cosine of the angles
  integer,intent(out)                  :: ns0
  !-------------------------------------------------
  integer                       :: kmostrefringent,k,nsk,i
  integer,dimension(1)          :: kmostrefringent_arr
  real*8                        :: relindex,relsin
  !-------------------------------------------------


  !----
  ! search the most refringent layer
  kmostrefringent_arr=maxloc(real(sqrt(eps)))
  kmostrefringent=kmostrefringent_arr(1)
  !----
  ! calculate the gaussian weights and nodes for the most refringent layer
  call gaussquad(n,mhu(:,kmostrefringent),weight(:,kmostrefringent))

  ns(kmostrefringent)=n

  !print *,'#kmostrefringent=',kmostrefringent

  !----
  ! calculate the nodes and weights of all the other layers
  do k=1,l
     if (k/=kmostrefringent) then
        relindex=real(sqrt(eps(kmostrefringent)/eps(k))) 
        !! real is an approximation. See for instance "Semiconductor Physics, Quantum Electronics & Optoelectronics. 2001. V. 4, N 3. P. 214-218."

        !! calculate the angles (=node)
        do i=1,n
           relsin=relindex*sqrt(1-mhu(i,kmostrefringent)**2)
           if (relsin<1.0) then
              mhu(i,k)=sqrt(1-relsin**2)
              nsk=i ! nsk is the max i
           else
	      mhu(i:n,k)=0.0
              exit
           endif
        enddo
        ns(k)=nsk
        !! calculate the weight ("a" in Y-Q Jin)
        weight(1,k)=1-0.5*(mhu(1,k)+mhu(2,k))
        weight(nsk,k)=0.5*(mhu(nsk-1,k)+mhu(nsk,k))
        weight(2:nsk-1,k)=0.5*(mhu(1:nsk-2,k)-mhu(3:nsk,k))
     endif
  enddo
  !print *,'ns=',ns

  relindex=real(sqrt(eps(1)/1.0))
  !! calculate the angles (=node) in the air
  do i=1,n
     relsin=relindex*sqrt(1-mhu(i,1)**2)
     if (relsin<1.0) then
        outmhu(i)=sqrt(1-relsin**2)
        ns0=i
     else
        outmhu(i:n)=0.0
	exit
     endif
  enddo
end subroutine compute_streams


subroutine gaussquad(n,mhu,weight)
  implicit none
  integer,intent(in)                  :: n
  real*8,dimension(n),intent(out)     :: mhu,weight
  !-------------------------
  real*8,dimension(2*n+1)             :: mhu0,weight0
  real*8,dimension(2*n+1)             :: b,b2,d  ! work space
  !-------------------------

  ! calculate the gaussian weights and nodes

  call czerg(mhu0,weight0,2*n,2*n+1,b,b2,d)

  ! reverse the order (needed by mldisort)
  mhu(:)=mhu0(2*n:n+1:-1) ! czerg calculate all the root even the negative ones... we only need the positive ones
  weight(:)=weight0(2*n:n+1:-1)

end subroutine gaussquad


!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine savemat(filename,mat)
  implicit none
  character(*),intent(in)          :: filename
  real*8,intent(in),dimension(:,:) :: mat
  !----------------------
  integer :: i,j,k
  !----------------------
  

  print *,'#size=',size(mat,1),size(mat,2)
  open(10,file=filename,form='unformatted',access='direct',status='unknown',recl=8)
  k=1

  do i=1,size(mat,1)
     do j=1,size(mat,2)
        write (10,rec=k) mat(i,j)
        k=k+1
     enddo
  enddo
  close(10)

  open(10,file=trim(filename)//'.hdr',status='unknown')
  write (10,'(A)') 'ENVI'
  write (10,'(A)') 'description = {'
  write (10,'(A)') '  File Imported into ENVI.}'
  write (10,'(A,I5)') 'samples = ',size(mat,2)
  write (10,'(A,I5)') 'lines   = ',size(mat,1)
  write (10,'(A)') 'bands   = 1'
  write (10,'(A)') 'header offset = 0'
  write (10,'(A)') 'file type = ENVI Standard'
  write (10,'(A)') 'data type = 5'
  write (10,'(A)') 'interleave = bsq'
  write (10,'(A)') 'sensor type = Unknown'
  write (10,'(A)') 'byte order = 0'
  write (10,'(A)') 'wavelength units = Unknown'
  close (10)


end subroutine savemat

end module mod_disort
