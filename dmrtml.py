#----------------------------------------------------------------------------#
#									     #
# DMRT-ML (Dense Media Radiative Transfer - Multi-Layer)		     #
# Copyright (c), all rights reserved, 2007-2012,  Ghislain Picard	     #
# email: Ghislain.Picard@ujf-grenoble.fr				     #
#									     #
# Main contributors: Ludovic Brucker, Alexandre Roy, Florent Dupont	     #
#									     #
#									     #
# Licensed under the terms of the GNU General Public License version 3:	     #
# http://www.opensource.org/licenses/gpl-3.0.html 			     #
#									     #
# dmrtml'url: http://lgge.osug.fr/~picard/dmrtml/                            #
#                                                                            #
# Recommended citation:							     #
#    G. Picard, L. Brucker, A. Roy, F. Dupont, M. Fily, and A. Royer,        #
#    Simulation of the microwave emission of multi-layered snowpacks using   #
#    the Dense Media Radiative Transfer theory: the DMRT-ML model,           #
#    Geoscientific Model Development, 6, 1061-1078, 2013                     #
#    doi:10.5194/gmd-6-1061-2013                                             #
#    http://www.geosci-model-dev.net/6/1061/2013/gmd-6-1061-2013.html        #
#									     #
#----------------------------------------------------------------------------#


import dmrtml_for
import numpy

NONSTICKY=1e6

#
# Main routine to call dmrtml
#
def dmrtml(freq,n,depth,density,radius,temp,tau=NONSTICKY,fwetness=0,medium='S',dist=False,soilp=None,tbatmodown=0,eps_ice=(0.0,0.0)):
    """run dmrtml and return a DMRTMLResult object"""
    if soilp is None:
        soilp=SoilParams()

    depth,density,radius,temp,tau,fwetness,medium=_clean_parameters(depth,density,radius,temp,tau,fwetness,medium)

    restuple=dmrtml_for.dmrtml_pywrapper(freq,n,depth,density,radius,temp,tau,fwetness,medium,dist,
                                       soilp.imodel,soilp.temp,soilp.eps[0],soilp.eps[1],
                                       soilp.sigma,soilp.SM,soilp.sand,soilp.clay,soilp.dm_rho,
                                       soilp.Q, soilp.N,
                                       tbatmodown,
                                       eps_ice[0],eps_ice[1])
    
    return DMRTMLResult(restuple)

#
# class returned by the dmrtml routine
#
class DMRTMLResult:
    """results from the DMRTML model"""
    def __init__(self,restuple):
        """create a DMRTMLResult from a tuple returned by dmrtml_wrapper"""
        self.TbVarr,self.TbHarr,self.mhu=restuple
        mask=self.mhu>0
        self.TbVarr=self.TbVarr[mask]
        self.TbHarr=self.TbHarr[mask]
        self.mhu=self.mhu[mask]

    def TbH(self,theta=None):
        """return TbH for the angle(s) theta)"""
        if theta is None:
            return self.TbHarr
        else:
            return numpy.interp(theta, self.thetas(), self.TbHarr)

    def TbV(self,theta=None):
        """return TbV for the angle(s) theta)"""
        if theta is None:
            return self.TbVarr
        else:
            return numpy.interp(theta, self.thetas(), self.TbVarr)

    def thetas(self):
        """return the actual zenith angles of the streams used for the calculation"""
        return numpy.arccos(self.mhu)*180/numpy.pi



def albedobeta(frequency,density,radius,temp,tau=NONSTICKY,fwetness=0,medium='S',dist=False,grodyapproach=True):
    """compute the albedo and extinction with the DMRT theory"""
    return dmrtml_for.albedobeta_pywrapper(frequency,density,radius,temp,
                                           tau,fwetness,medium,dist,grodyapproach)


def compute_streams(freq,n,density,radius=0,temp=273,
                    tau=NONSTICKY,fwetness=0,medium='S',dist=False):
    """compute the streams angles used by DISORT"""

    depth=numpy.ones_like(density)
    depth,density,radius,temp,tau,fwetness,medium=_clean_parameters(depth,density,radius,temp,tau,fwetness,medium)

    restuple=dmrtml_for.compute_streams_pywrapper(freq,n,density,radius,temp,tau,fwetness,medium,dist)
    return Streams(restuple) #! return mhu,weight,ns,outmhu,ns0

class Streams:
    def __init__(self,restuple):
        """create a Streams from a tuple returned by compute_streams_pywrapper"""
        self.mhu,self.weight,self.ns,self.outmhu,self.ns0=restuple

def ice_dielectric_constant(frequency,temperature):
    return dmrtml_for.icedielectric_pywrapper(frequency,temperature)





class SoilParams:
    """Soil parameters"""
    def __init__(self):
        self.imodel = 0 # no soil
        self.temp   = 273
        self.eps    = 1,0
        self.sigma  = 0.005
        self.SM     = 0.2
        self.sand   = 0.4
        self.clay   = 0.3
        self.dm_rho = 1100.0

        self.Q  = 0
        self.N  = 0

    def set_temperature(self,temp):
        self.temp = temp

    def set_soilmoisture(self,SM):
        if not SM is None: self.SM = SM

    def set_soiltexture(self,sand,clay,drymatter_density):
        if not sand is None: self.sand = sand
        if not clay is None: self.clay = clay
        if not drymatter_density is None: self.dm_rho = drymatter_density

    def set_roughness(self,sigma):
        if not sigma is None: self.sigma = sigma

    def set_QNH(self, q, n, h):
        if not q is None: self.Q = q
        if not n is None: self.N = n
        if not h is None: self.sigma = h

    def set_dielectricconstant(self, epsreal, epsimag=None):
        try:
            self.eps = epsreal.real, epsreal.imag
        except:
            if epsimag is not None:
                self.eps = epsreal,epsimag
            else:
                raise Exception("invalid value type for epsreal and/or epsimag")



class NoSoilParams(SoilParams):
    """No soil"""
    def __init__(self , temperature=273):
        SoilParams.__init__(self)
        self.imodel = 0        # imodel=0 no soil (rh=rv=0)
        self.set_temperature(temperature)


class FlatSoilParams(SoilParams):
    """Flat soil"""
    def __init__(self, temperature=273, eps="epspulliainen",sand=None,clay=None,drymatter_density=None, SM=None):
        SoilParams.__init__(self)
        self.set_temperature(temperature)

        if eps=="epspulliainen" or eps=="epsmodel":
            self.imodel = 2 # imodel=2 flat surface, fresnel coefficient with ESS epsilon
            self.set_soilmoisture(SM)
            self.set_soiltexture(sand,clay,drymatter_density)
            self.eps = 0,0

        elif eps=="epsdobson":
            self.imodel = 3 # imodel=3 flat surface, fresnel coefficient with Dobson epsilon
            self.set_soilmoisture(SM)
            self.set_soiltexture(sand,clay,drymatter_density)
            self.eps = 0,0

        elif eps=="epsmironov":
            self.imodel = 4 # imodel=4 flat surface, fresnel coefficient with Mironov epsilon
            self.set_soilmoisture(SM)
            self.set_soiltexture(sand,clay,drymatter_density)
            self.eps = 0,0

        else:
            self.imodel = 1 # imodel=1 flat surface, fresnel coefficient with prescribed eps
            self.eps = eps

class HUTRoughSoilParams(FlatSoilParams):
    """Rough soil"""
    def __init__(self, temperature=273, eps="epspulliainen",sand=None,clay=None,drymatter_density=None, SM=None,sigma=None):
        FlatSoilParams.__init__(self,temperature,eps,sand,clay,drymatter_density,SM)
        #        print "HUTRoughSoilParams",self.eps
        # imodel=101 HUT roughsoil reflectivity for rough surface with prescribed eps
        # imodel=102 HUT roughsoil reflectivity for rough surface with EPSS eps
        # imodel=103 HUT roughsoil reflectivity for rough surface with Dobson eps
        # imodel=104 HUT roughsoil reflectivity for rough surface with Mironov eps
        self.imodel += 100
        self.set_roughness(sigma)


class QNHRoughSoilParams(FlatSoilParams):
    """Rough soil"""
    def __init__(self, temperature=273, eps="epspulliainen",sand=None,clay=None,drymatter_density=None, SM=None,Q=None, N=None, H=None):
        FlatSoilParams.__init__(self,temperature, eps,sand,clay,drymatter_density, SM)

        # imodel=301 HUT roughsoil reflectivity for rough surface with prescribed eps
        # imodel=302 HUT roughsoil reflectivity for rough surface with EPSS eps

        # imodel=303 HUT roughsoil reflectivity for rough surface with Dobson eps
        # imodel=304 HUT roughsoil reflectivity for rough surface with Mironov eps
        self.set_QNH(Q,N,H)

        self.imodel += 300



class FlatIceParams(SoilParams):
    """Flat ice"""
    def __init__(self, temperature=273, eps="epsmodel"):

        SoilParams.__init__(self)
        if eps=="epsmodel":
            self.imodel = 202 # imodel=202 ice flat surface, fresnel coefficient with eps model
        else:
            self.imodel = 201 # imodel=201 ice flat surface, fresnel coefficient with prescribed eps (same as imodel=1)
            self.eps = eps




def _clean_parameters(depth,density,radius,temp,tau,fwetness,medium):

    depth=numpy.asfarray(depth)

    density=numpy.asfarray(density)
    if density.size==1:
        density=density*numpy.ones(depth.size)    

    radius=numpy.asfarray(radius)
    if radius.size==1:
        radius=radius*numpy.ones(depth.size)    

    temp=numpy.asfarray(temp)
    if temp.size==1:
        temp=temp*numpy.ones(depth.size)    

    tau=numpy.asfarray(tau)
    if tau.size==1:
        tau=tau*numpy.ones(depth.size)

    fwetness=numpy.asfarray(fwetness)
    if fwetness.size==1:
        fwetness=fwetness*numpy.ones(depth.size)

    medium=numpy.asarray(medium)
    if medium.size==1:
        m = numpy.chararray(depth.size)
        m[:]=medium
        medium=m

    return depth,density,radius,temp,tau,fwetness,medium
