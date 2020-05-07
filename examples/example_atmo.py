#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import sys
import dmrtml
from pylab import *


freq = 37.0e9                 # frequency unity: Hz 
height = array([10.0])        # height unity: m
temp = array([250.0])         # temperature unity: K
radius = array([800.0e-6])    # radius unity: meter
density = array([200.0])      # density unity: kg.m-3
dist = False                  # if True => use RAYLEIGH distribution of particles


        
        
tbatmodown=1
        
soilp = None
res_atmo = dmrtml.dmrtml(freq,64,height,density,radius,temp,
                         tau=dmrtml.NONSTICKY,dist=dist,soilp=soilp,tbatmodown=tbatmodown)

theta = range(0,75)

res = dmrtml.dmrtml(freq,64,height,density,radius,temp,
                    tau=dmrtml.NONSTICKY,dist=dist,soilp=soilp,tbatmodown=0)

rv= (res_atmo.TbV(theta) - res.TbV(theta))/tbatmodown
plot (theta,rv,'bo-',label='V-pol')
plot (theta,1-res.TbV(theta)/mean(temp))

rh= (res_atmo.TbH(theta) - res.TbH(theta))/tbatmodown
plot (theta,rh,'ro-',label='H-pol')
plot (theta,1-res.TbH(theta)/mean(temp))


xlabel('Incidence angle (Degree)')
ylabel('Brightness temperature (K)')
legend(loc=2)
show()



