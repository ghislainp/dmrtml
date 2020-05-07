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


freq = 19.0e9                       # frequency: Hz 
height = array([0.10, 100.0])       # height: m
temp = array([273.0, 273])          # temperature: K
density = array([350.0, 350.0])     # density: kg.m-3
dist = False                        # if True => use RAYLEIGH distribution of particles

soilp = None

for r in arange(100, 2000, 400):
    radius = array([r, r])*1e-6
    x, y = list(),list()
    for w in arange(0, 10, 0.5):
        fwetness = (array([w, 0]) / height) /density # m3/m3
        res = dmrtml.dmrtml(freq, 64, height, density, radius, temp,
                            tau=dmrtml.NONSTICKY, fwetness=fwetness, dist=dist, soilp=soilp)
        
        x.append(w)
        y.append(res.TbV(53)) #/temp[0]) for emissivity calculation
    
    plot(array(x), y, label='radius %i microns' % r)


xlabel('Liquid Water Content (m3/m3)')
ylabel('Brightness temperature (K)')
legend()
show()



