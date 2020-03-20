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
sys.path.append('..')      # add where dmrtml.py is located
import dmrtml
from pylab import *


fig = figure(figsize=(7, 5))
fig.set_facecolor('w')

freq = 19.0e9                       # frequency: Hz 
height = array([0.10, 100.0])       # height: m
temp = array([272.0, 272.0])        # temperature: K
density = array([300.0, 300.0])     # density: kg.m-3
dist = False                        # if True => use RAYLEIGH distribution of particles
soilp = None

for r, s in [(500, 'o-'), (1000, 's-')]:
    radius = array([r, r])*1e-6
    x,y = list(), list()
    for w in arange(0, 2, 0.02):
        #lwc = array([w/100.0, 0])
        fwetness = array([w/height[0], 0])/density
        res = dmrtml.dmrtml(freq, 64, height, density, radius, temp,
                            tau=dmrtml.NONSTICKY, fwetness=fwetness, dist=dist, soilp=soilp)
        
        #x.append(sum(lwc*height)*1e3)
        x.append(sum(fwetness*height*density))
        y.append(res.TbH(53)) #/temp[0]) for emissivity calculation
    
    plot(array(x), y, s, label='radius %3.1f mm' % (r/1000.0), markevery=10)


xlabel('Amount of Liquid Water (kg m$^{-2}$)')
ylabel('Brightness temperature (K)')
legend(loc='lower right')
show()



