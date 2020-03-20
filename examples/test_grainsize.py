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
sys.path.append('..')      # add where dmrtml.py is located
import dmrtml
from pylab import *


freq = 37.0e9                 # frequency : Hz 
height = array([100.0])       # height : m
temp = array([250.0])         # temperature : K
density = array([300.0])      # density : kg.m-3
dist = False                  # if True => use RAYLEIGH distribution of particles
soilp = None


x,y=list(),list()
for r in arange(100,1000,50):
    radius=array([r])*1e-6
    res = dmrtml.dmrtml(freq,64,height,density,radius,temp,
                        tau=dmrtml.NONSTICKY,dist=dist,soilp=soilp)
    
    x.append(r)
    y.append(res.TbV(53))   #/temp[0]) # divide for emissivity instead of Tb

plot(x,y)
xlabel('Grain size (microns)')
ylabel('Brightness temperature (K)')
legend()
show()



