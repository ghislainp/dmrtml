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


import dmrtml
from pylab import *


freq = 37.0e9                 # frequency unity: Hz 
height = array([10.0])        # height unity: m
temp = array([255.0])         # temperature unity: K
radius = array([100.0e-6])    # radius unity: meter
density = array([400.0])      # density unity: kg.m-3
dist = False                  # if True => use RAYLEIGH distribution of particles

medium='S' # SNOW
               
soilp = dmrtml.HUTRoughSoilParams(273) # other parameters have their default value
res = dmrtml.dmrtml(freq,64,height,density,radius,temp,tau=dmrtml.NONSTICKY,medium=medium,dist=dist,soilp=soilp)


theta = range(0,75)
plot(theta,res.TbV(theta),'bo-',label='Vertical polarization')
plot(theta,res.TbH(theta),'go-',label='Horizontal polarization')
xlabel('Incidence angle (Degree)')
ylabel ('Brightness temperature (K)')
legend(loc=8)
show()



