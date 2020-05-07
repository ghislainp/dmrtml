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
import dmrtml
from pylab import *

fig=figure(figsize=(7,5))
fig.set_facecolor('w')


freq = 37.0e9                 # frequency unity: Hz 
height = array([100.0])        # height unity: m
temp = array([260.0])         # temperature unity: K
radius = array([300.0])       # radius unity: micrometer
dist = False                  # if True => use RAYLEIGH distribution of particles


               
soilp = None
hold(True)

def coefficients(medium):
    x=[]
    scatt=[]
    abso=[]
    for d in arange(1,917,10):
        density=array([d])
        albedo,beta=dmrtml.albedobeta(freq,density,radius*1e-6,temp,grodyapproach=False,medium=medium,tau=tau)
        x.append(d)
        scatt.append(beta*albedo)  #/temp[0]) for emissivity calculation
        abso.append(beta*(1-albedo))
    return array(x),array(scatt),array(abso)



for (tau,color,label) in [(0.3,'go-','$\\tau$=0.3'),(1,'ks-','$\\tau$=1.0'),(1000,'b^-','non-sticky')]:
    x,scatt,abso=coefficients('S')
    mask=x < 917/2
    plot(x[mask],scatt[mask],color,alpha=0.7,label='%s, $a$=%g mm' % (label,radius[0]/1000),markevery=5)


tau=100
radius=array([400.0])
x,scatt,abso=coefficients('S')
mask=x < 917/2
plot(x[mask],scatt[mask],'rd-',alpha=0.7,label='%s, $a$=%g mm' % (label,radius[0]/1000),markevery=5)
ylim((0,1))


xlabel('Density (kg m$^{-3}$)')
ylabel('Scattering coefficient (m$^{-1}$)')

legend(loc='upper right')
show()



