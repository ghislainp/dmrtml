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
    for d in arange(0,917,10,dtype=float32):
        density=array([d])
        #albedo,beta=dmrtml.albedobeta(freq,density,radius*1e-6,temp,grodyapproach=False,medium=medium)
        albedo,beta=dmrtml.albedobeta(freq,density,radius*1e-6,temp,grodyapproach=True,medium=medium)
        x.append(d)
        scatt.append(beta*albedo)  #   /temp[0]) for emissivity calculation
        abso.append(beta*(1-albedo))
    return array(x),array(scatt),array(abso)

### SNOW
x,scatt,abso=coefficients('S')
mask=x < 917/2
plot(x[mask],scatt[mask],'b-',alpha=0.7)
plot(x[mask],abso[mask],'r-',alpha=0.7)
mask=x < 0.3 * 917
plot(x[mask][::3],scatt[mask][::3],'bo')
plot(x[mask][::3],abso[mask][::3],'ro',markerfacecolor='w',markeredgecolor='r')

# tangente
newx=arange(0,200)
newy=newx*(scatt[2]-scatt[1])/(x[2]-x[1])
plot(newx,newy,'b--',alpha=0.7)

### ICE
x,scatt,abso=coefficients('I')
mask=x > 917/2
plot(x[mask],scatt[mask],'b-',alpha=0.7)
plot(x[mask],abso[mask],'r-',alpha=0.7)
mask=x > (1 - 0.3) * 917
plot(x[mask][::3],scatt[mask][::3],'bs')
plot(x[mask][::3],abso[mask][::3],'rs',markerfacecolor='w',markeredgecolor='r')


### Auto
x,scatt,abso=coefficients('A')
plot(x,scatt,'b--',alpha=0.7)
plot(x,abso,'r--',alpha=0.7)

### BRIDGING
x,scatt,abso=coefficients('B')
plot(x,scatt,'b-.',lw=4,alpha=0.7)
plot(x,abso,'r-.',lw=4,alpha=0.7)


xlabel('Density (kg m$^{-3}$)')
ylabel('Scattering and absorption coefficients (m$^{-1}$)')

gca().xaxis.set_major_locator(MaxNLocator(10))

legend(loc=4)
show()



