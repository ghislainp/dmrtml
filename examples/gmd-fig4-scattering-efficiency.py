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


freq=89e9      # frequency in Hz
radius=0.25e-3 # radius in m
temp=260       # temperature in m

fig=figure(figsize=(9,5))
fig.set_facecolor('w')

shift=0
for (density,color,label) in [(200,'r','DMRT 200 kg m$^{-3}$'),
                              (300,'k','DMRT 300 kg m$^{-3}$'),
                              (1,'b','DMRT 0 kg m$^{-3}$ (indep.)')]:
    f=density/917.0

    radii=arange(0,3500.0,10.0)

    x=[]
    y=[]
    for radius in radii:
        a=radius*1e-6
        albedo,beta=dmrtml.albedobeta(freq,density,a,temp,grodyapproach=False)
        x.append(radius*1e-3)

        # compute efficiency from coefficient
        Qe=(beta*a*4) / (3*f)
        Qs=albedo*Qe
        y.append(Qs)

    radii=arange(0,3500.0,100.0)+shift
    shift+=33

    xgrody=[]
    ygrody=[]
    yabsgrody=[]
    for radius in radii:
        a=radius*1e-6
        albedo,beta=dmrtml.albedobeta(freq,density,a,temp) #,tau=0.3)
        xgrody.append(radius*1e-3)

        Qe=(beta*a*4) / (3*f)
        Qs=albedo*Qe
        ygrody.append(Qs)
        yabsgrody.append(Qe-Qs)

    x=array(x)
    y=array(y)

    xo=400.0e-3

    plot(x,y,'-',linewidth=2,color=color,label=label)
    if density>50:
        plot(xgrody,ygrody,':o',lw=1,color=color,alpha=0.6)


data=loadtxt("gmd-fig4-mie.out")
plot(data[:,0]*1e-3,data[:,1],'--',linewidth=2,label='Mie (indep.)')
xlim((0,3.500))
ylim((0,6))
xlabel('Grain radius (mm)')
ylabel('Scattering efficiency')
legend(loc=1,frameon=False)
show()
