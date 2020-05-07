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
from scipy.optimize.minpack import leastsq


temp=260

fig=figure(figsize=(7,5))
fig.set_facecolor('w')

shift=0
for (density,color,label) in [(200,'r','DMRT 200 kg m$^{-3}$'),
                              (300,'k','DMRT 300 kg m$^{-3}$')]:

    x=[]
    y=[]

    def efficiency(params):
        a=params[0]*1e-6
        albedo,beta=dmrtml.albedobeta(freq*1e9,density,a,temp,grodyapproach=False)
        f=density/917.0
        Qe=(beta*a*4) / (3*f)
        Qs=albedo*Qe
        return Qs-2

    for freq in range(6,200,5):
        params=[100.0]
        params,it=leastsq(efficiency,params)
        x.append(freq)
        y.append(params[0]*1e-3)

    plot(x,y,c=color, label=label,lw=2)
    fill_between(array(x),zeros_like(y),array(y),alpha=0.15,color='k')

xlabel('Frequency (GHz) ')
ylabel('Maximum grain radius (mm)')
ylim((0,2))
gca().xaxis.set_major_locator(MaxNLocator(10))

legend(loc=1)
show()
