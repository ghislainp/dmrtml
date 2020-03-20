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


freq=89e9
radius=0.25e-3
temp=260

fig=figure(figsize=(9,5))
fig.set_facecolor('w')


ax2=subplot(1,1,1)

shift=0
for (density,color,label) in [(200,'r',' 200 kg m$^{-3}$'),
                              (300,'k',' 300 kg m$^{-3}$')]:
    f=density/917.0

    radii=arange(0,1000.0,30.0)

    x=[]
    y2v=[]
    y2h=[]

    for radius in radii:
        a=radius*1e-6
        x.append(radius*1e-3)

        res=dmrtml.dmrtml(freq,64,[1000],density,a,temp)

        y2v.append(res.TbV(53))
        y2h.append(res.TbH(53))

    if density==200:
        s=''
    else:
        s='o'
    ax2.plot(x,y2v,s+'-',linewidth=2,color=color,markevery=5,label='V-pol '+label)
    ax2.plot(x,y2h,s+'--',linewidth=2,color=color,markevery=5,label='H-pol '+label)

    
ax2.set_xlabel('Grain radius (mm)')
ax2.set_ylabel('Brightness temperature (K)')

ax2.legend(loc=1)


show()
