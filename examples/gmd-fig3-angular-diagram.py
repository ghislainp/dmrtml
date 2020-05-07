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


freq = 18e9                 # frequency unity: Hz 
height = array([100.0])     # height unity: m
temp = array([272.0])       # temperature unity: K
radius = array([1750e-6])   # radius unity: meter
density = array([275.1])    # density unity: kg.m-3

def original():
    files=[('gmd-fig3-data-v.csv','bo',None),
           ('gmd-fig3-data-h.csv','go',None),
           ('gmd-fig3-simul-v.csv','b:',None),
           ('gmd-fig3-simul-h.csv','g:',None)]
    for f in files:
        filename,style,label=f
        data=loadtxt(filename,delimiter=',',skiprows=1)
        plot(data[:,0],data[:,1],style,label=label)

fig=figure(figsize=(11,5),dpi=150)
fig.set_facecolor('w')
subplot(1,2,1)

original()

res = dmrtml.dmrtml(freq,64,height,density,radius,temp,eps_ice=(3.2,0.016))
theta = range(0,85)
plot(theta,res.TbV(theta),'b-',label='V-pol')
plot(theta,res.TbH(theta),'g-',label='H-pol')


xlabel('Incidence angle (Degree)')
ylabel ('Brightness temperature (K)')
ylim((120,250))
legend(loc=6,frameon=False)

graphletter={ 'xy': (0,1), 'xytext': (0.03,0.95), 'xycoords': 'axes fraction', 'fontsize': 'x-large'}
annotate('a',**graphletter)

ap=dict(arrowstyle="->")
annotate('DMRT-ML', xytext=(30,135), xy=(78,140),arrowprops=ap )
annotate('Tsang and\nKong, 2001', xytext=(30,157), xy=(72,165),arrowprops=ap)

subplot(1,2,2)
original()

res = dmrtml.dmrtml(freq,64,height,density,radius,temp)
theta = range(0,81)
plot(theta,res.TbV(theta),'b--')
plot(theta,res.TbH(theta),'g--')

radius=array([830e-6])    # radius unity: meter

res = dmrtml.dmrtml(freq,64,height,density,radius,temp)
theta = range(0,81)
plot(theta,res.TbV(theta),'b-',label='V-pol')
plot(theta,res.TbH(theta),'g-',label='H-pol')


res_upper = dmrtml.dmrtml(freq,64,height,density,radius,temp,eps_ice=(3.18,0.0016*0.8))
res_lower = dmrtml.dmrtml(freq,64,height,density,radius,temp,eps_ice=(3.18,0.0016*1.2))
theta = arange(0,81,0.1)

y=res.TbV(theta)
yerr=[res_lower.TbV(theta)-y, y-res_upper.TbV(theta)]
errorbar(theta,y,yerr=yerr,color='0.8',zorder=-2,elinewidth=1,capsize=0)

y=res.TbH(theta)
yerr=[res_lower.TbH(theta)-y, y-res_upper.TbH(theta)]
errorbar(theta,y,yerr=yerr,color='0.7',zorder=-3,elinewidth=1,capsize=0)


theta=53.0

annotate('b',**graphletter)

xlabel('Incidence angle (Degree)')
#ylabel ('Brightness temperature (K)') # already on the leftmost y axis
ylim((120,250))
legend(loc=6,frameon=False)
show()
