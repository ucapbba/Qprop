#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:25:58 2022

@author: thomas
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.special

pi = np.pi 

### Insert desired path here ###
qprop_path = 'C:/GIT/AttoPhysics/Qprop/dat/'

##---------------------------------##
#    Data  type  and  plot  type    #
##---------------------------------##
log_scale = 1              # True or False
min_order = -10          # Data below 10^(min_order) not shown

##---------------------------------##
#         Colormap for plot         #
##---------------------------------##
FONTSIZE =10         # General fontsize on all axes and labels
# 2D #
cmap  = 'cividis'       
TICKS_COLOR = 'black'
##---------------------------------##


units_label=' (a.u.)'
Ntheta = 301 # should be ODD, should coincide with tsurff.param. If 1, theta=0.5pi taken
Nphi   = 1   # should coincide with tsurff.param

tsurff_method = 'isurfv' # 'tsurff' or 'isurfv'
file_name = tsurff_method + '_polar.dat'

##---------------------------------##
#    Data  type  and  plot  type    #
##---------------------------------##
data_type = 'complex'      # 'abs' for [ k|a|^2 ] or 'complex' [ Re(a), Im(a) ] data
multiply_by_k = 0          # True or False
wrt_energy = 0             # With respect to 0.5k^2 or to k (True/False)
##---------------------------------##

if Ntheta>1:
    if Ntheta%2 == 0:
        Ntheta += 1
        print('Number of theta angles increased by one, as it should be ODD (as in the QPROP code)')

x=np.loadtxt(qprop_path+file_name)
print('Making Qprop plot.')
x = np.transpose(x)
Nk = int(len(x[0])/Ntheta/Nphi)
ik_max = Nk # Reduce to truncate data
Energy = x[0].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
K      = x[1].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
Theta  = x[2].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
Phi    = x[3].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
kmax=K[-1,-1,-1]
if data_type == 'abs':
    W = x[4].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::] # k |a|^2
if data_type == 'complex':
    W = (x[5]+1.j*x[6]).reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::] # real(a)+j imag(a)
    W = abs(W*W)
Nk = ik_max


if multiply_by_k == True:
    W *= K

phi_ind  = 0
Angles2D =  Theta[::,::,phi_ind]
E2D      = Energy[::,::,phi_ind]
K2D      =      K[::,::,phi_ind]
W2D      =      W[::,::,phi_ind]
       
fig = plt.figure() 
if wrt_energy==True:
    R = E2D    
if wrt_energy==False:
    R = K2D

# W2D = W2D/W2D.max()+np.max(10**(min_order-1)-W2D/W2D.max(), 0)
vmax=W2D.max()
ax = plt.subplot()

X = R*np.cos(Angles2D)
Y = R*np.sin(Angles2D)
extra=2
if log_scale==True:
    cplt = ax.pcolor( X, Y, (W2D), norm = LogNorm(), vmin=10**min_order, vmax=1, cmap=cmap)
if log_scale==False:
    cplt = ax.pcolor( X, Y, (W2D+1e-100),                   vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
if wrt_energy==True:
    plt.xlabel(r'${\bf k}_{\parallel} k_{\|}/2$'+units_label,fontsize=FONTSIZE+extra)
    plt.ylabel(r'${\bf k}_{\perp} k_{\perp}/2$'+units_label,fontsize=FONTSIZE+extra)
if wrt_energy==False:
    plt.xlabel(r'$p_{\parallel}$'+units_label,fontsize=FONTSIZE+2*extra)
    plt.ylabel(r'$p_{\perp}$'+units_label,fontsize=FONTSIZE+2*extra)
                    
plt.axis('scaled')
plt.xticks(fontsize=FONTSIZE+extra)
plt.yticks(fontsize=FONTSIZE+extra)
cbar = plt.colorbar(cplt)

cbar.ax.tick_params(labelsize=FONTSIZE-4)
cbar.set_label(label="Normalised Differential ionization probability ",fontsize=FONTSIZE)
plt.show()