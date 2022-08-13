##--------------------------##
# V.A. Tulsky, Rostock, 2021 #
##--------------------------##
from __future__ import print_function
from math import *
import numpy as np
#import matplotlib
#matplotlib.use('Agg') #not todisplay plots unless plt.show() used
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.special

pi = np.pi # Why not?


def Ylm(l, m, theta, phi):
    # Pay attention to a reverse order of arguments!
    result = 0.0
    if abs(m)<=l:
        result = scipy.special.sph_harm(m, l, phi, theta)
    return result


##---------------------------------##
# Choose a folder to take data from #
##---------------------------------##
folder_name = 'long-wavelength'
folder_name = 'excited'
Ntheta = 301 # should be ODD, should coincide with tsurff.param. If 1, theta=0.5pi taken
Nphi   = 1   # should coincide with tsurff.param

# folder_name = 'attoclock'
# folder_name = 'vortex'
# folder_name = 'lissajous'
# Ntheta = 1    # should be ODD, should coincide with tsurff.param. If 1, theta=0.5pi taken
# Nphi   = 120  # should coincide with tsurff.param
##---------------------------------##
path = '../'+folder_name+'/dat/'
tsurff_method = 'isurfv' # 'tsurff' or 'isurfv'

##---------------------------------##
#    Data  type  and  plot  type    #
##---------------------------------##
plot_type = '2D'           # '2D', '1D' or 'angular'
polarization = 'xz'        # 'xz' or 'xy'
m_number = 0               # Only used for 'xz' polarization
expansion_method = 'polar' # 'Ylm' or 'polar' # Ylm if data is recorded in lm-amplitudes
data_type = 'complex'      # 'abs' for [ k|a|^2 ] or 'complex' [ Re(a), Im(a) ] data
log_scale = 1              # True or False
min_order = -4             # Data below 10^(min_order) not shown
multiply_by_k = 1          # True or False
wrt_energy = 0             # With respect to 0.5k^2 or to k (True/False)
eV_units   = 0             # Energy in eV (True) or a.u. (False)
##---------------------------------##

##---------------------------------##
# Linestyle for 1D, colormap for 2D #
##---------------------------------##
FONTSIZE = 16         # General fontsize on all axes and labels
# 2D #
polar_canvas = 0      # True - round frame, False - square frame
cartesian = 1         # If True and if polar_canvas=False, plots in xz/xy coordinates
cmap  = 'hot_r'       # pylab.cm.get_cmap('hot_r'6)
Custom_ticks = 0      # True or False
TICKS_COLOR = 'black'
# 1D #
COLOR = 'black'
LINE='-'
##---------------------------------##







units_label=' (a.u.)'
if wrt_energy==True:
    if eV_units==True:
        units_label=' (eV)'

if data_type == 'abs':
    multiply_by_k = False # Data already in k|a|^2 format

if polarization == 'xz':
    Nphi = 1 # Spectrum is axial symmetric

if Ntheta>1:
    if Ntheta%2 == 0:
        Ntheta += 1
        print('Number of theta angles increased by one, as it should be ODD (as in the QPROP code)')

if expansion_method == 'polar':
    fname = tsurff_method + '_polar.dat'
if expansion_method == 'Ylm':
    fname = tsurff_method + '_partial.dat'

print('Taking data from', path+fname)
print('Making', plot_type, 'plot.')

if expansion_method == 'polar': 
    x=np.loadtxt(path+fname)
    x = np.transpose(x)
    Nk = int(len(x[0])/Ntheta/Nphi)
    ik_max = Nk # Reduce to truncate data
    Energy = x[0].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
    K      = x[1].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
    Theta  = x[2].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
    Phi    = x[3].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::]
    kmax=K[-1,-1,-1]
    if wrt_energy==True:
        if eV_units==True:
            Energy*=27.211
            K*=(27.211*2)**0.5
    if data_type == 'abs':
        W = x[4].reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::] # k |a|^2
    if data_type == 'complex':
        W = (x[5]+1.j*x[6]).reshape(Nk,Ntheta,Nphi)[0:ik_max,::,::] # real(a)+j imag(a)
        W = abs(W*W)
    Nk = ik_max
if expansion_method == 'Ylm':
    x=np.loadtxt(path+fname)
    x = np.transpose(x)
    Nk= len(x[0])
    ik_max = Nk # Reduce to truncate data
    Lm_max = len(x)-2-1 # First 2 rows are Ek and k, last row is k*(total sum)
    if data_type == 'complex':
        Lm_max = int(Lm_max/2)
    Energy = x[0][0:ik_max]
    k      = x[1][0:ik_max]
    if wrt_energy==True:
        if eV_units==True:
            Energy*=27.211
            k*=(27.211*2)**0.5
    if Ntheta>1:
        theta  = pi*np.arange(Ntheta)/(Ntheta-1)
    else:
        theta = np.arange(1)+0.5*pi
    if polarization=='xz':
        phi = np.array([0])
    if polarization=='xy':
        phi = 2*pi*np.arange(Nphi)/(Nphi-1)
    kmax=k[-1]

    Theta,K,Phi = np.meshgrid(theta,k,phi)
    Theta,Energy,Phi = np.meshgrid(theta,Energy,phi)
    W = np.zeros_like(K)+0.j
    Wlm = [np.zeros_like(k)+0.j for i in range(Lm_max)]
    if data_type == 'abs':
        print('This data does not allow resolving in angles')
    if data_type == 'complex':  
        if polarization=='xz':
            print('l=',)
            for l in range(0,Lm_max):
                print(l,)
                Wlm[l] = (x[2+2*l]+1.j*x[2+2*l+1])[0:ik_max]
                for itheta in range(Ntheta):
                    YLM = Ylm(l, m_number, Theta[0,itheta,0], 0.)
                    W[::,itheta,0] += Wlm[l]*YLM
        if polarization=='xy':
            Lmax = int(Lm_max**0.5)
            print('l=',)
            for l in range(0,Lmax):
                print(l,' m=',)
                for m in range(-l,l+1):
                    print(m,)
                    lm = l*(l+1)+m
                    Wlm[lm] = (x[2+2*lm]+1.j*x[2+2*lm+1])[0:ik_max]
                    for itheta in range(Ntheta):
                        for iphi in range(Nphi):
                            YLM = Ylm(l, m, Theta[0,itheta,0], Phi[0,0,iphi])+0.j
                            W[::,itheta,iphi] += Wlm[lm]*YLM
                print()
        W = abs(W*W)
    Nk = ik_max

if multiply_by_k == True:
    W *= K
    

if plot_type == '2D':
    if polarization == 'xz':
        phi_ind  = 0
        Angles2D =  Theta[::,::,phi_ind]
        E2D      = Energy[::,::,phi_ind]
        K2D      =      K[::,::,phi_ind]
        W2D      =      W[::,::,phi_ind]
    elif polarization == 'xy':
        theta_ind = int(Ntheta/2)
        Angles2D  =    Phi[::,theta_ind,::]
        E2D       = Energy[::,theta_ind,::]
        K2D       =      K[::,theta_ind,::]
        W2D       =      W[::,theta_ind,::]     
    fig = plt.figure() 
    if wrt_energy==True:
        R = E2D
    if wrt_energy==False:
        R = K2D
    vmax=W2D.max()
    if polar_canvas==True:
        ax = plt.subplot(polar=True)
        if log_scale==True:
            cplt = ax.pcolor( Angles2D, R, (W2D+1e-100), norm = LogNorm(), vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
        if log_scale==False:
            cplt = ax.pcolor( Angles2D, R, (W2D+1e-100),                   vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
        if Custom_ticks==True:
            ticks_step = 0.2*ceil(kmax*10)*0.1
            TICKS = [ticks_step*i for i in range(1,5)]
            Angle_ticks = [45*i for i in range(8)]
            plt.plot([0,0],[R.min(),R.max()], linewidth=1, color=TICKS_COLOR)
            plt.plot(np.zeros_like(TICKS),TICKS,'.',color=TICKS_COLOR)
            ax.set_rgrids(TICKS, angle=0, color=TICKS_COLOR, horizontalalignment='center', verticalalignment='bottom', fontsize=FONTSIZE)
            ax.set_thetagrids(Angle_ticks,color='black', fontsize=FONTSIZE)
    if polar_canvas==False:
        ax = plt.subplot()
        if cartesian==True:
            X = R*np.cos(Angles2D)
            Y = R*np.sin(Angles2D)
        if cartesian==False:
            X = Angles2D
            Y = R
        if log_scale==True:
            cplt = ax.pcolor( X, Y, (W2D+1e-100), norm = LogNorm(), vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
        if log_scale==False:
            cplt = ax.pcolor( X, Y, (W2D+1e-100),                   vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
        if cartesian==True:
            if polarization == 'xz':
                if wrt_energy==True:
                    plt.xlabel(r'${\bf k}_{\|} k_{\|}/2$'+units_label,fontsize=FONTSIZE)
                    plt.ylabel(r'${\bf k}_{\perp} k_{\perp}/2$'+units_label,fontsize=FONTSIZE)
                if wrt_energy==False:
                    plt.xlabel(r'$k_{\|}$'+units_label,fontsize=FONTSIZE)
                    plt.ylabel(r'$k_{\perp}$'+units_label,fontsize=FONTSIZE)
            if polarization == 'xy':
                if wrt_energy==True:
                    plt.xlabel(r'${\bf k}_x k_x/2$'+units_label,fontsize=FONTSIZE)
                    plt.ylabel(r'${\bf k}_y k_y/2$'+units_label,fontsize=FONTSIZE)
                if wrt_energy==False:
                    plt.xlabel(r'$k_x$'+units_label,fontsize=FONTSIZE)
                    plt.ylabel(r'$k_y$'+units_label,fontsize=FONTSIZE)
        if cartesian==False:
            if polarization == 'xz':
                plt.xlabel(r'$\theta$',fontsize=FONTSIZE)
                if Custom_ticks==True:
                    plt.xticks([0.5*pi*i for i in range(-2,3,1)],fontsize=FONTSIZE)
                    ax.set_xticklabels([str(0.5*i)+r"$\pi$" for i in range(-2,3,1)],fontsize=FONTSIZE)
            if polarization == 'xy':
                plt.xlabel(r'$\varphi$',fontsize=FONTSIZE)
                if Custom_ticks==True:
                    plt.xticks([0.5*pi*i for i in range(0,5,1)],fontsize=FONTSIZE)
                    ax.set_xticklabels([str(0.5*i)+r"$\pi$" for i in range(0,5,1)],fontsize=FONTSIZE)
            plt.xlim(Angles2D.min(),Angles2D.max())
            if wrt_energy==True:
                plt.ylabel(r'$k^2/2$'+units_label,fontsize=FONTSIZE)
            if wrt_energy==False:
                plt.ylabel(r'$k$'+units_label,fontsize=FONTSIZE)
    plt.xticks(fontsize=FONTSIZE)
    plt.yticks(fontsize=FONTSIZE)
    cbar = plt.colorbar(cplt)

    cbar.ax.tick_params(labelsize=FONTSIZE)
    cbar.set_label(label="Diff. ioniz. probability (a.u.)",fontsize=FONTSIZE)
    plt.show()
    
if plot_type == '1D':
    phi_angles_ind = [0,Nphi/2]
    for i in range(2):
        if polarization == 'xz':
            itheta = -i # First & last
            iphi = 0
        if polarization == 'xy':
            itheta = Ntheta/2
            iphi = phi_angles_ind[i]
        W1D = W[1:Nk,itheta,iphi]
        if wrt_energy==True:
            X = Energy[1:Nk,itheta,iphi]
            xlabel = r'${\bf k}k/2$'+units_label
        if wrt_energy==False:
            X = K[1:Nk,itheta,iphi]
            xlabel = r'$k$'+units_label
        Y = W1D
        plt.plot((1-2*i)*X, Y, linestyle=LINE,color=COLOR,linewidth=2)
        plt.plot((1-2*i)*X, Y, 'o')
    if log_scale==True:
        plt.yscale('log')
        plt.ylim(1e-6,1e-0)
    plt.xticks(size=FONTSIZE)
    plt.yticks(size=FONTSIZE)
    plt.xlabel(xlabel,fontsize=FONTSIZE)
    plt.ylabel(r"Diff. ioniz. probability (a.u.)",fontsize=FONTSIZE)
    plt.show()
    
if plot_type == 'angular':
    fig = plt.figure()
    ax = plt.subplot()
    if polarization == 'xz':
        iphi = 0
        Wang = (W[::,::,iphi]).sum(0)
        Wang *= 1./Wang.max()
        angles = Theta[0,::,iphi]
        plt.plot(angles,Wang,color='black',linewidth=3)
        plt.plot(-angles,Wang,color='black',linewidth=3)
        plt.xlim(-pi,pi)
        plt.ylim(0,1.1)
        plt.xticks([0.5*pi*i for i in range(-2,3,1)],fontsize=FONTSIZE)
        ax.set_xticklabels([str(0.5*i)+r"$\pi$" for i in range(-2,3,1)],fontsize=FONTSIZE)
        plt.xlabel(r"$\theta$",fontsize=FONTSIZE)
    if polarization == 'xy':
        itheta = Ntheta/2
        Wang = (W[::,itheta,::]).sum(0)
        Wang *= 1./Wang.max()
        angles = Phi[0,itheta,::]
        plt.plot(angles,Wang,color='black',linewidth=2)
        plt.xlim(0,2*pi)
        plt.xticks([0.5*pi*i for i in range(0,5,1)],fontsize=FONTSIZE)
        ax.set_xticklabels([str(0.5*i)+r"$\pi$" for i in range(0,5,1)],fontsize=FONTSIZE)
        plt.xlabel(r"$\varphi$",fontsize=FONTSIZE)
    plt.ylabel("Angular distribution (arb.u.)",fontsize=FONTSIZE)
    plt.show()