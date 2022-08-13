##--------------------------##
# V.A. Tulsky, Rostock, 2021 #
##--------------------------##
from math import *
import numpy as np
#import matplotlib
#matplotlib.use('Agg') #not todisplay plots unless plt.show() used
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.special
pi = np.pi # why not?

def Ylm(l, m, theta, phi):
    return scipy.special.sph_harm(m, l, phi, theta) # pay attention to a reverse order of arguments!


# ================================================================== #

path = '../long-wavelength/'
path = '../excited/'

# path = '../attoclock/'
# path = '../hhg/'
# path = '../vortex/'
# path = '../lissajous/'

polarization = '34'     #  qprop-dim in "initial.param"

Nl     = 25              # According to "initial.param" or "propagate.param"
m      = 0              # According to "initial.param"
dr     = 0.2            # According to "initial.param"
Nr_max = 200            # Limit of radial points to save time

fname = 'dat/imag_prop_'+str(m)+'_wf_fin.dat'
fname = 'dat/real_prop_wf.dat'
# it = 2000                                     # To make a snapshot at timestep "it"
# fname = 'wf/real_prop_wf_'+str(it)+'.dat'   #  of a real-time propagated w.f.

plot_type  = '2D'       # '1D' or '2D'
plot_plane = 'xz'       # Plotting plane 'xz' or 'xy'
Nangle     = 101        # Ntheta for 'xz', Nphi for 'xy'
mirrored   = True       # If True then 'xz' plot is mirrored

min_order  = -20         # data below 10^(min_order) is not plotted
cmap       = 'hot_r'    # colormap for plot_type = '2D'

# ================================================================== #

if plot_plane == 'xz':
    angle = pi/max(1,Nangle-1)*np.arange(Nangle)
    fixed_phi = 0.0*pi
if plot_plane == 'xy':
    angle = 2.*pi/max(1,Nangle-1)*np.arange(Nangle)
    fixed_theta = 0.5*pi

print('Loading')  
x=np.loadtxt(path+fname)
r = np.arange(Nr_max)*dr

print('Constructing array')  
if polarization == '34':
    Nr = np.size(x)/2/Nl
    Psi = np.zeros((Nl,Nr_max))+0.j
    for il in range(m,Nl):
        for ir in range(Nr_max):
            Psi[il,ir]=x[int(Nr*il+ir)][0]+1.j*x[int(Nr*il+ir)][1]
    psi = np.zeros((Nr_max,Nangle))+0.j
    for il in range(m,Nl):
        if plot_plane == 'xz':
            for i_theta in range(Nangle):
                YLM = Ylm(il,m,angle[i_theta],fixed_phi)
                for ir in range(Nr_max):
                    psi[ir,i_theta] += Psi[il,ir]*YLM
        if plot_plane == 'xy':
            for i_phi in range(Nangle):
                YLM = Ylm(il,m,fixed_theta,angle[i_phi])
                for ir in range(Nr_max):
                    psi[ir,i_phi] += Psi[il,ir]*YLM

if polarization == '44':
    Nr = np.size(x)/2/(Nl*Nl)
    Psi = np.zeros((Nl*Nl,Nr_max))+0.j
    for il in range(0,Nl):
        for im in range(-il,il+1):
            ilm = il*(il+1)+im
            for ir in range(Nr_max):
                Psi[ilm,ir]=x[int(Nr*ilm+ir)][0]+1.j*x[int(Nr*ilm+ir)][1]
    psi = np.zeros((Nr_max,Nangle))+0.j
    for il in range(0,Nl):
        for im in range(-il,il+1):
            ilm = il*(il+1)+im
            if plot_plane == 'xz':
                for i_theta in range(Nangle):
                    YLM = Ylm(il,im,angle[i_theta],fixed_phi)
                    for ir in range(Nr_max):
                        psi[ir,i_theta] += Psi[ilm,ir]*YLM
            if plot_plane == 'xy':
                for i_phi in range(Nangle):
                    YLM = Ylm(il,im,fixed_theta,angle[i_phi])    
                    for ir in range(Nr_max):
                        psi[ir,i_phi] += Psi[ilm,ir]*YLM

Angle, R = np.meshgrid(angle,r)

print('Plotting')
if plot_type == '2D':
    ax = plt.subplot()
    Min_r_ind = 1
    Max_r_ind = Nr_max
    X = R*np.cos(Angle)
    Y = R*np.sin(Angle)
    X = X[Min_r_ind:Max_r_ind,::]
    Y = Y[Min_r_ind:Max_r_ind,::]
    w = (abs(psi**2))[Min_r_ind:Max_r_ind,::]
    R2 = (R**2)[Min_r_ind:Max_r_ind,::]
    Func = (w/R2+1e-40)
#        Func *= 1./Func.max()
    vmax = Func.max()
    cplt = plt.pcolor(X, Y, Func, norm = LogNorm(), vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
    if mirrored==True and plot_plane == 'xz':
        cplt = plt.pcolor(X, -Y, Func, norm = LogNorm(), vmin=vmax*10**min_order, vmax=vmax, cmap=cmap)
    if plot_plane == 'xz':
        plt.xlabel('z, (a.u.)')
        plt.ylabel('x, (a.u.)')
    if plot_plane == 'xy':
        plt.xlabel('x, (a.u.)')
        plt.ylabel('y, (a.u.)')
    plt.colorbar(cplt)
    ax.set_aspect('equal')
    
if plot_type == '1D':
    Min_r_ind = 1
    Max_r_ind = Nr_max
    X = r[Min_r_ind:Max_r_ind]
    for il in range(0,Nl):
        w = (abs(Psi)**2)[il][Min_r_ind:Max_r_ind]
        r2 = (X**2)
        Y = w/r2
        plt.plot(X, Y)
    plt.xlim(0,50)