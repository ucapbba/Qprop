##--------------------------##
# V.A. Tulsky, Rostock, 2021 #
##--------------------------##
from math import *
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi # Why not?





##---------------------------------##
# Choose a folder to take data from #
##---------------------------------##
fname  = '../attoclock/dat/real_prop_vpot_xy.dat'
fname  = '../long-wavelength/dat/real_prop_vpot_z.dat'
fname  = '../vortex/dat/real_prop_vpot_xy.dat'
fname  = '../lissajous/dat/real_prop_vpot_xy.dat'
fname  = '../excited/dat/real_prop_vpot_z.dat'
# fname  = '../hhg/dat/real_prop_vpot_z.dat'
##---------------------------------##

plot_E = 1 # True or False
plot_A = 1 # True or False

# --------------------------------- #
# If True, plots E(t) vs t, -A(t) vs t.
# If False, plots Ex vs Ey, -Ax vs -Ay (only for xy fields).
vs_time = 0
# --------------------------------- #





if fname[-5]=='y': # Check if it is "*xy.dat"
    x = np.genfromtxt(fname)
    Nt = np.size(x)/3
    x = np.transpose(x)
    T  = x[0]
    Ax = x[1]
    Ay = x[2]
    dt = T[1]-T[0]
    if plot_E == True:
        Ex = -(np.roll(Ax,-1)-np.roll(Ax,1))*0.5/dt
        Ey = -(np.roll(Ay,-1)-np.roll(Ay,1))*0.5/dt
        
        fig = plt.figure()
        ax = plt.subplot()
        if vs_time == True:
            ax.plot(T, Ex, linewidth=2, color='red')
            ax.plot(T, Ey, linewidth=2, color='blue')
            plt.xlabel(r"$t$", fontsize=18)
            plt.ylabel(r"$E$", fontsize=18)
        else:
            ax.plot(Ex, Ey, linewidth=2, color='black')
            ax.set_aspect('equal')
            plt.xlabel(r"$E_x$", fontsize=18)
            plt.ylabel(r"$E_y$", fontsize=18)
    if plot_A == True:
        fig = plt.figure()
        ax = plt.subplot()
        if vs_time == True:
            ax.plot(T, -Ax, linewidth=2, color='red')
            ax.plot(T, -Ay, linewidth=2, color='blue')
            plt.xlabel(r"$t$",  fontsize=18)
            plt.ylabel(r"$-A$", fontsize=18)
        else:
            ax.plot(-Ax, -Ay, linewidth=2, color='black')
            ax.set_aspect('equal')
            plt.xlabel(r"$-A_x$", fontsize=18)
            plt.ylabel(r"$-A_y$", fontsize=18)
    plt.show()

if fname[-5]=='z': # Check if it is "*z.dat"
    x = np.genfromtxt(fname)
    Nt = np.size(x)/2
    x = np.transpose(x)
    T  = x[0]
    Az = x[1]
    dt = T[1]-T[0]
    if plot_E == True:
        Ez = -(np.roll(Az,-1)-np.roll(Az,1))*0.5/dt
        
#        fig = plt.figure()
#        ax = plt.subplot()
        plt.plot(T, Ez,':', linewidth=4, color='blue')
        plt.xlabel(r"$t$", fontsize=18)
        plt.ylabel(r"$E_z$", fontsize=18)
    if plot_A == True:
        fig = plt.figure()
        ax = plt.subplot()
        ax.plot(T, -Az, linewidth=2, color='black')
        plt.xlabel(r"$t$", fontsize=18)
        plt.ylabel(r"$-A_z$", fontsize=18)
    plt.show()