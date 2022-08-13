##--------------------------##
# V.A. Tulsky, Rostock, 2021 #
##--------------------------##
import numpy as np
import matplotlib
#matplotlib.use('Agg') #not todisplay plots unless plt.show() used
import matplotlib.pyplot as plt

pi = np.pi # why not?



##---------------------------------##
# Choose a folder to take data from #
##---------------------------------##
fname  = 'attoclock'
fname  = 'long-wavelength'
fname  = 'vortex'
fname  = 'lissajous'
fname  = 'excited'
fname  = 'hhg'

fname = '../' + fname + '/dat/real_prop_obser.dat'
##---------------------------------##

plot_type = 'hhg' # 'hhg' or 'accel'

# Makes window(t<Ncut*dt)=0 and window(t>(Nt-Ncut)*dt)=0
Ncut = 201 



x = np.genfromtxt(fname) # read data from file
x = np.transpose(x)

Nt = len(x[0])       # Take less to truncate the time interval

T        = x[0][0:Nt] # Time
d2s1_dt2 = x[5][0:Nt] # = -d^2<S>/dt^2 = - ( -<dU/dS> - ES ), S = z or x
d2s2_dt2 = x[6][0:Nt] # = -d^2<S>/dt^2 = - ( -<dU/dS> - ES ), S = 0 or y

dt        = T[1]-T[0] # Timestep
omega_max = pi/dt     # Frequencies greater definitely cannot be resolved with such timestep
domega    = pi/T[-1]  # Same way from the time-frequency uncertainty: Dt*Domega > pi
Nomega    = int(omega_max/domega)+1

window = np.zeros_like(T)
T_win  = (T[Ncut:-Ncut-1]-T[Ncut])/T[-2*Ncut-1]
# This window is applied to the integrand in Fourier transrorm
window[Ncut:-Ncut-1] = np.sin(pi*T_win)**2 

if plot_type == 'hhg':
    Omega = np.arange(Nomega)*domega
    HHG1 = np.zeros_like(Omega)+0.j
    HHG2 = np.zeros_like(Omega)+0.j
#    for iom in range(Nomega):
#        HHG1[iom] = (d2s1_dt2*np.exp(-1.j*Omega[iom]*T)*window).sum()*dt
#        HHG2[iom] = (d2s2_dt2*np.exp(-1.j*Omega[iom]*T)*window).sum()*dt
    HHG1 = np.fft.rfft(d2s1_dt2*window, n = 2*Nomega-1)
    HHG2 = np.fft.rfft(d2s2_dt2*window, n = 2*Nomega-1)
    fig = plt.figure()
    X  = Omega
    Y1 = abs(HHG1**2).real
    Y2 = abs(HHG2**2).real
    plt.plot(X,  Y1+Y2,    linewidth=2, color='black')
    # plt.plot(X,  Y1,    linewidth=2, color='blue')
    # plt.plot(X,  Y2,    linewidth=2, color='red')
    plt.yscale('log')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(r"$ \Omega$ (a.u.)", fontsize=18)
    plt.ylabel(r"$ \vert \langle \ddot{\bf r}_\Omega \rangle \vert^2$ (a.u.)", fontsize=18)
    plt.show()

if plot_type == 'accel':
    X  = T*dt
    Y1 = d2s1_dt2#/abs(d2s1_dt2).max()
    Y2 = d2s2_dt2#/abs(d2s1_dt2).max()
    Y3 = window
    plt.plot(X,  Y3,'--',linewidth=2, color='grey')
    plt.plot(X,  Y1,     linewidth=2, color='black')
#    plt.plot(X,  Y2,     linewidth=2, color='blue')
    plt.xlabel(r"$t$ (a.u.)", fontsize=18)
    plt.ylabel(r"$\langle \ddot{x}(t) \rangle$ (arb.u.)", fontsize=18)
    plt.xlim(X[0],X[-1])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()