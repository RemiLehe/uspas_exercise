import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from scipy.constants import c

def get_exp_ikz( dz, dt, w ):
    return( 1.j*(dz/(c*dt))*np.sin( w*dt*0.5 ) + \
            np.sqrt( 1.+0.j - (dz/(c*dt))**2*np.sin( w*dt*0.5 )**2 ) )

def plot_reflection( dz1, dz2 ):

    dt = dz1/c
    beta = 2*c*dt/(dz1 + dz2)
    w = np.pi/(dt*200) * (np.arange(200) + 0.5 )

    eikz1 = get_exp_ikz( dz1, dt, w )
    eikz2 = get_exp_ikz( dz2, dt, w )

    R = abs( ( np.exp(-1.j*w*dt*0.5) - np.exp(1.j*w*dt*0.5) + \
             beta*( eikz2 - 1./eikz1 ) )/
             ( np.exp(-1.j*w*dt*0.5) - np.exp(1.j*w*dt*0.5) + \
             beta*( eikz2 + eikz1 ) ) )

    plt.subplot(121)
    plt.plot( w, R )
    plt.xlabel('$\omega$')
    plt.ylabel('R')

    plt.subplot(122)
    plt.semilogx( 2*np.pi*c/(w*dz2), R )
    plt.xlabel('$\lambda/\Delta z_2$')
    plt.ylabel('R')

    plt.show()
    
if __name__ == '__main__':
    dz1 = 1.
    dz2 = 5.
    plot_reflection( dz1, dz2 )
    

