"""
This script simulates the propagation of electromagnetic fields in 1D
Author: Remi Lehe

Usage
-----
- To run the code in non-interactive mode:
Type `python em_pic_1d.py`

- To run the code in interactive mode:
Type `ipython --matplotlib` and then:
`from em_pic_1d import step, plot_fields`
Then call e.g. `step(10)` and `plot_fields()` as many times as you want.
"""
import os
import shutil # os and shutil are utility packages, that help manipulate files
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c # Numerical value of the speed of light

# Numerical parameters of the simulation
# --------------------------------------
Lz = 100.        # Total length of the physical domain
Nz = 100         # Number of grid points
dz = Lz/Nz       # Spacing between gridpoints along the z axis
dt = dz/c  # Timestep of the simulation

# Initialize the field arrays
# ---------------------------
n = 0   # Iteration number
Ex = np.zeros( Nz+2 )    # At iteration n, Ex[k] represents the field
                         # at time n*dt and position (k + 1/2)*dz
By = np.zeros( Nz+2 )    # At iteration n, By[k] represents the field
                         # at time (n - 1/2)*dz and position k*dz
# The first and last element of Ex and By are ghost cell: they simply duplicate
# the last-but-one and second element respectively (periodic boundaries)

# Initialize the field values with an oscillating function
# --------------------------------------------------------
k = 2*np.pi/5. # Wavevector of the pulse
L = 20.  # Length of the pulse
def pulse_shape(z, t):
    """Return the profile of the pulse: gaussian envelope + oscillations"""
    return( np.cos(k*(z-c*t)) * np.exp(-(z-c*t-Lz/2)**2/L**2) )
z = np.arange(1,Nz+1)*dz
Ex[1:-1] = pulse_shape( z + 0.5*dz, 0 )
By[1:-1] = 1./c * pulse_shape( z, -0.5*dt )

# Definition of the step and plot_fields functions
# ------------------------------------------------
def step(N_steps):
    """
    Perform `N_steps` iterations of the field update
    """
    global n, Ex, By
    print( 'Performing %d iterations' %N_steps )
    
    # Loop over the timesteps
    for i_step in range(N_steps):

        # Update the By field from time (n-1/2)*dt to time (n+1/2)*dt
        for k in range(1,Nz+1):
            # ASSIGNMENT: Replace "+ 0" by the appropriate expression
            By[k] = By[k] + 0
        # Apply periodic boundary conditions (do not modify)
        By[0] = By[Nz-1]
        By[Nz] = By[1]
        
        # Update the Ex field from time n*dt to time (n+1)*dt
        # Loop over the gridpoints
        for k in range(1,Nz+1):
            # ASSIGNMENT: Replace "+ 0" by the appropriate expression
            Ex[k] = Ex[k] + 0
        # Apply periodic boundary conditions (do not modify)
        Ex[0] = Ex[Nz-1]
        Ex[Nz] = Ex[1]

        # Increment the iteration number
        n = n+1

def plot_fields( save_figure=False ):
    """
    Plot the Ex and By field using matplotlib
    If save_figure is True, the plots are saved as PNG files,
    in a folder named `diagnostics`
    """
    print( 'Plotting the fields at iteration %d' %n )
    
    plt.clf()
    plt.suptitle('Fields at iteration %d' %n)
    # Plot of Ex
    plt.subplot(211)
    plt.plot( z+0.5*dz, Ex[1:-1], 'o-' )
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, Lz)
    plt.ylabel('$E_x^n$')
    plt.xlabel('z')
    # Plot of By
    plt.subplot(212)
    plt.plot( z, By[1:-1], 'o-' )
    plt.ylim(-1.1/c, 1.1/c)
    plt.xlim(0, Lz)
    plt.ylabel('$B_y^{n-1/2}$')
    plt.xlabel('z')
    
    if save_figure is True:
        # Check that the diagnostics folder exists
        if os.path.exists('diagnostics') is False:
            os.mkdir('diagnostics')
        plt.savefig( "diagnostics/iteration%03d.png" %n) 


if __name__ == '__main__':

    # Remove the diagnostics folder if it exists
    if os.path.exists('diagnostics') is True:
        shutil.rmtree('diagnostics')

    # Run Nz iterations (in 10 batches, with plotting inbetween)
    for i in range(10):
        plot_fields( save_figure=True )
        step(Nz/10)
    plot_fields( save_figure=True )
