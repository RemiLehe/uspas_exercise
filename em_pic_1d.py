"""
This script simulates the propagation of electromagnetic fields in 1D
Author: Remi Lehe

Usage
-----
- To run the code in non-interactive mode:
Type `python em_pic_1d.py`
"""
import os
import shutil # os and shutil are utility packages, that help manipulate files
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c # Numerical value of the speed of light

class EM1DSolver(object):

    def __init__(self, Nz=100, Lz=100., dtcoef=1.):
        """
        Initialize the EM1DSolver object

        Parameters
        ----------
        Nz: int
           Number of grid points

        Lz: real
           Total length of the physical domain

        dtcoef: real
           The timestep will be:
           dt = dtcoef*dz/c  with dz = Lz/Nz
        """
        self.Nz = Nz
        self.dz = Lz/Nz  # Spacing between gridpoints along the z axis
        self.dt = dtcoef*self.dz/c  # Timestep of the simulation
        self.Lz = Lz

        self.n = 0  # Iteration number
        self.z = np.arange(1,Nz+1)*self.dz  # Position of the grid points
        self.Ex = np.zeros( Nz+2 )  # At iteration n, Ex[k] represents the 
                                # field at time n*dt and position k*dz
        self.By = np.zeros( Nz+2 )  # At iteration n, By[k] represents the
                         # field at time (n - 1/2)*dt and position (k+1/2)*dz
        # The first and last element of Ex and By are ghost cell: they simply
        # duplicate the last-but-one and second element respectively
        # (periodic boundaries)

        # Initialize the field values with an oscillating function
        k = 2*np.pi/5. # Wavevector of the pulse
        L = 20.  # Length of the pulse
        def pulse_shape(z, t):
            return( np.cos(k*(z-c*t)) * np.exp(-(z-c*t-Lz/2)**2/L**2) )

        self.Ex[1:-1] = pulse_shape( self.z, 0 )
        self.By[1:-1] = 1./c * pulse_shape( self.z + 0.5*self.dz, -0.5*self.dt )


    def step(self, N_steps):
        """
        Perform `N_steps` iterations of the field update
        """
        print( 'Performing %d iterations' %N_steps )
        
        # Loop over the timesteps
        for i_step in range(N_steps):

            # Update the By field from time (n-1/2)*dt to time (n+1/2)*dt
            for k in range(1,self.Nz+1):
                # ASSIGNMENT: Replace "+ 0" by the appropriate expression
                self.By[k] = self.By[k] + 0
            # Apply periodic boundary conditions (do not modify)
            self.By[0] = self.By[self.Nz]
            self.By[self.Nz+1] = self.By[1]
            
            # Update the Ex field from time n*dt to time (n+1)*dt
            # Loop over the gridpoints
            for k in range(1,self.Nz+1):
                # ASSIGNMENT: Replace "+ 0" by the appropriate expression
                self.Ex[k] = self.Ex[k] + 0
            # Apply periodic boundary conditions (do not modify)
            self.Ex[0] = self.Ex[self.Nz]
            self.Ex[self.Nz+1] = self.Ex[1]
    
            # Increment the iteration number
            self.n = self.n+1
    
    def plot_fields( self, save_figure=False ):
        """
        Plot the Ex and By field using matplotlib
        If save_figure is True, the plots are saved as PNG files,
        in a folder named `diagnostics`
        """
        print( 'Plotting the fields at iteration %d' %self.n )

        plt.clf()
        plt.suptitle('Fields at iteration %d' %self.n)
        # Plot of Ex
        plt.subplot(211)
        plt.plot( self.z, self.Ex[1:-1], 'o-' )
        plt.ylim(-1.1, 1.1)
        plt.xlim(0, self.Lz)
        plt.ylabel('$E_x^n$')
        plt.xlabel('z')
        # Plot of By
        plt.subplot(212)
        plt.plot( self.z+0.5*self.dz, self.By[1:-1], 'o-' )
        plt.ylim(-1.1/c, 1.1/c)
        plt.xlim(0, self.Lz)
        plt.ylabel('$B_y^{n-1/2}$')
        plt.xlabel('z')

        if save_figure is True:
            # Check that the diagnostics folder exists
            if os.path.exists('diagnostics') is False:
                os.mkdir('diagnostics')
            plt.savefig( "diagnostics/iteration%03d.png" %self.n)

if __name__ == '__main__':
    
    # Remove the diagnostics folder if it exists
    if os.path.exists('diagnostics') is True:
        shutil.rmtree('diagnostics')
    
    # Run Nz iterations (in 10 batches, with plotting inbetween)
    em = EM1DSolver( )
    for i in range( 10 ):
        em.plot_fields( save_figure=True )
        em.step( 10 )
    em.plot_fields( save_figure=True )
