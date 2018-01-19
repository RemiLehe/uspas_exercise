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
# Import MPI communicator
from mpi4py.MPI import COMM_WORLD as mpi_comm

class EM1DSolver(object):

    def __init__(self, Nz_global=200, Lz=200., dtcoef=1.):
        """
        Initialize the EM1DSolver object

        Parameters
        ----------
        Nz_global: int
           Total number of grid points (across all MPI domains)

        Lz: real
           Total length of the physical domain

        dtcoef: real
           The timestep will be:
           dt = dtcoef*dz/c  with dz = Lz/Nz
        """
        self.Nz_global = Nz_global
        self.dz = Lz/Nz_global  # Spacing between gridpoints along the z axis
        self.dt = dtcoef*self.dz/c  # Timestep of the simulation
        self.Lz = Lz
        self.n = 0  # Iteration number
    
        # INITIALIZATION: NEW LINES RELATED TO MPI
        self.Nz_local = int( self.Nz_global/mpi_comm.size )
        self.Ex = np.zeros( self.Nz_local+2 )
        self.By = np.zeros( self.Nz_local+2 )
        print('MPI rank %d initialized its local sub-domain.' %mpi_comm.rank) 
        # The first and last element of Ex and By are guard cell:
        # They duplicate elements from the neighboring MPI processes

        # Initialize the field values with an oscillating function
        k = 2*np.pi/5. # Wavevector of the pulse
        L = 20.  # Length of the pulse
        def pulse_shape(z, t):
            return( np.cos(k*(z-c*t)) * np.exp(-(z-c*t-Lz/2)**2/L**2) )
        z = self.dz*(mpi_comm.rank*self.Nz_local + np.arange(1,self.Nz_local+1))
        self.Ex[1:-1] = pulse_shape( z, 0 )
        self.By[1:-1] = 1./c * pulse_shape( z + 0.5*self.dz, -0.5*self.dt )


    def step(self, N_steps):
        """
        Perform `N_steps` iterations of the field update
        """
        print( 'Performing %d iterations' %N_steps )
        
        # Loop over the timesteps
        for i_step in range(N_steps):

            # Update the By field from time (n-1/2)*dt to time (n+1/2)*dt
            for k in range(1,self.Nz_local+1):
                self.By[k] = self.By[k] \
                    - self.dt*(self.Ex[k+1]-self.Ex[k])/self.dz
            # Apply periodic boundary conditions
            # Set the values of the guard cells
            self.By[0] = 0
            self.By[self.Nz_local+1] = 0

            # Update the Ex field from time n*dt to time (n+1)*dt
            # Loop over the gridpoints
            for k in range(1,self.Nz_local+1):
                self.Ex[k] = self.Ex[k] \
                   - c**2*self.dt*(self.By[k]-self.By[k-1])/self.dz
            # Apply periodic boundary conditions (do not modify)
            self.Ex[0] = 0
            self.Ex[self.Nz_local+1] = 0
    
            # Increment the iteration number
            self.n = self.n+1
    
    def plot_fields( self, save_figure=False ):
        """
        Plot the Ex and By field using matplotlib
        If save_figure is True, the plots are saved as PNG files,
        in a folder named `diagnostics`
        """
        print( 'Plotting the fields at iteration %d' %self.n )
        
        # PLOTTING: NEW LINES RELATED TO MPI
        global_Ex = np.concatenate( mpi_comm.gather( self.Ex ) )
        global_By = np.concatenate( mpi_comm.gather( self.By ) )

        plt.clf()
        plt.suptitle('Fields at iteration %d' %self.n)
        # Plot of Ex
        plt.subplot(211)
        plt.plot( global_Ex, 'o-' )
        plt.ylim(-1.1, 1.1)
        plt.xlim(0, self.Lz)
        plt.ylabel('$E_x^n$')
        plt.xlabel('z')
        # Plot of By
        plt.subplot(212)
        plt.plot( global_By, 'o-' )
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
    if mpi_comm.rank==0 and os.path.exists('diagnostics') is True:
        shutil.rmtree('diagnostics')
    
    # Run Nz iterations (in 10 batches, with plotting inbetween)
    em = EM1DSolver( )
    for i in range( 10 ):
        em.plot_fields( save_figure=True )
        em.step( 10 )
    em.plot_fields( save_figure=True )
