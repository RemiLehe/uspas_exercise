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

    def __init__(self, Nz_global=10000, Lz=200., dtcoef=1.):
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
        z = self.dz*(mpi_comm.rank*self.Nz_local + np.arange(0,self.Nz_local+2))
        self.Ex = pulse_shape( z, 0 )
        self.By = 1./c * pulse_shape( z + 0.5*self.dz, -0.5*self.dt )


    def step(self, N_steps):
        """
        Perform `N_steps` iterations of the field update
        """
        if mpi_comm.rank == 0:
            print( 'Performing %d iterations' %N_steps )

        # Loop over the timesteps
        for i_step in range(N_steps):

            # Update the By field from time (n-1/2)*dt to time (n+1/2)*dt
            for k in range(1,self.Nz_local+1):
                self.By[k] = self.By[k] \
                    - self.dt*(self.Ex[k+1]-self.Ex[k])/self.dz
            # MPI exchange: send values of physical cells to neighbors
            # receive values to be put into guard cells
            By_from_left_proc, By_from_right_proc = \
                exchange_guard_cells( self.By[1], self.By[self.Nz_local] )
            # ASSIGNEMENT: Set the guard cells with the right value
            self.By[0] = 0
            self.By[self.Nz_local+1] = 0

            # Update the Ex field from time n*dt to time (n+1)*dt
            # Loop over the gridpoints
            for k in range(1,self.Nz_local+1):
                self.Ex[k] = self.Ex[k] \
                   - c**2*self.dt*(self.By[k]-self.By[k-1])/self.dz
            # MPI exchange: send values of physical cells to neighbors
            # receive values to be put into guard cells
            Ex_from_left_proc, Ex_from_right_proc = \
                exchange_guard_cells( self.Ex[1], self.Ex[self.Nz_local] )
            # ASSIGNEMENT: Set the guard cells with the right value
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
        # PLOTTING: NEW LINES RELATED TO MPI
        Ex_from_all_procs = mpi_comm.gather( self.Ex[1:-1] )
        By_from_all_procs = mpi_comm.gather( self.By[1:-1] )

        if mpi_comm.rank == 0:
            print( 'Plotting the fields at iteration %d' %self.n )
            
            global_Ex = np.concatenate( Ex_from_all_procs )
            global_By = np.concatenate( By_from_all_procs )

            plt.clf()
            plt.suptitle('Fields at iteration %d' %self.n)
            # Plot of Ex
            plt.subplot(211)
            z = self.dz*np.arange( self.Nz_global )
            plt.plot( z, global_Ex, '-' )
            plt.ylim(-1.1, 1.1)
            plt.xlim(0, self.Lz)
            plt.ylabel('$E_x^n$')
            plt.xlabel('z')
            # Plot of By
            plt.subplot(212)
            z = self.dz*np.arange( self.Nz_global ) + 0.5*self.dz
            plt.plot( z, global_By, '-' )
            plt.ylim(-1.1/c, 1.1/c)
            plt.xlim(0, self.Lz)
            plt.ylabel('$B_y^{n-1/2}$')
            plt.xlabel('z')

            if save_figure is True:
                # Check that the diagnostics folder exists
                if os.path.exists('diagnostics') is False:
                    os.mkdir('diagnostics')
                plt.savefig( "diagnostics/iteration%03d.png" %self.n)


def exchange_guard_cells( physical_F_left, physical_F_right ):
    # MPI exchanges of guard cells
    # Send physical cell to left proc
    req1 = mpi_comm.isend( physical_F_left,
                dest=(mpi_comm.rank-1)%mpi_comm.size )
    # Send physical cell to right proc
    req2 = mpi_comm.isend( physical_F_right,
                dest=(mpi_comm.rank+1)%mpi_comm.size )
    # Receive value from right proc
    req3 = mpi_comm.irecv( source=(mpi_comm.rank+1)%mpi_comm.size )
    # Receive value from left proc
    req4 = mpi_comm.irecv( source=(mpi_comm.rank-1)%mpi_comm.size )
    # Wait for the processors to finish sending/receiving
    req1.wait()
    req2.wait()
    F_from_right = req3.wait()
    F_from_left = req4.wait()

    return F_from_left, F_from_right

if __name__ == '__main__':

    # Remove the diagnostics folder if it exists
    if mpi_comm.rank==0 and os.path.exists('diagnostics'):
        shutil.rmtree('diagnostics')

    # Run Nz iterations (in 10 batches, with plotting inbetween)
    em = EM1DSolver( )
    for i in range( 5 ):
        em.plot_fields( save_figure=True )
        em.step( 200 )
    em.plot_fields( save_figure=True )
