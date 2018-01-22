"""
This is a script that integrates the equations of motion for a particle
in constant electric and magnetic fields.

Usage:
------
python particle_pusher.py
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, c

# Set the default mass and charge to that of an electron
m = m_e
q = -e

class ParticleIntegrator( object ):
    """
    Class that integrates the motion of a single particle in constant
    E and B field
    """

    def __init__(self, Ex, Ey, Ez, Bx, By, Bz, px0, py0, pz0, dt):
        """
        Initialize the particle integrator

        Parameters:
        -----------
        Ex, Ey, Ez: values of the constant electric field
        Bx, By, Bz: values of the constant magnetic field
        px0, py0, pz0: values of the initial momenta
        dt: timestep
        """
        # Fields
        self.Ex = Ex
        self.Ey = Ey
        self.Ez = Ez
        self.Bx = Bx
        self.By = By
        self.Bz = Bz

        # Iteration counter and timestep
        self.n = 0
        self.dt = dt
        
        # Initial positions and momenta
        # At a given iteration n, the positions correspond to time n*dt
        # and the momenta correspond to time (n-1/2)*dt.
        self.x = 0
        self.y = 0
        self.z = 0
        self.px = px0
        self.py = py0
        self.pz = pz0

        # Corresponding history
        self.x_history = [ ]
        self.y_history = [ ]
        self.z_history = [ ]
        self.px_history = [ ]
        self.py_history = [ ]
        self.pz_history = [ ]
        self.t_history = [ ]

    def step( self, N=1 ):
        """Perform N timesteps of the particle pusher"""

        # Loop over timesteps
        for i in range(N):
            
            # Push px, py, pz from (n-1/2)*dt to (n+1/2)*dt
            self.update_momenta()
            
            # Push x, y, z from n*dt to (n+1)*dt
            self.update_positions()

            # Increment iteration number
            self.n = self.n + 1
            
            # Register the new position in the history (for diagnostics)
            self.x_history.append( self.x )
            self.y_history.append( self.y )
            self.z_history.append( self.z )
            self.px_history.append( self.px )
            self.py_history.append( self.py )
            self.pz_history.append( self.pz )
            self.t_history.append( self.n*self.dt )

    def update_positions( self ):
        """Update x, y, z over one timestep"""
        
        # Compute Lorentz factor
        gamma = np.sqrt(1 + (self.px**2 + self.py**2 + self.pz**2)/(m*c)**2 )

        # Update the particle positions
        self.x = self.x + self.px*dt/(gamma*m)
        self.y = self.y + self.py*dt/(gamma*m)
        self.z = self.z + self.pz*dt/(gamma*m)
            
    def update_momenta( self ):
        """Update px, py, pz over one timestep using the alternative pusher"""

        # ASSIGNEMENT: COMPLETE THIS METHOD

        # Calculate tau and s
        tau_x = q*self.Bx*dt/(2*m)
        tau_y = q*self.By*dt/(2*m)
        tau_z = q*self.Bz*dt/(2*m)
        gamma_old = np.sqrt(1 + (self.px**2 + self.py**2 + self.pz**2)/(m*c)**2)
        s_x = q*self.Bx*dt/(2*m*gamma_old)
        s_y = q*self.By*dt/(2*m*gamma_old)
        s_z = q*self.Bz*dt/(2*m*gamma_old)

        # Calculte p' (represented in the code by pp_x, pp_y, pp_z)
        pp_x = self.px + q*self.Ex*dt + (self.py*s_z - self.pz*s_y)
        ## ASSIGNEMENT: COMPLETE THE LINES BELOW
        pp_y = self.py + 0
        pp_z = self.pz + 0

        # Calculate gamma_new (i.e. gamma at (n+1/2)*dt)
        tau2 = tau_x**2 + tau_y**2 + tau_z**2
        ## ASSIGNEMENT: COMPLETE THE EXPRESSION OF gamma_prime and u
        gamma_prime = 0
        u = 0
        ## ASSIGNEMENT: REPLACE THE EXPRESSION OF gamma_new BY THE CORRECT ONE
        gamma_new = 1

        # Calculate t
        ## ASSIGNEMENT: REPLACE THE EXPRESSION OF t_x, t_y, t_z
        ## BELOW BY THE CORRECT ONES
        t_x = tau_x
        t_y = tau_y
        t_z = tau_z
        
        # Calculate the new px, py, pz
        t2 = t_x**2 + t_y**2 + t_z**2
        p_dot_t = pp_x*t_x + pp_y*t_y + pp_z*t_z
        self.px = ( pp_x + (pp_y*t_z-pp_z*t_y) + p_dot_t*t_x )/(1 + t2)
        ## ASSIGNEMENT: REPLACE THE EXPRESSION BELOW BY THE CORRECT ONES
        self.py = 0
        self.pz = 0

    def plot_history( self ):
        """ Plot the history of the positions and momenta as a function of time"""
        plt.subplot(211)
        plt.title('Positions')
        plt.plot( self.t_history, self.x_history, '-o', label='x' )
        plt.plot( self.t_history, self.y_history, '-o', label='y' )
        plt.plot( self.t_history, self.z_history, '-o', label='z' )
        plt.legend(loc=0)

        plt.subplot(212)
        plt.title('Momenta')
        plt.plot( self.t_history, self.px_history, '-o', label='px' )
        plt.plot( self.t_history, self.py_history, '-o', label='py' )
        plt.plot( self.t_history, self.pz_history, '-o', label='pz' )
        plt.legend(loc=0)
        plt.xlabel('Time')
        
if __name__ == '__main__':

    Bz0 = 1.
    py0 = 100*m*c
    dt = 1.e-10

    # Integrate the equations of motion with a timestep dt
    pusher = ParticleIntegrator( Ex=0, Ey=0, Ez=0, Bx=0, By=0, Bz=Bz0,
                                        px0=0, py0=py0, pz0=0, dt=dt )
    pusher.step( 100 )
    pusher.plot_history()
    plt.show()
    
    
    

        
        
        
