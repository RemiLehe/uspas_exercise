"""
This script simulates the propagation of electromagnetic fields in 1D 
on a succession of 2 grids that can be set a different resolutions 
and are linked by an algorithm selected by the user.
The code prints the coefficients of reflection and transmission, and can perform 
scans on the wavelength and the method used to connect the grids.
Author: J.-L. Vay and Remi Lehe

Usage
-----
- To run the code in non-interactive mode:
Type `python em_pic_1d_mr.py`
"""

from warp import *

l_scan = 0           # performs scan if True, run one case otherwise

method = 1           # method for linking the two grids
                     # 1 => dE/dt = (B2-B1)*2./(dx1+dx2)
                     # 2 => dE/dt = B2/dx2-B1/dx1
                     # 3 => dE/dt = B2/dx2-B1/dx2 with jumping in grid 1 to get B1

Nz     = 300         # number of grid cells
MRcoef = 3           # ratio between grid cells size of grid 2 vs grid 1
lw     = 25.         # wavenumber
k      = 2.*pi/lw    # wave vector (normalized by dx)
dtcoef = 1.          # factor to multiply dz/c for setting the time step on grid 1 and the reference grid
c      = clight      # shortcut for the speed of light

if l_scan:
    methods = [1,2,3]
    lws = 10.**(arange(log10(4.5),log10(200.),log10(1.5)))
    ks = 2.*pi/lws
    
def getw_yee(k,dx,dt):
    # return omega (normalized by dx/c) as a function of normalized k
    return 2/dt*arcsin(dt/dx*sin(k*dx/2))

def getvgroup(k,dx,dt):
    # return group velocity (normalized by c) as a function of normalized k
    return cos(k*dx/2)/sqrt(1.-(dt/dx)**2*(sin(k*dx/2)**2))

class EM1DSolver(object):

    def __init__(self, Nz=100, Lz=100., dtcoef=1., omegadx=None, l_init=0, vgroup=1.):
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
        self.omegadx = omegadx
        self.t = 0.
        self.L = 100./vgroup

        self.n = 0  # Iteration number
        self.z = arange(Nz)*self.dz  # Position of the grid points
        self.Ex = zeros( Nz )  # At iteration n, Ex[k] represents the 
                                # field at time n*dt and position k*dz
        self.By = zeros( Nz )  # At iteration n, By[k] represents the
                         # field at time (n - 1/2)*dz and position (k+1/2)*dz

        self.Exh = AppendableArray()
        
    def step(self, N_steps):
        """
        Perform `N_steps` iterations of the field update
        """
        
        # Loop over the timesteps
        for i_step in range(N_steps):

            
            # Update the By field from time (n-1/2)*dt to time (n+1/2)*dt
            for k in range(0,self.Nz-1):
                # ASSIGNMENT: Replace "+ 0" by the appropriate expression
                self.By[k] = self.By[k] \
                    - self.dt*(self.Ex[k+1]-self.Ex[k])/self.dz

            # impose field at at left end
            if self.omegadx is not None:
                if (c*self.t)<self.L:
                    coef = 2.*pi*c*self.t/self.L
                    self.Ex[0] =  (10.-15.*cos(coef)+6.*cos(2*coef)-cos(3*coef))/32. \
                               *  cos(self.omegadx*self.t/(self.dz/c)) 
                else:
                    self.Ex[0] = 0.
                self.Exh.append(self.Ex[0])
                               
            # Update the Ex field from time n*dt to time (n+1)*dt
            # Loop over the gridpoints
            for k in range(1,self.Nz):
                # ASSIGNMENT: Replace "+ 0" by the appropriate expression
                self.Ex[k] = self.Ex[k] \
                   - c**2*self.dt*(self.By[k]-self.By[k-1])/self.dz

            # Increment the iteration number and time
            self.n = self.n+1
            self.t = self.t+self.dt
    
    def get_energy(self):
        return self.dz*(self.Ex*self.Ex+self.By*self.By*c*c)
        
    def get_totenergy(self):
        return sum(self.get_energy())

def myrun(k,method,l_plot=False): 
    # --- get omega of injected wave for the specified wave vector 
    omegadx = getw_yee(k,1.,dtcoef)

    # --- get group velocity for the specified wave vector 
    vgroup  = getvgroup(k,1.,dtcoef)
    
    # --- initializes grid 1, 2 and reference grid
    emgrid1 = EM1DSolver(Nz=Nz, Lz=100., dtcoef=dtcoef, omegadx = omegadx, vgroup=vgroup)
    emgrid2 = EM1DSolver(Nz=int(Nz/MRcoef), Lz=100., dtcoef=dtcoef/MRcoef)
    emgridref = EM1DSolver(Nz=2*Nz, Lz=200., dtcoef=dtcoef, omegadx = omegadx, vgroup=vgroup)

    # --- setup diagnostics
    if l_plot:
        setup()
        winon()

        def myplot():
            # --- plots electric fields 
            plsys(9)
            plg(emgridref.Ex,arange(emgridref.Nz)*emgrid1.dz,marks=0,color='blue',width=3)
            plg(emgrid1.Ex,arange(emgrid1.Nz)*emgrid1.dz,marks=0)
            plg(emgrid2.Ex,arange(emgrid2.Nz)*emgrid2.dz+emgrid1.Lz-emgrid1.dz,marks=0,color='red')
            limits(0.,200.,-1.5,1.5)
            ptitles('','Electric field','Z')
            # --- plots energy
            plsys(10)
            plg(emgridref.get_energy(),arange(emgridref.Nz)*emgrid1.dz,marks=0,color='blue',width=3)
            plg(emgrid1.get_energy(),arange(emgrid1.Nz)*emgrid1.dz,marks=0)
            plg(emgrid2.get_energy(),arange(emgrid2.Nz)*emgrid2.dz+emgrid1.Lz-emgrid1.dz,marks=0,color='red')
            limits(0.,200.,-1.5,1.5)
            ptitles('','Field energy','Z')
        
    def step(n=1):
        for i in range(n):
            if l_plot and i%10==0:
                fma();
                myplot()
                pyg_pending()

            # --- advance magnetic field on last node of first grid using first node of second grid
            emgrid1.By[-1] = emgrid1.By[-1] \
                           - emgrid1.dt*(emgrid2.Ex[0]-emgrid1.Ex[-1])/emgrid1.dz
            # --- advance electric and magnetic field in the core of grid 1
            emgrid1.step(1)

            # --- advance electric field on last node of first grid using first node of second grid
            # --- with selected method
            if method==1:
                dzMR = 0.5*(emgrid1.dz+emgrid2.dz)
                emgrid2.Ex[0] = emgrid2.Ex[0] \
                              - c**2*emgrid2.dt*(emgrid2.By[0]-emgrid1.By[-1])/dzMR

            if method==2:
                emgrid2.Ex[0] = emgrid2.Ex[0] \
                              - c**2*emgrid2.dt*(emgrid2.By[0]/emgrid2.dz-emgrid1.By[-1]/emgrid1.dz)

            if method==3:
                emgrid2.Ex[0] = emgrid2.Ex[0] \
                              - c**2*emgrid2.dt*(emgrid2.By[0]-emgrid1.By[-MRcoef+1])/emgrid2.dz

            # --- advance electric and magnetic field in the core of grid 2
            emgrid2.step(1)

            # --- advance electric and magnetic field in the reference grid
            emgridref.step(1)

    # --- run until the pulse has reflected from interface between grids 1 and 2 
    # --- and propagated back to the middle of grid 1
    step(int(2*Nz/(dtcoef*vgroup)))
        
    print '    total energy 1, 2, ref = ',emgrid1.get_totenergy(),emgrid2.get_totenergy(),emgridref.get_totenergy()
    print '    Incident wave number     = ',2.*pi/k   

    R = sqrt(emgrid1.get_totenergy()/emgridref.get_totenergy())
    T = sqrt(emgrid2.get_totenergy()/emgridref.get_totenergy())
    print '    Reflection coefficient   = ', R
    print '    Transmission coefficient = ', T
    print '    R+T = ',R+T
    
    return R

if l_scan:
    R = []           
    for method in methods:
        R.append(AppendableArray())
        print 'Start scan for method ',method
        for k in ks:
            print '  k = ',k
            R[-1].append(myrun(k,method))
else:
    R = myrun(k,method,l_plot=True)

if l_scan:
    n = 1000    
    k1 = arange(n)*pi*0.99999/(n-1)
    k1 = k1[1:]
    dt  = dtcoef
    dx1 = 1.
    dx2 = MRcoef*dx1

    def getw_yee(k,dt,dx):
        return 2/dt*arcsin(dt/dx*sin(k*dx/2))
 
    def getk_yee(w,dt,dx):
        a=dt/dx
        c=-dt/dx
        b=exp(1j*w*dt/2)-exp(-1j*w*dt/2)
        k = 2.*1j*log((-b+sqrt(b*b-4*a*c))/(2*a))/dx
        return k

    omega = getw_yee(k1,dt,dx1)
    k2 = getk_yee(omega,dt,dx2)
    
    def refl(omega,k1,k2,dt,dx1,dx2,beta1,beta2):
        j=1j
        num   = exp(j*omega*dt/2)-exp(-j*omega*dt/2)+beta2*exp(-j*k2*dx2/2)-beta1*exp(j*k1*dx1/2)
        denom = exp(j*omega*dt/2)-exp(-j*omega*dt/2)+beta2*exp(-j*k2*dx2/2)+beta1*exp(-j*k1*dx1/2)
        return abs(num/denom)
    
    beta1 = beta2 = dt*2./(dx1+dx2)
    Rscheme1 = refl(omega,k1,k2,dt,dx1,dx2,beta1,beta2)

    beta1 = dt/dx1
    beta2 = dt/dx2
    Rscheme2 = refl(omega,k1,k2,dt,dx1,dx2,beta1,beta2)

    pla(Rscheme1,2*pi/k1,color=blue,width=4)
    pla(Rscheme2,2*pi/k1,color=red,width=4)
    limits(2.,100.,1.e-4,10)
    logxy(1,1)
    ptitles('','!l/!dx','|R|')
    ppgeneric(R[0][...],lws,msize=1,marker='o',color=blue)
    ppgeneric(R[1][...],lws,msize=1,marker='+',color=red)
    ppgeneric(R[2][...],lws,msize=1,marker='x',color=black)
    pdf('coef_refl')
    
