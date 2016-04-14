import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from copy import deepcopy
from numpy.linalg import eig

import energy
import kdtree as K
from energy import get_particles_E
from misc import crossmatch
import cluster_density
from binary_motion import rotate

AU_to_parsec = 4.84e-6
RSun_parsec  = 2.2506e-8

def frmt(x):
    return "%19.11e" % x
    
def SphVol(radius):
    return np.pi*4./3 * radius**3



class PySnap(object):
    """
    This class is a dynamical state of an nbody simulation.

    It is initialized with the corresponding time, names of
    stars, positions and velocities.

    S = PySnap(t,n,m,x,y,z,vx,vy,vz)
    """
    def __init__(self,t,n,m,x,y,z,vx,vy,vz):
        self.t=t
        self.N=len(n)
        self.n=np.array(n)
        self.m=np.array(m)
        self.x=np.array(x)
        self.y=np.array(y)
        self.z=np.array(z)
        self.vx=np.array(vx)
        self.vy=np.array(vy)
        self.vz=np.array(vz)
        self.units=["length","time","mass"]
        self.UptoDate = False
        self.Silent = False
        self.G = 1.0

    def MakeSilent(self):
        """
        Right now this only stops the printing of "computing energies".
        Was and can be used for more.
        """
        self.Silent=True

    #=====================================================================
    # Mass, radius, velocity and scaling
    #=====================================================================

    def Mt(self):
        """Total mass"""
        return sum(self.m)

    def scale_r(self, factor):
        """Apply a multiplicative factor to x,y,z."""
        self.x = self.x * factor
        self.y = self.y * factor
        self.z = self.z * factor
        if hasattr(self,"spanning_tree"):
            del self.spanning_tree
        self.UptoDate = False
    def Scale_r(self,factor):
        """
        Apply a multiplicative factor to x,y,z, then correct
        velocities and time to preserve the dynamical state.
        """
        self.t = self.t * factor**(3./2)
        self.scale_r( factor )
        self.scale_v(1/np.sqrt(factor))
        self.UptoDate = False

    def scale_v(self, factor):
        """Apply a multiplicative factor to vx,vy,vz"""
        self.vx = self.vx * factor
        self.vy = self.vy * factor
        self.vz = self.vz * factor
        self.UptoDate = False

    def scale_m(self, factor, CorrectVelocities=False):
        """Apply a multiplicative factor to masses"""
        self.m = factor * self.m
        self.UptoDate = False

    def r(self):
        """ 2D array of positions. """
        r=np.concatenate(([self.x],[self.y],[self.z]),axis=0)
        r=np.swapaxes(r,0,1)
        return r
    def v(self):
        """ 2D array of velocities. """
        v=np.concatenate(([self.vx],[self.vy],[self.vz]),axis=0)
        v=np.swapaxes(v,0,1)
        return v
    def V(self):
        """ Array of velocity : sqrt( vx**2 + vy**2 + vz**2 )"""
        return  np.sqrt(self.vx**2 + self.vy**2+ self.vz**2)
    def R(self):    
        """ Array of radius : sqrt( x**2 + y**2 + z**2 ) """
        return  np.sqrt(self.x**2 + self.y**2+ self.z**2)

    def LagrangianRadius(self, fraction, Return_n = False):
        """
        Return the Lagrangian radius of the snap for a given fraction of mass. 
        For example to get the radius of 25% of the mass:
           Rhm = snap.LagrangianRadius(0.25)
        If Return_n is True, also returns the names of stars within this radius
        """
        S = self.SelectParticles( indStars = np.argsort(self.R()) )
        R = S.R()
        cum_m = np.cumsum(S.m)
        ind = np.max(np.nonzero( cum_m < fraction*cum_m[-1] ) )
        if not Return_n:
            return R[ind]
        else:
            return R[ind], S.n[ind]

    def sigmaV(self):
        """Velocity dispersion"""
        return np.sqrt(np.mean(self.V()**2))
        
    def CrossingTime(self):
        """Return the crossing time of the system: Rhm/sigmaV"""
        return self.LagrangianRadius(0.5)/self.sigmaV()

    def Rhm(self):
        """Half-mass radius (after correction for center of mass)"""
        TempCopy = deepcopy(self)
        TempCopy.CorrectionCenterOfMass()
        return TempCopy.LagrangianRadius(0.5)

    def Rv(self):
        """Virial radius :  G * Mt**2 / ( 3 * Ep ) """ 
        return - self.G * self.Mt()**2 / ( 3 * self.tEp ) 

    def AngularMomentum(self):
        """Total angular momentum of the system. Should be conserved."""
        return np.sum( np.transpose(np.tile(self.m, (3,1)))*np.cross(self.r(),self.v()),
                       axis=0 )

    def UniformFreeFallTime(self,G=1):
        """Returns the free fall time, only for Uniform clusters"""
        rho = self.Mt() / (SphVol(max(self.R())))
        return 0.25*np.sqrt(3*np.pi/(2.*rho))

    def CenterOfMass(self):
        """Return coordinates of the center of mass"""
        com_x,com_y,com_z= sum(self.m*self.x)/self.Mt(),\
                           sum(self.m*self.y)/self.Mt(),\
                           sum(self.m*self.z)/self.Mt()
        return [com_x,com_y,com_z]

    def CenterOfVelocity(self):    
        """Return coordinates of the center of velocities"""
        cov_x,cov_y,cov_z= sum(self.m*self.vx),\
                           sum(self.m*self.vy),\
                           sum(self.m*self.vz)
        return [cov_x,cov_y,cov_z]
        

    def CorrectionCenterOfMass(self,com=None):
        """
        Set the center of mass to provided coordinates in com,
        if com not set, set the center of mass to 0,0,0.
        """
        if com is not None:
            self.x=  self.x-com[0]
            self.y=  self.y-com[1]
            self.z=  self.z-com[2]
        else:
            self.x=  self.x-sum(self.m*self.x)/self.Mt()
            self.y=  self.y-sum(self.m*self.y)/self.Mt()
            self.z=  self.z-sum(self.m*self.z)/self.Mt()


    def CorrectionGeometricalCenter(self):
        """Set the geometrical center (not mass weighted) to 0,0,0"""
        com = [ np.mean(self.x), np.mean(self.y), np.mean(self.z) ]
        self.CorrectionCenterOfMass( com = com )

    def CorrectionCenterOfVelocity(self):
        """Set the center of velocity to 0,0,0"""
        self.vx=  self.vx-sum(self.m*self.vx)/self.Mt()
        self.vy=  self.vy-sum(self.m*self.vy)/self.Mt()
        self.vz=  self.vz-sum(self.m*self.vz)/self.Mt()
        self.UptoDate=False



    #=================================================================================
    # Energy stuff
    #=================================================================================

    def compute_energies(self, Ek=True, Ep=True, E=True, G=1.):
        """
        Default behaviour returns the array of each particles total, 
        kinetic and potential energy:
            E, Ek, Ep = snap.compute_energies()
        Ek and Ep keywords allow to choose not to compute one of these.
        """
        Ek_arr,Ep_arr = energy.get_particles_E(self.x,self.y,self.z,
                                               self.V(),self.m,G)
        ret=[]
        if E:  ret.append(Ep_arr + Ek_arr)
        if Ek:  ret.append(Ek_arr)
        if Ep:  ret.append(Ep_arr)
        if len(ret)==1: ret=ret[0]
        return ret

    def _get_energy(self,attr):
        """
        This function shouldn't be explicitely used. It's part of the
        energy system of PySnap:

        To get the array of particles energies, simply use "snap.E".
        The @property below associated with this function allows the
        instance to check if energies were already calculated  when the 
        energy array is requested, and compute it if needed. In this 
        latter case, a message "Computing energies" will be printed.

        The same works with kinetic and potential energies: "snap.Ek" 
        and "snap.Ep". To get total system energies: "snap.tE", "snap.tEk"
        and "snap.tEp". 
        """
        if hasattr(self,attr) and self.UptoDate:
            return getattr(self,attr)
        else:
            if self.UptoDate:
                if "t" in attr and hasattr(self,"_"+attr[2:]):
                    self._tEk = self._Ek.sum()
                    self._tEp = 0.5*self._Ep.sum()
                    self._tE = self._tEk + self._tEp
                    return getattr(self,attr)
            else:
                if not self.Silent:
                    print "Computing energies"
                self._E, self._Ek, self._Ep = \
                    self.compute_energies(Ek=True,Ep=True,E=True,G=self.G)
                self._tEk = self._Ek.sum()
                self._tEp = 0.5 * self._Ep.sum() # Not counting each pair twice !
                self._tE = self._tEk + self._tEp
                self.UptoDate = True
                return getattr(self,attr)
    
    @property
    def E(self): return self._get_energy("_E")
    @E.setter
    def E(self,E_array): self._E = E_array
    @property
    def Ek(self): return self._get_energy("_Ek")
    @Ek.setter
    def Ek(self,E_array): self._Ek = E_array
    @property
    def Ep(self): return self._get_energy("_Ep")
    @Ep.setter
    def Ep(self,E_array): self._Ep = E_array
    @property
    def tE(self): return self._get_energy("_tE")
    @property
    def tEk(self): return self._get_energy("_tEk")
    @property
    def tEp(self): return self._get_energy("_tEp")

    def Q(self, G=1):
        """Return virial ratio: -Ek/Ep"""
        return abs(self.tEk/self.tEp)

    def virialise(self,Q,Standard=True,G=1):
        """
        Scales velocities to specified virial ratio, computes it 
        and returns it for verification. 
        0.5 is equilibrium, 0 is no kinetic energy. 
        """
        tol=1e-8
        if Q==0.:
            self.scale_v(0.)
            self.UptoDate=False
            return 0.
        else:
            x,y,z,vx,vy,vz = energy.virialise(Q, self.x,self.y,self.z,
                                              self.vx,self.vy,self.vz,self.m,
                                              Standard=Standard, G=G)
            [self.x,self.y,self.z,self.vx,self.vy,self.vz]=[x,y,z,vx,vy,vz]
            self.UptoDate=False
            return self.Q()



    #=====================================================================================
    # Hubble field
    #=====================================================================================

    def HubbleField(self, Mean=True):
        """ 
        Gives out the Hubble parameter of the relation V = H*R
        (from the radial projection). If Mean is set to False, 
        one H parameter is returned by particle, as an array.
        By default, the mean H will be returned.
        """
        r = self.r()
        v = self.v()
        prod = r*v
        arr =  np.sum( r*v, axis=1) / np.sum( r*r, axis=1) 
        if Mean:
            return np.mean(arr)
        else:
            return arr

    def ApplyHubbleField(self, value):
        """
        Impose an Hubble field V=H*R by adding appropriate velocity
        """
        currentfield = np.mean(self.HubbleField())
        newv = (value - currentfield) * self.r()
        self.vx = self.vx + newv[:,0]
        self.vy = self.vy + newv[:,1]
        self.vz = self.vz + newv[:,2]



    #=====================================================================================
    # Densities
    #=====================================================================================


    def GetDensities(self, Nnb=10, Number=False, Recompute=False):
        """
        Returns an array of local density for each particles, obtained through the 
        evaluation of mass in the sphere defined by the Nnb-th neighbour. This is from
        Casertano & Hut 1985 and uses the kdtree.
        Nnb       :  Number of neighbours to use to determine density
        Number    :  Boolean. If True, number density is returned instead of mass density
        Recompute :  Boolean. If True, pre-existing self.rho or rho_number is ignored
        """
        if Number:
            if Recompute or not hasattr(self,"rho_number"):
                self.rho_number = cluster_density.GetDensity(self.m, 
                                                             self.x, self.y, self.z, 
                                                             Nnb, Number)
            return self.rho_number
        if not Number:
            if Recompute or not hasattr(self,"rho"):
                self.rho = cluster_density.GetDensity(self.m, self.x, self.y, self.z, 
                                                       Nnb, Number)
            return self.rho

    def DensityCenter(self,Number=False):
        """Returns the density-weighted barycenter (Casertano&Hut,85)"""
        self.densities = self.GetDensities(6,Number=Number)
        return np.sum( self.densities * np.transpose(self.r()), axis=1 ) / \
            np.sum(self.densities)

    def CorrectionDensityCenter(self,Number=False):
        """Set the center to the density-weighted barycenter (Casertano&Hut,85)"""
        self.densities = self.GetDensities(6,Number=Number)
        r_center = np.sum( self.densities * np.transpose(self.r()), axis=1 ) / \
                   np.sum(self.densities)
        self.CorrectionCenterOfMass( com = r_center )

    def DensityRadius(self):
        """Returns the density-weighted mean radius (Casertano&Hut,85)"""
        self.CorrectionDensityCenter()
        return np.sum( self.densities * self.R() ) / np.sum(self.densities)

    def RandomRotation(self, NewSnap=False):
        """
        Apply a rotation of a random angle along all 3 directions. 
        If you do not want to modify the current snap, you can generate a 
        rotated model by setting NewSnap=True
        """
        angles = 2*np.pi*np.random.random(3)
        r = rotate( self.r(), *angles)
        v = rotate( self.v(), *angles)
        x = r[:,0]; y = r[:,1]; z = r[:,2]
        vx = v[:,0]; vy = v[:,1]; vz = v[:,2]
        if NewSnap:
            return PySnap(self.t,self.n,self.m,x,y,z,vx,vy,vz)
        else:
            self.x = x; self.y = y; self.z = z
            self.vx = vx; self.vy = vy; self.vz = vz


    #=====================================================================================
    # Various file printing and plots
    #=====================================================================================
            
    def WriteSnapshot(self, filename, NB6_InputFormat=False, GLnemoFormat=False):
        """
        filename = snap.Writesnapshot(filename)

        Write an ascii snapshot:
             N
             DIM
             t
             n1 m1 x1 y1 z1 vx1 vy1 vz1
             n2 m2 x2 y2 z2 vx2 vy2 vz2
             ...
             nN mN xN yN zN vxN vyN vzN

        If NB6_InputFormat is set, the 3 variables at the top are
        omitted and particles identities (n) are not printed. This is the
        kind of file Nbody6 expects as an input ("fort.10").
        If GLnemoFormat is set, the following format applies:
            N
            DIM
            t
            m1
            m2
            ....
            x1 y1 z1
            x2 y2 z2
            .....
            vx1 vy1 vz1
            vx2 vy2 vz2
        """
        data=[self.n,self.m,self.x,self.y,self.z,self.vx,self.vy,self.vz]
        with open(filename,"w") as f:
            # Nbody6 doesn't want a header or the identity (n) of 
            # stars when the snapshot is taken as initial state
            if NB6_InputFormat:
                for i in range(self.N):
                    for a in data[1:]:
                        f.write(frmt(a[i]))
                    f.write("\n")
                return filename
            else:
                f.write(str(self.N)+"\n")
                f.write(str(3)+"\n")
                f.write(str(self.t)+"\n")
            if GLnemoFormat:
                for i in range(self.N):
                    f.write(frmt(self.m[i])+"\n")
                for i in range(self.N):
                    f.write(frmt(self.x[i])+frmt(self.y[i])+frmt(self.z[i])+"\n")
                for i in range(self.N):
                    f.write(frmt(self.vx[i])+frmt(self.vy[i])+frmt(self.vz[i])+"\n")
                return filename
            else:
                for i in range(self.N):
                    for a in data:
                        f.write(frmt(a[i]))
                    f.write("\n")
                return filename


    def Plot(self,symbol='o',markersize=2, 
             subsets=None, ax_object = None, DontTouchAxis=False,**kwargs):
        """
        ax = snap.Plot()

        Open an interactive 3d representation of the snapshot and return
        the axis object. You can choose to plot only a part of the cluster, or
        to successively plot several subsets of particles, which
        results in them being displayed in different colors. For 
        example, if you want to only plot the particles within a 0.4 
        radius from origin:

               snap.Plot(subsets=[snap.R() < 0.4])

        If you want the particles with x > 1 in a different color:

               snap.Plot(subsets =[[snap.x < 0.4], [snap.x > 0.4] ])

        You can provide the axis object to plot the PySnap on an
        existing plot, but it should be created with a 3d projection:
        
               import mpl_toolkits.mplot3d.axes3d as p3
               fig = matplotlib.pyplot.figure() 
               ax = fig.add_subplot( 111, projection='3d' )
               ax.plot(range(5),range(5),range(5))
               snap.Plot(ax_object=ax)
        
        By default, the axis range of an existing plot is remodified be 
        centered on the origin and contain every data point. To turn this
        off and maintain existing range, set DontTouchAxis True.

        All keywords for the normal plot function can be passed.
        """
        if ax_object is None:
            fig=plt.figure(figsize=(10,10))
            ax=fig.add_subplot(111,projection='3d',adjustable='box')
        else:
            ax = ax_object
        ax.set_aspect('equal')
        plt.tight_layout()
        if subsets is None:
            subsets = [[range(self.N)]]
        for sub in subsets:
            ax.plot(self.x[sub],self.y[sub],self.z[sub],
                    symbol, markersize=markersize, **kwargs)
        Max_coord=np.max(np.abs(np.array([self.x,self.y,self.z])))
        if not DontTouchAxis:
            ax.set_xlim3d(-1.1*Max_coord,1.1*Max_coord)
            ax.set_ylim3d(-1.1*Max_coord,1.1*Max_coord)
            ax.set_zlim3d(-1.1*Max_coord,1.1*Max_coord)
        plt.show()
        return ax


    def DensityHistogram(self, bins=20, Plot = True, Log = True):
        """
        radiusbins, density = snap.DensityHistogram()
        
        Separate the system in radial bins, then compute the volumic stellar 
        density inside each. 
           bins : if integer, sets the number of bins
                  if array, sets bins edges
           Plot : True will display density=f(radius) in linear scale
           Log  : True will change the plot to log-log scale.
        """
        R=self.R()
        if type(bins) is int:
            radiusbins=np.logspace( np.log10(np.min(R)) ,
                                    np.log10(np.max(R)) , bins)
        else:
            radiusbins = bins
        volumebins= SphVol(radiusbins[1:]) - SphVol(radiusbins[:-1])
        Mbins=[]
        for i in range(len(radiusbins)-1):
            ind= np.logical_and(R>radiusbins[i], R<radiusbins[i+1])
            Mbins.append(sum(self.m[ind]))
        density = np.array(Mbins) / volumebins
        if Plot:
            if Log:
                plt.plot(np.log10(radiusbins[0:len(radiusbins)-1]), 
                         np.log10(density))
            else:
                plt.plot(radiusbins[0:len(radiusbins)-1], density)
            plt.show()
            plt.set_xlabel("Radius")
            plt.set_xlabel("Density")
        return radiusbins[0:len(radiusbins)-1], density


    def SurfaceDensity(self, bins=20, Plot = True, Log = True ):
        """
        radiusbins, sigma = snap.SurfaceDensity()

        Return the stellar projected surface density in logarithmic radial bins
           bins : if integer, sets the number of bins
                  if array, sets bins edges
           Plot : True will display density=f(radius) in linear scale
           Log  : True will change the plot to log-log scale.
        """
        self.CorrectionCenterOfMass()
        Rproj=np.sqrt(self.x**2+self.y**2)
        if type(bins) is int:
            radiusbins=np.logspace( np.log10(np.min(Rproj)) ,
                                    np.log10(np.max(Rproj)) , bins)
        else:
            radiusbins = bins
        nbins = len(radiusbins)
        sigma=[]
        for i in range(nbins-1):
            ind= np.logical_and(Rproj>radiusbins[i], 
                                Rproj<radiusbins[i+1])
            sigma.append(sum(self.m[ind]) /   \
                ( np.pi * (radiusbins[i+1]**2-radiusbins[i]**2) ))
        if Plot:
            if Log:
                plt.plot(np.log10(radiusbins[0:len(radiusbins)-1]), 
                         np.log10(np.array(sigma)))
            else:
                plt.plot(radiusbins[0:len(radiusbins)-1], sigma)
            plt.show()
            plt.set_xlabel("Radius")
            plt.set_xlabel("Surface Density")
        return radiusbins[0:len(radiusbins)-1], np.array(sigma)


    def CumulatedMass(self, bins=25, Plot=False, Logr= False):
        """
        radiusbins, cumM = snap.SurfaceDensity()

        Return the cumulated mass starting from the center in logarithmic radial bins
           bins : if integer, sets the number of bins
                  if array, sets bins edges
           Plot : True will display density=f(radius) in linear scale
           Logr  : True will create logarithmic radius bins
        """
        R=self.R()
        if type(bins) is int:
            if Logr:
                radiusbins=np.logspace( np.log10(np.min(R)) , 
                                        np.log10(np.max(R)) , bins)
            else:
                radiusbins=np.linspace( np.min(R) , 
                                        np.max(R) , bins)
        else:
            radiusbins = bins
        nbins = len(radiusbins)
        massbins = np.zeros(nbins)
        for i in range(nbins):
            massbins[i] = sum(self.m[R < radiusbins[i+1]])
        if Plot:
            plt.plot(radiusbins[1:], massbins)
            plt.show()
            plt.set_xlabel("Radius")
            plt.set_xlabel("Cumulated mass")
        return radiusbins[1:], massbins


    def MassFunction(self, bins=25, Plot=False, **kwargs):
        """
        bins, N_arr = snap.MassFunction( bins = 15, linewidth = 2  )

        Return the mass function.
           bins : if integer, sets the number of bins from min to max
                  if array, sets bins edges
           Plot : True will display the mass function histogram
        All normal plot keywords can be passed.
        """
        m = self.m
        if type(bins) is int:
            massbins=np.logspace( np.log10(np.min(m)) , 
                                  np.log10(np.max(m)) , bins)
        else:
            massbins = bins
        Nstar = np.zeros(len(massbins)-1)
        for i in range(len(massbins)-1):
            Nstar[i] = len(m[ (massbins[i] < m) & ( m < massbins[i+1] ) ])
        massbins = (massbins[1:] + massbins[:-1])/2
        if Plot:
            plt.plot(np.log10(massbins), np.log10(Nstar),**kwargs)
            plt.show()
            plt.set_xlabel("Mass")
            plt.set_xlabel("N")
        return massbins, Nstar


    #=====================================================================================
    # Select, delete, add particles
    #=====================================================================================


    def SelectParticles(self, bool_index=None, nStars=None, indStars=None,
                        TakeEnergies=False):
        """
        subsnap =  snap.SelectParticles( snap.y > 0.5 )
        subsnap =  snap.SelectParticles( nStars = StarsNames )
        subsnap =  snap.SelectParticles( indStars = np.array([0,12,14,85,114]) )

        Return a subset of the current instance as a new PySnap. You
        can use it in three ways:
            bool_index :   For an array a, the statement "a > 0.5" returns a 
                         boolean numpy array of the size of a, with True at each 
                         position where the condition is met. So if you want to 
                         select particles with distance to origin below 0.8, just 
                         use: 
                           snap.SelectParticles( snap.R() < 0.8)
            nStars     :   A numpy array containing the indices of the particles
                         you want to extract.
            indStars    :  A numpy array of the actual indices of stars in
                         the snap arrays. Use at your own risk.
            TakeEnergies: If True, the stars in the new system bring their energies 
                          from the full system with them. If energies were not previously
                          computed when you use this, they will be, which can take a long
                          time, hence the default "False".
        """
        if bool_index is not None:
            S =  PySnap(self.t, self.n[bool_index], self.m[bool_index],
                          self.x[bool_index], self.y[bool_index], 
                          self.z[bool_index], self.vx[bool_index], 
                          self.vy[bool_index], self.vz[bool_index])
            S.Silent=self.Silent
            if TakeEnergies:
                S.E = self.E[bool_index]
                S.Ep = self.Ep[bool_index]
                S.Ek = self.Ek[bool_index]
                S.UptoDate=True
            if hasattr(self,"rho"):
                S.rho =  self.rho[bool_index]
            return S
        else:
            ref_array = nStars if nStars is not None else indStars
            n,m,x,y,z,vx,vy,vz=[np.zeros(len(ref_array)) for k in range(8)]
            ind_arr = []
            if nStars is None and indStars is not None:
                if max(indStars) > self.N-1:
                    raise Exception("Bad star index.")
                ind_arr = indStars
            elif nStars is not None and indStars is None:
                n = list(self.n)
                ind_arr = np.nonzero(np.in1d(n,nStars))[0]
                if len(ind_arr) != len(nStars):
                    raise Exception("At least one n from nStars was not found.")
            else:
                raise Exception("Please provide bool_index, nStars or indStars")
            S =  PySnap(self.t, self.n[ind_arr], self.m[ind_arr],
                        self.x[ind_arr], self.y[ind_arr], self.z[ind_arr],
                        self.vx[ind_arr], self.vy[ind_arr], self.vz[ind_arr])
            S.Silent=self.Silent
            if TakeEnergies:
                S.E = self.E[ind_arr]
                S.Ep = self.Ep[ind_arr]
                S.Ek = self.Ek[ind_arr]
                S.UptoDate=True
            if hasattr(self,"rho"):
                S.rho =  self.rho[ind_arr]
            return S


    def DeleteParticles(self, bool_index=None, nStars=None, indStars=None ):
        """
        Delete a subset of stars from the snap. Can be specified through bool_index,
        nStars, indStars. See the documentation of PySnap.SelectParticles.
        """
        if bool_index is not None:
            indices = np.nonzero( bool_index )[0]
        if nStars is not None:
            n = list(self.n)
            indices= np.nonzero(np.in1d(n,nStars))[0]
            if len(ind_arr) != len(nStars):
                raise Exception("At least one n from nStars was not found.")
        if indStars is not None:
            indices = indStars
        arrays = [self.x, self.y, self.z, self.vx, self.vy, self.vz, 
                  self.m, self.n]
        [ self.x, self.y, self.z, self.vx, self.vy, self.vz, 
          self.m, self.n] = [ np.delete(arr,indices) 
                              for arr in arrays ]
        self.UptoDate=False
        self.N = len(self.x)

    def AddParticles(self, rp, vp, mp, nstars=None):
        """
        Add particles. rp and vp are 2d numpy array [N,DIM], of
        position and velocities, mp is a 1d numpy array for masses.
        You can provide the new names you want for these stars in the
        resulting snap.
        """
        pos_arr = [self.x,self.y,self.z]
        self.x, self.y, self.z = [ np.hstack((arr,np.array(rp)[:,i])) 
                                   for i,arr in enumerate(pos_arr) ]
        vel_arr = [self.vx,self.vy,self.vz]
        self.vx,self.vy,self.vz = [ np.hstack((arr,np.array(vp)[:,i])) 
                                    for i,arr in enumerate(vel_arr) ]
        self.m = np.hstack((self.m, mp))
        if nstars is None:
            new_n = range(self.n.max()+1, self.n.max()+ 1 + len(mp) )
        else:
            if len(np.nonzero(np.in1d(self.n, nstars))[0])!=0:
                raise Exception("Requested identity for added particle already taken")
            new_n = nstars
        self.n = np.hstack((self.n, new_n))
        self.N = self.N + len(mp)
        self.UptoDate=False

    def FuseWith(self,pysnap):
        """
        Fuse the PySnap with another one.
        """
        if not isinstance(pysnap,PySnap):
            raise Exception("Error, you tried to fuse a PySnap with something"
                            " that was not a PySnap. Go to your room and think"
                            " about what you've done.")
        self.AddParticles(pysnap.r(), pysnap.v(), pysnap.m)


    def SubSnap(self,path):
        """
        This method is some fractal magic. I'm too lazy to write a proper docstring.
        """
        if not hasattr(self,"Tree"):
            raise Exception("Nope. This method only works on fractal models.")
        T = self.Tree
        ind = T.Leaf.FillIndices(T.nlevel)
        for p in path:
            ind = ind[p]
        for j in range(T.nlevel-1-len(path)):
            ind = [ item for sublist in ind for item in sublist]
        return self.SelectParticles(nStars = ind)


    #=====================================================================================
    # Various interesting stuff
    #=====================================================================================


    def kdTree(self):
        """
        Return a KDtree (see kdtree module)
        """
        return K.Tree(self.x, self.y, self.z)

    
    def PCA(self,SpatialCorrection=False):
        """
        Perform a Principal Component Analysis matrix for the snap.
        Return the eigen values of the matrix.
        """
        if SpatialCorrection:
            s= deepcopy(self)
            s.CorrectionCenterOfMass()
        else:
            s = self
        m = np.matrix(s.r())
        M = m.transpose() * m
        return eig(M)[0]
        
       
    def ConvertToPhysical(self, Mfactor, Lfactor):
        """
        Convert all units into solar mass, pc and pc/Myr. Convert self.t in Myr as well. 
        The input is the total stellar mass of the system and the length conversion 
        (how much parsec per unit length). 
        Example:
        15000 stars with a mean stellar mass of 0.5 and 1 length u = 1 parsec
        S.ConvertToPhysical(7500, 1.)
        """
        newG = 4.49e-3
        self.Tfactor = np.sqrt( Lfactor**3 / ( newG * Mfactor ) )
        self.scale_m( Mfactor )
        self.scale_r( Lfactor )
        self.scale_v ( Lfactor / self.Tfactor )
        self.t = self.Tfactor * self.t
        self.G = newG
