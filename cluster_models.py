"""
This module can create PySnap objects for several cluster structures.
Here are a few examples showcasing the different structures and parameters
that can be tuned.

S = Uniform(mass_range=[0.01, 100] )
 | This create a N=1000 (default) uniform sphere model. Masses are
 | chosen from a Salpeter mass function between 0.01 and 100 solar 
 | masses.

S = King(N=10000, virial_ratio = 0.4)
 | The King model is a widely used centrally concentrated model for 
 | star clusters. Here, its virial_ratio is set to 0.4, which is
 | slightly colder than equilibrium (0.5). The cluster, when left
 | to evolve, will contract a little. If virial_ratio had been set 
 | to 0.6, the cluster would have been hotter and would have expanded.

S = Hubble( Hub=1.5 )
 | The Hubble model was created by Julien Dorval and is a way to obtain
 | self-consistent substructures. See Dorval et al, 2016. Here, the Hub
 | parameter tunes the strength of the initial expansion. 1.5 is higher
 | than the critical value 1.4, thus this model will never stop expanding.

S = Plummer( virial_ratio=0.01, truncature=3. )
 | Plummer model is another famous cluster model. It is simpler to build
 | than King. Here it is made very cold (virial_ratio=0.01) and will 
 | collapse when left to evolve. truncature is special to plummer models, 
 | it specifies the level of outliers. Here, stars with less than 3% chance
 | to be spawned are rejected, much stricter than the default 0.1, meaning
 | the cluster outer regions are severely depleted.


If you want to take a look at a cluster you created, just use:
  S.Plot()

"""

import os
import numpy as np
import random
from copy import deepcopy
import subprocess as sub
from scipy.interpolate import interp1d
import cPickle as pk
import ctypes as C
import inspect

import energy
import binary_motion as bm
from nb6.pysnap import PySnap
from nb6.miscnb6 import snapname,read_snapshot

# c library for king model
try:
    lib_path=os.path.join( os.environ["STARFIDDLELIB"], "libking.so")
except KeyError:
    lib_path="/usr/local/lib/libking.so"
king_lib = C.CDLL(lib_path)
double_pointer = C.POINTER(C.c_double)


Pi=np.pi
AU_to_parsec=4.84e-6  # 1 Au in parsec
RSun_parsec=2.2506e-8     # 1 solar radius in parsec



def salpeter(N,alpha=2.37,mass_range=[0.35, 20]):
    """
    masses = salpeter(N, alpha=2.37, mass_range=[0.35, 20])
    ---------------------------------------------------------
    N          : Number of requested masses
    alpha      : Slope of mass funtion
    mass_range : Minimum and maximum mass in solar mass.

    Create a sample of stellar masses from a Salpeter mass function.
    """
    N = int(N)
    m1,m2 = mass_range
    p1 = m1**(1-alpha)
    p2 = m2**(1-alpha)
    P = p1+ np.random.random(N)*(p2-p1)
    return P**(1. / (1-alpha) )

# Default arguments for cluster creation
def make_dict(**kwargs):
    return kwargs
standard_arguments = make_dict(
    N=1000, 
    mass_range=[0.35,20], 
    alpha=2.37,
    virial_ratio=0.5,
    Silent=False,
# King
    W0 = 4,
# Plummer 
    A0=0.2,
    truncature=0.1,
# Hubble,
    Hub=1.00)
    

def AttributeDirections(C):
    """
    Take an array of vector norms, return three arrays of 
    vector coordinates of these norms in random directions.
    """
    X1,X2 = [ np.random.random(len(C)) for i in [1,2] ]
    thetap=np.arccos(1-2*X1)
    phip=2*Pi*X2
    cx = C*np.sin(thetap)*np.cos(phip)
    cy = C*np.sin(thetap)*np.sin(phip)
    cz = C*np.cos(thetap)
    return cx,cy,cz



class Dummie(object):
    """
    Basic structure for the creation of a model
    """
    def __init__(self,**kwargs):
        self.kwargs = dict(standard_arguments.items() + kwargs.items())
        for key in self.kwargs:
            setattr(self,key,self.kwargs[key])
        self.get_masses()
        self.get_positions()
        self.get_velocities()
        self.build_snap()
        self.virialise_snap()
        self.S.o_Mt = self.o_Mt
    def get_masses(self):
        self.m = salpeter(self.N,mass_range=self.mass_range,alpha=self.alpha)
        self.o_Mt = np.sum(self.m)
        self.m = self.m / self.m.sum()
    def get_positions(self):
        pass
    def get_velocities(self):
        pass
    def build_snap(self):
        self.S=PySnap(0,range(1,self.N+1),
                      self.m,self.x,self.y,self.z,
                      self.vx,self.vy,self.vz)
        self.S.Silent= self.Silent
        self.S.o_Mt = self.o_Mt
    def virialise_snap(self):
        if self.virial_ratio >= 0:
            if not self.Silent:
                print "virialising to ",self.virial_ratio
            self.S.CorrectionCenterOfMass()
            self.S.CorrectionCenterOfVelocity()
            Q = self.S.virialise(self.virial_ratio)
            if np.isnan(Q):
                raise Exception("Error during the virialization.")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # #        UNIFORM        # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def Uniform_GetR(N):
    """
    Return a radius between 0 and 1 according to a R^2 distribution.
    """
    result = []
    for i in range(N):
        R,P=0,1
        while P>R**2:  R,P=np.random.random(2)
        result.append(R)
    return result

def get_uniform_coordinates(self):
    """
    Replace get_positions method in Dummie class
    """
    X1,X2 = [ np.random.random(self.N) for i in [1,2] ]
    theta=np.arccos(1-2*X1)
    phi=2*Pi*X2
    self.R= Uniform_GetR(self.N)
    self.x,self.y,self.z = AttributeDirections(self.R)

class CreateUniform(Dummie):
    get_positions = get_uniform_coordinates
    def get_velocities(self):
        # We populate the velocity space
        V = self.R * np.random.random(self.N)
        self.vx, self.vy, self.vz = AttributeDirections(V)

def Uniform(**kwargs):  
    """
    S = Uniform( N=1000, mass_range=[0.35,20], alpha=2.37,
                 virial_ratio=0.5, Silent=False)
    ---------------------------------------------------
    N            : Number of particles
    mass_range   : Min and max  mass for the salpeter mass function
    alpha        : Slope for the salpeter mass function
    virial_ratio : Initial virial state of the cluster. 0.5
                      is equilibrium, 0 is no velocity. In this model
                      velocities are taken from a uniform probability law

    Create an Uniform model, as a PySnap instance.
    """
    return CreateUniform(**kwargs).S


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # #         KING          # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def King_fE(W):
    return np.exp(W)-1

def King_GetEnergy(Elim):
    """
    Get an energy from a King energy distribution from 0 to specified Elim.
    """
    result = []
    for i in range(len(Elim)):
        Ws,P=0,1000
        while(P > King_fE(Ws)):
            Ws=np.random.random()*Elim[i]
            P=1.05*np.random.random()*King_fE(Elim[i])
        result.append(Ws)
    return result

def King_GetPotential(W0):
    """ 
    radius, potential, mass = King_GetPotential(W0)
    Calls the c program to compute the potential and cumulated 
    mass along the radius. W0 is the King concentration parameter.
    """   
    king_lib.King_model.argtypes = [C.c_double, C.c_double, 
                                    C.c_double, C.c_int,
                                    double_pointer, double_pointer, 
                                    double_pointer]
    king_lib.King_model.restype = C.c_int
    size = 1000; # should be enough
    radius = (size*C.c_double)()
    mass = (size*C.c_double)()
    potential = (size*C.c_double)()
    n = king_lib.King_model(W0, 1e-5, 0.039414, size, radius, potential, mass)
    radius = np.asarray( radius[:n], dtype=np.double)
    potential = np.asarray( potential[:n], dtype=np.double)
    mass = np.asarray( mass[:n], dtype=np.double)
    return radius, potential, mass

class CreateKing(Dummie):
    def get_positions(self):
        self.rad,self.pot,self.mass=King_GetPotential(self.W0)
        R_m=interp1d(self.mass,self.rad,kind="slinear")
        Psi=interp1d(self.rad,self.pot,kind="slinear")
        mass_min,mass_max=self.mass[0],self.mass[-1]
        Dm=mass_max-mass_min
        self.R = R_m(mass_min+Dm*np.random.random(self.N))
        self.x,self.y,self.z = AttributeDirections(self.R)
    def get_velocities(self):
        Psi=interp1d(self.rad,self.pot,kind="slinear")
        Psi_local=Psi(self.R)
        W = King_GetEnergy(Psi_local)
        self.V =  np.sqrt(2*(Psi_local - W))
        self.vx, self.vy, self.vz = AttributeDirections(self.V)

def King(**kwargs): 
    """
    S = King( N=1000, mass_range=[0.35,20], alpha=2.37,
                 W0 = 4, virial_ratio=0.5, Silent=False)
    ------------------------------------------------------
       N            : Number of particles
       mass_range   : Min and max  mass for the salpeter mass function
       alpha        : Slope for the salpeter mass function
       W0           : King concentration parameter   
       virial_ratio : Initial virial state of the cluster. 0.5
                      is equilibrium, 0 is no velocity. King model
                      have their own specific velocity distribution.

    Create a King model, as a PySnap instance. Calls a c code to compute 
    the King potential, based on a Fortran code by Gerry Quinlan, modified 
    by Christian Boily and converted to C by Julien Dorval.
    """
    return CreateKing(**kwargs).S



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # #       HUBBLE          # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class CreateHubble(Dummie):
    get_positions = get_uniform_coordinates
    def get_velocities(self):
        self.vx = self.Hub * self.x
        self.vy = self.Hub * self.y
        self.vz = self.Hub * self.z
    def virialise_snap(self):
        pass

def Hubble(**kwargs): 
    """
    S = King( N=1000, mass_range=[0.35,20], alpha=2.37,
              Hub = 1., Silent=False)
    ------------------------------------------------------
       N            : Number of particles
       mass_range   : Min and max  mass for the salpeter mass function
       alpha        : Slope for the salpeter mass function
       Hub          : Hubble parameter v = Hub * r

    Create a Hubble model, as a PySnap instance. Hubble model are uniform 
    sphere with radial velocities following a Hubble  velocity field: 
        v = Hub * r 
    There is no virial ratio argument for Hubble models.
    """
    return CreateHubble(**kwargs).S



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # #       PLUMMER         # # # # # # # # # # # # # 
# # # # # # # # # # #                       # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def PlummerR(N,A0,truncature):
    #We reject stars with a radius less than truncature% likely
    result = []
    for i in range(N):
        R,P=1,1
        while P>(1 + (R/A0)**2 )**(-5./2):
            R=np.random.random()*A0*np.sqrt(1+(100./truncature)**(2./5))
            P=np.random.random()
        result.append(R)
    return np.array(result)

def PlummerV(N,R,A0):
    vmax=np.sqrt(2)*(R+A0)**(-1./4)
    V = []
    for i in range(N):
        q,P=1,0
        while P>(1-q)**(7./2)*q**2:
            q,P=np.random.random(2)
        V.append(q*vmax[i])
    return np.array(V)

class CreatePlummer(Dummie):
    def get_positions(self):
        self.R = PlummerR(self.N,self.A0,self.truncature)
        self.x,self.y,self.z = AttributeDirections(self.R)
    def get_velocities(self):
        self.V = PlummerV(self.N, self.R, self.A0)
        self.vx, self.vy, self.vz = AttributeDirections(self.V)

def Plummer(**kwargs): 
    """
    S = Plummer( N=1000, mass_range=[0.35,20], alpha=2.37,
              A0 = 2,  truncature = 0.1,
              virial_ratio=0.5, Silent=False)
    -----------------------------------------------------------------------
       N            : Number of particles
       mass_range   : Min and max  mass for the salpeter mass function
       alpha        : Slope for the salpeter mass function
       A0           : Plummer parameter
       truncature   : Particles with less than truncature % chance to spawn
                      are rejected, to avoid outliers.
       virial_ratio : Initial virial state of the cluster. 0.5
                      is equilibrium, 0 is no velocity.

    Create a Plummer model, as a PySnap instance.
    """
    return CreatePlummer(**kwargs).S
