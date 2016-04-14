"""
This program generates positions of stars born through a fractal tree 
(see Goodwin, Whitworth 2004)
"""

import numpy as np
from random import *
from matplotlib import pyplot as plt
import itertools
from copy import deepcopy

from nb6.pysnap import PySnap
import cluster_models as CM
from cluster_models import salpeter
from energy import virialise



class Leaf(object):
    """
    Basis unit of the fractal tree. 
    Contains a recursive self-replicating method (GrowNextLeaves)
    """
    coords = np.array([[ 1, 1,-1], [-1, 1, 1], [-1,-1, 1], [ 1,-1,-1],
                       [ 1,-1, 1], [-1, 1,-1], [ 1, 1, 1], [-1,-1,-1]])
    def __init__(self,position,level):
        self.position = position
        self.level = level
        self.LastLeaf = True
        self.mass = None
        self.velocity = None
    def GenerateNewLeafPosition(self,length,noise,C):
        return [self.position[0] + C[0] * length/4. + np.random.normal(0,noise*length),
                self.position[1] + C[1] * length/4. + np.random.normal(0,noise*length),
                self.position[2] + C[2] * length/4. + np.random.normal(0,noise*length)]
    def GenerateNewLeaves(self,length,dimension,noise):
        """
        Computes 8 random numbers, compares to the creation probability that
        is computed from the fractal domension, then create child leaves to the
        current leaves.
        """
        P = np.random.random(8) < 2**(dimension-3)
        self.ChildLeaves = []
        for C in Leaf.coords[P]:
            NewPositions = self.GenerateNewLeafPosition(length,noise,C)
            self.ChildLeaves.append( Leaf(NewPositions, self.level+1))
        if len(self.ChildLeaves) is not 0:
            self.LastLeaf = False
    def GrowNextLeaves(self,nlevel,length,noise,dimension):
        """
        Create child leaves and go into child leaves to create
        "grand-child" leaves, and so on util nlevel layers of 
        leaves are created.
        """
        self.GenerateNewLeaves(length,dimension,noise)
        for C in self.ChildLeaves:
            if C.level < self.level + nlevel:
                C.GrowNextLeaves(nlevel-1,length/2.,noise,dimension)
    def PrintLeaf(self,prefix=""):
        x,y,z =  self.position[0], self.position[1], self.position[2]
        print prefix+" "+str(x)+" "+str(y)+" "+str(z)
        if not self.LastLeaf:
            for L in self.ChildLeaves:
                L.PrintLeaf(prefix=prefix+"   |   ")
    def CollectLeaves(self,level,RequestedLeaves=None):
        """
        Grabs all leaves objects at a specified level
        """
        if RequestedLeaves is None:
            RequestedLeaves =[]
        for C in self.ChildLeaves:
            if C.level == level:
                RequestedLeaves.append(C)
            else:
                RequestedLeaves = C.CollectLeaves(level,RequestedLeaves=RequestedLeaves)
        return RequestedLeaves
    def GetDescendingMass(self,TotalMass=None):
        """
        Go down the fractal tree, adding any mass present in leaves at each level.
        """
        if TotalMass is None:
            TotalMass = 0
        LocalMass = 0
        for C in self.ChildLeaves:
            if C.mass is not None:
                LocalMass += C.mass
            else:
                LocalMass += C.GetDescendingMass(TotalMass=LocalMass)
        self.mass = LocalMass
        return LocalMass
    def GetChildPositions(self):
        r = []
        for C in self.ChildLeaves:
            r.append(C.position)
        return np.array(r)
    def GetChildMasses(self):
        m = []
        for C in self.ChildLeaves:
            m.append(C.mass)
        return np.array(m)
    def GetChildVelocities(self):
        v = []
        for C in self.ChildLeaves:
            v.append(C.velocity)
        return np.array(v)
    def VirialiseChild(self,Q):
        """Virialise the child leaves velocities to match the resquested Q"""
        if self.level == 0:
            self.velocity = np.array([0.,0.,0.])
        if not self.LastLeaf :
            if len(self.ChildLeaves) is not 1:
                m        = self.GetChildMasses()
                x,y,z    = np.transpose( self.GetChildPositions() )
                vx,vy,vz = np.transpose(np.random.normal(0.,1.0,(len(self.ChildLeaves),3)))
                v = np.array([ [VX,VY,VZ]  for (VX,VY,VZ) in zip(vx,vy,vz)])
                for (vv,C) in zip(v,self.ChildLeaves):
                    C.velocity = vv
                    C.VirialiseChild(Q)
            else:
                self.ChildLeaves[0].velocity = self.velocity
                self.ChildLeaves[0].VirialiseChild(Q)
    def TransmitVelocity(self):
        """Apply inheritance of velocity from parent to child"""
        if not self.LastLeaf:
            for C in self.ChildLeaves:
                C.velocity += self.velocity
                C.TransmitVelocity()
    def FillIndices(self,level,inc=0):
        indices = []
        for C in self.ChildLeaves:
            if C.level == level:
                inc +=1
                indices.append(inc)
            else:
                if len(C.ChildLeaves) is not 0:
                    ind, inc = C.FillIndices(level,inc=inc)
                    indices.append(ind)
        if self.level == 0:
            return indices
        else:
            return indices, inc


class FractalTree(object):
    """
    Tree = FractalTree( nlevel, dimension, noise, length=1.0)
    ---------------------------------------------
      nlevel     : how many layers of leaves
      dimension  : fractal dimension, must be <3. Leaf spawning depends on
                   the probability 2^(dimension-3). Can be seen as a 
                   reversed filling factor: 3 is full filling, no fractality.
      noise      : To avoid grid aspect, some noise is applied to positions
                   at each generation.
      length     : Side length of total system.
      mass_range : Stellar mass function. Same than standard cluster creation
      alpha      : Stellar mass function. Same than standard cluster creation

    Generate a fractal tree made of Leaf objects. The final number of particles
    is uncertain, trial and error is advised to see how fractal tree building work.
    """
    def __init__(self,nlevel,dimension,noise,length = 1.0,
                 mass_range=[0.35,20], alpha=2.37, AttributeVelocities=True):
        self.nlevel = nlevel
        self.N = 0
        while self.N is 0:
            self.Leaf = Leaf([0,0,0],0)
            self.Leaf.GrowNextLeaves(nlevel,length,noise,dimension)
            self.particles = self.Leaf.CollectLeaves(self.nlevel)
            self.N = len(self.particles)
        self.AttributeMasses(mass_range,alpha)
    def GetVelocities(self):
        self.velocities = [[L.velocity[0],L.velocity[1],L.velocity[2]] 
                           for L in self.particles ]
        return np.array(self.velocities)
    def GetPositions(self):
        self.positions = [ [L.position[0],L.position[1],L.position[2]] 
                           for L in self.particles ]
        return np.array(self.positions)
    def Plot(self,**kwargs):
        P = self.GetPositions()
        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(111,projection="3d")
        ax.plot(P[:,0],P[:,1],P[:,2],"o",markersize=2,**kwargs)
        plt.tight_layout()
        plt.show()
    def AttributeMasses(self,mass_range=[0.2,20],alpha=2.37):
        masses = CM.salpeter(self.N,alpha,mass_range)
        self.O_masses = deepcopy(masses)
        self.masses = masses/masses.sum()
        for m,p in zip(masses,self.particles):
            p.mass = m
        self.Mt = masses.sum()
        self.Leaf.GetDescendingMass()
    def Snap(self,AttributeVelocities=True):
        """
        Convert Fractal tree to PySnap.
        """
        if AttributeVelocities:
            self.Leaf.VirialiseChild(0.5)
            self.Leaf.TransmitVelocity()
            v = self.GetVelocities()
        else:
            v = np.zeros((self.N,3)) 
        r = self.GetPositions()
        S = PySnap(0,range(1,self.N+1),self.masses,r[:,0],r[:,1],r[:,2],
                      v[:,0],v[:,1],v[:,2])
        S.Leaf = self.Leaf
        return S




def FractalModel(nlevel=5, dimension=2.0, noise=0.5, 
                 length=1.0, mass_range=[0.2,20],Velocities=True):
    """
    Sf = FractalTree( nlevel, dimension, noise, length=1.0)
    ---------------------------------------------
      nlevel     : how many layers of leaves
      dimension  : fractal dimension, must be <3. Leaf spawning depends on
                   the probability 2^(dimension-3). Can be seen as a 
                   reversed filling factor: 3 is full filling, no fractality.
      noise      : To avoid grid aspect, some noise is applied to positions
                   at each generation.
      length     : Side length of total system.
      mass_range : Stellar mass function. Same than standard cluster creation
      alpha      : Stellar mass function. Same than standard cluster creation

    Generate PySnap of a fractal model, then virialise it. 
    The final number of particles is uncertain, trial and error is advised to see
    how fractal tree building work.
    """
    T = FractalTree(nlevel,dimension,noise,AttributeVelocities=Velocities,mass_range=mass_range)
    print T.N, " particles created"
    S = T.Snap()
    S.Tree = T
    S.virialise(0.5)
    return S














