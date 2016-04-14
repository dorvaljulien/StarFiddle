import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import cPickle as pk
import collections
import os
import subprocess
import datetime
import distutils
from copy import deepcopy
import h5py

from nb6.input_list import input_variables
from nb6.pysnap import PySnap
import animation3d as anim3d
import anim_binary_histogram as animbin
import anim_binary_parameters as animparam
import cluster_models as CM
from nb6.miscnb6 import read_run, snapname
from binary_motion import RandomBinary
import binaries

AU_to_parsec=4.84e-6 


class ReadRun:
    """
    Contains all Nbody6 snapshots from a given run. It can be 
    initialized in two ways:

    - From the path at which the snapshots are found:
    R = ReadRun( path = "path/to/run/" )

    - From an actual Pysnap array:
    R = ReadRun.FromData( S = snaparray )
    """

    def __init__(self, path=None,Scale=True, Silent=False):
        if path is not None:
            self.S = read_run(path,Silent=Silent)
            if os.path.exists(os.path.join(path,"scaling")):
                with open(os.path.join(path,"scaling"), "r") as f:
                    self.scaling = float(f.read())
            else:
                self.scaling = 1
            self.snaprange = [0,len(self.S)]
    
    @classmethod
    def FromData(self, S, Time=False):
        R = ReadRun(None)
        R.S=S
        if Time:
            for i,s in enumerate(R.S):
                s.t = i
        R.scaling = None
        R.snaprange = [0, len(S)-1]
        return R

    # The following makes the readrun instance behave like a list
    # of PySnap: R[4] is actually R.S[4] and len(R) is len(R.S)
    def __getitem__(self,i):
        return self.S[i]
    def __len__(self):
        return len(self.S)

    #=================================================================================
    # Time, Scaling.
    #=================================================================================

    def t(self):
        """Array of all snapshot times"""
        return [s.t for s in self.S]

    def Scale_r(self,factor):
        """
        Apply a factor to all positions, then scale time and velocities to
        preserve the dynamical state.
        """
        for s in self.S:
            s.t = s.t * factor**(3./2)
            s.scale_r( factor )
            s.scale_v(1/np.sqrt(factor))
            s.UptoDate = False

    def Scale_t(self,factor):
        """
        Apply a factor to time, then scale positions and velocities to
        preserve the dynamical state.
        """
        for s in self.S:
            s.t = s.t * factor
            s.scale_r( factor**(-3./2) )
            s.scale_v(1/np.sqrt( factor**(-3./2)) )
            s.UptoDate = False

    #=================================================================================
    # Animation
    #=================================================================================
    
    def Animation(self,symbol='bo',markersize=2, delay=20,
                  SliderTime=True,
                  FromAnim=None,
                  **kwargs):
        """
        Launch an interactive 3d plot of the cluster evolution.
        
           R1.Animation()

        The FromAnim keyword allows to add another run, on the
        same animation:
        
           R2.Animation(symbol="rx", FromAnim = R1.Anim)

        If SliderTime is False, the slider of the animation will
        display the number of snaps.
        """        
        if FromAnim is None:
            self.Anim=anim3d.Animation(self.S, **kwargs)
            self.Anim.symbol=symbol
            self.Anim.markersize=markersize
            self.Anim.delay=delay
            self.Anim.BackgroundColor='white'
            self.Anim.launch()
        else:
            self.Anim=anim3d.AddAnimation(FromAnim,self.S, **kwargs)
            self.Anim.symbol=symbol
            self.Anim.markersize=markersize
            self.Anim.launch()
        plt.show()


    #=================================================================================
    # Selection
    #=================================================================================

    def SelectParticles(self, nStars, DontCorrect=True, TakeEnergies=True):
        """
        Take a list of names of stars, then return a Run instance of only these stars
        
        DontCorrect : Boolean, if false, correct center of mass and velocity
        TakeEnergies: Boolean, if false, leave out stars energies. Useful as
                      if energies are not computed yet in the source snaps, 
                      they will be as the function tries to access them-> can 
                      slow down the process quite a bit.
        """
        newS=[]
        for i,snap in enumerate(self.S):
            ind_arr = np.nonzero(np.in1d(snap.n,nStars))[0]
            if len(ind_arr) != len(nStars):
                raise Exception("At least one n from nStars was not found.")
            newS.append(PySnap(snap.t,snap.n[ind_arr], snap.m[ind_arr], 
                               snap.x[ind_arr], snap.y[ind_arr], snap.z[ind_arr], 
                               snap.vx[ind_arr], snap.vy[ind_arr], snap.vz[ind_arr]))
            if TakeEnergies:
                newS[-1].E = snap.E[ind_arr]
                newS[-1].Ep = snap.Ep[ind_arr]
                newS[-1].Ek = snap.Ek[ind_arr]
                newS[-1].UptoDate=True
        if not DontCorrect:
            for S in newS:
                S.CorrectionCenterOfMass()
                S.CorrectionCenterOfVelocity()
        return ReadRun.FromData(newS)

    #=================================================================================
    # Data writing
    #=================================================================================    

    def WriteHDF5(self, filename, precision="d", IncludeEnergies=True, Densities =False):
        """
        Writes the run as a HDF5 file.

           R.WriteHDF5( "path/to/file") -> writes path/to/file.hdf5 (automatic extension)

        precision       : default is double precision "d". Single precision is "f"
        IncludeEnergies : Boolean. Particles energies should be included in the file. If
                          energies were not computed yet, they will be as the function
                          tries to access them. Can be very slow.
        Densities       : Boolean. If True, local densities for each particles are included
                          (computed for 6 neighbours, this is hardcoded for now).
        """
        f = h5py.File(filename+".hdf5", "w")
        f.attrs["nstep"]=len(self)
        f.attrs["precision"]=precision
        for i,s in enumerate(self.S):
            grp = f.create_group(str(i))
            grp.attrs["N"] = s.N
            grp.attrs["t"] = s.t
            n_set = grp.create_dataset("n", (s.N,), dtype="i");  n_set[...] = s.n
            m_set = grp.create_dataset("m", (s.N,), dtype=precision);  m_set[...] = s.m
            x_set = grp.create_dataset("x", (s.N,), dtype=precision); x_set[...] = s.x
            y_set = grp.create_dataset("y", (s.N,), dtype=precision); y_set[...] = s.y
            z_set = grp.create_dataset("z", (s.N,), dtype=precision); z_set[...] = s.z
            vx_set = grp.create_dataset("vx", (s.N,), dtype=precision); vx_set[...] = s.vx
            vy_set = grp.create_dataset("vy", (s.N,), dtype=precision); vy_set[...] = s.vy
            vz_set = grp.create_dataset("vz", (s.N,), dtype=precision); vz_set[...] = s.vz
            if IncludeEnergies:
                E_set =grp.create_dataset("E" ,(s.N,), dtype=precision); E_set[...]  = s.E
                Ep_set=grp.create_dataset("Ep",(s.N,), dtype=precision); Ep_set[...] = s.Ep
                Ek_set=grp.create_dataset("Ek",(s.N,), dtype=precision); Ek_set[...] = s.Ek
            if Densities:
                rho_set = grp.create_dataset("rho", (s.N,), dtype=precision)
                rho_set[...] = s.GetDensities(6,Silent=True)
        f.close()
        return filename+".hdf5"




def ReadHDF5(filename):
    """
    Reads a hdf5 file and returns a ReadRun instance:
       R = ReadHDF5( "path/to/file.hdf5" )
    """
    path= "/"+os.path.join( *filename.split("/")[:-1] )
    f =  h5py.File(filename, "r")
    nstep = f.attrs["nstep"]; precision = f.attrs["precision"]
    S = []
    try:
        if precision is "f":
            print "Warning, precision is 32 bits. Binaries business is risky."
    except: pass
    print "Reading "+filename
    for i in range(nstep):
        i = str(i)
        t = f[i].attrs["t"]; n = f[i]["n"];  m = f[i]["m"]
        x,y,z    = [ f[i][name] for name in ["x","y","z"] ]
        vx,vy,vz = [ f[i][name] for name in ["vx","vy","vz"] ]
        Rewrite = False
        S.append( PySnap(t,n,m,x,y,z,vx,vy,vz) )
        if "E" in f[i]:
            S[-1].E  = np.array(f[i]["E"])
            S[-1].Ep = np.array(f[i]["Ep"])
            S[-1].Ek = np.array(f[i]["Ek"])
            S[-1].UptoDate = True
        if "rho" in f[i]:
            S[-1].rho = np.array(f[i]["rho"])
    f.close()
    R = ReadRun.FromData(S=S)
    return R

