"""
This module is about computing star kinetic and potential energies. It
uses an algorithm coded and c, accessible by python through Ctypes and
a shared library.
"""

import numpy as np
import ctypes as C
import os
import copy

try:
    lib_path=os.path.join( os.environ["STARFIDDLELIB"],"libenergy.so")
except KeyError:
    lib_path="/usr/local/lib/libenergy.so"

_lib = C.CDLL(lib_path)
double_pointer = C.POINTER(C.c_double)
int_pointer = C.POINTER(C.c_int)


def virialise(Q,x,y,z,vx,vy,vz,M, Standard=True, G=1):
    """
    x,y,z,vx,vy,vz = virialise(Q,x,y,z,vx,vy,vz,M)
    ----------------------------------------------
    Q        : Virial ratio to apply to the system -Ek/Ep
    x,y,z    : Position arrays
    vx,vy,vz : Velocity arrays
    M        : mass array
    
    Take an existing system and scale positions and velocities to
    reach the required virial ratio, while fitting the nbody units:
       Ep = -0.5
    """
    
    _lib.virialise.argtypes = [C.c_int, C.c_double, C.c_int,
                                     double_pointer, 
                                     double_pointer, double_pointer, double_pointer, 
                                     double_pointer, double_pointer,double_pointer,
                                     C.c_double] 
    _lib.virialise.restype = C.c_int
    N = len(M)
    Standard = 1 if Standard else 0
    #We have to convert the numpy arrays of x, y, etc to ctypes arrays:
    xc = (N*C.c_double)(); xc[:] = x; x = xc
    yc = (N*C.c_double)(); yc[:] = y; y = yc
    zc = (N*C.c_double)(); zc[:] = z; z = zc
    vxc = (N*C.c_double)(); vxc[:] = vx; vx = vxc
    vyc = (N*C.c_double)(); vyc[:] = vy; vy = vyc
    vzc = (N*C.c_double)(); vzc[:] = vz; vz = vzc
    Mc = (N*C.c_double)(); Mc[:] = M; M = Mc
    _lib.virialise(N,Q, Standard, M,x,y,z,vx,vy,vz,G)
    x = np.asarray(x[:N], dtype=np.float)
    y = np.asarray(y[:N], dtype=np.float)
    z = np.asarray(z[:N], dtype=np.float)
    vx = np.asarray(vx[:N], dtype=np.float)
    vy = np.asarray(vy[:N], dtype=np.float)
    vz = np.asarray(vz[:N], dtype=np.float)
    return x,y,z,vx,vy,vz



def get_particles_E(x,y,z,V,M,G=1):
    """
    Ek,Ep = get_particles_E(x,y,z,V,M,G=1)
    -----------------------------------------------
    x,y,z : Position arrays
    V     : Velocities norm
    M     : masses

    Compute Kinetic and potential energy of a given system.
    """
    _lib.get_particles_E.argtypes = [C.c_int, double_pointer, 
                                     double_pointer, double_pointer, double_pointer, 
                                     double_pointer, double_pointer,double_pointer,
                                     C.c_double] 
    _lib.get_particles_E.restype = C.c_void_p
    N = len(M)
    #We have to convert the numpy arrays of x, y, etc to ctypes arrays:
    xc = (N*C.c_double)();  xc[:]  = x; x = xc
    yc = (N*C.c_double)();  yc[:]  = y; y = yc
    zc = (N*C.c_double)();  zc[:]  = z;    z = zc
    V2c = (N*C.c_double)(); V2c[:] = V**2; V2 = V2c
    Mc = (N*C.c_double)();  Mc[:]  = M; M = Mc
    Ek, Ep =  [ (N*C.c_double)() for i in range(2) ]
    _lib.get_particles_E(N,M,x,y,z,V2,Ek,Ep,G)
    Ek = np.asarray(Ek[:N], dtype=np.float)
    Ep = np.asarray(Ep[:N], dtype=np.float)
    return Ek,Ep
