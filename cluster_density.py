import numpy as np
import ctypes as C
import os

try:
    lib_path = os.path.join( os.environ["STARFIDDLELIB"], "libdensity.so")
except KeyError:
    lib_path="/usr/local/lib/libdensity.so"
_lib = C.CDLL(lib_path)

double_pointer = C.POINTER(C.c_double)
int_pointer = C.POINTER(C.c_int)



def GetDensity(m,x,y,z, Nnb, Number=False):
    """
    rho = GetDensity(m, x, y, z, Nnb)
    m     : Masses array
    x,y,z : Positions array
    Nnb   : Number of neighbours used to compute the density for each star

    Use a kdtree to get Nnb neighbours of each star, and compute the corresponding
    local density. Return a density array the same length as m,x,y,z.
    """
    _lib.density_distribution.argtypes = [ C.c_int, C.c_int, C.c_int,
                                           double_pointer, double_pointer, 
                                           double_pointer, double_pointer,
                                           double_pointer ]
    _lib.density_distribution.restype = C.c_int
    Number = 1 if Number else 0
    N = len(m)
    mc = (N*C.c_double)(); mc[:] = m; m = mc
    xc = (N*C.c_double)(); xc[:] = x; x = xc
    yc = (N*C.c_double)(); yc[:] = y; y = yc
    zc = (N*C.c_double)(); zc[:] = z; z = zc
    density = ( N*C.c_double )()
    n = _lib.density_distribution( N, Nnb, Number, m, x, y, z, density  )
    density = np.asarray(density[:n] , dtype=np.double) 
    return density
