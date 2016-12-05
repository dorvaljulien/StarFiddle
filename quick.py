"""
Quick launch for the environment:
   from quick import *
"""

from nb6.readrun import ReadRun,ReadHDF5
from nb6.tools import LaunchNbody6
from nb6.miscnb6 import read_snapshot, snapname
import cluster_models as CM
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from fractals import FractalModel
from animation3d import plot3d
import h5py
from density_scatter_plot import DensityScatter
Colors = ["b","r","g","c","m","y","k","w"]


def NBODY6(*args, **kwargs):
    I = LaunchNbody6(*args,**kwargs)
    return ReadRun(I.path, Silent=True)

class Container(object):
    def __init__(self):
        pass

log = np.log10

def SphVol(radius):
    return np.pi*4./3 * radius**3







