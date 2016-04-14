"""
A "density scatter plot" is a scatter plot with points displaying a color
corresponding to the local "point density".
"""

import matplotlib.pyplot as plt, numpy as np, numpy.random, scipy
import cubehelix
from kdtree import Tree
log=np.log10


def DensityScatter(xdat,ydat,ax=None,NNb=15,Nbins=20,logx=False,logy=False,**kwargs):
    """
    ax = DensityScatter( xdat, ydat,ax=None, NNb=15, Nbins=20,
                         logx=False, logy=False, **kwargs)
    ------------------------------------------------------------
    xdat      : data array for x 
    ydat      : data array for y 
    ax        : If needed, previously existing matplotlib axis object
    Nnb       : Number of neighbour points to compute local density
    Nbins     : Number of density(colour) bins
    logx      : Boolean, do you want xdata to be displayed on a logscale ?
    logy      : Boolean, do you want ydata to be displayed on a logscale ?
    **kwargs  : Means any additionnal keyword will be passed to plt.plot

    Display a scatter plot of xdat, ydat and attribute colors to points 
    according to the local point density. Allows to visualize the distribution
    of points in high density regions without doing an histogram2d.
    """
    N=len(xdat)
    xdat = np.array(xdat); ydat = np.array(ydat)
    X = (xdat - min(xdat))/(max(xdat) - min(xdat))
    Y = (ydat - min(ydat))/(max(ydat) - min(ydat))
    if logx:
        X = log(xdat/max(xdat))
    if logy:
        Y = log(ydat/max(ydat))
    T = Tree(X, Y)
    density = np.zeros(N)
    def ComputeDensity(nb,d):
        return nb/( np.pi*d**2)
    for i in range(N):
        _, dist = T.nnearest(i, NNb)
        density[i] = ComputeDensity(NNb,dist[-1])
    density_bins = np.logspace( 0.5*(log(min(density))+log(max(density))), log(max(density)), 
                                Nbins)
    density_bins = np.array( [0] + list(density_bins) )
    SelectionIndices = []
    for i in range(Nbins):
        ind_arr = np.nonzero(( density_bins[i] <= density ) * ( density < density_bins[i+1]))[0]
        SelectionIndices.append(ind_arr)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    cm = plt.get_cmap("rainbow")
    for i,(ind,alph) in enumerate(zip(SelectionIndices,np.linspace(1.,0.,Nbins))):
        color = cm(1.*(i)/Nbins)
        ax.plot( xdat[ind], ydat[ind], "o", color=color, alpha=alph, 
                 markeredgecolor="none",**kwargs)
    if logy:
        ax.set_yscale("log")
    if logx:
        ax.set_xscale("log")
    return ax



def DensityScatter3D(xdat,ydat,zdat,ax=None,NNb=15,Nbins=20,**kwargs):
    """
    ax = DensityScatter3D( xdat, ydat, zdat, ax=None, NNb=15, Nbins=20, **kwargs)
    ------------------------------------------------------------
    xdat      : data array for x 
    ydat      : data array for y 
    zdat      : data array for z 
    ax        : If needed, previously existing matplotlib axis object
    Nnb       : Number of neighbour points to compute local density
    Nbins     : Number of density(colour) bins
    **kwargs  : Means any additionnal keyword will be passed to plt.plot

    Display a 3d scatter plot of xdat, ydat, zdat and attribute colors to points 
    according to the local point density. It's kind of experimental, I played with
    transparency and order or display to be able to see what's going on in high density
    regions.
    """
    N=len(xdat)
    xdat = np.array(xdat); ydat = np.array(ydat)
    X = (xdat - min(xdat))/(max(xdat) - min(xdat))
    Y = (ydat - min(ydat))/(max(ydat) - min(ydat))
    Z = (zdat - min(zdat))/(max(zdat) - min(zdat))
    T = Tree(X, Y, Z)
    density = np.zeros(N)
    def ComputeDensity(nb,d):
        return nb/( 4./3 *np.pi*d**3)
    for i in range(N):
        _, dist = T.nnearest(i, NNb)
        density[i] = ComputeDensity(NNb,dist[-1])
    density_bins = np.logspace( 0.5*(log(min(density))+log(max(density))), log(max(density)), 
                                Nbins)
    density_bins = np.array( [0] + list(density_bins) )
    SelectionIndices = []
    for i in range(Nbins):
        ind_arr = np.nonzero(( density_bins[i] <= density ) * ( density < density_bins[i+1]))[0]
        SelectionIndices.append(ind_arr)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    cm = plt.get_cmap("rainbow")
    for i,(ind,alph) in enumerate(zip(SelectionIndices,np.linspace(1.,0.,Nbins))):
        color = cm(1.*(i)/Nbins)
        ax.plot( xdat[ind], ydat[ind], zdat[ind], "o", color=color, alpha=alph, 
                 markeredgecolor="none",**kwargs)
    return ax


















