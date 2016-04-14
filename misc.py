"""
This contains several useful functions.
"""

import numpy as np


def crossmatch(A,b):
    """
    Returns an array of indices such as:
      ind = whereare(a,b)
    means that a[ind] == b
    /!\ a should have unique occurences of b elements.
    """
    inda = np.arange(len(A))
    ind = np.argsort(A)
    a = A[ind]; inda=inda[ind]
    ind = np.nonzero(np.in1d(a,b))[0]
    if len(ind)<len(b):
        raise Exception("Some element of second array was not found"
                        " in first array")
    elif len(ind)>len(b):
        raise Exception("First array may have not unique elements.")
    inda= inda[ind]; a = a[ind]
    indb = b.argsort().argsort()
    return inda[indb]


def Log3DPlot(x,y,z, *args, **kwargs):
    """
    Produces a 3d scatter plot in log scale.
    INPUT: 
    - 3 arrays x, y, z of same length
    - 'Xrange', 'Yrange' and 'Zrange' keywords to specify the plot ranges
      ( they will be automatically computed if not specified )
    - Any other simple argument or keyword argument will be passed to the plot function
    OUTPUT:
      Returns an axis object which you can use to modify further the plot.
    EXAMPLE:
       x,y,z = [ 10e4*np.random.random(500) for i in range(3)]
       Log3DPlot( x,y,z, "bo", Xrange=[100,1e6],Yrange=[100,1e6],Zrange=[100,1e6])
    WARNING:
       This is a rather ugly hack:
       - The ticks won't change when zoom in or out
       - The actual values in the plot are log(x), log(y) and log(z). The ticklabels
         were changed by hand to display the corresponding x,y,z values. Bear that in
         mind when adding data to the plot.
    """
    d = {}
    for key,val in zip( ["Xrange", "Yrange",  "Zrange" ],  [x,y,z] ):
        if kwargs.has_key(key):
            d[key] = kwargs.pop(key)
        else:
            d[key] = [ min(val), max(val) ]
    Xrange, Yrange, Zrange = d["Xrange"], d["Yrange"], d["Zrange"]
    X,Y,Z = [ np.log10(a) for a in [x,y,z] ]
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")
    ax.plot(X, Y, Z,  *args, **kwargs)
    ax.set_xlim( np.log10(Xrange[0]), np.log10(Xrange[1]) )
    ax.set_ylim( np.log10(Yrange[0]), np.log10(Yrange[1]) )
    ax.set_zlim( np.log10(Zrange[0]), np.log10(Zrange[1]) )
    Xlogs = range(int(math.ceil(np.log10(Xrange[0]))),int(math.floor(np.log10(Xrange[1]))+1))
    Ylogs = range(int(math.ceil(np.log10(Yrange[0]))),int(math.floor(np.log10(Yrange[1]))+1))
    Zlogs = range(int(math.ceil(np.log10(Zrange[0]))),int(math.floor(np.log10(Zrange[1]))+1))
    ax.set_xticks( Xlogs )
    ax.set_xticklabels( [ r'$10^{{{val}}}$'.format(val=str(int(l))) for l in Xlogs ] )
    ax.set_yticks( Ylogs )
    ax.set_yticklabels( [ r'$10^{{{val}}}$'.format(val=str(int(l))) for l in Ylogs ] )
    ax.set_zticks( Zlogs )
    ax.set_zticklabels( [ r'$10^{{{val}}}$'.format(val=str(int(l))) for l in Zlogs ] )
    return ax


import matplotlib.pyplot as plt
import numpy as np

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

