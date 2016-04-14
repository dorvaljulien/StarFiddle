import numpy as np
import os
import ctypes as C
import time

from nb6.pysnap import PySnap



def frmt(x):
    return "%15.7e" % x

def dist(r1,r2):
    return np.sqrt( (r1[0]-r2[0])**2 + 
                    (r1[1]-r2[1])**2 + 
                    (r1[2]-r2[2])**2   )

def snapname(nsnap):
    """snapname(17) -> snap_0000017 """
    prefix="snap_"
    suffix="%07d" % nsnap
    name = prefix + suffix
    return name


def get_nfiles(path):
    n_start=1
    nfiles=0
    if os.path.isfile(path+snapname(n_start)):
        while os.path.isfile(path+snapname(n_start+nfiles)):
            nfiles+=1
    else:
        raise Exception("No snap_0000001 file found.")
    return nfiles


def read_header(filename):
    """N, dim, t = read_header(filename)"""
    with open(filename) as f:

        N=int(f.readline().strip())
        dim=int(f.readline().strip())
        t=float(f.readline().strip())

    return N, dim, t




def read_snapshot(filename, NB6_InputFormat=False):
    """
    This function gets the name of a snap_******* file 
    ( ex: snap_0000017), which is a homemade output file for nbody6 
    and returns the data as a PySnap instance.
    P = rs.read_snapshot("snap_0000001")
    """
    if NB6_InputFormat:
        with open(filename) as f:
            N=len(f.readlines())
        m,x,y,z,vx,vy,vz=np.loadtxt(filename,unpack=True)
        dynstate=PySnap(0,range(1,N+1),m,x,y,z,vx,vy,vz)
    else:
        with open(filename) as f:
            try:
                N=int(f.readline().strip())
                dim=int(f.readline().strip())
                t=float(f.readline().strip())
                InputFormat = False
            except ValueError:
                N=len(f.readlines())
                InputFormat = True
        if InputFormat:
            m,x,y,z,vx,vy,vz=np.loadtxt(filename,unpack=True)
            dynstate=PySnap(0,range(1,N+1),m,x,y,z,vx,vy,vz)
        else:
            n,m,x,y,z,vx,vy,vz=np.loadtxt(filename,unpack=True,skiprows=3)
            n = np.asarray(n,dtype=np.int)
            dynstate=PySnap(t,n,m,x,y,z,vx,vy,vz)
    return dynstate




def read_run(path, startat=None, upto=None, Silent=False):
    """
    This function gets the path of a NBODY6 run, reads the snapshots
    and return a list of PySnap objects. Used by ReadRun.
    If you're not using NBODY6, just put your snapshots in a given
    directory "example_dir" as ascii files in the format described in
    PySnap.WriteSnapshot and name them
    snap_0000001
    snap_0000002
    snap_0000003
    snap_0000004
    ...

    Then you can create your run as R = ReadRun("example_dir")

    You can choose to start reading the run at the snap number "startat"
    and end it at "upto".
    """
    firstfile=0
    if path[-1]!="/":
        path=path+"/"
    while not os.path.exists(path+snapname(firstfile)) and firstfile < 1000:
        firstfile+=1
    if firstfile == 1000:
        raise Exception("No snapshot file found at "+path)

    lastfile=firstfile
    while os.path.exists(path+snapname(lastfile+1)):
        lastfile+=1

    if startat is not None:
        if os.path.exists(path+snapname(startat)):
            firstfile = startat
    if upto is not None:
        if os.path.exists(path+snapname(upto)):
            lastfile = upto

    if not Silent:
        print "Reading from file", firstfile, " to file", lastfile
    snaplist=[]
    for n in range(firstfile, lastfile+1):
        nfiles = lastfile - firstfile
        currentfile = n - firstfile
        if nfiles <= 10:
            if not Silent:
                print "  Reading file "+str(currentfile)+" of "+str(nfiles)
            snaplist.append(read_snapshot(path+snapname(n)))
        else:
            if currentfile % ( (nfiles-nfiles%10) / 10 )==0:
                if not Silent:
                    print "  Reading file "+str(currentfile)+" of "+str(nfiles)
            snaplist.append(read_snapshot(path+snapname(n)))

    return snaplist
    

