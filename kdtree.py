"""
Module to create and use KD-trees: 
all you ever ever need for your neighbours searches.

Python code: Julien Dorval
C code: Adapted from c++ to c by Julien Dorval
        original c++ code from Numerical Recipees(3rd edition)

Any questions: dorvaljulien@gmail.com

"""

import numpy as np
import random
import ctypes as C
import matplotlib.pyplot as plt
import time
import struct
import sys
import os
from distutils.util import strtobool

try:
    lib_path=os.path.join( os.environ["STARFIDDLELIB"],"libkdtree.so")
except KeyError:
    lib_path="/usr/local/lib/libkdtree.so"
_lib = C.CDLL(lib_path)



#float_choice = "float"
float_choice = "double"

FLOAT = C.c_double
numpy_float = np.double

FLOAT_pointer=C.POINTER(C.c_double)
int_pointer=C.POINTER(C.c_int)



# -------------------------------------------------------------------
#                          C  Structures 
# -------------------------------------------------------------------
class C_Point(C.Structure):
    pass
C_Point._fields_=[("DIM",C.c_int),("x",FLOAT_pointer)]


class C_Box(C.Structure):
    pass
C_Box._fields_=[("lo",C_Point),("hi",C_Point),
              ("mom",C.c_int),("dau1",C.c_int),("dau2",C.c_int),
              ("ptlo",C.c_int),("pthi",C.c_int)]

class C_Tree(C.Structure):
    pass
C_Tree._fields_=[("N",C.c_int),("DIM",C.c_int),
                 ("pts",C.POINTER(C_Point)),
                 ("ptindx", int_pointer), ("rptindx",int_pointer),
                 ("nboxes",C.c_int), ("boxes",C.POINTER(C_Box))  ]


# -------------------------------------------------------------------
#                          Tree structure 
# ------------------------------------------------------------------
class Tree:
    """

    User manual:
    Create the tree with 
         tree = Tree(x,y,z,..)
    ( works with as many dimensions as you like )

    Or read the tree from a kdtree file (that was written by
    the WriteToFile method) :
         tree = Tree.FromFile("tree.kd")

    A Tree structure contains:
    - N     :   number of particles
    - DIM   :   dimension
    - nboxes:   number of boxes
    - c_tree:   A C_Tree instance mimicking the c Tree structure


    List of methods:
       - Return the coordinates as an array of vectors:
             tree.vectors()

       - Visualize a 2d tree
             tree.draw()

       - Find the closest neighbour of an arbitrary point x:
             nb = tree.nearest(x)

       - Find the n neighbours (identity and distances) of point j:
             indexes, distances = tree.nnearest(j,n)
         the indexes and distances are sorted by decreasing distance

       - Find all points within r of an arbitrary point [x,y,..]:
             nblist, distances = tree.locatenear([x,y,..],r)
         the indexes and distances are sorted by decreasing distance

       - Writes the tree data into the specified binary file.
             tree.WriteToFile("tree.kd")
       
       - Free the memory used by the tree
             tree.FreeMemory()
             -> The memory should be cleared as soon as you quit
                python anyway. However, python can sometimes have
                a hard time with ctypes and things can get funky.
                This method is here as an emergency tool.
 """

    def __init__(self,*vectors):
        if len(vectors)!=0:
            # We explain to python what types c expects
            _lib.KDtree.argtype = [ FLOAT_pointer, C.c_int, C.POINTER(C_Tree) ]
            _lib.KDtree.restype =  C.c_void_p

            DIM=len(vectors)
            N=len(vectors[0])
            coord=(DIM*N*FLOAT)() # This is a c_types pointer
            for i in range(N):
                for ndim,x in enumerate(vectors):
                    coord[N*ndim+i]=x[i]
            c_tree = C_Tree()
            c_tree_p = C.pointer(c_tree)
            _lib.KDtree(coord,N,DIM,c_tree_p)
            self.c_tree=c_tree
            self.N=c_tree.N
            self.DIM=c_tree.DIM
            self.nboxes=c_tree.nboxes
            del coord
        else: # no vectors were provided: empty structure
            pass
        
    @classmethod
    def FromFile(cls,filename):
        """
        Use the c function read_tree to read a kdtree file:
             tree = Tree.FromFile("tree.kd")
        """
        _lib.read_tree.argtypes = [C.c_char_p, C.POINTER(C_Tree)]
        _lib.read_tree.restype = C.c_void_p
        
        # Let's first get useful quantities from the file:
        with open(filename,"rb") as f:
            N = struct.unpack('i',f.read(4))[0]
            DIM = struct.unpack('i',f.read(4))[0]
        c_tree = C_Tree()
        c_tree_p = C.pointer(c_tree)
        _lib.read_tree(filename, c_tree_p)
        tree = Tree()
        for a in ["N","DIM","nboxes","boxes"]:
            setattr(tree,a,getattr(c_tree,a))
        tree.c_tree = c_tree
        return tree
    
    def vectors(self):
        coords = np.zeros((self.N,self.DIM))
        for i in range(self.N):
            coords[i,:]=[ self.c_tree.pts[i].x[j] 
                          for j in range(self.DIM)]
        return [coords[:,i] for i in range(self.DIM)]


    def draw(self, **kwargs):
        """Visualize a 2d tree. Accept standard plot keyword arguments."""
        if self.DIM != 2: 
            raise Exception("Cannot draw non 2D tree." 
                            "Sorry about that.")
        if self.N > 1000:
            R = user_yes_no_query("Are you sure you want to draw a " 
                                  "tree with "+str(self.N)+" points ?")
            if not R:
                return
        x,y = self.vectors()    
        plt.plot(x,y,'bo',**kwargs)
        for i in range(self.c_tree.nboxes): 
            draw_box(self.c_tree.boxes[i])
        M=1.05*np.max(abs(self.coords))
        plt.xlim(-M,M)
        plt.ylim(-M,M)
        plt.show()

    def nearest( self , point ):
        """
        Find the closest neighbour of an arbitrary point x:
            nb = tree.nearest(x)
        """
        cpoint=C_Point()
        cpoint.DIM = len(point)
        cpoint.x= (len(point)*FLOAT)(*point)
        _lib.nearest.argtype = [ C_Tree, C_Point  ]
        _lib.nearest.restype = C.c_int
        return _lib.nearest( self.c_tree, cpoint )

    def nnearest( self, n_point, n_neighbours ):
        """
        Find the n neighbours (identity and distances) of point j:
            indexes, distances = tree.nnearest(j,n)
        """
        _lib.nearest.argtype = [ C_Tree, C.c_int, int_pointer, 
                                 FLOAT_pointer, C.c_int ]
        _lib.nnearest.restype = C.c_void_p
        ind_neighbours=(n_neighbours*C.c_int)()
        distances=(n_neighbours*FLOAT)()
        _lib.nnearest(self.c_tree, n_point, 
                      ind_neighbours, distances, n_neighbours)
        ind = np.asarray(ind_neighbours[:int(n_neighbours)], 
                         dtype=np.int)
        dist = np.asarray(distances[:int(n_neighbours)], 
                          dtype=numpy_float)
        return ind, dist


    def locatenear( self, point, r, nmax=100000):
        """
        Find all points within r of an arbitrary point x:
             nblist, dist = tree.locatenear(x,r)
        """
        cpoint=C_Point()
        cpoint.DIM = len(point)
        cpoint.x= (len(point)*FLOAT)(*point)

        _lib.locatenear.argtype = [ C_Tree, C_Point, 
                                    FLOAT, 
                                    int_pointer, FLOAT_pointer,
                                    C.c_int  ]
        _lib.locatenear.restype = C.c_int

        nb_list=(self.N * C.c_int)()
        dist_list=(self.N * FLOAT)()
        r = (FLOAT)(r) 
        n_nb=_lib.locatenear(self.c_tree, 
                             cpoint, r, nb_list, dist_list, nmax )
        nb_list= np.asarray( nb_list[:n_nb], dtype=np.int)
        dist_list= np.asarray( dist_list[:n_nb], dtype=numpy_float)
        return nb_list, dist_list

    def WriteToFile(self,filename):
        """
        Writes the tree data into the specified binary file.
        """
        _lib.write_tree.argtype = [  C_Tree, C.c_char_p  ]
        _lib.write_tree.restype = C.c_void_p
        _lib.write_tree(self.c_tree,filename)

    def FreeMemory(self):
        """
        Free the memory used by the tree. 
        """
        _lib.free_tree.argtype = [  C_Tree  ]
        _lib.free_tree.restype = C.c_void_p
        _lib.free_tree(self.c_tree)


def draw_box(box):
    xlo,ylo,xhi,yhi=box.lo.x[0],box.lo.x[1],box.hi.x[0],box.hi.x[1],
    l = 1 # linewidth
    plt.plot([xlo,xlo],[ylo,yhi],'k',linewidth=l)
    plt.plot([xlo,xhi],[yhi,yhi],'k',linewidth=l)
    plt.plot([xhi,xhi],[yhi,ylo],'k',linewidth=l)
    plt.plot([xhi,xlo],[ylo,ylo],'k',linewidth=l)

def user_yes_no_query(question):
    sys.stdout.write('%s [y/n]\n' % question)
    while True:
        try:
            return strtobool(raw_input().lower())
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')



def CompanionSurfaceDensity(vectors,NNb):
    _lib.companion_surface_density.argtype = [ C.c_int, FLOAT_pointer, C.c_int, 
                                               C.c_int, FLOAT_pointer]
    _lib.companion_surface_density.restype = C.c_voidp
    DIM=len(vectors)
    N=len(vectors[0])
    coord=(DIM*N*FLOAT)() # This is a c_types pointer
    for i in range(N):
        for ndim,x in enumerate(vectors):
            coord[N*ndim+i]=x[i]
    density = (N*FLOAT)() # This is a c_types pointer
    _lib.companion_surface_density(NNb,coord,N,DIM,density)
    return np.asarray(density[:N],dtype=np.double)




if __name__ == "__main__":
    """
    Testing the various functions for a 2d tree.
    """

    N=1000000

    print "Generating coordinates"
    x = 2*np.random.random(N) -1
    y = 2*np.random.random(N) -1
    # z = 2*np.random.random(N) -1


    t0 = time.time()
    # Creating the initial tree
    tree0 = Tree(x,y)
    t1 = time.time()



    # Writing it to a file
    print "Writing to a file..."
    tree0.WriteToFile("test_tree.kd")
    t2 = time.time()

    # Reading it back into another tree
    print "Reading from the file..."
    tree = Tree.FromFile("test_tree.kd")
    t3 = time.time()

    # print t1-t0, "to compute the tree"
    # print t2-t1, "to write it to a file"
    # print t3-t2, "to read from the file"


    #  Chose what you want to test
    #testing = "draw"
    testing = "locatenear"
    #testing = "nnearest"
    #testing = "nearest"
    #testing = "none"


    # Testing the draw function
    if (testing=="draw"):
        tree.draw(markeredgecolor="b",markersize=4)
        plt.show()

    # Testing locatenear:
    if (testing=="locatenear"):
        print "Testing locatenear:"
        print "   We pick coordinates, then use the locatenear function"
        print "   to get all tree members within a chosen radius of these"
        print "   coordinates.\n.... "
        
        point=[0.2,0.2]     # Central point
        radius=0.01        # Maximum distance to the point
        ind , dist = tree.locatenear(point,radius)

        print "Plotting..."
        fig=plt.figure(1)
        ax=fig.add_subplot(111)
        x,y = tree.vectors()
        ax.plot(x,y,'bo')
        # Rescaling the plot to have a good view of the circle
        zoom=0.5
        plt.xlim((point[0]-(1/zoom)*radius, point[0]+(1/zoom)*radius))
        plt.ylim((point[1]-(1/zoom)*radius, point[1]+(1/zoom)*radius))
        # Plotting the circle around the point:
        th=np.linspace(0,2*np.pi,100)
        xcircle=np.array(point[0]+radius*np.cos(th))
        ycircle=np.array(point[1]+radius*np.sin(th))
        ax.plot(xcircle, ycircle)
        ax.set_aspect('equal')

        # Plotting central point as a green dot
        ax.plot([point[0],point[0]],[point[1],point[1]],'g.')
        # Getting the indices of all points lying in the circle:

        # Plotting these in a different colour:
        for i in ind:
            ax.plot([x[i],x[i]]
                    ,[y[i],y[i]],'ro')
        print "Done."
        plt.show()


    # Testing nnearest:
    if (testing=="nnearest"):
        print "Testing nnearest"
        fig=plt.figure(1)
        ax=fig.add_subplot(111)
        x,y = tree.vectors()
        ax.plot(x,y,'bo')

        index = random.randint(0,N-1)  # getting a particle at random
        Nnb = 50 # How many neighbours we want

        # Getting the indexes and distance to the Nnb closest particles
        ind, dist =tree.nnearest(index, Nnb)
        print "ind :", ind
        print "dist :", dist
        # Rescaling the plot to have a good view of the circle
        zoom=0.5
        maxdist = np.max(dist)
        plt.xlim((x[index]-(1/zoom)*maxdist, x[index]+(1/zoom)*maxdist))
        plt.ylim((y[index]-(1/zoom)*maxdist, y[index]+(1/zoom)*maxdist))
        # Plotting a circle enclosing all the neighbours
        th=np.linspace(0,2*np.pi,100)
        xcircle=np.array(x[index]+maxdist*np.cos(th))
        ycircle=np.array(y[index]+maxdist*np.sin(th))
        ax.plot(xcircle, ycircle)
        ax.set_aspect('equal')

        # Plotting all neighbours in red
        for i in ind:
            ax.plot([x[i],x[i]],[y[i],y[i]],'ro')
        # Plotting central particle in green
        ax.plot([x[index],x[index]],[y[index],y[index]],'go')
        
        plt.show()
    


    # Testing nearest:
    if (testing=="nearest"):
        print "Testing nearest"
        fig=plt.figure(1)
        ax=fig.add_subplot(111)
        x,y = tree.vectors()
        ax.plot(x,y,'bo')

        point=[0.2,0.2]
        ind = tree.nearest(point)
        dist= np.sqrt( (point[0]-x[ind])**2 + (point[1]-y[ind])**2  )
        zoom=0.001
        plt.xlim((point[0]-(1/zoom)*dist, point[0]+(1/zoom)*dist))
        plt.ylim((point[1]-(1/zoom)*dist, point[1]+(1/zoom)*dist))
        # Plotting central point as a green dot
        ax.plot([point[0],point[0]],[point[1],point[1]],'g.')
        # Plotting closest particle in red
        ax.plot([x[ind],x[ind]],[y[ind],y[ind]],'ro')
        plt.show()




