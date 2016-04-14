"""
This module deals with binary stars. The two main functions are:

- a , e = get_parameters(r,v,mt,G):
    This retrieves orbital parameters of a binary system when 
    only relative position, relative velocity and masses are known.

- r1,r2,v1,v2 = RandomBinary( a,e,m1,m2,G )
    This retrieves position and velocity as seen from center of mass 
    given orbital parameters. The phase and orientation of the binary
     are randomly chosen.

In many functions, the default gravitationnal constant is 1. Positions
masses and velocities should be in nbody units for this to make sense.
"""

from numpy import cos,sin,sqrt
from numpy.linalg import norm
import numpy as np
Pi=np.pi

#********************************************************************
#***  From position and velocity, get the orbital parameters  *******
#********************************************************************

def scalar(a,b):
    """Returns the scalar product of 2 3d vectors (list/array)"""
    s = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    return s

def vectorial(a,b):
    """Returns the scalar product of 2 3d vectors (list/array)"""
    v1= a[1]*b[2]-a[2]*b[1]
    v2= a[2]*b[0]-a[0]*b[2]
    v3= a[0]*b[1]-a[1]*b[0]
    return np.array([v1,v2,v3])

def get_a(r,v,mt,G=1):
    """
    This function takes the relative radius and velocity vector of 
    a system, the total mass, and returns the semi-major axis.
    """
    r=norm(r)
    v=norm(v)
    a= r / ( 2 - r*v**2 / (G*mt) )
    return a 
    
def get_e(r,v,a,mt,G=1):
    """
    This function takes the relative radius, velocity, as well as the 
    previously obtained semi-major axis of the system and the total 
    mass to return the computed eccentricity.
    """
    h=vectorial(r,v)
    h=norm(h)
    e2= 1 - h**2 / ( a * G * mt )
    e=sqrt(e2)
    return e

def get_parameters(r,v,mt,G=1):
    """
    This functions takes relative distance (3d vector) relative 
    velocity (3d vector) and total mass to return the semi major axis 
    and eccentricity as a,e.
    """
    V2 = norm(v)**2
    R = norm(r)
    ebin = V2/2. - G * mt / R
    if ebin > 0:
        raise Exception("get_parameters was fed with data from" 
                        "two unbound stars. Cannot compute a and e.")
    a=get_a(r,v,mt,G)
    e=get_e(r,v,a,mt,G)
    return a,e


#*********************************************************************
#*****  From the orbital parameters, get position and velocity  ******
#*********************************************************************

def rotate(vec,phi1,phi2,phi3):
    """
    Rotate the given 3d vector along the three axis with given angles. 
    Useful for rotated_binary_sample
    """
    matrot1=np.matrix([ [1,0,0] , 
                        [0, cos(phi1) , -sin(phi1)] , 
                        [0, sin(phi1) , cos(phi1)] , ])
    matrot2=np.matrix([ [cos(phi2), 0, sin(phi2)] , 
                        [0,1,0] , 
                        [-sin(phi2), 0, cos(phi2)]  ])
    matrot3=np.matrix([ [cos(phi3) , -sin(phi3) , 0] , 
                        [sin(phi3) , cos(phi3) , 0] , 
                        [0,0,1] ])
    vec=np.matrix(vec)
    vec=vec*matrot1*matrot2*matrot3
    return np.array(vec)

def ang_momentum(a,e,m1,m2,G=1):
    """returns the angular momentum  in the appropriate units."""
    return sqrt(a*(1-e**2)*G*(m1+m2))

def r_orbit(a,e,m1,m2,theta=0,G=1):
    """
    For two orbiting particles and a given theta, 
    returns relative radius. 
    """
    r = a* ( 1 - e**2 ) / ( 1 - e*cos(theta))
    return r

def r_dot(a,e,m1,m2,theta=0,G=1):
    """
    For two orbiting particles and a given theta, 
    returns relative radius derivative. 
    """
    h=ang_momentum(a,e,m1,m2,G)
    rd= - h  *  e*sin(theta) / ( a * (1-e**2) )
    return rd

def theta_dot(a,e,m1,m2,theta=0,G=1):
    """
    For two orbiting particles and a given theta, 
    returns theta derivative. 
    """
    r=r_orbit(a,e,m1,m2,theta,G)
    h=ang_momentum(a,e,m1,m2,G)
    return h /r**2

def elliptical_orbit(e,a,n):
    """returns a [0:2Pi] array and the coordinates of a 2d ellipsis"""
    th=np.arange(float(n))/n*2*3.1415
    r=a*(1-e**2)/(1+e*np.cos(th))
    return th,r

def velocity_binary(a,e,m1,m2,theta=0,G=1):
    """
    Returns x1,y1,x2,y2,vx1,vy1,vx2,vy2 
    for two gravitating particles (planar motion)
    """
    if e<0 or e>1:
        raise Exception("Error, bad eccentricity fed to the function.")
    r=r_orbit(a,e,m1,m2,theta,G)
    rd=r_dot(a,e,m1,m2,theta,G)
    thd=theta_dot(a,e,m1,m2,theta,G)
    mu1,mu2= m2/(m1+m2) , m1/(m1+m2)

    x1=-mu1*r*cos(theta)  
    y1=-mu1*r*sin(theta)
    x2=mu2*r*cos(theta)
    y2=mu2*r*sin(theta)
    vx1= - mu1 * (rd*cos(theta)-r*thd*sin(theta))
    vy1= - mu1 * (rd*sin(theta)+r*thd*cos(theta))
    vx2=  mu2 * (rd*cos(theta)-r*thd*sin(theta))
    vy2= mu2 * (rd*sin(theta)+r*thd*cos(theta))

    return x1,y1,x2,y2,vx1,vy1,vx2,vy2


def rotated_binary_sample(a,e,m1,m2,theta,phi1,phi2,phi3,G=1):
    """
    Returns r1,r2,v1,v2: numpy arrays of positions and velocities of a 
    binary rotated along three axes with phi1,phi2,phi3, 
    at an internal angle theta
    """
    x1,y1,x2,y2,vx1,vy1,vx2,vy2 = velocity_binary( a,e,m1,m2,theta,G)
    r1=rotate([x1,y1,0],phi1,phi2,phi3)[0]
    r2=rotate([x2,y2,0],phi1,phi2,phi3)[0]
    v1=rotate([vx1,vy1,0],phi1,phi2,phi3)[0]
    v2=rotate([vx2,vy2,0],phi1,phi2,phi3)[0]
    return r1,r2,v1,v2

def RandomBinary(a,e,m1,m2,G=1):
    """
    Return position and velocity of generated binary system.

    r1,r2,v1,v2 = RandomBinary(a,e,m1,m2,G=1)
    
    The binary system is created, rotated with random angles along 3 
    axis, then its internal angle is also set to a random angle.
    """
    theta,phi1,phi2,phi3=[np.random.random()*2*Pi for i in [1,2,3,4]]
    return rotated_binary_sample(a,e,m1,m2,theta,phi1,phi2,phi3,G=G)

  
