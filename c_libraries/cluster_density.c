#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <kdtree.h>
#include <argsort.h>
 
void salpeter(int N, double *m, double alpha, double m1, double m2){
    /*
      Generate N masses according to a salpeter mass function of slope
      alpha and a mass range of [m1, m2].
    */
    int i;
    srand(time(NULL));
    double lm,prob, mass;
    double lmin = log10(m1);
    double lmax = log10(m2);
    for( i=0; i<N; i++ ){
	mass = m2;
	prob = 1.;
	while(prob > pow(mass,1-alpha)){
	    mass = pow(10, lmin + ( (double)rand()/(double)RAND_MAX )*(lmax-lmin) );
	    prob = 1.05 * ((double)rand()/(double)RAND_MAX) * pow(m1,1-alpha);
	}
	m[i] = mass;
    }
}



void generate_uniform(int N, double *x, double*y, double*z, double *m, double R){
    /*
      Takes allocated (size N) arrays to build an uniform model of bounding radius R.
      Uses a standard salpeter mass function between 0.2 and 20 solar masses.
     */
    salpeter(N, m, 2.37, 0.2, 20.);
    int i;
    double r, p, x1, x2, phi, theta;
    srand(time(NULL));
    for( i=0; i<N; i++ ){
	r = 0;
	p = 1;
	while(p>r*r){
	    r = (double)rand()/(double)RAND_MAX ;
	    p = (double)rand()/(double)RAND_MAX ;
	}
	x1 = (double)rand()/(double)RAND_MAX ;
	x2 = (double)rand()/(double)RAND_MAX ; 
	theta = acos(1-2*x1);
	phi = 2*M_PI*x2;
	x[i] = r*R*sin(theta)*cos(phi);
        y[i] = r*R*sin(theta)*sin(phi);
        z[i] = r*R*cos(theta);
    }
}



int density_distribution(int N, int Nnb, int Number, double *m, 
			 double *x, double *y, double *z,
			  double *density){
    /*
        Takes masses (m) and coordinates (x,y,z) of a system of N
      particles, then builds a kdtree and finds local density
      for each point by evaluating Nnb neighbours around it.
        If Number is set to 1, the returned density is a number 
      density. If set to 0, it is a normal mass density.
    */
    int i,k;
    double m_local, loc_r, V_local, rho_local;
    int *n_neighbours=malloc(Nnb*sizeof(int));
    double *d_neighbours=malloc(Nnb*sizeof(double));
    /*To build the KDtree, we need a 1d array of coordinates*/

    double *coords=malloc(3*N*sizeof(double));
    for(i=0;i<N;i++){
        coords[i]=x[i];
        coords[N+i]=y[i];
        coords[2*N+i]=z[i]; 
    }
    Tree tree;
    KDtree(coords,N,3,&tree);

    if(Number==0){
	for( i=0; i<N; i++ ){
	    m_local = 0;
	    (void)nnearest(tree,i,n_neighbours, d_neighbours, Nnb);
	    for(k=0;k<Nnb-1;k++) m_local += m[n_neighbours[k]]; 
	    loc_r = d_neighbours[Nnb-1];
	    V_local = (4./3)*M_PI*loc_r*loc_r*loc_r;
	    density[i] = m_local/V_local;
	}
    }
    else{
	for( i=0; i<N; i++ ){
	    (void)nnearest(tree,i,n_neighbours, d_neighbours, Nnb);
	    loc_r = d_neighbours[Nnb-1];
	    V_local = (4./3)*M_PI*loc_r*loc_r*loc_r;
	    density[i] = Nnb/V_local;
	}
    }
    return N;
}




