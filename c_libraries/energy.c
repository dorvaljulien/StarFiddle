#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>


double inv_distance(double x1, double x2, double y1, double y2, double z1, double z2)
{
    double distance2;
    distance2= (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    return 1/ sqrt(distance2);
}

void scale(int array_size, double factor, int count, ...)
/* 
   count is the number of arrays you want to scale with the supplied factor,
   for example to multiply arrays bananas, apple and oranges by 8.5:
   scale(N, 8.5 ,3 ,bananas, oranges, apple); 
*/
{
    va_list ap;
    int i,j;
    va_start (ap, count);                     /* Initialize the argument list. */
    for (i = 0; i < count; i++)
    {
        double *array = va_arg (ap, double*); /* Get the next argument value. */
	for( j=0; j<array_size; j++ )
	    array[j] = factor * array[j];
    }   
    va_end (ap);                              /* Clean up. */
}

void get_particles_E(int N, double *m, double *x, double* y, double *z, 
		     double *V2, double *Ek, double *Ep, double G){
    /*
      Take number of particles, mass, position and squared velocities array,
      return kinetic energy and potential energy array.
    */
    int i=0,j=0,k=0;
    double Ep_tmp, invd;
    for( i=0 ; i < N ; i++ ){
        Ek[i] = 0.5*m[i]*V2[i];
	Ep_tmp = 0.0;
        for( j=0 ; j < N ; j++ ){
	    k=0;
            if(i!=j) {
                invd=inv_distance(x[i],x[j],y[i],y[j],z[i],z[j]);
                Ep_tmp += - G*m[i]*m[j]*invd;
            }
	}
	Ep[i] =  Ep_tmp;
    }
}

int virialise(int N, double Q, int Standard,
	      double *m, double *rx, double *ry, double *rz, 
	      double *vx, double *vy, double *vz, double G ){
    /* 
          Take the number of particles N, the requested virial ratio Q, 
       masses, positions and velocities array, and scale these to
       fit Q.
          If Standard is set to 0, only velocities are scaled
       to reach Q, only modifying Ek. If set to 1, both 
       velocities and positions are modified to reach Q and 
       have a total energy of -0.25 (Hénon units) 
    */
    double tEk=0, tEp=0, target_tEk, target_tEp;
    double tol = 1e-8;
    int i,inc=0;
    double *V2 = malloc( N* sizeof(double) );
    double *Ek = malloc( N* sizeof(double) );
    double *Ep = malloc( N* sizeof(double) );
    for( i=0; i<N; i++ )
	V2[i] = vx[i]*vx[i] + vy[i]*vy[i]+ vz[i]*vz[i];
    get_particles_E(N,m,rx,ry,rz,V2,Ek,Ep,G);
    for( i=0; i<N; i++ ){
	tEk += Ek[i];
	tEp += 0.5*Ep[i];
    }
    if (Standard == 1){
	target_tEp = - 0.25 / ( 1 - Q );
	target_tEk = - Q * target_tEp;
	scale(N,tEp/target_tEp,3,rx,ry,rz);
    }else{
	target_tEk = - Q * tEp;
    }
    if (Q==0){
	scale(N, 0., 1, V2);
	scale(N, sqrt(target_tEk/tEk), 3, vx, vy, vz);
    }else{
	scale(N, sqrt(target_tEk/tEk), 3, vx, vy, vz);
    }
    free(V2); free(Ek); free(Ep);
    return 0;
}
 

int locate(int N, int *array, int target)
/*
  Find the indice of a particular element, target, inside array.
  If it's not found, returns -1.
*/
{
    int i;
    for(i=0;i<N;i++){
        if(array[i]==target) return i;
    }
    return -1;
}






