double inv_distance(double x1, double x2, double y1, double y2, double z1, double z2);
void scale(int array_size, double factor, int count, ...);
void get_particles_E(int N, double *m, double *x, double* y, double *z, 
		     double *V2, double *Ek, double *Ep, double G);
int virialise(int N, double Q,
	      double *m, double *rx, double *ry, double *rz, 
	      double *vx, double *vy, double *vz, double G );
int locate(int N, int *array, int target);
