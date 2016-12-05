#include<stdlib.h>
#include<stdio.h>
#include<math.h>



#define PrI(_i) {printf("%d\n",_i);}

// Defining clean memory allocation.
#define ALLOC(_n, _type)                                                \
    ({ void *x = malloc(_n * sizeof(_type));                            \
	if(x==NULL){							\
	    fprintf(stderr,                                             \
                   "Error - %s (line %d) : Memory allocation failed\n",	\
		    __FILE__, __LINE__);                                \
	    exit(EXIT_FAILURE);						\
	}								\
	x;})

int King_model(double w0, double rmin, double dlogr, int Nmax,
	       double *radius, double *potential, double *mass);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void potdir(double r, double *pot, double *dpotdr, double sigma2, double den1);
void odeint(double *ystart, int nvar, double x1, double x2, double eps, double h1,
	    double hmin, int *nok, int *nbad, double sigma2, double den1,
	    void (*derivs)(double, double* , double*, double, double ),
	    void (*rkqs)(double*, double*, int, double*, double, 
			 double, double*,  double*, double*, double, double,
			 void (*)(double, double*, double*, double, double)
	          	)
           );
void rkqs(double *y, double *dydx, int n, double *x, 
	  double htry, double eps, double *yscal, 
	  double *hdid, double *hnext, double sigma2, double den1,   
	  void (*derivs)(double, double *, double *, double, double)); 
void rkck(double *y, double *dydx, int n, double x, 
	  double h, double *yout, double *yerr, double sigma2, double den1,
	  void (*derivs)(double, double*, double*, double, double));
void polint(double *xa, double *ya, int n, double x, double *y, double *dy);
double kgvfnc(double v, int nk, double ppot, double sigma2);
double qromo( double a, double b, double eps, 
	      double (*kgvfnc)(double,int,double,double),
	      double (*choose)(double (*)(double, int, double, double), 
			       double, double, int, int, double, double),
	      int nk, double ppot, double sigma2); 
double midpoint( double (*func)(double,int,double,double), 
		 double a, double b, int n,
		 int nk, double ppot, double sigma2);

