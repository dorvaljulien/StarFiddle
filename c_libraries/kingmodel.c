/* 
This code allows the user to build King model for nbody dynamics of star clusters. 
It is based on a Fortran code by Gerry Quinlan, modified by Christian Boily and 
converted to C by Julien Dorval.
*/


#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"kingmodel.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
double sqre(double a) {return a*a;}
double cub(double a) {return a*a*a;}
double SIGN(double a, double b)
    {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}


int main(){
    int i;
    
    double w0 = 2;
    double rmin = 1e-5;
    double dlogr = 0.039414;
    int Nmax = 300;
    
    double *radius = ALLOC(Nmax,double);
    double *mass = ALLOC(Nmax,double);
    double *potential = ALLOC(Nmax,double);
    
    int N = King_model(w0,rmin,dlogr,Nmax, radius, potential, mass);
    if (N == Nmax){
	fprintf(stderr,"Warning, didn't reach tidal radius.");
	fprintf(stderr," This means the returned model is\n ");
	fprintf(stderr,"truncated, you may want to try again");
	fprintf(stderr," with a higher Nmax or higher dlogr. \n");
    }

    FILE *f = fopen("c_function.out","w");
    for( i=0; i<N; i++ ){
	fprintf(f,"%le %le %le\n",radius[i], potential[i], mass[i]);
    }
    fclose(f);
    
    free(radius); free(potential); free(mass);
}


int King_model(double w0, double rmin, double dlogr, int Nmax,
	       double *radius, double *potential, double *mass)
/*
Purpose:
     This function will compute a king model of parameter w0 from rmin
     up to the tidal radius of the cluster, logarithmically stepping 
     by dlogr.
Arguments:
     double w0      - King parameter, the higher the deeper.
     double rmin    - starting radius, will be radius[0],
     double dlogr   - log radius bin
     int Nmax       - Maximum number of points we want
     double *radius - Malloc'd (Nmax) array
     double *potential - Malloc'd (Nmax) array
     double *mass   - Malloc'd (Nmax) array
Returns:
     nrg - actual number of data points stored in the arrays
*/
{
/* ERROR CONTROL */
    double eps = 1e-6;

/* MISC VARIABLES */
    double sigma2=1.0;
    double den1=1.0;
    int i,j,nok,nbad,nk,nr,nrg;
    int nvar = 2, Nrmax = Nmax;

    double h1=1e-3, hmin = 1e-8;
    double r1, r2, rt, ppot, a,b;
    double oldpot1, oldpot2;
    double p,tote,totm,h,r0,p0,ss,d0;

    double pot[2];
    double *den = ALLOC(Nrmax, double);
    double *dene = ALLOC(Nrmax, double);
    double *denm = ALLOC(Nrmax, double);
    double *dene2 = ALLOC(Nrmax, double);
    double *denm2 = ALLOC(Nrmax, double);
    double *der = ALLOC(Nrmax, double);
    double *v2 = ALLOC(Nrmax, double);
   
    pot[0] = w0 * sigma2;
    pot[1] = 0;
    r2 = rmin;

/*--------------------------------------------------------------------
       Integrate outward to tidal boundary. Store potential and
       density, and mass at every step.
--------------------------------------------------------------------*/
    for( nr=0; nr<Nrmax; nr++ ){
	r1 = r2;
	r2 = pow(10, log10(r1)+dlogr);
	oldpot1 = pot[0];
	oldpot2 = pot[1];
	odeint( pot, nvar, r1, r2, eps, h1, hmin, &nok, &nbad, 
	        sigma2 , den1, potdir, rkqs);
	//printf("nr %d pot %le %le \n", nr, pot[0], pot[1]);

	if(pot[0] > 0){
	    /*  Have not reached rt yet; store data and continue.*/
	    radius[nr] = r2;
	    potential[nr] = pot[0];
	    mass[nr] = -pot[1]*r2*r2;
	    p = pot[0]/sigma2;
	    den[nr] = den1 * ( exp(p) * erf(sqrt(p)) - 
	                      sqrt(4.*p/M_PI)*(1. +2.*p/3) );
	}else{
	    if (pot[0] == 0){ 
		break; 
	    }
	    /* Have passed rt. Use bisection to find it.*/
	    for( i=0; i<20; i++ ){
		pot[0] = oldpot1;
		pot[1] = oldpot2;
		rt = 0.5 * (r1+r2);
		odeint( pot, nvar, r1, rt, eps, h1, hmin, &nok, 
	               &nbad, sigma2 , den1, potdir, rkqs);
		if(pot[0] > 0){
		    r1 = rt;
		    oldpot1 = pot[0];
		    oldpot2 = pot[1];
		}else if(pot[0] < 0){ 
		    r2 = rt; 
		}else{ 
		    break; 
		}
	    }
	    radius[nr] = rt;
	    potential[nr] = 0;
	    mass[nr] = -pot[1]*r2*r2;
	    den[nr] = 0;
	    break;
	}
	if(den[nr] < 0) break;

    }
    nrg = nr;

/*-------------------------------------------------------------------
           Compute total mass and energy.
--------------------------------------------------------------------*/
    for( nr=0; nr<nrg; nr++ ){
	dene[nr] = radius[nr]*radius[nr] * den[nr] * potential[nr];
	denm[nr] = radius[nr]*radius[nr] * den[nr];
    }
    spline(radius, dene, nrg, 0, 2e30, dene2);
    spline(radius, denm, nrg, 0, 2e30, denm2);
    tote = 0;
    totm = 0;

    for( nr=0; nr<nrg-1; nr++ ){
	h = radius[nr+1]-radius[nr];
	tote = tote + 0.5*h*( dene[nr+1] + dene[nr] ) -
	    cub(h)/24.*( dene2[nr+1] + dene2[nr] );
	totm = totm + 0.5*h*( denm[nr+1] + denm[nr] ) -
	    cub(h)/24.*( denm2[nr+1] + denm2[nr] );
    }

    totm = 4.*M_PI*totm;
    tote = 2.*M_PI*tote + 0.5*totm*totm/radius[nrg-1];
    tote = 0.5*tote;

    r0 = 4.*tote/sqre(totm);
    p0 = totm/(4.*tote);
    d0 = pow(totm,5)/cub(4.*tote);

    for( nr=0; nr<nrg; nr++ ){
        radius[nr] = r0*radius[nr];
        potential[nr] = p0*potential[nr];
        mass[nr] = mass[nr]/totm;
	den[nr] = d0*den[nr];
    }

/*--------------------------------------------------------------------
       Get the mean square velocity v2 as a function of radius.
--------------------------------------------------------------------*/
    /* for( nr=0; nr<nrg; nr++ ){ */
    /* 	nk = 2; */
    /* 	ppot = potential[nr]; */
    /* 	a = 0; b = sqrt(2*ppot); */
    /* 	ss = qromo(a, b, eps, kgvfnc, midpoint, nk, ppot, sigma2); */
    /* 	v2[nr] = ss; */
    /* 	nk = 0; */
    /* 	ss = qromo(a, b, eps, kgvfnc, midpoint, nk, ppot, sigma2); */
    /* 	v2[nr] = v2[nr]/ss; */
    /* } */

    free(v2);
    free(den); free(dene); free(denm);
    free(dene2); free(denm2); free(der);

    return nrg;
}








void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
/* Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e.,
 yi = f (xi ), with x1 < x2 < . . . < xN , and given values yp1 and ypn 
for the first derivative of the interpolating function at points 1 and n, 
respectively, this routine returns an array y2[1..n] that contains the 
second derivatives of the interpolating function at the tabulated points 
xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is 
signaled to set the corresponding boundary condition for a natural spline, 
with zero second derivative on that boundary. */
{

    int i,k;
    double p,qn,sig,un;
    double *u = ALLOC(n, double);

    if (yp1 > 0.99e30)	y2[0]=u[0]=0.0;
    else {
	    y2[0] = -0.5;
	    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
    for (i=1; i<n-1; i++) {
	sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
	p=sig*y2[i-1]+2.0;
	y2[i]=(sig-1.0)/p;
	u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
	u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
	qn=un=0.0;
    else {
	qn=0.5;
	un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k=n-2;k>=0;k--)
	y2[k]=y2[k]*y2[k+1]+u[k];
    free(u);
}



void potdir(double r, double *pot, double *dpotdr, 
	    double sigma2, double den1){
    dpotdr[0] = pot[1];
    double den, p = pot[0]/sigma2;
    if(p>0)
	den = den1 * ( exp(p)*erf(sqrt(p)) - 
		       sqrt(4.*p/M_PI)*(1. +2.*p/3)    );
    else
	den = 0;
    dpotdr[1] = -(2./r)*dpotdr[0] - 4.*M_PI*den;
}




#define MAXSTP 10000
#define TINY 1.0e-30
void odeint(double *ystart, int nvar, double x1, double x2, 
	    double eps, double h1, double hmin, int *nok, int *nbad, 
	    double sigma2, double den1,
	    void (*derivs)(double, double* , double*, double, double ),
	    void (*rkqs)(double*, double*, int, double*, double, 
			 double, double*,  double*, double*, double, double,
			 void (*)(double, double*, double*, double, double)
	          	)
            )
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	// kmax is for storage of intermediate data (for inspection)
	// These parts were deleted.
	int kmax=0,kount=0, Nmax =2;
	double *xp,**yp,dxsav;
	double *y = ALLOC(Nmax, double);
	double *dydx = ALLOC(Nmax, double);
	double *yscal = ALLOC(Nmax, double);
	int sw=0;
	

	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=0;i<nvar;i++) y[i]=ystart[i];
	for (nstp=0;nstp<MAXSTP;nstp++) {
	    (*derivs)(x,y,dydx,sigma2,den1);

	    for (i=0;i<nvar;i++)
		yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,sigma2,den1,derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=0;i<nvar;i++) ystart[i]=y[i];
			free(y); free(dydx); free(yscal);
			return;
		}
		if (fabs(hnext) <= hmin) {
		    printf("Step size too small in odeint");
		    exit(EXIT_FAILURE);
		}
		h=hnext;
	}
	printf("Too many steps in routine odeint");
	exit(EXIT_FAILURE);
}
#undef MAXSTP
#undef TINY

 
#define SAFETY 0.9   
#define PGROW -0.2   
#define PSHRNK -0.25   
#define ERRCON 1.89e-4 
void rkqs(double *y, double *dydx, int n, double *x, 
	      double htry, double eps, double *yscal, 
	      double *hdid, double *hnext, double sigma2, double den1,   
	      void (*derivs)(double, double *, double *, double, double))   
{   
    int i;   
    double errmax,h,htemp,xnew;   
    double *yerr=ALLOC(n,double);   
    double *ytemp=ALLOC(n,double);   
    
    h=htry;   
    for (;;) {
	rkck(y,dydx,n,*x,h,ytemp,yerr,sigma2,den1,derivs);   
        errmax=0.0;   
        for (i=0;i<n;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));   
        errmax /= eps;   
        if (errmax <= 1.0) break;   
        htemp=SAFETY*h*pow(errmax,PSHRNK);   
        h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));   
        xnew=(*x)+h;   
    }   
    if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);   
    else *hnext=5.0*h;   
    *x += (*hdid=h);   
    for (i=0;i<n;i++) y[i]=ytemp[i]; 
    free(ytemp); free(yerr);
}   
#undef SAFETY   
#undef PGROW   
#undef PSHRNK   
#undef ERRCON   
  


  
   
void rkck(double *y, double *dydx, int n, double x, 
	  double h, double *yout, double *yerr, 
	  double sigma2, double den1,
	  void (*derivs)(double, double*, double*, double, double))
{   
    int i;   
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,   
        b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,   
        b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,   
        b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,   
        b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,   
        c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,   
        dc5 = -277.00/14336.0;   
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,   
        dc4=c4-13525.0/55296.0,dc6=c6-0.25;   
    double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;   
    int sw =0;
    if (y[1] == 0) sw =1;
    
    ak2=ALLOC(n,double);   
    ak3=ALLOC(n,double);   
    ak4=ALLOC(n,double);   
    ak5=ALLOC(n,double);   
    ak6=ALLOC(n,double);   
    ytemp = ALLOC(n,double);
   
    for (i=0;i<n;i++)   
        ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2, sigma2, den1);
    for (i=0;i<n;i++){
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]); 
    }
    (*derivs)(x+a3*h,ytemp,ak3, sigma2, den1);   
    for (i=0;i<n;i++)   
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);   
    (*derivs)(x+a4*h,ytemp,ak4, sigma2, den1);   
    for (i=0;i<n;i++)   
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);   
    (*derivs)(x+a5*h,ytemp,ak5, sigma2, den1);   
    for (i=0;i<n;i++)   
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);   
    (*derivs)(x+a6*h,ytemp,ak6, sigma2, den1);   
    for (i=0;i<n;i++){
	yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);   
    }
    for (i=0;i<n;i++)   
        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);   
    free(ytemp);   
    free(ak6);   
    free(ak5);   
    free(ak4);   
    free(ak3);   
    free(ak2);   
}   
   




void polint(double *xa, double *ya, int n, double x, 
	    double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;
	dif=fabs(x-xa[1]);
	c=ALLOC(n,double);
	d=ALLOC(n,double);

	for (i=0;i<n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
			    printf("Error in routine polint\n");
			    exit(EXIT_FAILURE);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free(d); free(c);
}



double kgvfnc(double v, int nk, double ppot, double sigma2){
    return pow(v,(nk+2)) * ( exp( (ppot-0.5d*v*v)/sigma2 ) - 1. );
}



double qromo( double a, double b, double eps, 
	    double (*kgvfnc)(double,int,double,double),
	    double (*choose)(double (*)(double, int, double, double), 
			     double, double, int, int, double, double),
	    int nk, double ppot, double sigma2)
{
    int j,Jmax = 14, Jmaxp = Jmax +1, K=5, Km=K-1;
    double ss,dss;
    double *h = ALLOC(Jmaxp, double);
    double *s = ALLOC(Jmaxp, double);

    h[0] = 1;
    for( j=0; j<Jmax; j++ ){
	s[j]=(*choose)(kgvfnc,a,b,j, nk, ppot, sigma2);
	if (j >= K) {
	    polint( &h[j-K], &s[j-K], K, 0., &ss, &dss);
	    if (fabs(dss) <= eps*fabs(ss)){
		free(h); free(s);
		return ss;
	    }
	}
	h[j+1]=h[j]/9.0;
    }
    printf("Too many steps in routing qromo\n");
    exit(EXIT_FAILURE);
}



double midpoint( double (*func)(double,int,double,double), 
		 double a, double b, int n,
		 int nk, double ppot, double sigma2){
    int it,i;
    double tnm, del, ddel, sum, x;
    static double s;
    if(n==1) return (s=(b-a) * func(0.5*(a+b), nk, ppot, sigma2));
    else{
	it = pow(3,n-2);
	tnm = it;
	del = (b-a) / (3.*tnm);
	ddel = del + del;
	x = a + 0.5 * del;
	sum = 0;
	for( i=0; i<it; i++ ){
	    sum += func(x, nk, ppot, sigma2);
	    x += ddel;
	    sum += func(x, nk, ppot, sigma2);
	    x += del;
	}
	s=(s+(b-a)*sum/tnm)/3.;
	return s;
    }
}
    






