/* 
This is a C implementation of a KDtree, used in neighbour searches.
It is translated from c++, the original code can be found in the
Numerical Recipees, third edition.
It was converted to c by Julien Dorval - dorvaljulien@gmail.com

An example of usage can be found in the main at the end.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <time.h>  
#include "argsort.h" 
#define BIG 1e38
#define FLOAT double


// Defining clean memory allocation.
#define ALLOC(_n, _type, _label)				            \
    ({void *x = malloc(_n * sizeof(_type));				    \
	if(x==NULL){						            \
	    fprintf(stderr,"%s (line %d) : Error, Could not allocate %s\n", \
		    __FILE__, __LINE__, _label);			    \
	    exit(EXIT_FAILURE);						    \
	}								    \
	x;})


//FLOAT inline square(FLOAT a) {return a*a;}
FLOAT square(FLOAT a) {return a*a;}

/* Define the Point, Box and Tree structures */

typedef struct Point 
{
    int DIM;
    FLOAT *x;
}Point;

typedef struct Box
// Base element of the tree
{
    Point hi,lo;
    int mom, dau1, dau2, ptlo, pthi;
}Box;

typedef struct Tree
{
    int N;               // Number of points in the tree
    int DIM;             // Dimension of the system
    Point *pts;          // Array of Points with coordinates
    int *ptindx;         // Keep track of the position of each point
    int *rptindx;        // reverse index: rptindx[ptindx[i]]=i
    int nboxes;          // Total number of boxes in the tree
    Box *boxes;          // Array of boxes
    FLOAT *global_coords; // Array containing all coordinates for all points
}Tree;



//--------------------------------------------------------------------
//-----   Structure management: creation, initialization       -------
//--------------------------------------------------------------------


void modify_Point(Point *pt, 
		  FLOAT *coords){
    int i;
    for( i=0; i<(pt->DIM); i++ ) pt->x[i] = coords[i];
}


void init_and_fill_Points(FLOAT *coords, 
			  int N, 
			  int DIM, 
			  Point *pts, 
			  FLOAT *global_coords, 
			  long *ind)
/*
   This takes a previously malloced array of Point (pts)
   and initializes it:
   - setting the DIM integer
   - Attributing to .x a location in the global_coords
   array that holds all coordinates for all points. 
*/

{
    int i,j;
    FLOAT *tmp_coords=ALLOC(DIM,FLOAT,"tmp_coords");
    for(i=0;i<N;i++){
	pts[i].DIM = DIM;
	for (j=0;j<DIM;j++) tmp_coords[j]=coords[N*j+i];   
	pts[i].x = &global_coords[*ind]; *ind+=DIM;
	for( j=0; j<DIM; j++ ){
	    pts[i].x[j] = tmp_coords[j];
	}
	//modify_Point(&pts[i],tmp_coords);
    }
    free(tmp_coords);
}




void modify_Box(Point *lo, Point *hi, 
		int mom, int dau1, int dau2, 
		int ptlo, int pthi, 
		Box *b)
// Fills a previously empty box.
{
    int DIM=lo->DIM;
    modify_Point(&(b->lo),lo->x);
    modify_Point(&(b->hi),hi->x);
    b->mom=mom; b->dau1=dau1; b->dau2=dau2; b->ptlo=ptlo; b->pthi=pthi;
}


//--------------------------------------------------------------------
//------------------   Distances functions       ---------------------
//-------------------------------------------------------------------



FLOAT dist(Point a, Point b)
// Gives distance between 2 points
{
    int i, DIM=a.DIM;
    FLOAT dist=0;
    for(i=0;i<DIM;i++)  dist+=square(a.x[i]-b.x[i]);
    return sqrt(dist);
}

FLOAT distBox(Box b, Point p)
// Gives distance of a point from a box
{
    int DIM=p.DIM;
    FLOAT dd = 0;
    int i;
    for (i=0; i<DIM; i++){
        if (p.x[i]<b.lo.x[i]) dd += square(p.x[i]-b.lo.x[i]);
        if (p.x[i]>b.hi.x[i]) dd += square(p.x[i]-b.hi.x[i]);
    }
    return sqrt(dd);
}

FLOAT disti(Tree tree, int j, int k) 
//Return distance between point j and k in the tree
{
    if (j == k) return BIG;
    else return dist(tree.pts[j], tree.pts[k]);
}



//--------------------------------------------------------------------
//-------------   Array manipulation: swap, sort      ----------------
//--------------------------------------------------------------------


void swap(int *array, int i, int j)
// swap elements i and j in array.
{
    int tmp;
    tmp=array[i];
    array[i]=array[j];
    array[j]=tmp;
}

int selecti(const int k, int *indx, int n, FLOAT *arr)
/* All the work of the tree is done here: the array arr is provided, then
   partially sorted through an index array indx (arr itself is untouched).
   At the end of the routine:  
   arr[indx[0..k-1]] < arr[indx[k]] < arr[indx[k+1..n-1]]  */
{
    int i,ia,ir,j,l,mid;
    FLOAT a;

    l=0;
    ir=n-1;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[indx[ir]] < arr[indx[l]])
                swap(indx,l,ir);
            return indx[k];
        } else {
            mid=(l+ir) >> 1;
            swap(indx,mid,l+1);
            if (arr[indx[l]] > arr[indx[ir]]) swap(indx,l,ir);
            if (arr[indx[l+1]] > arr[indx[ir]]) swap(indx,l+1,ir);
            if (arr[indx[l]] > arr[indx[l+1]]) swap(indx,l,l+1);
            i=l+1;
            j=ir;
            ia = indx[l+1];
            a=arr[ia];
            for (;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if (j < i) break;
                swap(indx,i,j);
            }
            indx[l+1]=indx[j];
            indx[j]=ia;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
        }
    }
}

void sift_down(FLOAT *heap, int *ndx, int nn) 
// This is a sorting routine, related to the heap-sort algorithm.
// It is used by the nnearest routine
{
    int n = nn - 1;
    int j,jold,ia;
    FLOAT a;
    a = heap[0];
    ia = ndx[0];
    jold = 0;
    j = 1;
    while (j <= n) {
        if (j < n && heap[j] < heap[j+1]) j++;
        if (a >= heap[j]) break;
        heap[jold] = heap[j];
        ndx[jold] = ndx[j];
        jold = j;
        j = 2*j + 1;
    }
    heap[jold] = a;
    ndx[jold] = ia;
}



//---------------------------------------------------------------------
//------         MAIN FUNCTION: Kdtree building          -------------
//---------------------------------------------------------------------
 

void KDtree(FLOAT *coords, int N, int DIM, Tree *tree)
/*
 Build a KD tree from a 1D FLOAT array containing all the coordinates: 
 [x[...],y[...],..etc]
 The dimension and total number of points have to be provided. The tree
 pointer should be malloced first.
*/
{
    if (coords == NULL) printf("Provided coords pointer non allocated !\n");
    if (tree == NULL) printf("Provided tree pointer non allocated !\n");

    // Defining variables.
    int jtmp,i,k,kk,j,nowtask,jbox,np,tmom,tdim,ptlo,pthi;
    int nboxes, m, ntmp;
    int *hp;
    FLOAT *cp;
    int taskmom[50], taskdim[50];
    long Ncoords;


    int *ptindx=ALLOC(N, int, "ptindx");
    int *rptindx=ALLOC(N, int, "rptindx");
    // Computing number of boxes from N
    for (k=0; k<N; k++) ptindx[k] = k; m = 1;
    for (ntmp = N; ntmp; ntmp >>= 1)  m <<= 1;
    nboxes = 2*N - (m >> 1);
    if (m < nboxes) nboxes = m;
    nboxes--;

    // Any Point in the tree will have its .x actually stored in here.
    // This avoid memory overhead when mallocing the points one by one.
    long ind = 0;
    Ncoords = DIM * ( 2 + (long)N + 2*(long)nboxes);
    FLOAT *global_coords = ALLOC(Ncoords,FLOAT, "global_coords");
    Point *pts = ALLOC( N, Point, "pts");
    init_and_fill_Points(coords, N, DIM, pts, &global_coords[ind], &ind);

    Box *boxes = ALLOC( nboxes, Box, "boxes");
    for( i=0; i<nboxes; i++ ){
	boxes[i].lo.DIM = DIM;
    	boxes[i].hi.DIM = DIM;
    	boxes[i].hi.x = &global_coords[ind]; ind += DIM;
    	boxes[i].lo.x = &global_coords[ind]; ind += DIM;
	boxes[i].dau1 = 0;
	boxes[i].dau2 = 0;
    }

    // Preparing the lo and hi reference Points
    FLOAT *xlo = ALLOC( DIM, FLOAT, "xlo");
    FLOAT *xhi = ALLOC( DIM, FLOAT, "xhi");
    for(i=0;i<DIM;i++) { xlo[i]=-BIG; xhi[i]=BIG; }
    Point *lo = ALLOC(1, Point, "lo");
    Point *hi = ALLOC(1, Point, "hi");
    lo->x = &global_coords[ind]; ind+=DIM;
    hi->x = &global_coords[ind]; ind+=DIM;
    lo->DIM = hi ->DIM = DIM;
    modify_Point(lo,xlo); modify_Point(hi,xhi);
    // First box:
    modify_Box(lo,hi,0,0,0,0,N-1,&boxes[0]);
    jbox = 0;
    taskmom[1] = 0;
    taskdim[1] = 0;
    nowtask = 1;
    
    // Getting into main loop 
    int step = nboxes / 20;

    while (nowtask) {
        tmom = taskmom[nowtask];
        tdim = taskdim[nowtask--];
        ptlo = boxes[tmom].ptlo;
        pthi = boxes[tmom].pthi;
        hp = &ptindx[ptlo];
        cp = &coords[tdim*N];
        np = pthi - ptlo + 1;
        kk = (np-1)/2;
        (void) selecti(kk,hp,np,cp);
        modify_Point(hi,boxes[tmom].hi.x);
        modify_Point(lo,boxes[tmom].lo.x);
        hi->x[tdim] =  coords[tdim*N + hp[kk]];
        lo->x[tdim] =  coords[tdim*N + hp[kk]];
        jbox++;
        modify_Box(&(boxes[tmom].lo), hi, tmom, 0, 0, ptlo, ptlo+kk, &boxes[jbox]);
        jbox++;
        modify_Box(lo, &(boxes[tmom].hi), tmom, 0, 0, ptlo+kk+1,  pthi, &boxes[jbox]);
        boxes[tmom].dau1 = jbox-1;
        boxes[tmom].dau2 = jbox;
        if (kk > 1) {
            taskmom[++nowtask] = jbox-1;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
        if (np - kk > 3) {
            taskmom[++nowtask] = jbox;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
    }

    for (j=0; j<N; j++) rptindx[ptindx[j]] = j;
    tree->N = N;
    tree->DIM = DIM;
    tree->pts = pts;
    tree->ptindx = ptindx;
    tree->rptindx = rptindx;
    tree->nboxes = nboxes;
    tree->boxes = boxes;
    tree->global_coords = global_coords;

    free(lo); free(xlo); free(hi); free(xhi);
}



//---------------------------------------------------------------------
//------------        reading and writing         --------------------
//---------------------------------------------------------------------


void write_tree(Tree tree, char *name)
/*
  Allows to write the tree as a binary file named "name".
*/
{
    int i,j;
    FILE *f=fopen(name,"wb");
    
    fwrite(&tree.N,sizeof(int),1,f);
    fwrite(&tree.DIM,sizeof(int),1,f);
    fwrite(&tree.nboxes,sizeof(int),1,f);

    /* Write the coordinates first*/
    for( i=0; i<tree.N; i++ ){
	fwrite(tree.pts[i].x,tree.DIM*sizeof(FLOAT),1,f);
    }
    
    /* Write the boxes */
    for( i=0; i<tree.nboxes; i++ ){
	fwrite(tree.boxes[i].hi.x,tree.DIM*sizeof(FLOAT),1,f);
	fwrite(tree.boxes[i].lo.x,tree.DIM*sizeof(FLOAT),1,f);
	fwrite(&tree.boxes[i].mom,sizeof(int),1,f);
	fwrite(&tree.boxes[i].dau1,sizeof(int),1,f);
	fwrite(&tree.boxes[i].dau2,sizeof(int),1,f);
	fwrite(&tree.boxes[i].ptlo,sizeof(int),1,f);
	fwrite(&tree.boxes[i].pthi,sizeof(int),1,f);
    }

    fwrite(tree.ptindx,tree.N*sizeof(int),1,f);
    
    fclose(f);
}

void read_tree(char *name, Tree *tree)
/*
  Allows to read a tree which was previously written as a binary file through write_tree.
*/
{
    int i,ii,j;
    FILE *f=fopen(name,"rb");
    if(f==NULL){
	printf("couldn't open %s\n",name);
	return;
    }
    int N,DIM,nboxes;

    fread(&N,sizeof(int),1,f);
    fread(&DIM,sizeof(int),1,f);
    fread(&nboxes,sizeof(int),1,f);

    long ind = 0;
    long Ncoords = DIM * (2 + (long)N + 2*(long)nboxes );
    FLOAT *global_coords = ALLOC(Ncoords, FLOAT, "global_coords");

    FLOAT *coords = ALLOC(DIM * N, FLOAT, "coords");
    for( i=0; i<N; i++ ){
	for( ii=0; ii<DIM; ii++ ){
	    fread(&coords[N*ii+i],sizeof(FLOAT),1,f);
	}
    }

    Point *pts = ALLOC(N, Point, "pts");
    init_and_fill_Points(coords, N, DIM, pts, &global_coords[ind], &ind);

    Box *boxes = ALLOC(nboxes, Box, "boxes");
    for( i=0; i<nboxes; i++ ){
	boxes[i].lo.DIM = DIM;
    	boxes[i].hi.DIM = DIM;
    	boxes[i].hi.x = &global_coords[ind]; ind += DIM;
    	boxes[i].lo.x = &global_coords[ind]; ind += DIM;
    }

    FLOAT *xlo = ALLOC(DIM, FLOAT, "xlo");
    FLOAT *xhi = ALLOC(DIM, FLOAT, "xhi");
    int mom,dau1,dau2,ptlo,pthi;
    Point *lo = ALLOC(1,Point,"lo");
    Point *hi = ALLOC(1,Point,"hi");
    lo->x = &global_coords[ind]; ind+=DIM;
    hi->x = &global_coords[ind]; ind+=DIM;
    lo->DIM = hi->DIM = DIM;
     /* Read the boxes */

    for( i=0; i<nboxes; i++ ){
	fread(xhi,DIM*sizeof(FLOAT),1,f);
	fread(xlo,DIM*sizeof(FLOAT),1,f);
	modify_Point(lo,xlo);
	modify_Point(hi,xhi);
	fread(&mom,sizeof(int),1,f);
	fread(&dau1,sizeof(int),1,f);
	fread(&dau2,sizeof(int),1,f);
	fread(&ptlo,sizeof(int),1,f);
	fread(&pthi,sizeof(int),1,f);
	modify_Box(lo,hi,mom,dau1,dau2,ptlo,pthi,&boxes[i]);
    }

    int *ptindx = ALLOC(N, int, "ptindx");
    int *rptindx = ALLOC(N, int, "rptindx");

    fread(ptindx, N * sizeof(int), 1, f);
    fclose(f);

    /* Reconstitute rptindx */
    for (j=0; j<N; j++) rptindx[ptindx[j]] = j;

    tree->N = N;
    tree->DIM = DIM;
    tree->pts = pts;
    tree->ptindx = ptindx;
    tree->rptindx = rptindx;
    tree->nboxes = nboxes;
    tree->boxes = boxes;
    tree->global_coords = global_coords;

    free(coords); free(xlo); free(xhi); free(hi); free(lo);
}


//---------------------------------------------------------------------
//--------------------   Locate functions      -----------------------
//---------------------------------------------------------------------


int locatePoint(Tree tree, Point pt) 
// " In which box of the tree is the arbitrary pt lying ?"
{
    int nb,d1,jdim;
    nb = jdim = 0;
    while (tree.boxes[nb].dau1) {
        d1 = tree.boxes[nb].dau1;
        if (pt.x[jdim] <= tree.boxes[d1].hi.x[jdim]) nb=d1;
        else nb=tree.boxes[nb].dau2;
        jdim = ++jdim % tree.DIM;
    }
    return nb;
}

int locateMember(Tree tree, int jpt) 
// " In which box of the tree is the member point lying ?"
{
    int nb,d1,jh;
    jh = tree.rptindx[jpt];
    nb = 0;
    while (tree.boxes[nb].dau1) {
        d1 = tree.boxes[nb].dau1;
        if (jh <= tree.boxes[d1].pthi) nb=d1;
        else nb = tree.boxes[nb].dau2;
    }
    return nb;
}



//--------------------------------------------------------------------
//--------------        KDtree applications         ------------------
//--------------------------------------------------------------------



int nearest(Tree tree, Point pt) 
// Returns the index of the closest neighbour of the arbitrary point x
{
    int i,k,nrst,ntask;
    int task[50];
    FLOAT dnrst = BIG, d;
    k = locatePoint(tree,pt);
    for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
        d = dist(tree.pts[tree.ptindx[i]],pt);
        if (d < dnrst) {
            nrst = tree.ptindx[i];
            dnrst = d;
        }
    }
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (distBox(tree.boxes[k],pt) < dnrst) {
            if (tree.boxes[k].dau1) {
                task[++ntask] = tree.boxes[k].dau1;
                task[++ntask] = tree.boxes[k].dau2;
            } else {
                for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
                    d = dist(tree.pts[tree.ptindx[i]],pt);
                    if (d < dnrst) {
                        nrst = tree.ptindx[i];
                        dnrst = d;
                    }
                }
            }
        }
    }
    return nrst;
}



void nnearest(Tree tree, int jpt, int *nn, FLOAT *dn, int n)
// finds the n neighbours of the point of index jpt, returns their
// indexes in nn and their distances to jpt in dn
{
    int j,i,k,ntask,kp;
    int task[50];
    FLOAT d;
    if (jpt > tree.N -1){
	printf("Error, index of point out of bounds.\n"); 
	return;
    }
    if (n > tree.N-1){ 
	printf("Error: too many neighbors requested\n"); 
	return;
    }
    for (i=0; i<n; i++) dn[i] = BIG;
    kp = tree.boxes[locateMember(tree,jpt)].mom;
    while (tree.boxes[kp].pthi - tree.boxes[kp].ptlo < n){
	kp = tree.boxes[kp].mom;
    }
    for (i=tree.boxes[kp].ptlo; i<=tree.boxes[kp].pthi; i++) {
        if (jpt == tree.ptindx[i]) continue;
        d = disti(tree,tree.ptindx[i],jpt);
        if (d < dn[0]) {
            dn[0] = d;
            nn[0] = tree.ptindx[i];
            if (n>1) sift_down(dn,nn,n);
        }
    }
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (k == kp) continue;
        if (distBox(tree.boxes[k],tree.pts[jpt]) < dn[0]) {
            if (tree.boxes[k].dau1) {
                task[++ntask] = tree.boxes[k].dau1;
                task[++ntask] = tree.boxes[k].dau2;
            } else {
                for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
                    d = disti(tree,tree.ptindx[i],jpt);
                    if (d < dn[0]) {
                        dn[0] = d;
                        nn[0] = tree.ptindx[i];
                        if (n>1) sift_down(dn,nn,n);
                    }
                }
            }
        }
    }
    //    We now sort the neighbours from closest to furthest:
    int *idx=malloc(n*sizeof(int));
    sort(n,dn,idx);
    permutation_int(n,idx,1,nn);
    free(idx);
    return;
}


int locatenear(Tree tree, Point pt, FLOAT r, int *list, FLOAT *dn, int nmax) 
/* 
   Given an arbitrary point of coordinates x, returns in list the 
   indexes of all members of the tree lying within r of it. nmax 
   is the maximum number of possible neighbours.
*/
{
    int n,k,i,nb,nbold,nret,ntask,jdim,d1,d2;
    int task[50];
    FLOAT d;
    int *idx=malloc(nmax*sizeof(int));
    nb = jdim = nret = 0;
    if (r < 0.0) { printf("Error: radius must be nonnegative\n"); return 1;}
    while (tree.boxes[nb].dau1) {
        nbold = nb;
        d1 = tree.boxes[nb].dau1;
        d2 = tree.boxes[nb].dau2;
        if (pt.x[jdim] + r <= tree.boxes[d1].hi.x[jdim]) nb = d1;
        else if (pt.x[jdim] - r >= tree.boxes[d2].lo.x[jdim]) nb = d2;
        jdim = ++jdim % tree.DIM;
        if (nb == nbold) break;
    }
    task[1] = nb;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (distBox(tree.boxes[k],pt) > r) continue;
        if (tree.boxes[k].dau1) {
            task[++ntask] = tree.boxes[k].dau1;
            task[++ntask] = tree.boxes[k].dau2;
        } else {
	    for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
		d = dist(tree.pts[tree.ptindx[i]],pt); 
		if (d <= r && nret < nmax){
		    dn[nret] = d;
		    list[nret] = tree.ptindx[i];
		    nret++;
		}
		if (nret == nmax){
		    sort(nret,dn,idx);
		    permutation_int(n,idx,1,list);
		    free(idx);
		    return nret;
		} 
	    }
        }
    }
    //    We now sort the neighbours from closest to furthest:
    sort(nret,dn,idx);
    permutation_int(nret,idx,1,list);
    free(idx);
    return nret;
}


//--------------------------------------------------------------------
//------------------        Memory cleaning         ------------------
//--------------------------------------------------------------------

void free_tree(Tree *tree){
    free(tree->global_coords);
    free(tree->boxes);
    free(tree->ptindx);
    free(tree->rptindx);
    free(tree->pts);
    free(tree);
}

FLOAT volume(FLOAT r, int DIM){
    if(DIM==1){ return 2*r;}
    if(DIM==2){ return M_PI*r*r; }
    if(DIM==3){ return 4*M_PI*r*r*r/3.; }
    fprintf(stderr,"Could not compute volume for DIM>3\n");
    return 0;
}




void companion_surface_density(int NNb, FLOAT *x, int N, int DIM, FLOAT *density){
    Tree *tree = ALLOC(1, Tree, "tree");
    FLOAT *nb_dist = ALLOC(N,FLOAT,"nb_dist");
    int *nb_list = ALLOC(N,int,"nb_list");
    int i;
    KDtree(x,N,DIM,tree);
    for( i=0; i<N; i++ ){
	nnearest( *tree, i, nb_list, nb_dist, NNb);
	density[i] = (NNb-1)/volume(nb_dist[NNb-1],DIM);
    }
    free_tree(tree);
    free(nb_dist);
    free(nb_list);
    return;
}



int main()
{

    /* Prepare variables */
    int N=100;
    int i,j;
    int a,DIM=3;
    FLOAT *x=ALLOC(DIM*N, FLOAT, "x"); // This is actually x,y,z in a single array
    Tree *tree = ALLOC(1, Tree, "tree");
    srand(time(NULL));

    /* Fill coords array with random values */
    for(i=0;i<DIM*N;i++)
    {
        x[i]=(FLOAT)rand()/(FLOAT)RAND_MAX;
    }
    printf("Coordinates generated...\n");

    /* Build the tree */
    KDtree(x,N,DIM,tree);
    free(x);
    printf("Tree built...\n");
   
    /* Create a point */
    Point pt;
    pt.DIM = DIM;
    pt.x = ALLOC(DIM, FLOAT, "pt.x");
    for(i=0;i<DIM;i++)
    {
        pt.x[i]=(FLOAT)rand()/(FLOAT)RAND_MAX;
    }
    printf("Random point created...\n");
    
    /* Prepare the call to locatenear */
    int *ind = ALLOC(N, int, "ind");
    FLOAT *dist = ALLOC(N, FLOAT, "dist");

    /* Ask the tree: which points of the tree are within 0.5 of the
       point we created earlier ? */
    printf("Looking for neighbours of the random point...\n");
    int n = locatenear(*tree, pt, 0.5, ind, dist, 1000 );

    /* Print a part of the result */
    if (n>10){
	printf("%d neighbours of the random point were found. First 10 results:\n",n);
	n=10;
    }else{
	printf("%d neighbours of random point were found:\n",n);
    }
    printf(" n_neighbour  neighbour-indice    neighbour-distance\n");
    for( i=0; i<n; i++ ){
	printf("%6.1d %15.1d %20.3lf\n",i,ind[i],dist[i]);
    }


    free(ind);
    free(dist);
    free_tree(tree);
}








