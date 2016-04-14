#define FLOAT double

//--------------------------------------------------------------------
//------------------        STRUCTURES             -------------------
//--------------------------------------------------------------------

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
    int N;         // Number of points in the tree
    int DIM;       // Dimension of the system
    Point *pts;    // Array of Points with coordinates
    int *ptindx;   // Keep track of the position of each point
    int *rptindx;  // reverse index: rptindx[ptindx[i]]=i
    int nboxes;    // Total number of boxes in the tree
    Box *boxes;    // Array of boxes
}Tree;



//--------------------------------------------------------------------
//-----   Structure management: creation, initialization       -------
//--------------------------------------------------------------------
void modify_Point(Point *pt, FLOAT *coords);
void init_and_fill_Points(FLOAT *coords,  int N,   int DIM, 
			  Point *pts, 
			  FLOAT *global_coords,  long *ind);
void modify_Box(Point *lo, Point *hi, 
		int mom, int dau1, int dau2, 
		int ptlo, int pthi, 
		Box *b);

//--------------------------------------------------------------------
//------------------   Distances functions       ---------------------
//-------------------------------------------------------------------
void swap(int *array, int i, int j);
FLOAT dist(Point a, Point b);
FLOAT distBox(Box b, Point p);
FLOAT disti(Tree tree, int j, int k);
int selecti(const int k, int *indx, int n, FLOAT *arr);
void sift_down(FLOAT *heap, int *ndx, int nn);

//---------------------------------------------------------------------
//--------        MAIN FUNCTION: Kdtree building          -------------
//---------------------------------------------------------------------
void  KDtree(FLOAT *coords, int N, int DIM, Tree *tree);

//---------------------------------------------------------------------
//------------        reading and writing         --------------------
//---------------------------------------------------------------------
void write_tree(Tree tree, char *name);
void read_tree(char *name, Tree *tree);

//---------------------------------------------------------------------
//--------------------   Locate functions      -----------------------
//---------------------------------------------------------------------
int locatePoint(Tree tree, Point pt);
int locateMember(Tree tree, int jpt);

//--------------------------------------------------------------------
//--------------        KDtree applications         ------------------
//--------------------------------------------------------------------
int nearest(Tree tree, Point pt);
void nnearest(Tree tree, int jpt, int *nn, FLOAT *dn, int n);
int locatenear(Tree tree, Point pt, FLOAT r, 
	       int *list, FLOAT *dist, int nmax);

//--------------------------------------------------------------------
//------------------        Memory cleaning         ------------------
//--------------------------------------------------------------------
void free_tree(Tree *tree);

