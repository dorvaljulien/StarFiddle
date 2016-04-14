/*
This module contains sort functions that keep track of the old
indices of sorted elements in the array.

An example can be found in the main at the end.

All functions are available for double, float, int. float and int versions
are explicit. For example:
"permutation" is for double
"permutation_float" is for float
"permutation_int" is for int
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
 
static int compar (const void *a, const void *b);
void sort(int n, double *arr, int *idx);
void sort_float(int n, float *arr, int *idx);
void permutation(int array_size, int *idx, int count, ...);
void permutation_float(int array_size, int *idx, int count, ...);
void permutation_int(int array_size, int *idx, int count, ...);

double *base_arr;
static int compar (const void *a, const void *b)
{
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])
    return -1;
  if (base_arr[aa] == base_arr[bb])
    return 0;
  if (base_arr[aa] > base_arr[bb])
    return 1;
}

int *ibase_arr;
static int icompar (const void *a, const void *b)
{
  int aa = *((int *) a), bb = *((int *) b);
  if (ibase_arr[aa] < ibase_arr[bb])
    return -1;
  if (ibase_arr[aa] == ibase_arr[bb])
    return 0;
  if (ibase_arr[aa] > ibase_arr[bb])
    return 1;
}


float *fbase_arr;
static int fcompar (const void *a, const void *b)
{
  int aa = *((int *) a), bb = *((int *) b);
  if (fbase_arr[aa] < fbase_arr[bb])
    return -1;
  if (fbase_arr[aa] == fbase_arr[bb])
    return 0;
  if (fbase_arr[aa] > fbase_arr[bb])
    return 1;
}
 
void sort(int n, double *arr, int *idx)
{
  int i,j;
  double *copy = malloc (sizeof (double) * n);
  for (i = 0; i < n; i++)
      idx[i] = i;
  base_arr = arr;
  qsort (idx, n, sizeof (int), compar);
  for(j=0;j<n;j++)
      copy[j]=arr[j];
  for(j=0;j<n;j++)
      arr[j]=copy[idx[j]];
  free(copy);
}

void sort_int(int n, int *arr, int *idx)
{
  int i,j;
  int *copy = malloc (sizeof (int) * n);
  for (i = 0; i < n; i++)
      idx[i] = i;
  ibase_arr = arr;
  qsort (idx, n, sizeof (int), icompar);
  for(j=0;j<n;j++)
      copy[j]=arr[j];
  for(j=0;j<n;j++)
      arr[j]=copy[idx[j]];
  free(copy);
}

void sort_float(int n, float *arr, int *idx)
{
  int i,j;
  float *copy = malloc (sizeof (float) * n);
  for (i = 0; i < n; i++)
      idx[i] = i;
  fbase_arr = arr;
  qsort (idx, n, sizeof (int), fcompar);
  for(j=0;j<n;j++)
      copy[j]=arr[j];
  for(j=0;j<n;j++)
      arr[j]=copy[idx[j]];
  free(copy);
}


void permutation(int array_size, int *idx, int count, ...)
/* 
   count is the number of arrays you want to arrange along the index array idx,
   for example to rearrange arrays bananas, apple and oranges along idx:
   permutation_FLOAT(N,idx,3,bananas,oranges,apple); 
*/
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    double *copy = malloc(array_size * sizeof(double));
    for (i = 0; i < count; i++)
    {
        double *array = va_arg (ap, double*);    /* Get the next argument value. */
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);     /* Clean up. */
    free(copy);
}

void permutation_float(int array_size, int *idx, int count, ...)
/* 
   count is the number of arrays you want to arrange along the index array idx,
   for example to rearrange arrays bananas, apple and oranges along idx:
   permutation_FLOAT(N,idx,3,bananas,oranges,apple); 
*/
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    float *copy = malloc(array_size * sizeof(float));
    for (i = 0; i < count; i++)
    {
        float *array = va_arg (ap, float*);    /* Get the next argument value. */
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);     /* Clean up. */
    free(copy);
}

void permutation_int(int array_size, int *idx, int count, ...)
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    int *copy = malloc(array_size * sizeof(int));
    for (i = 0; i < count; i++)
    {
        int *array = va_arg (ap, int*);    /* Get the next argument value. */
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);                  /* Clean up. */
    free(copy);
}



int main()
{

    int i,N=10;
    int *x = malloc(N * sizeof(int));
    int *y = malloc(N * sizeof(int));
    int *z = malloc(N * sizeof(int));
    int *idx = malloc(N * sizeof(int));
    
    srand(time(NULL));
    for (i=0;i<N;i++)
    {   x[i] = rand() % 100;
        y[i] = i;
        z[i]=i*i;
    }
    
    printf("Random |  i  |  i*i\n");
    for (i=0;i<N;i++)
    {
        printf("%5.1d %5.1d %5.1d\n",x[i],y[i],z[i]);
    }

    sort_int(N,x,idx);
    permutation_int(N,idx,2,y,z);
    printf("We sort the random column, then arrange i and i*i with the same permutation.\n");


    for (i=0;i<N;i++)
    {
        printf("%5.1d %5.1d %5.1d\n",x[i],y[i],z[i]);
    }

    free(x); free(y); free(z); free(idx);


}



