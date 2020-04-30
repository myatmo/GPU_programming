#include <stdio.h>
#include <string.h>
#include <stdlib.h> 

void merge (  int n, int m, int *a, int *b, int *c)
{
  // Given two arrays a, and b, that are sorted
  // sort the union of a and b into the array c.
  // n = length of a, 
  // m = length of b
  // c must be of length n+m

  int j, ia = 0, ib = 0, ic= 0;
  while (ia < n && ib < m) {
    if (a[ia] < b[ib]) 
      c[ic++] = a[ia++]; 
    else c[ic++] = b[ib++];
  }
//-------------------- a not emptied: copy remainder to c[]
  if (ia < n) 
    for (j=ia; j<n;j++) c[ic++] = a[j]; 
//-------------------- b not emptied: copy remainder to c[]  
  if (ib < m) 
      for (j=ib; j<m;j++) c[ic++] = b[j]; 
}

void MergeSort (int n, int *a, int *b)
{
/*---------- Standard Mergesort algorithm 
 array a is to be sorted. array b is used as
 temp space. both are assumed to be of length n
 on return a is sorted increasingly.
 Algorithm is recursive. 
 ---------- */
  int n1, n2; 
  if (n > 1)      {
    n1 = (n+1)/2;
    n2 = n-n1;
    MergeSort(n1,a,b);
    MergeSort(n2,&a[n1],b);
    merge (n1,n2,a,&a[n1],b); 
    memcpy(a,b,n*sizeof(int));
  }
}

int makList (int *list, int *a, int size, int hsize){
  /*---------- given a sorted array a of length 'size' 
   * we want to set up a binary tree of `medians'. 
   * NOTE: hsize = 2^k-1 for some k. List is of length hsize.
   ---------- */
  int i0, h0=0, j, stra;
  i0 = (size-1)/2;  
  stra = size/2;
  list[h0++] = a[i0];
//-------------------- main loop
  while (h0<hsize-1) {
    i0 = i0/2;
    for (j=i0; j<size; j+=stra) 
      list[h0++] = a[j];
    stra /= 2; 
  } 
  return(0); 
}

