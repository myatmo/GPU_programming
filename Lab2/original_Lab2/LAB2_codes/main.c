#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> 
#define MAXSIZE 100000
#define MAXPES 32 
#define DEBUG_MODE 0

int main(int argc, char *argv[]){
/*---------------------------------------------------------------------
 *
 * Main. Each process generates n data integers  [semi randomly] which
 * it puts  in array Loc_dat. It  then sorts this array  and obtains a
 * list[] of pivots  to be used by the  HQuicksort function. This list
 * is broadcast  to all  PEs. HQuicksort is  then invoked to  sort the
 * data. On return the root PE gathers all data [using an MPI_Gatherv]
 * and  checks that the  array is  indeed sorted.
 * ------------------------------ 
 * YS. 03/06/2020 for csci 5451 /
 *-----------------------------------------------------------------------*/
  int ntot, n, naa, i, err, myid, seed, j, count, nprocs;
  int Loc_dat[MAXSIZE],  Glob_dat[MAXSIZE];
  int *sav_dat;
  /*-------------------- aa and list have max size of p*log p */
  int aa[20*MAXPES]; 
  int list[20*MAXPES]; 
  double etim;
  /*---------- file for seeing results on each pe (debugging) */
  FILE *fd;
  char strn[9];  
/*---------- prototype functions   */
  int HQuicksort(int *Abuf, int *lenA,  int *list, int N, MPI_Comm comm) ;
  void MergeSort (int n, int *a, int* b);
  int makList (int *list, int *a, int size, int hsize);
/*---------- begin */
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  sprintf(strn,"DAT%1d.dat",myid);
  if (DEBUG_MODE) 
    fd =fopen(strn,"w");
/*---------- define ntot and n (n=actual size of Loc_dat) */
  // just a test
  n = 1000 ;
  ntot = n*nprocs;
/*---------- generate some data */
  seed = 163 + myid*211;
  for (j=0; j<n; j++) 
    Loc_dat[j] = ((j+1)*157 + 397*seed)  %   278459;  
  /* this is for checking correctness at end */
  sav_dat = (int *) malloc(ntot*sizeof(int));
  /*---------- save data into an array and gather it into
               into root for later use */
  memcpy(sav_dat,Loc_dat,n*sizeof(int));	   
  MPI_Gather(sav_dat,n,MPI_INT,sav_dat,n,MPI_INT,0,MPI_COMM_WORLD);  
/*---------- initial sort local  data */
/*           here Glob_dat used as temp. buffer */
  MergeSort (n,Loc_dat,Glob_dat);
/*---------- sampling to get pivots. The global final list 
             is built from local lists of the same length*/
  count = nprocs-1;
/*---------- get my local list */
  makList (list, Loc_dat, n, count);
/*---------- gather  these lists */
  MPI_Gather(list,count,MPI_INT,aa,count,MPI_INT,0, MPI_COMM_WORLD) ;
/*---------- sort them and get the pivots by calling makList */
  if (myid == 0) {
    naa = nprocs*count;
    MergeSort (naa, aa, &aa[count*nprocs]);
    makList(list, aa, naa, nprocs) ;
    etim = MPI_Wtime();
  } 
  MPI_Bcast(list, count, MPI_INT, 0, MPI_COMM_WORLD) ; 
  MPI_Barrier(MPI_COMM_WORLD);
/*---------- call hypercube sort */  
  err = HQuicksort(Loc_dat, &n, list,  n*nprocs, MPI_COMM_WORLD);
  if (err) {
    printf(" ERR : error in HQquicksort err= %d\n",err);
    MPI_Finalize(); 
    return(1);
  }
/*---------- use aa to gather sizes of arrays */
  MPI_Gather(&n,1,MPI_INT,aa,1,MPI_INT,0,MPI_COMM_WORLD) ;
  if (myid==0) {
    for (i=0; i< nprocs ; i++) 
      printf (" Number of items in Pe %2d = %4d\n",i,aa[i]);
/*---------- use list for pointers to global array */
    list[0] =0;
    for (i=0; i< nprocs ; i++) 
/*---------- transform list into pointer to recv buffer */
      list[i+1] = list[i] + aa[i];
/*---------- now ready for gatgerv to gather all data */    
  }
  MPI_Gatherv (Loc_dat,n,MPI_INT,Glob_dat,aa,list,MPI_INT,0,MPI_COMM_WORLD);
  /*also get saved data for comparison */
  if (myid == 0) {
    etim = MPI_Wtime()-etim;
    
    MergeSort(ntot, sav_dat, Loc_dat);
    for (i=0; i<ntot; i++){
      if (Glob_dat[i] != sav_dat[i]){
	printf("ERROR: SORT FAILED at i = %d\n",i);
	MPI_Finalize(); 
	return(1);
      }
    }
    printf("SORT SUCCESSFULL! \n");
    printf("Elaps. time %e sec.  \n",etim);
  }
/*---DEBUG Mode------- print processed list from each PE into a file */
  if (DEBUG_MODE) {
    fprintf (fd," Number of items in Pe %2d = %4d\n",myid,n);
    for (i=0; i<n; i++)
      fprintf(fd," i %4d a[i] %4d\n",i,Loc_dat[i]);
    fprintf (fd," \n ------  end list in PE %d  --------------\n",myid);
    fclose(fd);  
  }
  free(sav_dat);
  MPI_Finalize(); 
  return(0);
}

