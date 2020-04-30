#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#define NMAX 100
int pipe_ge(double *AA, int n,  MPI_Comm comm);
int back_solve(double *A, int n, int nloc, double *x, MPI_Comm comm);
/*--------------- Gaussian elimination -- Pipelined version */
int main(int argc, char *argv[]){
  int n, myid, err, nprocs;
  double AA[NMAX*NMAX+NMAX];
  double x[NMAX];
  int i, j, nloc, np, iglob;
  double t;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  n = 6;
  nloc = n/nprocs;
  np = n+1;
/*--------------- generate local matrix and rhs */
  for (i=0; i<nloc; i++) {
    iglob = myid + i*nprocs;
    t = 0;
    for (j=0; j<n; j++) {
      AA[i*np+j]= 1.0 /(double)(iglob+j+1);
      if (j != iglob) t += AA[i*np+j]* (double)j;
    }
    /*--------------- reset diagonal entry       */
    AA[i*np+iglob] = (double) n ;
    /*--------------- set rhs[i]                */
    //    AA[i*np+n] = (double) n + t;
    AA[i*np+n] = (double) (n*iglob) + t;
  }
  //for (int i=0; i<10; i++)
    //printf("myid: %d, In main AA %.3f\n", myid, AA[i]);
/*--------------- call pipelined GE                */
  err = pipe_ge(AA,  n,  MPI_COMM_WORLD);
  if (err) {
    printf(" ERR : error in pipe_ge --  err= %d\n",err);
    MPI_Finalize();
    return(1);
  }
/*-------------------- done with GE */
/*-------------------- back-solve. last column of A contains RHS */
  err = back_solve(AA, n, nloc,  x, MPI_COMM_WORLD);
  if (err) {
    printf(" ERR : error in backsolv err= %d\n",err);
    MPI_Finalize();
    return(2);
  }
/*--------------- printing sol */
  MPI_Gather(&x,nloc,MPI_DOUBLE,&x,nloc,MPI_DOUBLE,0,MPI_COMM_WORLD);
  /* -------------------- solution is scattered -- need to write it
                          in  correct order */
  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0) {
    t = 0.0;
    printf(" Solution on Pe %d \n",myid);
    for (i=0; i<nloc; i++) {
      for (j=0; j<nprocs; j++){
	t+= pow(x[j*nloc+i]- (double) (i*nprocs+j),2);
	/*-------------------- print only when n is not large*/
	if (n<=64)
	  printf(" %9.3e", x[j*nloc+i]) ;
      }
      if (n <= 64)
	printf(" \n");
    }
    printf(" \n \n Err in solution = %8.3e \n",sqrt(t));
    printf(" \n");
  }
  MPI_Finalize();
}
