#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#define NMAX 100
/*--------------- Gaussian elimination -- Pipelined version */
int pipe_ge(double *AA, int n,  MPI_Comm comm){
  /* Pipelined Gaussian Elimination
     AA = pointer to matrix and right-hand side. AA is a one-dimensional
	  array holding a matrix of size n*(n+1)-- row-wise. The
          right-hand side occupies the last column.
     n  = dimension of problem
     On return Gaussian elimination has been applied and only upper triangular
     part is relevant -- back-solve should be called to get solution x.
     *----------------------------------------------------------------------*/
  int nprocs, myid, nloc, start, South, North;
  double *x, *piv;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  South = (myid+1)%nprocs;  // South and North values are between 0 and nprocs
  North = (nprocs+(myid-1))%nprocs;
  nloc = n/nprocs;
  x = (double *) malloc((NMAX*NMAX+NMAX)*sizeof(double)); // buffer to receive row k
  piv = (double *) malloc((NMAX*NMAX+NMAX)*sizeof(double)); // buffer to save pivots

  int row;
  int j = 0;
  for (int i=0; i<(n+1)*nloc; i++)
    printf("In proc %d, AA %.3f\n", myid, AA[i]);
  for (int k=0; k<n-1; k++){
    // if row k is in proc myid
    if (k%nprocs == myid){
      row = j*(n+1);  // which row to send
      //  send row k to proc South
      MPI_Send(&AA[row], n+1, MPI_DOUBLE, South, 0, MPI_COMM_WORLD) ;
      //printf("At k = %d, Sending row j %d (%d) from proc %d to proc %d (South)\n",k, j, row, myid, South);
      j += 1;  // counter to keep track of the current row to send in each proc
      // since A is 1D array, the first row is from A[0:n+1], next row is from A[n+1:2n+1] etc.
    }
    else{
      //row = j*(n+1);  // which row to receive
      //  receive row k from proc North
      MPI_Recv(&x[0], n+1, MPI_DOUBLE, North, 0, MPI_COMM_WORLD, &status);
      //printf("At k = %d, Receiving row j %d (%d) on proc %d from proc %d (North)\n",k, j, row, myid, North);

      // for the rest of the proc, send row k to South only if South is not in proc myid
      // no need to send to its own proc; for eg, row 0 from proc 0 doen't need to receive row 0
      if (South != k%nprocs){
        MPI_Send(&x[0], n+1, MPI_DOUBLE, South, 0, MPI_COMM_WORLD) ;
        //printf("At k = %d, extra Sending row j %d (%d) from proc %d to proc %d (South)\n", k, j, row, myid, South);
      }

    }

    // get starting point in each proc
    start = (myid <= k%nprocs) ? (int) k/nprocs+1 : (int) k/nprocs;

    // loop through rows in each proc to compute pivot
    for (int i=start; i<nloc; i++){
      int count = 0;
      int r = i*(n+1) + k;  // index of A[i,k] in 1D
      if(k%nprocs==myid){ // for the same proc that sent A[k]
        // compute piv = A[i][k]/A[k][k]
        piv[i] = AA[r]/AA[k];
        //printf("at k = %d, In proc %d, AA[%d] %.3f, AA[%d] %3f, piv[%d] = %3f\n",k, myid, r, AA[r], k, AA[k], i, piv[i]);

        // process elimination for relevant row in proc myid; A[i][j] -= piv*A[k][j]
        for (int j=r+1; j<(i+1)*(n+1); j++){
          count += 1;
          AA[j] -= piv[i]*AA[row+k+count];
          //printf("elimination: at k = %d, myid = %d, A[%d]: %3f -= piv[%d] * AA[%d] = %3f\n",k, myid, j, AA[j], i, row+k+count, x[row+k+count]);
        }
      }
      else{
        // we received A[k] from North and saved in x[k] so use x[k] instead of A[k]
        piv[i] = AA[r]/x[k];
        //printf("at k = %d, In proc %d, AA[%d] %.3f, x[%d] %3f, piv[%d] = %3f\n",k, myid, r, AA[r], k, x[k], i, piv[i]);

        // process elimination for relevant row in proc myid; A[i][j] -= piv*A[k][j]
        for (int j=r+1; j<(i+1)*(n+1); j++){
          count += 1;
          AA[j] -= piv[i]*x[k+count];
          //printf("elimination: at k = %d, myid = %d, A[%d]: %3f -= piv[%d] * x[%d] = %3f\n",k, myid, j, AA[j], i, k+count, x[k+count]);
        }
      }
      AA[r] = 0;

    }
    // at this point, we want every proc to be finished so wait for all of them before starting next k loop
    MPI_Barrier(MPI_COMM_WORLD);

  }
  // After all k loop, AA is a upper
  //for (int i=0; i<(n+1)*nloc; i++)
    //printf("After In proc %d, AA %.3f\n", myid, AA[i]);



  return(0);
}
