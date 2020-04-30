#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
int back_solve(double *A, int n, int nloc, double *x, MPI_Comm comm){
  /*--------------------
  / upper triangular solution
  / A = pointer to an n * (n+1) matrix. Only the upper triangular part of
        A is used. The right-hand side occupies the last column of A.
    The column version of the algorithm is implemented.
  */
  /*-------------------- NOTE: THIS JUST SETS X = RHS AND RETURNS */
  /* it is provided so that you can see how to access cerain elements of
     the matrix stored as a one-dimensional array                 */
  int k, id,  nprocs, np=n+1, myid, kloc;
  double t;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /*--------------- back-solve loop          */
  for (k=n-1; k>=0; k--) {
    id = (k % nprocs) ;
    if (id == myid) {
      kloc = (k-id)/nprocs;
      int r = kloc*np+n;  // the index of last column of A
      t = A[r]/A[k];

      /*-------------------- set x=rhs*/
      x[kloc] = t;
      printf("k = %d, myid = %d, x[%d] = %3f, A[%d]=%3f\n", k, myid, kloc, x[kloc], r, A[r]);
      // broadcase t to all rows
      MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
      
      // subtract multple of (relevant) part of column k from rhs
      // A[:][n] = A[:][n] - t*A[:][k]
      //x[kloc] -= t*A[r-1];

    }

    printf("at k = %d, myid = %d, recieved t = %3f\n", k, myid, t);

  }
  printf (" ------  end back-solve  in PE %d   ---------\n \n",myid);
  return(0);
}
