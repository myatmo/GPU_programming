#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#define MAXPES 32
#define MAXSIZE 100000

int HQuicksort(int *Abuf, int *lenA,  int *list, int N, MPI_Comm comm) {
  // Abuf = local array
  //list = list of pivots  list from top to bottom of tree.
  //size limit for local buffers
  int myid, nprocs, n, ncub, q, name_len;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  int A_low[MAXSIZE], A_high[MAXSIZE], temp[MAXSIZE];
  int *recvcounts, *displs;
  MPI_Status status;
  MPI_Request req;
  const int root = 0;
  lenA = Abuf;

  // prototype functions
  int check_bit(int num, int n);
  void merge(int n, int m, int *a, int *b, int *c);

  //get myid (rank)
  MPI_Comm_rank(comm, &myid);
  //get size
  MPI_Comm_size(comm, &nprocs);
  // get the processor name
  MPI_Get_processor_name(proc_name, &name_len);
  // get the number of data in each proc
  n = N/nprocs;
  //get ncub so that 2^ncub = size
  ncub = log2(nprocs);
  //printf("from processor %d out of %d, executing on %s\n", myid, nprocs, proc_name);

  // check if the nprocs is the power of 2, if not exit with error
  if (myid == root){
    if(ceil(log2(nprocs)) != floor(log2(nprocs))){
      MPI_Finalize();
      printf("nprocs is not a power of 2!");
      exit(0);
      return(1);
    }
  }

  // initialize recvcounts which will collect the different lengths in each proc
  // after Abuf is sorted in each proc
  if (myid == root){
    recvcounts = (int *)malloc(nprocs*sizeof(int));
    displs = (int *)malloc(nprocs*sizeof(int));
  }

  // loop through each dimension in hypercube starting from highest to lowest
  for(int j=ncub; j>0; j--){
    // myid+q is the number of porc to exchange with proc myid; myid < q
    q = pow(2, j-1);

    // t_idx is the strating index of list in each level of tree
    int t_idx = pow(2, ncub-j)-1;
    int nsub = nprocs/(pow(2, ncub-j)); // number of nodes in sub-cubes at current level
    int t = t_idx + floor(myid/nsub); // pivot idex to map on hypercube
    //printf("j = %d, nsub = %d, myid = %d maps t_p= %d \n", j, nsub, myid, t);

    int piv = list[t];
    int s = 0;// from A[0] to A[s-1] is A_low and from A[s] to A[n] is A_high

    while (Abuf[s] <= piv)
      s += 1;

    // check the len of msg that is to be split
    // if the len is equal to n or 0 then one proc is trying to send all the data
    // and there is no data left to be sent, which cause an error count in MPI_Send and Recv
    // due to the load-imbalance in data, the program will quit with an err for now
    if (!(s < n && s > 0)){
      printf("Load-imbalance occurs!");
      exit(0);
      return(1);
    }

    int bit = check_bit(myid, j); // get the jth bit
    int send, recv, total;  // total is the number in each proc after exchanging messages

    if (bit == 0){
      // at first total is n in each proc
      if (j==ncub)
        send = n-s;
      else
        send = total-s;

      // send higher part of array to higher dimension
      MPI_Isend(&send, 1, MPI_INT, myid+q, 0, comm, &req); // send the size of array first
      MPI_Isend(&Abuf[s], send, MPI_INT, myid+q, 0, comm, &req);

      // receive lower part of array from higher dimension and store it in A_low
      MPI_Recv(&recv, 1, MPI_INT, myid+q, 0, comm, &status);  // receive the size of array
      MPI_Recv(A_low, recv, MPI_INT, myid+q, 0, comm, &status);

      // merge the arrays
      merge(s, recv, &Abuf[0], A_low, temp);
      total = s+recv;
    }
    else{
      if (j==ncub)
        send = n-s;
      else
        send = total-s;

      // receive higher part of array from lower dimension and store it in A_high
      MPI_Recv(&recv, 1, MPI_INT, myid-q, 0, comm, &status);
      MPI_Recv(A_high, recv, MPI_INT, myid-q, 0, comm, &status);

      // send lower part of array to lower dimension
      MPI_Isend(&s, 1, MPI_INT, myid-q, 0, comm, &req);
      MPI_Isend(&Abuf[0], s, MPI_INT, myid-q, 0, comm, &req);

      // merge the arrays
      merge(send, recv, &Abuf[s], A_high, temp);
      total = send+recv;
    }
    MPI_Barrier(comm);
    memset(Abuf, 0, N*sizeof(Abuf));
    memcpy(Abuf, temp, total*sizeof(int));
    memset(temp, 0, N*sizeof(int));

    if (j == 1){
      // gather the total length in each proc
      MPI_Gather(&total, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);

      if (myid == root){
        // count displacements
        for (int i=0; i<nprocs; i++) {
          if (i==0)
            displs[i] = 0;
          displs[i] = displs[i-1] + recvcounts[i-1];
        }
      }
      // gather all sorted values from all proc in root
      MPI_Gatherv(Abuf, total, MPI_INT, temp, recvcounts, displs, MPI_INT, root, comm);

      // scatter collected sorted values to all proc; each proc gets n numbers
      MPI_Scatter(temp, n, MPI_INT, Abuf, n, MPI_INT, root, comm) ;
    }
  }
  return(0);
}

int check_bit(int num, int n){
  // right shift num n times and bitwise AND with 1 to get the nth bit in num
  return (num >> (n-1)) & 1;
}
