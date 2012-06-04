#include "ki.h"
#include "comm.h"

/* MPI wrapper */
void kiBarrier() {
  KI_Barrier(ki_cmm_domain);
}

int KI_Barrier(MPI_Comm comm) {
  if (kiIsParallel()) 
    return MPI_Barrier(comm);
  else 
    return 0;
}

int KI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  if (kiIsParallel()) 
    return MPI_Bcast(buffer, count, datatype, root, comm);
  else 
    return 0;
}

int KI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  if (kiIsParallel()) 
    return MPI_Send(buf, count, datatype, dest, tag, comm);
  else
    return 0;
}

int KI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status) {
  if (kiIsParallel())
    return MPI_Recv(buf, count, datatype, source, tag, comm, status);
  else 
    return 0;
}

int KI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  if (kiIsParallel())
    return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  else {
    int sz;
    MPI_Type_size(datatype, &sz);
    memcpy(recvbuf, sendbuf, count * sz);
    return 0;
  }
}


int KI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  if (kiIsParallel())
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  else {
    int sz;
    MPI_Type_size(datatype, &sz);
    memcpy(recvbuf, sendbuf, count * sz);
    return 0;
  }
}


int KI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  /* should work if !kiIsParallel() as well */
  return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

int KI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  /* should work if !kiIsParallel() as well */
  return MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm);
}

int KI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  /* should work if !kiIsParallel() as well */
  return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

int KI_Scatterv(void *sendbuf, int *sendcnts, int *displs, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  /* should work if !kiIsParallel() as well */
  return MPI_Scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}


/* Enhanced MPI calls */
void KI_Bcast_buffer(void* buf, int count, int root, MPI_Comm comm) {
  int bufSize = count;
  KI_Bcast(&bufSize, 1, MPI_INT, root, comm);
  KI_Bcast(buf, bufSize, MPI_CHAR, root, comm);
}

void KI_Bcast_string(char* str, int root, MPI_Comm comm) {
  KI_Bcast_buffer(str, strlen(str)+1, root, comm);
}


void KI_Gatherv_buffer(void* sendbuf, int sendcnt, void* recvbuf, int recvbuflimit, int root, MPI_Comm comm, /*OUT*/int* recvbuflen) {
  /* simply patch together buffers into a contiguous buffer with no delimiter */
  
  if (!kiIsParallel()) {
    assert(sendcnt <= recvbuflimit);
    memcpy(recvbuf, sendbuf, sendcnt);
    *recvbuflen = sendcnt;
    return;
  }

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int *recvcnts = 0, *displs = 0;
  if (rank == root) {
    recvcnts = (int*)kiArenaMalloc(KI_ARENA_FARMER, sizeof(int) * size);
    displs   = (int*)kiArenaMalloc(KI_ARENA_FARMER, sizeof(int) * size);
  }

  KI_Gather(&sendcnt, 1, MPI_INT, recvcnts, 1, MPI_INT, root, comm);
  if (rank == root) {
    int j;
    for (j = 0; j < size; ++j) {
      kipmsg(7, "recvcnts[%d] = %d\n", j, recvcnts[j]);
    }
  }
  int i, *pr, *pd, displ_sum = 0;
  if (rank == root) {
    for (i = 0, pr = recvcnts, pd = displs; i < size; ++i, ++pr, ++pd) {
      if (displ_sum + *pr > recvbuflimit) {
        kipm("displ_sum + *pr = %d > recvbuflimit\n", displ_sum + *pr);
        *pr = 0;
      }
      *pd = displ_sum;
      displ_sum += *pr;
      kipmsg(7, "displ[%d] = %d\n", i, *pd);
    }
    /* assert(displ_sum <= recvbuflimit); */
  }
  
  KI_Gatherv(sendbuf, sendcnt, MPI_CHAR, recvbuf, recvcnts, displs, MPI_CHAR, root, comm);
  *recvbuflen = displ_sum;

  if (!KI_USE_ARENA) {
    if (rank == root) {
      kifree(recvcnts, (sizeof(int) * size));
      kifree(displs, (sizeof(int) * size));
    }
  } else {
    kiArenaClear(KI_ARENA_FARMER);
  }
}


