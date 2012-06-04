#ifndef _COMM_H_
#define _COMM_H_

#include "mpi.h"

/* MPI wrapper */
void kiBarrier();

int  KI_Barrier(MPI_Comm comm);
int  KI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int  KI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int  KI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status);
int  KI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int  KI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int  KI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int  KI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm);
int  KI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int  KI_Scatterv(void *sendbuf, int *sendcnts, int *displs, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);


/* Enhanced MPI calls */
void KI_Bcast_buffer(void* buf, int count, int root, MPI_Comm comm);
void KI_Bcast_string(char* str, int root, MPI_Comm comm);
void KI_Gatherv_buffer(void* sendbuf, int sendcnt, void* recvbuf, int recvbuflimit, int root, MPI_Comm comm, /*OUT*/int* recvbuflen);

#endif /* _COMM_H_ */
