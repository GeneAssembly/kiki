#ifndef _IO_H_
#define _IO_H_

/* File wrapper */
typedef struct {
  FILE*    fp;
  MPI_File fh;
} KI_File;

int KI_File_open(char *filename, int amode, MPI_Info info, KI_File *fh);
int KI_File_close(KI_File *fh);
int KI_File_get_size(KI_File fh, MPI_Offset *size);
int KI_File_seek(KI_File fh, MPI_Offset offset, int whence);
int KI_File_read(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements);
int KI_File_write(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements);
int KI_File_read_shared(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements);
int KI_File_write_shared(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements);

int KI_File_copy(char* filename1, char* filename2, bool replace);
int KI_File_backup(char* filename);
int KI_File_remove(char* filename);

bool KI_File_exists(char* filename);

void kifprintf(KI_File fh, char* format, ...);


#endif /* _IO_H_ */
