#include "extern.h"
#include "io.h"

/* File wrapper */

int KI_File_open(char *filename, int amode, MPI_Info info, KI_File *fh) {
  if (kiIsParallel() && ki_domain_size > 1) {
    if ((amode & MPI_MODE_WRONLY) && !(amode & MPI_MODE_APPEND)) {
      MPI_File_open(ki_cmm_domain, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE, info, &(fh->fh));
      MPI_File_close(&(fh->fh));
      amode |= MPI_MODE_CREATE;
    }
    return MPI_File_open(ki_cmm_domain, filename, amode, info, &(fh->fh));
  } else {
    char mode[10];
    if (amode & MPI_MODE_RDONLY) {
      sprintf(mode, "r");
    } else if ((amode & MPI_MODE_WRONLY) && !(amode & MPI_MODE_APPEND)) {
      sprintf(mode, "w");
    } else if (amode & MPI_MODE_APPEND) {
      sprintf(mode, "a");
    }
    fh->fp = fopen(filename, mode);
    if (fh->fp == NULL) return !MPI_SUCCESS;
    else return MPI_SUCCESS;
  }
}

int KI_File_close(KI_File *fh) {
  if (kiIsParallel() && ki_domain_size > 1) 
    return MPI_File_close(&(fh->fh));
  else {
    fclose(fh->fp);
    return MPI_SUCCESS;
  }
}

int KI_File_get_size(KI_File fh, MPI_Offset *size) {
  if (kiIsParallel() && ki_domain_size > 1) 
    return MPI_File_get_size(fh.fh, size);
  else {
    fseek(fh.fp, 0L, SEEK_END);
    *size = ftell(fh.fp);
    fseek(fh.fp, 0L, SEEK_SET);
    return MPI_SUCCESS;
  }
}

int KI_File_seek(KI_File fh, MPI_Offset offset, int whence) {
  if (kiIsParallel() && ki_domain_size > 1) 
    return MPI_File_seek(fh.fh, offset, whence);
  else {
    int w = whence;
    if (whence == MPI_SEEK_SET) w = SEEK_SET;
    else if (whence == MPI_SEEK_END) w = SEEK_END;
    return fseek(fh.fp, offset, w);
  }
}

int KI_File_read(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements) {
  *elements = 0;
  int ret = 0;
  if (kiIsParallel() && ki_domain_size > 1)  {
    ret = MPI_File_read(fh.fh, buf, count, datatype, status);
    if (ret != MPI_SUCCESS) {
      char estr[MPI_MAX_ERROR_STRING];
      int len;
      MPI_Error_string(ret, estr, &len);
      fprintf(stderr, "MPI_File_read failed: %s\n", estr);
      return ret;
    }
    MPI_Get_count(status, datatype, (int *)elements);
  } else {
    int sz;
    MPI_Type_size(datatype, &sz);
    *elements = fread(buf, sz, count, fh.fp);
    return ferror(fh.fp) ? !MPI_SUCCESS : MPI_SUCCESS;
  }
  return ret;
}

int KI_File_write(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements) {
  *elements = 0;
  if (kiIsParallel() && ki_domain_size > 1) {
	int ret = MPI_File_write(fh.fh, buf, count, datatype, status);
  	if (ret == MPI_SUCCESS) {
  	  MPI_Get_count(status, datatype, (int *)elements);
  	}
    return ret;
  } else {
    int sz;
    MPI_Type_size(datatype, &sz);
    *elements = fwrite(buf, sz, count, fh.fp);
    return ferror(fh.fp) ? !MPI_SUCCESS : MPI_SUCCESS;
  }
}

int KI_File_read_shared(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements) {
  if (kiIsParallel() && ki_domain_size > 1) {
    int ret = MPI_File_read_shared(fh.fh, buf, count, datatype, status);
    if (ret == MPI_SUCCESS) {
      MPI_Get_count(status, datatype, (int *)elements);
    }
    return ret;
  } else {
    int sz;
    MPI_Type_size(datatype, &sz);
    *elements = fread(buf, sz, count, fh.fp);
    fprintf(stderr, "File_read_shared elements: %u\n", *elements);
    fflush(fh.fp);
    return ferror(fh.fp) ? !MPI_SUCCESS : MPI_SUCCESS;
  }
}

int KI_File_write_shared(KI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status, size_t *elements) {
  if (kiIsParallel() && ki_domain_size > 1) {
	int ret = MPI_File_write_shared(fh.fh, buf, count, datatype, status);
	if (ret == MPI_SUCCESS) {
		MPI_Get_count(status, datatype, (int *)elements);
	}
    return ret;
  } else {
    int sz;
    MPI_Type_size(datatype, &sz);
    *elements = fwrite(buf, sz, count, fh.fp);
    fprintf(stderr, "File_write_shared elements: %u\n", *elements);
    fflush(fh.fp);
    return ferror(fh.fp) ? !MPI_SUCCESS : MPI_SUCCESS;
  }
}

extern char ki_tmp_buf[KI_BUF_SIZE];
int KI_File_copy(char* filename1, char* filename2, bool replace) {
  if (!KI_File_exists(filename1)) return -1;
  if (KI_File_exists(filename2) && !replace) return !MPI_SUCCESS;

  FILE *fp1, *fp2;
  if ((fp1 = fopen(filename1, "rb")) == NULL) {
    fprintf(stderr, "Cannot read file: %s\n", filename1);
    return !MPI_SUCCESS;
  }
  if ((fp2 = fopen(filename2, "wb")) == NULL) {
    fprintf(stderr, "Cannot write to file: %s\n", filename2);
    return !MPI_SUCCESS;
  }

  int nRead;
  while (feof(fp1) == 0) {
    if ((nRead = fread(ki_tmp_buf, 1, KI_BUF_SIZE, fp1)) != KI_BUF_SIZE){
      if (ferror(fp1) != 0 || feof(fp1) == 0) {
        fprintf(stderr,"Read file error.\n");
        return !MPI_SUCCESS;
      }
    }
    if ((fwrite(ki_tmp_buf, 1, nRead, fp2)) != nRead) {
      fprintf(stderr, "Write file error.\n");
      return !MPI_SUCCESS;
    }
  }
  return 0;
}

int KI_File_backup(char* filename) {
  char backup[KI_NAME_BUF_SIZE];
  sprintf(backup, ".%s.bak", filename);
  return KI_File_copy(filename, backup, /*replace*/true);
}

int KI_File_remove(char* filename) {
  if (!KI_File_exists(filename)) return MPI_SUCCESS;
  /* MPI_File fh; */
  /* MPI_File_open(ki_cmm_domain, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &fh); */
  /* MPI_File_close(&fh); */
  int ret = remove(filename);
  return ret;
}

bool KI_File_exists(char* filename) {
  FILE* fp = fopen(filename, "r");
  if (fp != 0) {
    fclose(fp);
    return true;
  }
  return false;
}


void kifprintf(KI_File fh, char* format, ...) {
  va_list args;
  char buf[1024*1024];
  
  va_start(args, format);
  vsprintf(buf, format, args);
  va_end(args);

  int len = strlen(buf);
  MPI_Status status;
  size_t elements;
  KI_File_write_shared(fh, buf, len, MPI_CHAR, &status, &elements);
}


