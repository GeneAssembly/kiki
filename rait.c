#include "extern.h"
#include <math.h>
#include "raiphy.h"



void trainFastaList(char* fileName, char* dbName, int kmerSize) {
  FILE* fp = fopen(fileName, "r");
  if (fp == NULL) {
    fprintf(stderr, "Could not open trainning input file %s.\n", fileName);
    kiAbort(1);
  }
  
  FILE* fd = fopen(dbName, "w");
  if (fd == NULL) {
    fprintf(stderr, "Could not write to database file %s.\n", dbName);
    kiAbort(1);
  }

  kmer_freq_t* freq = kiAllocKmerFreq(kmerSize);
  int dim = 1 << (2*kmerSize);
  double vector[dim];
  rai_db_t* db = kiAllocRaiDb(kmerSize);

  char buf[KI_NAME_BUF_SIZE] = "";

  char oldClass[KI_NAME_BUF_SIZE] = "";
  char className[KI_NAME_BUF_SIZE] = "";
  char fastaName[KI_NAME_BUF_SIZE] = "";

  /* FILE* fa; */
  
  while (fgets(buf, sizeof(buf), fp) != NULL) {
    char *p, *q;

    for (p = buf; *p != '\0' && *p != '\n' &&  *p != '\t' && *p != ' ';  p++);  *p = '\0';
    for (p = p+1; *p != '\0' && *p != '\n' && (*p == '\t' || *p == ' '); p++);
    for (q = p;   *q != '\0' && *q != '\n' &&  *q != '\t' && *q != ' ';  q++);  *q = '\0';

    strcpy(className, buf);
    strcpy(fastaName, p);
    fprintf(stderr, "class=%s \t fa=%s\n", className, fastaName);
    
    if (strcmp(className, oldClass) != 0) {
      if (strlen(oldClass) > 0) {
        /* add old profile vector */
        kiFreqToRaiVector(freq, vector);
        kiRaiDbAdd(db, oldClass, vector);
      }
      
      /* start new vector */
      kiClearKmerFreq(freq);
      strcpy(oldClass, className);
    }

    kiScanFastaForFreq(fastaName, freq);
  }

  kiFreqToRaiVector(freq, vector);
  kiRaiDbAdd(db, className, vector);

  kiPrintRaiDb(fd, db);

  freq = kiFreeKmerFreq(freq);
  db = kiFreeRaiDb(db);

  fclose(fp);
  fclose(fd);
}

void trainFastaListOriginal(char* fileName, char* dbName, int kmerSize) {
  FILE* fp = fopen(fileName, "r");
  if (fp == NULL) {
    fprintf(stderr, "Could not open trainning input file %s.\n", fileName);
    kiAbort(1);
  }
  
  FILE* fd = fopen(dbName, "w");
  if (fd == NULL) {
    fprintf(stderr, "Could not write to database file %s.\n", dbName);
    kiAbort(1);
  }

  kmer_freq_t* freq = kiAllocKmerFreq(kmerSize);
  int dim = 1 << (2*kmerSize);
  double vector[dim];
  rai_db_t* db = kiAllocRaiDb(kmerSize);

  char buf[KI_NAME_BUF_SIZE] = "";

  char oldClass[KI_NAME_BUF_SIZE] = "";
  char className[KI_NAME_BUF_SIZE] = "";
  char fastaName[KI_NAME_BUF_SIZE] = "";

  /* FILE* fa; */
  
  while (fgets(buf, sizeof(buf), fp) != NULL) {
    char *p, *q;

    for (p = buf; *p != '\0' && *p != '\n' &&  *p != '\t' && *p != ' ';  p++);  *p = '\0';
    for (p = p+1; *p != '\0' && *p != '\n' && (*p == '\t' || *p == ' '); p++);
    for (q = p;   *q != '\0' && *q != '\n' &&  *q != '\t' && *q != ' ';  q++);  *q = '\0';

    strcpy(className, buf);
    strcpy(fastaName, p);
    fprintf(stderr, "class=%s \t fa=%s\n", className, fastaName);
    
    if (strcmp(className, oldClass) != 0) {
      if (strlen(oldClass) > 0) {
        /* add old profile vector */
        kiFreqToRaiVectorOriginal(freq, vector);
        kiRaiDbAdd(db, oldClass, vector);
      }
      
      /* start new vector */
      kiClearKmerFreq(freq);
      strcpy(oldClass, className);
    }

    kiScanFastaForFreq(fastaName, freq);
  }

  kiFreqToRaiVectorOriginal(freq, vector);
  kiRaiDbAdd(db, className, vector);

  kiPrintRaiDb(fd, db);

  freq = kiFreeKmerFreq(freq);
  db = kiFreeRaiDb(db);

  fclose(fp);
  fclose(fd);
}

void userMain(int argc, char *argv[]) {
  char usage[] = "rait [-k 7] -o database -i otu.fasta.table \n";

  char* db       = NULL;
  char* fileName = NULL;
  int   kmerSize = 7;
  bool  newImpl  = false;
  
  int i;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-i") == 0) { 
      fileName = argv[++i];
    } else if (strcmp(argv[i], "-o") == 0) { 
      db = argv[++i];
    } else if (strcmp(argv[i], "-k") == 0) { 
      kmerSize = atoi(argv[++i]);
      if (kmerSize < 3 || kmerSize > 10) {
        fprintf(stderr, "K value out of range [3, 10] (%d); use 7 instead\n", kmerSize);
        kmerSize = 7;
      }
    } else if (strcmp(argv[i], "-new") == 0) { 
        newImpl = true;
    }
  }

  if (!fileName) {
    fprintf(stderr, "Usage: %s\n", usage);
    pmSetLevel(0);
    return;
  }
  

  char dbName[KI_NAME_BUF_SIZE];
  if (db) {
    strcpy(dbName, db);
  } else {
    char* slash = strrchr(fileName, '/');
    char* name  = slash ? slash+1 : fileName;
    sprintf(dbName, "%s.db", name);
  }


  if (newImpl) 
      trainFastaList(fileName, dbName, kmerSize);
  else 
      trainFastaListOriginal(fileName, dbName, kmerSize);
}

int main(int argc, char *argv[])
{
  kiInit(&argc, &argv);
  kiStart();

  if (kiIsUserDomain())
    userMain(argc, argv);
  
  kiStop();
  kiFinalize();
  return 0;
}
