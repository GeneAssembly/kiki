#include "extern.h"
#include <math.h>
#include "raiphy.h"


void classifySequential(char* fileName, rai_db_t* db) { /* online */

  FILE* fp = fopen(fileName, "r");
  char  buf[KI_SEQ_BUF_SIZE] = "";
  char *seqSkip  = " \t\r\n";
  int   nSaved = 10;
  char* seq = (char*)kimalloc(nSaved);
  double margin = 0.;
  int class;
  
  *seq = '\0';
  /* int nDone = 0; */  

  while(fgets(buf,sizeof(buf),fp) != NULL) {
    if (buf[0] == '>') {

      if (*seq != '\0') {
        class = classifySequence(seq, db, &margin);
        printf("%s\n", db->names[class]);
      }
    
      printf("%s", buf);
      if (nSaved > 0) *seq = '\0';

    } else {
      /* count non-space characters and append to sequence */
      int nKeep = 0;
      char *p, *q;
      for (p=buf; *p != '\0'; p++) {
        for (q=seqSkip; *q != '\0'; q++) {
          if (*p == *q)
            break;
        }
        if (*p != *q)
          nKeep++;
      }
      int nOld = (seq == NULL) ? 0 : strlen(seq);
      if (nOld + nKeep + 1 > nSaved) {
        int newSaved = nOld + nKeep + 1;
        seq = (char*)kirealloc(seq, nSaved, newSaved, /*copy*/false);
        /* seq = (char*)kirealloc(seq, nOld, nSaved, /\*copy*\/false); */
        nSaved = newSaved;
      }
      char *out = seq + nOld;
      for (p = buf; *p != '\0'; p++) {
        for (q = seqSkip; *q != '\0'; q++) {
          if (*p == *q)
            break;
        }
        if (*p != *q) {
          *out = *p;
          out++;
        }
      }
      assert(out - seq == nKeep + nOld);
      *out = '\0';
      
      /* kipm("nDone = %d    \r", ++nDone); */
      
    }
  }

  if (*seq != '\0') {
    class = classifySequence(seq, db, &margin);
    printf("%s\n", db->names[class]);
  }
  
  fclose(fp);
  kifree(seq, nSaved);

}

void raiTest(int argc, char *argv[]) {
  /* char usage[] = "rai -i database.to.train.fasta"; */

  char* fileName = NULL;
  int k = 7;
  int i;
  rai_db_t* db = NULL;
  
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-k") == 0) {
      if (i < argc) k = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-d") == 0) { 
      fileName = argv[++i];
      db = loadDatabase(fileName);
    } else if (strcmp(argv[i], "-i") == 0) { 
      while (i+1 < argc && argv[i+1][0] != '-') {
        fileName = argv[++i];
        classifySequential(fileName, db);
      }
    }
  }
  kiFreeRaiDb(db);
}

void raiParallel(int argc, char *argv[]) {
  char usage[] = "rai -d database -i fasta1 fasta2 ...  \n   or: rai -d database -I fasta.file.list\n";
  
  int nSeqs      = 0;
  int nTotalSeqs = 0;
  char* fileName = NULL;
  char* dbName   = NULL;

  int i;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-i") == 0) { 
      while (i+1 < argc && argv[i+1][0] != '-') {
        fileName = argv[++i];
        if (kiIsDomainRoot()) {
          kiUserReadFasta(fileName, &nSeqs);
        }
      }
    } else if (strcmp(argv[i], "-I") == 0) { 
      fileName = argv[++i];
      char faFileName[KI_NAME_BUF_SIZE];
      if (kiIsDomainRoot()) {
        FILE* fp = fopen(fileName, "r");
        if (fp != NULL) {
          while (fscanf(fp, "%s", faFileName) != EOF) {
            kipmsg(3, "faFileName = %s\n", faFileName);
            kiUserReadFasta(faFileName, &nSeqs);
          }
        }
        fclose(fp);
      }
    } else if (strcmp(argv[i], "-d") == 0 && i+1 < argc) {
      dbName = argv[++i];
    } 
  }
  if (nTotalSeqs > nSeqs && kiIsDomainRoot())
    printf("Total sequences read = %d\n", nTotalSeqs);

  if (!fileName || !dbName) {
    fprintf(stderr, "Usage: %s\n", usage);
    pmSetLevel(0);
    return;
  }

  char binName[200];
  char* slash = strrchr(fileName, '/');
  if (slash != NULL) fileName = slash + 1;
  sprintf(binName, "%s.bin", fileName);

  kiUserRAIphyOriginal(dbName, binName);

}

void userMain(int argc, char *argv[]) {

  kiInit(&argc, &argv);
  kiStart();

  if (kiIsUserDomain())
    raiParallel(argc, argv);
  
  kiStop();
  kiFinalize();
  
}

int main(int argc, char *argv[])
{
  int i;
  int old = 0;
  
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-serial") == 0) {
      old = 1;
    }
  }

  if (old)
    raiTest(argc, argv);
  else
    userMain(argc, argv);

  return 0;
}
