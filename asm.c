
#include "extern.h"

int global_start_time = 0;
void startTimer() {
  global_start_time = (int)MPI_Wtime();
}

bool stopTimer(int seconds) {
  
  return ((int)MPI_Wtime() - global_start_time >= seconds);
}

double global_time = 0.;
void checkTime(char* info) {
  if (global_time < 0.1) {
    global_time = MPI_Wtime();
    return;
  }
  double time = MPI_Wtime();
  /* double elapsed = time - global_time; */
  global_time = time;

  /* C timer */
  char buffer[30];
  getCurrentTime(buffer);
  
  /* struct timeval tv; */
  /* time_t curtime; */
  /* gettimeofday(&tv, NULL); */
  /* curtime=tv.tv_sec; */
  /* strftime(buffer, 30, "%m-%d-%Y %T.", localtime(&curtime)); */
  
  /* kipm("Time: %s (%.2f s elapsed) (MPI clock = %.2f s) (system clock = %s%ld)\n", info, elapsed, MPI_Wtime() - global_start_time, buffer, tv.tv_usec); */
  /* kipm("Time: %s (%.2f s elapsed) (MPI clock = %.2f s) (system clock = %s)\n", info, elapsed, MPI_Wtime() - global_start_time, buffer); */
}

void hybridAssemble(int argc, char *argv[]) {
  startTimer();

  if (argc < 2) {
    fprintf(stderr, "No input fasta files specified.\n");
    exit(1);
  }
  
  int kmerLen = 25;
  int i;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-k") == 0) { 
      ++i;
      if ((i < argc) && argv[i][0] >= '0' && argv[i][0] <= '9') {
        kmerLen = atoi(argv[i]);
      } else --i;
      kipmsg0(2, "kmerLen = %d\n", kmerLen);
    }
  }
  if (kiIsDomainRoot()) {
    kiUserSetK(kmerLen);
    kiUserSetHashMode(KI_HASH_ENDS);
  }
  int minOverlap = kmerLen;

  int optDot = 0;               /* dot option: 0: none, 1: small, 2: big */
  int optProgress = 0;
  int optPersist = 0;
  /* int optMinReadLen = 0;  */
  int optMinReadLen = 75; 
  float optDepletion = -1.0;
  /* int abruptStop = 0; */
  
  /* Load sequences from multiple files */
  char inputFiles[200][200];
  int nInputFiles = 0;
  long long nSeqs = 0, nTotalSeqs = 0;
  char* fileName = NULL;
  char* outputName = NULL;

  for (i = 1; i < argc; ++i) {

    if (strcmp(argv[i], "-h") == 0) {
      kipm0("Usage: ki (-i single.file | -I file.list) [options] \n");
      kipm0("\t -k int        kmer length, also used as minimum overlap (D = 25) \n");
      kipm0("\t -o outname    output prefix \n");
      kipm0("\t -u float|int  fraction of nodes used as users or the number of user nodes \n");
      kipm0("\t -v [level]    verbosity level (D = 3)  \n");
      kipm0("\t -persist [t]  checkpointing interval (eg, 30m, 6h, D = 1h) \n");
      kipm0("\t -dep float    early termination when a fraction of reads are depleted (D = 1.0) \n");
      kipm0("\t -dot file     generate graphviz output \n");
      kipm0("\t -h            show this usage information \n\n");

      return;
      
    } else if (strcmp(argv[i], "-i") == 0) { 
      while (i+1 < argc && argv[i+1][0] != '-') {
        fileName = argv[++i];
        strcpy(inputFiles[nInputFiles++], fileName);
      }
    } else if (strcmp(argv[i], "-I") == 0) { 
      fileName = argv[++i];
      char faFileName[KI_NAME_BUF_SIZE];
      if (kiIsDomainRoot()) {
        FILE* fp = fopen(fileName, "r");
        if (fp != NULL) {
          while (fscanf(fp, "%s", faFileName) != EOF) {
            kipmsg(3, "faFileName = %s\n", faFileName);
            strcpy(inputFiles[nInputFiles++], faFileName);
          }
        }
        fclose(fp);
      } 
    } else if (strcmp(argv[i], "-o") == 0) { 
      outputName = argv[++i];
    } else if (strcmp(argv[i], "-dep") == 0) { 
      if (i+1 < argc && argv[i+1][0] != '-') {        
        optDepletion = atof(argv[++i]);
        if (optDepletion <= 0 || optDepletion >= 1) {
          optDepletion = -1.0;
        } else {
          kipmsg0(2, "optDepletion = %.3f\n", optDepletion);
        }
      }
    } else if (strcmp(argv[i], "-dot") == 0) {
      if (i+1 < argc && argv[i+1][0] != '-')
        optDot = atoi(argv[++i]);
      else 
        optDot = 1;
    } else if (strcmp(argv[i], "-split") == 0) { 
      optMinReadLen = 75;
      if (i+1 < argc && argv[i+1][0] != '-') {
        char* lenStr = argv[++i];
        optMinReadLen = atoi(lenStr);
      }
    } else if (strcmp(argv[i], "-progress") == 0) {
      optProgress = 1;
    } else if (strcmp(argv[i], "-persist") == 0) {
      optPersist = 50 * 60;      
      if (i+1 < argc && argv[i+1][0] != '-') {
        char* timeStr = argv[++i];
        optPersist = atoi(timeStr);
        int unit = 60;
        switch (timeStr[strlen(timeStr)-1]) {
          case 's':
          case 'S': unit = 1;     break;
          case 'h':
          case 'H': unit = 3600;  break;
          case 'd':
          case 'D': unit = 86400; break;
          default:  unit = 60;
        }
        optPersist *= unit;
      }
      kipmsg0(2, "optPersist = %d\n", optPersist);
    }
    kiBarrier();
  }

  /* read fasta/fastq files */
  if (kiIsDomainRoot()) {
    if (nInputFiles > 0) {
      for (i = 0; i < nInputFiles; ++i) {
        kiUserReadFastaOrFastq(inputFiles[i], &nSeqs);
        printf("Read %lld sequences from '%s'.\n", nSeqs, inputFiles[i]);
        nTotalSeqs += nSeqs;
      }
    } else {
      fprintf(stderr, "Please use -i or -I to specify input files.\n");
      kiAbort(-1);
    }
  }
  
  if (nTotalSeqs > nSeqs && kiIsDomainRoot())
    printf("Total sequences = %lld\n", nTotalSeqs);

  if (kiIsDomainRoot()) {
    int minLen, maxLen;
    kiUserGetMinMaxReadLen(&minLen, &maxLen);
    kipmsg(2, "Read length range = [%d, %d]\n", minLen, maxLen);

    if (minLen < maxLen && optMinReadLen > 0) {
      long long nSeqsLeft;
      kiUserProcessVarLenReads(optMinReadLen, minLen, maxLen, minOverlap, &nSeqsLeft);
      kipmsg(2, "Total sequences after length processing = %ld\n", nSeqsLeft);
    }
  }

  if (optPersist > 0) {
    int timeLeft = optPersist - ((int)MPI_Wtime() - global_start_time);
    if (timeLeft <= 1) {       
      fprintf(stderr, "Used up specified execution time before the start of assembly; increase \"-persist\" time!\n");
      kiAbort(-1);
    }
  }
  
  if (kiIsDomainRoot()) {
    kiUserSetTermination(optPersist, optDepletion);
  }

  kiBarrier();

  char contigFile[KI_NAME_BUF_SIZE];
  char outputKey[KI_NAME_BUF_SIZE];

  if (outputName != NULL) {
    strcpy(outputKey, outputName);
  } else {
    /* char* slash = strrchr(fileName, '/'); */
    /* char* name  = slash ? slash+1 : fileName; */
    sprintf(outputKey, "%s", fileName);
  }
  
  if (ki_domain_size > 1) {
    sprintf(contigFile, "%s.contig.%d", outputKey, ki_domain_rank);
  }
  else {
    sprintf(contigFile, "%s.contig", outputKey);
  }
  
  if (optPersist > 0) {
    KI_File_backup(contigFile);
  }
    
  /* File name prefix for graphviz */
  char dotPrefix[200];
  char* dot = strrchr(outputKey, '.');
  if (dot != NULL) *dot = '\0';

  long long nProcessed;
  if (optPersist > 0) {    
    if (kiIsDomainRoot()) {
      kiUserRestoreState(outputKey, &nProcessed);
      kipm("Restored nProcessed = %ld\n", nProcessed);
    }
    KI_Bcast(&nProcessed, 1, MPI_LONG_LONG, 0, ki_cmm_domain);
  }

  FILE *fp;
  char fmode[2] = "w";
  if (optPersist > 0 && nProcessed > 0) {
    fmode[0] = 'a';
    fmode[1] = '+';
  }

  fp = fopen(contigFile, fmode);
  if (fp == NULL) {
    kipmsg(0, "Could not open output file '%s'\n", contigFile);
    kiAbort(-1);
  }

  /* return; */
  if (kiIsDomainRoot()) kiUserLoadReadSeqs();

  char query[KI_SEQ_BUF_SIZE]; 
  char seed[KI_SEQ_BUF_SIZE];
  char seedName[KI_SEQ_BUF_SIZE];
  int  iSeed;
  int  iContig = 0;
  int  d;

  char* contig = kimalloc(KI_CONTIG_SIZE);
  char* elongation[2];
  for (i = 0; i < 2; ++i) {
    elongation[i] = kimalloc(KI_CONTIG_SIZE);
  }

  bool bErase = false;
  
  alignment_t*     matches = kiAllocAlignment();
  overlap_graph_t* graph   = kiAllocOverlapGraph(minOverlap);

  /* int seedCount = 0; */
  int termStatus = 0;
  
  kiBarrier();

  kiUserGetSeedSeq(seed, &termStatus);
  
  while (termStatus == KI_TERMINATION_NONE && seed[0] !='\0') {
    
    kipmsg(4, "seed = '%s'\n", seed);
    /* if (kiIsDomainRoot() && seedCount++ % 10000 == 0) { */
      /* char ithSeed[100]; */
      /* sprintf(ithSeed, "New seed %d", ++seedCount); */
      /* checkTime(ithSeed); */
    /* } */

    for (d = 0; d < 2; ++d) {
      elongation[d][0] = '\0';

      if (d == 0) {             /* extending forward */
        strcpy(query, seed);

        kiUserGetOverlappingSeqs(query, minOverlap, .0, /*erase*/bErase, /*OUT*/matches);
        kipmsg(4, "Overlapping: matches->nSeq = %d\n", matches->nSeq);
        if (!kiAlignmentHasZeroOffset(matches)) { 
          /* Abort: between the time we got the seed and now,
             some other user already used that seq */
          kiClearAlignment(matches);
          kipmsg(5, "No base match, preempted\n");
          break;
        }
        kiOverlapGraphAddMatches(graph, matches, 0, minOverlap);
        kiOverlapGraphReviseSeed(graph, query, /*OUT*/&iSeed);
        if (kiOgvUnextended(graph, iSeed)) {   /* in case we need to start from a different seed */
          kiUserGetOverlappingSeqs(query, minOverlap, .0, /*erase*/bErase, /*OUT*/matches);
          kiOverlapGraphAddMatches(graph, matches, 0, minOverlap);
          kiOverlapGraphReviseSeed(graph, query, /*OUT*/&iSeed);
        }
        strcpy(seedName, graph->vertices[iSeed]->name);
        
      } else {                /* reverse direction */
        kiReverseComplementSeq(seed, query);
        kipmsg(5, "reverse seed = %s\n", query);

        kiAlignmentAdd(matches, seedName, query, true);
        kiOverlapGraphAddMatches(graph, matches, 0, minOverlap);
        kiUserGetOverlappingSeqs(query, minOverlap, .0, /*erase*/bErase, /*OUT*/matches);
        kiOverlapGraphAddMatches(graph, matches, 0, minOverlap);
        kiOverlapGraphReviseSeed(graph, query, /*OUT*/&iSeed);
      }
      
      sprintf(dotPrefix, "%s-%s-%c", outputKey, graph->vertices[iSeed]->name, (d > 0 ? 'b' : 'f'));

      char adventure[KI_SEQ_BUF_SIZE];
      int baseOffset;
      int iSolid = iSeed;
      bool progress = true;
      int step = 0;
      bool dotBig = (optDot > 1) ? true : false;
      if (optDot > 0) kiOverlapGraphExportDot(graph, dotPrefix, step, "-0", dotBig);

      while (progress) {
        char dstr[20];
        sprintf(dstr, "progress d=%d", d);
        checkTime(dstr);

        kipmsg(4, "iSolid = %d\n", iSolid);
        kipmsg(4, "v = %s\n", graph->vertices[iSolid]->name);
        kiOverlapGraphExplore(graph, iSolid, adventure);
        if (optDot > 0) kiOverlapGraphExportDot(graph, dotPrefix, step, "-1", dotBig);

        if (adventure[0] == '\0') break;
        kipmsg(4, "adventure = '%s'\n", adventure);

        if (kiOverlapGraphTooBig(graph)) {
          
          /* kipm("Oops\n"); */
          /* reap current elongation */
          kiOverlapGraphExtend(graph, &iSolid);
          kiOverlapGraphAppendContig(graph, iSeed, elongation[d]);

          if (graph->maxTailOffset >= KI_CONTIG_SIZE / 2 - 1000) {
            kiClearAlignment(matches);
            kiClearOverlapGraph(graph);
            break;
          }
          
          /* slide window in graph */
          kiClearAlignment(matches);
          ogvertex_t* u = graph->vertices[iSolid];
          int tmp_occr = u->occr;
          strcpy(query, u->seq);
          kiClearOverlapGraph(graph);

          ogvertex_t* v = kiAllocOgVertex();
          v->offset = 0;
          v->occr = tmp_occr;
          v->seq  = kiArenaStrdup(KI_ARENA_DEFAULT, query);; 
          v->flag = KI_OGV_SOLID;
          iSolid = 0;

          kiOverlapGraphAppendVertex(graph, v);

          kiUserGetOverlappingSeqs(v->seq + 1, minOverlap, .0, /*erase*/bErase, /*OUT*/matches);
          kiOverlapGraphAddMatches(graph, matches, 1, minOverlap);
          
          kiOverlapGraphExplore(graph, iSolid, adventure);
        }

        strcpy(query, graph->vertices[iSolid]->seq + strlen(graph->vertices[iSolid]->seq) - minOverlap + 1);
        strcpy(query+strlen(query), adventure);
        kiUserGetPrefixSeqs(query, minOverlap, .0, /*erase*/bErase, /*OUT*/matches);
        kipmsg(4, "Prefix: matches->nSeq = %d\n", matches->nSeq);
        baseOffset = graph->vertices[iSolid]->offset + strlen(graph->vertices[iSolid]->seq) - minOverlap + 1;
        kipmsg(4, "baseOffset = %d\n", baseOffset);
        kiOverlapGraphAddMatches(graph, matches, baseOffset, minOverlap);
        if (optDot > 0) kiOverlapGraphExportDot(graph, dotPrefix, step, "-2", dotBig);
        
        progress = kiOverlapGraphExtend(graph, &iSolid);
        if (optDot > 0) kiOverlapGraphExportDot(graph, dotPrefix, step, "-3", dotBig);

        if (optProgress != 0) {
          if (progress) printf("%c", d > 0 ? '+' : '-');
          else printf("\n");
          /* fflush(stdout); */
        }
        
        step++;
      }
        
      kiOverlapGraphAppendContig(graph, iSeed, elongation[d]);
      kipmsg(4, "elongation %s %c = '%s'\n", seedName, d > 0 ? '-' : '+', elongation[d]);
    
      if (optDot > 0) kiOverlapGraphExportDot(graph, dotPrefix, step, NULL, dotBig);

      kiClearAlignment(matches);
      kiClearOverlapGraph(graph);

    }

    if (d > 0) {                /* not preempted */
      kiCombineElongation(seed, elongation[0], elongation[1], contig);
      kipmsg(4, "contig '%s' (%d) =  %s\n", seedName, strlen(contig), contig);

      int seqlen = strlen(contig);
      if (seqlen > 125) {
        fprintf(fp, ">Contig_%d_%d Len_%d_Seed_%s \n%s\n", ki_domain_rank, iContig++, seqlen, seedName, contig);
      }
    }
    
    kiUserGetSeedSeq(seed, &termStatus);
  }
  kipmsg(3, "termStatus = %d\n", termStatus);

  kiFreeAlignment(matches);
  kiFreeOverlapGraph(graph);

  kifree(contig, KI_CONTIG_SIZE);
  for (i = 0; i < 2; ++i) {
    kifree(elongation[i], KI_CONTIG_SIZE);
  }
  
  fclose(fp);

  if (optPersist > 0) {
    if (kiIsDomainRoot()) {
      if (termStatus == KI_TERMINATION_NONE) {
        char storeFile[KI_NAME_BUF_SIZE];  
        sprintf(storeFile, "%s.store", outputKey);
        KI_File_backup(storeFile);
        KI_File_remove(storeFile);
      } else {
        long long nProcessed;
        kiUserStoreState(outputKey, &nProcessed);
        if (nProcessed > 0) kipm("Stored nProcessed = %ld\n", nProcessed);
      }
    }
    kiBarrier();
  }
}

void testLoadFasta(int argc, char *argv[]) {
  if (!kiIsDomainRoot()) return;
  
  int kmerLen = 13;
  kiUserSetK(kmerLen);
  kiUserSetHashMode(KI_HASH_ENDS);

  long long nSeqs = 0, nTotalSeqs = 0;
  char* fileName = NULL;
  int i;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-i") == 0) { 
      ++i;
      while (i < argc && argv[i][0] != '-') {
        fileName = argv[i++];
        kiUserLoadFasta(fileName, &nSeqs);
        printf("Loaded %lld sequences from '%s'.\n", nSeqs, fileName);
        nTotalSeqs += nSeqs;
      }
    }
  }
  if (i > 2)
    printf("Total sequences loaded = %lld\n", nTotalSeqs);

  return;
  

  char kmer[KI_SEQ_BUF_SIZE] = "TATTTACCAACTA";
  int nHits = -1;
  kiUserGetKmer(kmer, &nHits);
  kipm("before nHits = %d\n", nHits);
  
  char seed[KI_SEQ_BUF_SIZE] = "TAATTCACTCTATTGAAGATATTTACCAACTACCG";

  alignment_t* matches = kiAllocAlignment();
  
  kiUserGetOverlappingSeqs(seed, kmerLen, .0, /*erase*/true, /*OUT*/matches);
  kipmsg(2, "matches->nSeq = %d\n", matches->nSeq);

  kiFreeAlignment(matches);

  nHits = -1;
  kiUserGetKmer(kmer, &nHits);
  kipm("before nHits = %d\n", nHits);

  /* test arena */
  char* smallStrings[100];
  for (i = 0; i < 100; ++i) {
    smallStrings[i] = (char*)kiArenaMalloc(KI_ARENA_USER, i*10 + 1);
  }
}


void testHash(int argc, char *argv[]) {
  char kmer1[] = "TGGTTTTCTTTCGGG";
  char kmer2[] = "TGCAACCGACCCGCT";
  char kmer3[] = "ATCAAAAAAAAAAAA";
  /* char kmer3[] = "ATTGCAACCGACCCG"; */
  char* strings[3];
  strings[0] = kmer1;
  strings[1] = kmer2;
  strings[2] = kmer3;

  hashtable_t* hash = kiMakeHashtable(strings, 3, 16, false, false);
  int i;
  for (i = 0; i < 3; ++i) {
    kiHashtableAddCopy(hash, strings[i], strlen(strings[i]), i+3, 0);
  }
  hash = kiFreeHashtable(hash);
  /* pm("nChains = %d\n", hash->nChains); */
}

void userMain(int argc, char *argv[]) {
  /* testLoadFasta(argc, argv); */
  /* testHash(argc, argv); */
  hybridAssemble(argc, argv);
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
