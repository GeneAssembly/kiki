#include "extern.h"
#include <math.h>

#include "raiphy.h"

int rai_base2int[117] = {
// 0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //   0-19
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //  20-39
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //  40-59
   0,0,0,0,0,0,0,2,0,0,0,1,0,0,0,0,0,0,0,0, //  60-79
   0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2, //  80-99
   0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,3        // 100-116 
};


kmer_freq_t* kiAllocKmerFreq(int kmerSize) {
  kmer_freq_t* freq = (kmer_freq_t*)kimalloc(sizeof(kmer_freq_t));
  freq->kmerSize = kmerSize;
  freq->count = (long long**)kimalloc(sizeof(long long*) * (kmerSize+1));
  int i = 0;
  for (i = 0; i < kmerSize + 1; ++i) {
    freq->count[i] = kimalloc(sizeof(long long) * (1 << (2*i))); /* 4^k */
    memset(freq->count[i], 0, sizeof(long long) * (1 << (2*i)));
  }
  return freq;
}

kmer_freq_t* kiFreeKmerFreq(kmer_freq_t* freq) {
  if (freq == NULL) return NULL;
  if (freq->count != NULL) {
    int i;
    for (i = 0; i < freq->kmerSize + 1; ++i) {
      freq->count[i] = kifree(freq->count[i], sizeof(long long) * (1 << (2*i))); /* 4^k */
    }
    freq->count = kifree(freq->count, sizeof(long long*) * (freq->kmerSize + 1));
  }
  return kifree(freq, sizeof(kmer_freq_t));
}

void kiClearKmerFreq(kmer_freq_t* freq) {
  int i = 0;
  for (i = 0; i < freq->kmerSize + 1; ++i) {
    memset(freq->count[i], 0, sizeof(long long) * (1 << (2*i)));
  }
}

void kiScanSeqsForFreq(alignment_t* seqs, kmer_freq_t* freq) {
  int i, j, k;
  int index1, index2;
  int x1, x2;
  char* p;
  int dim = 1 << (2*freq->kmerSize);
  int mask = dim - 1;

  for (i = 0; i < seqs->nSeq; ++i) {
    j = 1;
    index1 = index2 = 0;
    
    for (p = seqs->seqs[i]; *p != '\0'; ++p) {
      /* kipm("%c", *p); */
      index1 <<= 2; index1 |= rai_base2int[(int)*p]; index1 &= mask;
      index2 >>= 2; index2 |= (rai_base2int[base2complement[(int)*p]] << 12); 

      if (j < freq->kmerSize) { ++j; continue; }

      for (k = freq->kmerSize, x1 = index1, x2 = index2; k >= 0; k--, x1 >>= 2, x2 >>= 2) {
        freq->count[k][x1]++;
        freq->count[k][x2]++;
      }
    }
    /* kipm("\n"); */
  }
  /* kipm("freq->count[0][0] = %d\t", freq->count[0][0]); */
  /* kipm("freq->count[1][0] = %d\n", freq->count[1][0]); */
  
}

void kiScanFastaForFreq(char* fileName, kmer_freq_t* freq) {
  alignment_t* seqs;
  FILE* fa;
  if ((fa = fopen(fileName, "r")) == NULL) {
    fprintf(stderr, "Could not open fasta file %s.\n", fileName);
    kiAbort(1);
  }
  seqs = kiReadFasta(fa);
  fclose(fa);

  kipmsg(3, "Read %s contains %d sequences.\n", fileName, seqs->nSeq);
  kiScanSeqsForFreq(seqs, freq);

  int i;
  for (i = 0; i < seqs->nSeq; i++) {  /* in case Arena is in used and these do not get freed */
    seqs->names[i] = kifreestr(seqs->names[i]);
    seqs->seqs[i]  = kifreestr(seqs->seqs[i]);
  }
  seqs = kiFreeAlignment(seqs);

  /* double* v = kiFreqToRai(freq); */
  /* v = kifree(v, sizeof(double)*(1<<(2*freq->kmerSize))); */
}

void kiFreqToRaiVector(kmer_freq_t* freq, double* v) {
  int i, j, k = freq->kmerSize;
  int dim = 1 << (2*k);
  double r, fd, gd, dmin = 0.1;
  long long f, f1, g, g1;
  int mask[k];
  for (j = 1; j <= k - 2; ++j) {
    mask[j] = (1 << (j*2+2)) - 1;           /* 4 ^ (k-1) */
  }
  memset(v, 0, sizeof(double)*dim);
  for (i = 0; i < dim; ++i) {
    v[i] = 0.;
    for (j = 1; j <= k - 2; ++j) {
      r = 0.;
      f  = freq->count[k][i];                      /* f(x_1, x_2, ..., x_k)    */
      f1 = freq->count[k-1][i >> 2];               /* f(x_1, x_2, ..., x_{k-1} */
      g  = freq->count[j+1][i & mask[j]];          /* f(x_{k-j},  ..., x_k     */
      g1 = freq->count[j][(i & mask[j]) >> 2];     /* f(x_{k-j},  ..., x_{k-1} */
      if (f1 == 0 || g1 == 0) continue;
      fd = f > 0 ? (double)f : dmin;
      gd = g > 0 ? (double)g : dmin;
      /* kipm("%d, %d, %d, %d, %g, %g\n", f, f1, g, g1, fd, gd); */
      r = log((fd/(double)f1) / (gd/(double)g1)) / log(2.);
      /* kipm("  i=%d  rai_%d = %f\n", i, j, r); */
      v[i] += r;
      /* break; */
    }
    /* break; */
  }
}

rai_db_t* kiAllocRaiDb(int kmerSize) {
  rai_db_t* db = (rai_db_t*)kimalloc(sizeof(rai_db_t));
  db->nameLen  = 200;
  db->kmerSize = kmerSize;
  db->nDim     = 1 << (kmerSize * 2);
  db->nClass   = 0;
  db->nSaved   = 0;
  db->names    = NULL;
  db->vectors  = NULL;
  return db;
}

rai_db_t* kiFreeRaiDb(rai_db_t* db) {
  if (db == NULL)
    return NULL;
  
  int i;
  for (i = 0; i < db->nClass; ++i) {
    db->names[i]   = kifreestr(db->names[i]);
    db->vectors[i] = kifree(db->vectors[i], sizeof(double)*db->nDim);
  }
  
  db->names   = kifree(db->names, sizeof(char*) * db->nSaved);
  db->vectors = kifree(db->vectors, sizeof(double*) * db->nSaved);
  
  kifree(db, sizeof(rai_db_t));
  return NULL;
}

void kiRaiDbGrow(rai_db_t* db) {
  db->nClass++;
  
  if (db->nSaved == 0) {
    db->nSaved  = 100;
    db->names   = (char**)kimalloc(sizeof(char*) * db->nSaved);
    db->vectors = (double**)kimalloc(sizeof(double*) * db->nSaved);
  }
  
  if (db->nClass > db->nSaved) {
    int nNewSaved = db->nSaved * 2;
    db->names     = kirealloc(db->names,  sizeof(char*)*(db->nSaved),   sizeof(char*)*nNewSaved,   /*copy*/false);
    db->vectors   = kirealloc(db->vectors, sizeof(double*)*(db->nSaved), sizeof(double*)*nNewSaved, /*copy*/false);
    db->nSaved    = nNewSaved;
  }
}

void kiRaiDbAdd(rai_db_t* db, char* name, double* vector) { /* copy */
  int top = db->nClass;
  int dim = db->nDim;
  double* v = kimalloc(sizeof(double)*dim);
  memcpy(v, vector, sizeof(double)*dim);
  kiRaiDbGrow(db);
  db->names[top]  = kistrdup(name);
  db->vectors[top] = v;
}

void kiPrintRaiDb(FILE* fp, rai_db_t* db) {
  int i, j;
  fprintf(fp, "%sWORD_LENGTH:%d\n", "%", db->kmerSize); 
  for (i = 0; i < db->nClass; ++i) {
    fprintf(fp, "%s:", db->names[i]);
    for (j = 0; j < db->nDim; ++j) {
      double v = db->vectors[i][j];
      int round = (int)v;
      if (fabs(v-round) < 0.00001)
        fprintf(fp, "%.0f", v);
      else
        fprintf(fp, "%7.5f", v);
      if (j < db->nDim - 1)
        fprintf(fp, ":");
    }
    fprintf(fp, "\n");
  }
}

rai_db_t* loadDatabase(char* fileName) {

  char buf[1000000] = "";
  char *p, *q;
  FILE* fp = fopen(fileName, "r");
  if (fgets(buf, sizeof(buf), fp) == NULL) {
    fprintf(stderr, "Error reading first line\n");
  }
  int k = 7;
  int nDim = 1 << (k * 2);
  int i;
  char* name;
  double* vector;
  
  rai_db_t* db = kiAllocRaiDb(k);
  
  while(fgets(buf,sizeof(buf),fp) != NULL) {
    for (p = buf; *p != ':'; ++p); *p = '\0';
    kiRaiDbGrow(db);
    name = kistrdup(buf);
    vector = (double*)kimalloc(sizeof(double)*nDim);
    /* kipm("name = %s\n", name); */
    for (i = 0; i < nDim; ++i) {
      p++;
      for (q = p; *q != '\0' && *q != ':'; q++); *q = '\0';
      vector[i] = atof(p);
      p = q;
    }

    db->vectors[db->nClass-1] = vector;
    db->names[db->nClass-1]   = name;
    /* printf("name = %s\n", db->names[db->nClass-1]); */
    /* printf("v[5] = %.3f\n", db->vectors[db->nClass-1][16000]); */
    
  }
  
  return db;
}

int classifySequence(char* seq, rai_db_t* db, double* margin) {
  return 0;
}


/* Original classification and training algorithms implemented in the RAIphy source codes */

int classifySequenceOriginal(char* seq, rai_db_t* db, double* margin) {
  
  int k = db->kmerSize;
  int dim = db->nDim;
  int mask = dim - 1;

  int v[dim];
  /* int nz[dim]; */
  /* double nv[dim]; */
  /* double dv[dim]; */

  memset(v, 0, sizeof(int)*dim);
  /* memset(nz, 0, sizeof(int)*dim); */
  /* memset(nv, 0.0, sizeof(double)*dim); */
  /* memset(dv, 0.0, sizeof(double)*dim); */
  
  char *p;
  int i, j, index1, index2;

  j = 1;
  index1 = 0;
  index2 = 0;
    
  for (p = seq; *p != '\0'; ++p) {
    index1 <<= 2; index1 |= rai_base2int[(int)*p]; index1 &= mask;
    index2 >>= 2; index2 |= (rai_base2int[base2complement[(int)*p]] << 12); 

    if (j < k) { ++j; continue; }

    v[index1]++;
    v[index2]++;
  }

  int nz[dim];
  int nzi = 0;
  memset(nz, 0, sizeof(int)*dim);
  
  for (j = 0; j < dim; j++) {
    if (v[j] > 0) {
      nz[nzi++] = j;
    }
  }
  nz[nzi] = -1;
  
  /* RAI replicate code */
  /* int tmp; */
  /* for (i = 0; i < dim; i+=4) { */
  /*   tmp = 0; */
  /*   for (j=0; j<4; j++) tmp += v[i+j]; */
  /*   nv[i/4]=tmp; */

  /*   if (tmp != 0) */
  /*   { */
  /*     for (j=0; j<4; j++) */
  /*       dv[i+j] = 1.0*v[i+j] / tmp;                  // Normalize first v */
  /*   } */
  /* } */

  /* for (i=0; i < dim>>2; i+=4) */
  /* { */
  /*   tmp = 0; */
  /*   for (j=0; j<4; j++) */
  /*     tmp += nv[i+j]; */

  /*   if (tmp != 0)                                             // Normalize second v */
  /*   { */
  /*     for (j=0; j<4; j++) */
  /*       nv[i+j] = nv[i+j] / tmp; */
  /*   } */
  /* } */

  /* for (i=0; i<dim; i+=4) */
  /* { */
  /*   for (j=0; j<4; j++) */
  /*   { */
  /*     if (dv[i+j] != 0 && nv[i/4] != 0) { */
  /*       nz[i+j] = 1; */
  /*       dv[i+j] = 11.0*log(dv[i+j]) + log(nv[i/4]); */
  /*     } */
  /*     else  */
  /*       dv[i+j] = -20.0; */
  /*   } */
  /* } */
  /* end of RAI replicate code */
  

  int    best_index   = -1;
  double second_score = -1e10;
  double best_score   = -1e10;
  double score;
  
  for (i = 0; i < db->nClass; i++) {
    score = 0.;
    for (nzi = 0; (j = nz[nzi]) >= 0; nzi++) {
    /* for (j = 0; j < dim; j++) { */
      /* if (v[j] > 0) { */
        /* score += dv[j] * db->vectors[i][j]; */
        score += v[j] * db->vectors[i][j];
        /* kipm("v = %d, db = %.3f\n", v[j], db->vectors[i][j]); */
      /* } */
    }
    /* printf("db->vectors[i][0] = %.3f\n", db->vectors[i][j]); */
    /* kipm("score = %.3f, best_score = %.3f\n", score, best_score); */
    
    if (score > best_score) {
      second_score = best_score;
      best_score = score;
      best_index = i;
    }
  }

  if (margin != NULL) *margin = (best_score - second_score);

  /* kipm("best_score = %.3f\n", best_score); */
  /* kipm("best_index = %.3f\n", best_index); */
  return best_index;
}


void kiFreqToRaiVectorOriginal(kmer_freq_t* freq, double* vector) {
  
  int k   = freq->kmerSize;
  int dim = 1 << (2*k);

  int i;
  long long c, c1, c2;
  
  memset(vector, 0, sizeof(double)*dim);
  
  for (i = 0; i < dim; ++i) {
    c  = freq->count[k][i];
    c1 = freq->count[k-1][i>>2];
    c2 = freq->count[k-2][i>>4];
    
    if (c != 0 && c1 != 0 && c2 != 0)
      vector[i] = 11.0 * log(1.*c/c1) + log(1.*c1/c2);
    else 
      vector[i] = -20.0;
    
  }
}


void trainFromFastaOriginal(char* fileName, int kmerSize) {
  kipm("training %s...\n", fileName);

  int k = kmerSize;             /* default: 7 */
  int dim = 1 << 14;
  int mask = dim - 1;

  double dbv[dim];
  memset(dbv, 0.0, sizeof(double)*dim);

  /* return; */
  
  FILE* fp = fopen(fileName, "r");
  alignment_t* seqs = kiReadFasta(fp);
  
  kipm("Read %d sequences.\n", seqs->nSeq);

  int kmerCounts[dim], v[dim];
  double nv[dim];
  double dv[dim];
  memset(kmerCounts, 0, sizeof(int)*dim);
  memset(v, 0, sizeof(int)*dim);
  memset(nv, 0.0, sizeof(double)*dim);
  memset(dv, 0.0, sizeof(double)*dim);
  
  char *p;
  int i, j, index1, index2;
  for (i = 0; i < seqs->nSeq; ++i) {
    j = 1;
    index1 = 0;
    index2 = 0;
    
    for (p = seqs->seqs[i]; *p != '\0'; ++p) {
      /* kipm("%c", *p); */
      index1 <<= 2; index1 |= rai_base2int[(int)*p]; index1 &= mask;
      index2 >>= 2; index2 |= (rai_base2int[base2complement[(int)*p]] << 12); 

      if (j < k) { ++j; continue; }

      /* kmerCounts[index]++; */
      /* kmerCounts[mask^index]++; */
      v[index1]++;
      v[index2]++;
      
    }
    /* kipm("\n"); */
  }

  int tmp;
  
  for (i = 0; i < dim; i+=4) {
    tmp = 0;
    for (j=0; j<4; j++) tmp += v[i+j];
    nv[i/4]=tmp;

    if (tmp != 0)
    {
      for (j=0; j<4; j++)
        dv[i+j] = 1.0*v[i+j] / tmp;                  // Normalize first v
    }
  }

  for (j = 0; j < 10; ++j) {
    /* kipm("%.3f:", dv[j]); */
    /* kipm("%.3f:", nv[j]); */
  }
  /* kipm("\n"); */
  /* for (i=0; i < (1 << ((7-1)*2)); i+=4) */
  for (i=0; i < dim>>2; i+=4)
  {
    tmp = 0;
    for (j=0; j<4; j++)
      tmp += nv[i+j];

    if (tmp != 0)                                             // Normalize second v
    {
      for (j=0; j<4; j++)
        nv[i+j] = nv[i+j] / tmp;
    }
  }

  for (i=0; i<dim; i+=4)
  {
    for (j=0; j<4; j++)
    {
      if (dv[i+j] != 0 && nv[i/4] != 0)
        dv[i+j] = 11.0*log(dv[i+j]) + log(nv[i/4]);
      else
        dv[i+j] = -20.0;
    }
  }


  for (j = 0; j < 10; ++j) {
    int round = (int)dv[j];
    if (fabs(dv[j]-round) < 0.00001)
      kipm("%.0f:", dv[j])
    else
      kipm("%7.5f:", dv[j]);
  }
  kipm("\n");
}
