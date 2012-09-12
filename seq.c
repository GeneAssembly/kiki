
#include "extern.h"
#include "math.h"

#include "seq.h"

/* Sequence manipulation */
void kiUpcaseSeq(char* seq) {
  char* p;
  for (p = seq; *p !=  '\0'; p++) {
    switch (*p) {
      case 'A': break;
      case 'T': break;
      case 'G': break;
      case 'C': break;
      case 'a': *p = 'A'; break;
      case 't': *p = 'T'; break;
      case 'g': *p = 'G'; break;
      case 'c': *p = 'C'; break;
      default:  *p = 'N';
    }
  }
}

void kiTimeToString(double time, /*OUT*/char* str) {
  if (time > 180.) {
    sprintf(str, "%.2f m", time / 60.);
  } else {
    sprintf(str, "%.2f s", time);
  }
}


/* Sequence manipulation */
void kiReverseComplementSeq(char* seq, /*OUT*/char* comp) {
  int len = strlen(seq);
  char *p, *q;
  for (p = seq, q = comp+len-1; *p != '\0'; ++p, --q) 
    *q = kiComplementBase(*p);
  comp[len] = '\0';
}


void kiReverseComplementSeqN(char* seq, int n, /*OUT*/char* comp) {
  int len = MIN(strlen(seq), n);
  char *p, *q;
  for (p = seq, q = comp+len-1; *p != '\0'; ++p, --q) 
    *q = kiComplementBase(*p);
  comp[len] = '\0';
}

void kiPartialReverseComplementSeq(char* seq, int beg, /*OUT*/char* comp) {
  /* beg: index of the beginning char on the complement seq */
  int len = strlen(seq) - beg;
  char *p, *q;
  for (p = seq, q = comp+len-1; q >= comp; ++p, --q) 
    *q = kiComplementBase(*p);
  comp[len] = '\0';
}

char kiComplementBaseSlow(char c) {
  if (c == 'a' || c == 'A') return 'T';
  if (c == 't' || c == 'T') return 'A';
  if (c == 'g' || c == 'G') return 'C';
  if (c == 'c' || c == 'C') return 'G';
  return 'N';
}

int base2complement[117] = {
// 0,1,2, 3, 4, 5, 6, 7, 8,9,0, 1, 2,3,4,5, 6, 7, 8, 9
   0,0,0, 0, 0, 0, 0, 0, 0,0,0, 0, 0,0,0,0, 0, 0, 0, 0,  //   0-19
   0,0,0, 0, 0, 0, 0, 0, 0,0,0, 0, 0,0,0,0, 0, 0, 0, 0,  //  20-39
   0,0,0, 0, 0, 0, 0, 0, 0,0,0, 0, 0,0,0,0, 0, 0, 0, 0,  //  40-59
   0,0,0, 0, 0,'T',0,'G',0,0,0,'C',0,0,0,0, 0, 0, 0, 0,  //  60-79
   0,0,0, 0,'A',0, 0, 0, 0,0,0, 0, 0,0,0,0, 0,'t',0,'g', //  80-99
   0,0,0,'c',0, 0, 0, 0, 0,0,0, 0, 0,0,0,0,'a'           // 100-116
};
inline char kiComplementBase(char c) {
  return base2complement[(int)c];
}

int base2int[117] = {
// 0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //   0-19 
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //  20-39
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //  40-59
   0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0, //  60-79
   0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1, //  80-99
   0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3        // 100-116 
};
inline int kiBaseIndex(char c) {
  return base2int[(int)c];
}

/* int aa2int[] = { */
/* // 0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9 */
/*    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, */
/*    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, */
/*    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, */
/* }; */

char dna2aa[64] = {             /* standard genetic code */
  'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',
  'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
  'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',
  'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
  'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',
  'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
  '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S',
  '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
};

char tmpCompSeq[KI_CONTIG_SIZE];
inline void kiDna2Protein(char* dna, char* aa, int frame) {
  
  if (frame > 3) frame = -(frame-3);
  if (frame < 0) {
    kiReverseComplementSeq(dna, tmpCompSeq);
    dna = tmpCompSeq;
    frame = -frame;
  }
  
  int i = frame - 1, len = strlen(dna);
  char *p, *q;
  for (p = aa, q = dna+i; i + 3 <= len; i += 3, ++p, q += 3) {
    /* kipm("%d%d%d -> %d\n", base2int[*q], base2int[*(q+1)], base2int[*(q+2)], base2int[*q] * 16 + base2int[*(q+1)] * 4 + base2int[*(q+2)]); */
    if (*q == 'N' || *(q+1) == 'N' || *(q+2) == 'N') break;
    *p = dna2aa[ base2int[(int)*q] * 16 + base2int[(int)*(q+1)] * 4 + base2int[(int)*(q+2)] ];
  }
  *p = '\0';
}

/* Sequence discrepancy rate */
float kiSeqCmp(char* s1, char* s2) {
  assert(s1 && s2);
  int len = 0, nDiff = 0;
  char *p, *q;
  for (p = s1, q = s2; *p != '\0' && *q != '\0'; ++p, ++q) {
    if (*p != *q) ++nDiff;
    ++len;
  }
  return (1.0 * nDiff/len);
}

float kiSeqNCmp(char* s1, char* s2, int n) {
  assert(s1 && s2);
  int len = 0, nDiff = 0;
  char *p, *q;
  for (p = s1, q = s2; *p != '\0' && *q != '\0' && len < n; ++p, ++q) {
    if (*p != *q) ++nDiff;
    ++len;
  }
  return (1.0 * nDiff/len);
}

float kiSeqNCmpX(char* s1, char* s2, int n, /*OUT*/int* iDiff) {
  assert(s1 && s2 && n>0 && iDiff);
  int len = 0, nDiff = 0;
  char *p, *q;
  *iDiff = 0;
  for (p = s1, q = s2; *p != '\0' && *q != '\0' && len < n; ++p, ++q) {
    if (*p != *q) {
      ++nDiff;
      if (*iDiff == 0) *iDiff = len;
    }
    ++len;
  }
  return (1.0 * nDiff/len);
}

bool kiIsProtein(char* s) {
  for (; *s; ++s) {
    if (*s == 'B' ||
        (*s > 'C' && *s < 'G') ||
        (*s > 'G' && *s < 'N') ||
        (*s > 'N' && *s < 'T') ||
        (*s > 'T' && *s <= 'Z')) return true;
  }
  return false;
}

float kiCalcGCContent(char* dna) {
  int ngc = 0;
  int n = 0;
  char* p;
  for (p = dna; *p != '\0'; ++p, ++n) {
    if (*p == 'G' || *p == 'C' || *p == 'g' || *p == 'c')
      ++ngc;
  }
  return 100.*ngc/n;            /* in percentage */
}


/* Alignment */

char blosum62indices[24] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'};

int blosum62compact[24][24] = {
  /* #  Matrix made by matblas from blosum62.iij */
  /* #  * column uses minimum score */
  /* #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units */
  /* #  Blocks Database = /data/blocks_5.0/blocks.dat */
  /* #  Cluster Percentage: >= 62 */
  /* #  Entropy =   0.6979, Expected =  -0.5209 */
  /*       A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
  /*A*/ {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4}, 
  /*R*/ { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
  /*N*/ { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
  /*D*/ { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
  /*C*/ {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
  /*Q*/ { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
  /*E*/ { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
  /*G*/ {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
  /*H*/ { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
  /*I*/ { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
  /*L*/ { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
  /*K*/ { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
  /*M*/ { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
  /*F*/ { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
  /*P*/ { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
  /*S*/ {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
  /*T*/ {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
  /*W*/ { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
  /*Y*/ { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
  /*V*/ {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
  /*B*/ { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
  /*Z*/ { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
  /*X*/ {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
  /***/ { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};

int blosum62[100][100];

#define KI_MAX_ALI_LEN 5000
int dp[KI_MAX_ALI_LEN][KI_MAX_ALI_LEN];
int bt[KI_MAX_ALI_LEN][KI_MAX_ALI_LEN];
char as1[KI_MAX_ALI_LEN*2];
char as2[KI_MAX_ALI_LEN*2];

bool loadedBlosum62 = false;

void loadBlosum62() {
  int i, j, n = 24;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      blosum62[(int)(blosum62indices[i])][(int)(blosum62indices[j])] = blosum62compact[i][j];
    }
  }
  loadedBlosum62 = true;
  /* for (i = 0; i < n; ++i) { */
  /*   pm("%c", blosum62indices[i]); */
  /*   for (j = 0; j < n; ++j) { */
  /*     pm("%3d ", blosum62[(int)(blosum62indices[i])][(int)(blosum62indices[j])]); */
  /*   } */
  /*   pm("\n"); */
  /* } */
}


float kiPairwiseAlignmentLocal(char* s1, char* s2) {
  if (!loadedBlosum62) loadBlosum62();

  /* float lambda   = 0.340; */
  /* float K        = 0.151; */
  /* float H        = 0.470; */

  float lambda_g = 0.267;
  float K_g      = 0.0410;
  /* float H_g      = 0.140; */

  int len1 = strlen(s1);
  int len2 = strlen(s2);

  assert(len1 < KI_MAX_ALI_LEN && len2 < KI_MAX_ALI_LEN);
  
  int i,j, k;
  for (i = 0; i <= len1; ++i) { dp[i][0] = 0; bt[i][0] = 1; }
  for (j = 1; j <= len2; ++j) { dp[0][j] = 0; bt[0][j] = 2;}

  /* pm("  %s\n", s2); */
  pmsg(8, "   ");
  for (i = 0; i < len2; ++i) pmsg(8, "%7c ", s2[i]); pmsg(8, "\n");

  int gpi = 11, gpe = 1;        /* gap penalties: initiation and extension */
  int max_dp = 0, best_i = 0, best_j = 0;
  int v0, v, bt0;
  for (i = 1; i <= len1; ++i) {
    pmsg(8, "%c ", s1[i-1]);
    for (j = 1; j <= len2; ++j) {
      v0 = dp[i-1][j-1] + blosum62[(int)s1[i-1]][(int)s2[j-1]];
      bt0 = 0;
      /* if (v0 < 0) { */
      /*   v0 = 0; */
      /*   bt0 = bts; */
      /* } */
      for (k = 1; k < i; ++k) {
        if (dp[k][j] > v0) {
          v = dp[k][j] - gpi - gpe * (i - k - 1);
          if (v > v0) {
            v0 = v;
            bt0 = (i-k);
          }
        }
      }
      for (k = 1; k < j; ++k) {
        if (dp[i][k] > v0) {
          v = dp[i][k] - gpi - gpe * (j - k - 1);
          if (v > v0) {
            v0 = v;
            bt0 = -(j-k);
          }
        }
      }
      dp[i][j] = v0;
      bt[i][j] = bt0;

      if (dp[i][j] > max_dp) {
        max_dp = dp[i][j];
        best_i = i; best_j = j; 
      }
      /* pm("%d", bt[i][j]); */
      pmsg(8, " %3d:%+3d", dp[i][j],bt[i][j]);
      /* pm("%d", dp[i][j]%10); */
      /* pm(" %3d", dp[i][j]); */
    }
    pmsg(8, "\n");
  }

  int b1, b2, e1, e2;
  b1 = b2 = e1 = e2 = KI_MAX_ALI_LEN - 1;
  as1[e1] = as2[e2] = '\0';

  /* for (i = len1, j = len2; i > 0 || j > 0; ) { */
  for (i = best_i, j = best_j; i > 0 || j > 0; ) {
    if (bt[i][j] == 0) {
      as1[--b1] = s1[--i];
      as2[--b2] = s2[--j];
    } else if (bt[i][j] > 0) {
      for (k = bt[i][j]; k > 0; --k) {
        as1[--b1] = s1[--i];
        as2[--b2] = '-';
      }
    } else {
      for (k = -bt[i][j]; k > 0; --k) {
        as1[--b1] = '-';
        as2[--b2] = s2[--j];
      }
    }
  }
  
  kipmsg(6, "s1  = %s\ns2  = %s\n", s1, s2);
  kipmsg(6, "as1 = %s\nas2 = %s\n", as1+b1, as2+b2);

  double score = (double)dp[best_i][best_j];
  double bs    = (lambda_g * score - log(K_g)) / log(2.0);
  double e     = best_i * best_j * pow(2, -bs);
  
  kipmsg(7, "Local aligner: score = %.1f bits (%.0f), e-value = %g\n", bs, score, e);
  kipmsg(7, "best i, j = %d/%d, %d/%d\n", best_i, len1, best_j, len2);

  return bs;
}


float kiPairwiseAlignmentGlobal(char* s1, char* s2) {
  if (!loadedBlosum62) loadBlosum62();

  /* float lambda   = 0.340; */
  /* float K        = 0.151; */
  /* float H        = 0.470; */

  float lambda_g = 0.267;
  float K_g      = 0.0410;
  /* float H_g      = 0.140; */
  
  int len1 = strlen(s1);
  int len2 = strlen(s2);

  assert(len1 < KI_MAX_ALI_LEN && len2 < KI_MAX_ALI_LEN);

  int gpi = 11, gpe = 1;        /* gap penalties: initiation and extension */
  
  int i,j, k;
  for (i = 0; i <= len1; ++i) { dp[i][0] = -gpe * i; bt[i][0] = 1; }
  for (j = 1; j <= len2; ++j) { dp[0][j] = -gpe * j; bt[0][j] = 2; }

  /* pm("  %s\n", s2); */
  /* pm("   "); */
  /* for (i = 0; i < len2; ++i) pm("%7c ", s2[i]); pm("\n"); */

  /* int max_dp = 0, best_i = 0, best_j = 0; */
  int v0, v, bt0;
  for (i = 1; i <= len1; ++i) {
    /* pm("%c ", s1[i-1]); */
    for (j = 1; j <= len2; ++j) {
      v0 = dp[i-1][j-1] + blosum62[(int)s1[i-1]][(int)s2[j-1]];
      bt0 = 0;
      for (k = 1; k < i; ++k) {
        if (dp[k][j] > v0) {
          v = dp[k][j] - gpi - gpe * (i - k - 1);
          if (v > v0) {
            v0 = v;
            bt0 = (i-k);
          }
        }
      }
      for (k = 1; k < j; ++k) {
        if (dp[i][k] > v0) {
          v = dp[i][k] - gpi - gpe * (j - k - 1);
          if (v > v0) {
            v0 = v;
            bt0 = -(j-k);
          }
        }
      }
      dp[i][j] = v0;
      bt[i][j] = bt0;

      /* if (dp[i][j] > max_dp) { */
      /*   max_dp = dp[i][j]; */
      /*   best_i = i; best_j = j;  */
      /* } */

      /* pm("%d", bt[i][j]); */
      /* pm(" %3d:%+3d", dp[i][j],bt[i][j]); */
      /* pm("%d", dp[i][j]%10); */
      /* pm(" %3d", dp[i][j]); */
    }
    /* pm("\n"); */
  }

  int b1, b2, e1, e2;
  b1 = b2 = e1 = e2 = KI_MAX_ALI_LEN - 1;
  as1[e1] = as2[e2] = '\0';

  for (i = len1, j = len2; i > 0 || j > 0; ) {
  /* for (i = best_i, j = best_j; i > 0 || j > 0; ) { */
    if (bt[i][j] == 0) {
      as1[--b1] = s1[--i];
      as2[--b2] = s2[--j];
    } else if (bt[i][j] > 0) {
      for (k = bt[i][j]; k > 0; --k) {
        as1[--b1] = s1[--i];
        as2[--b2] = '-';
      }
    } else {
      for (k = -bt[i][j]; k > 0; --k) {
        as1[--b1] = '-';
        as2[--b2] = s2[--j];
      }
    }
  }
  
  kipmsg(6, "s1  = %s\ns2  = %s\n", s1, s2);
  kipmsg(6, "as1 = %s\nas2 = %s\n", as1+b1, as2+b2);

  double score = (double)dp[len1][len2];
  double bs    = (lambda_g * score - log(K_g)) / log(2.0);
  double e     = len1 * len2 * pow(2, -bs);
  
  kipmsg(7, "Global aligner: score = %.1f bits (%.0f), e-value = %g\n", bs, score, e);
  return bs;
}

float kiPairwiseAlignmentSimpleGap(char* s1, char* s2) {
  if (!loadedBlosum62) loadBlosum62();

  /* float lambda   = 0.340; */
  /* float K        = 0.151; */
  /* float H        = 0.470; */

  float lambda_g = 0.267;
  float K_g      = 0.0410;
  /* float H_g      = 0.140; */

  int len1 = strlen(s1);
  int len2 = strlen(s2);

  assert(len1 < KI_MAX_ALI_LEN && len2 < KI_MAX_ALI_LEN);
  
  int i,j;
  for (i = 0; i <= len1; ++i) { dp[i][0] = 0; bt[i][0] =  1; }
  for (j = 1; j <= len2; ++j) { dp[0][j] = 0; bt[0][j] = -1; }

  /* pm("  %s\n", s2); */
  /* pm("   "); */
  /* for (i = 0; i < len2; ++i) pm("%5c ", s2[i]); pm("\n"); */
  
  int max_dp = 0, best_i = 0, best_j = 0;
  int v0, v1, v2;
  for (i = 1; i <= len1; ++i) {
    /* pm("%c ", s1[i-1]); */
    for (j = 1; j <= len2; ++j) {
      v0 = dp[i-1][j-1] + blosum62[(int)s1[i-1]][(int)s2[j-1]];
      v1 = dp[i-1][j] - ((i == 1 || bt[i-1][j] ==  1) ? 1 : 11);
      v2 = dp[i][j-1] - ((j == 1 || bt[i][j-1] == -1) ? 1 : 11);
      /* v0 = dp[i-1][j-1] + (s1[i-1] == s2[j-1] ? 2 : -1); */
      /* v1 = dp[i-1][j] + ((i == 1 || bt[i-1][j] == 1) ? 0 : -1); */
      /* v2 = dp[i][j-1] + ((j == 1 || bt[i][j-1] == 2) ? 0 : -1); */
      if (v0 > v1 && v0 > v2) { /* strictly greater to always prefer extending gaps */
        dp[i][j] = v0; bt[i][j] = 0;
      } else if (v1 >= v2) {
        dp[i][j] = v1; bt[i][j] = 1;
      } else {
        dp[i][j] = v2; bt[i][j] = -1;
      }
      if (dp[i][j] > max_dp) {
        max_dp = dp[i][j];
        best_i = i; best_j = j;
      }
      /* pm("%d", bt[i][j]); */
      /* pm(" %3d:%d", dp[i][j],bt[i][j]); */
      /* pm("%d", dp[i][j]%10); */
      /* pm(" %3d", dp[i][j]); */
    }
    /* pm("\n"); */
  }

  int b1, b2, e1, e2;
  b1 = b2 = e1 = e2 = KI_MAX_ALI_LEN - 1;
  as1[e1] = as2[e2] = '\0';

  /* for (i = len1, j = len2; i > 0 || j > 0; ) { */
  for (i = best_i, j = best_j; i > 0 || j > 0; ) {
    if (bt[i][j] == 0) {
      as1[--b1] = s1[--i];
      as2[--b2] = s2[--j];
    } else if (bt[i][j] == 1) {
      as1[--b1] = s1[--i];
      as2[--b2] = '-';
    } else {
      as1[--b1] = '-';
      as2[--b2] = s2[--j];
    }
  }
  
  kipmsg(6, "s1  = %s\ns2  = %s\n", s1, s2);
  kipmsg(6, "as1 = %s\nas2 = %s\n", as1+b1, as2+b2);

  double score = (double)dp[best_i][best_j];
  double bs    = (lambda_g * score - log(K_g)) / log(2.0);
  double e     = best_i * best_j * pow(2, -bs);

  kipmsg(6, "simple aligner: score = %.1f bits (%.0f), e-value = %g\n", bs, score, e);
  
  return bs;
}




/* Sequence container */
alignment_t* kiAllocAlignment() {
  alignment_t* aln = (alignment_t*)kimalloc(sizeof(alignment_t));
  aln->nPos = 0;
  aln->nSeq = 0;
  aln->names = NULL;
  aln->seqs  = NULL;
  aln->flags = NULL;
  aln->nSaved = 0;
  return aln;
}

alignment_t *kiFreeAlignment(alignment_t* aln) { /* for non-arena use only */
  if (aln == NULL)
    return NULL;

  if (!KI_USE_ARENA) {
    int i;
    for (i = 0; i < aln->nSeq; i++) {
      aln->names[i] = kifreestr(aln->names[i]);
      aln->seqs[i]  = kifreestr(aln->seqs[i]);
    }
  }
  
  aln->names = kifree(aln->names, sizeof(char*)*aln->nSaved);
  aln->seqs  = kifree(aln->seqs,  sizeof(char*)*aln->nSaved);
  aln->flags = kifree(aln->flags, sizeof(int)*aln->nSaved);
  kifree(aln, sizeof(alignment_t));
  return NULL;
}

void kiClearAlignment(/*IN/OUT*/alignment_t* aln) {  
  assert(aln != NULL);
  if (KI_USE_ARENA) {
    memset(aln->names, 0, sizeof(char*) * aln->nSeq);
    memset(aln->seqs,  0, sizeof(char*) * aln->nSeq);
    memset(aln->flags, 0, sizeof(int)   * aln->nSeq);
  } else {
    int i;
    char** p = aln->names;
    char** q = aln->seqs;
    int* f = aln->flags;
    
    for (i = 0; i < aln->nSeq; ++i, ++p, ++q, ++f) {
      if (*p) *p = kifreestr(*p);
      if (*q) *q = kifreestr(*q);
      *f = 0;
    }
  }
  aln->nSeq = 0;
  /* leave aln->nSaved and the memory allocated unchanged */
}

void kiAlignmentGrow(alignment_t* aln) {
  aln->nSeq++;
  
  if (aln->nSaved == 0) {
    aln->nSaved = 100;
    aln->seqs   = (char**)kimalloc(sizeof(char*) * aln->nSaved);
    aln->names  = (char**)kimalloc(sizeof(char*) * aln->nSaved);
    aln->flags  = (int*)kimalloc(sizeof(int) * aln->nSaved);
  }
  if (aln->nSeq > aln->nSaved) {
    int nNewSaved = aln->nSaved * 2;
    aln->seqs   = kirealloc(aln->seqs,  sizeof(char*)*(aln->nSaved), sizeof(char*)*nNewSaved, /*copy*/false);
    aln->names  = kirealloc(aln->names, sizeof(char*)*(aln->nSaved), sizeof(char*)*nNewSaved, /*copy*/false);
    aln->flags  = kirealloc(aln->flags, sizeof(int)*(aln->nSaved),   sizeof(int)*nNewSaved,   /*copy*/false);
    aln->nSaved = nNewSaved;
  }
}

seq_id_t kiAlignmentParseAndAdd(alignment_t* aln, char* name, char* seq) {
  kiAlignmentGrow(aln);
  
  int currIndex = aln->nSeq - 1;
  aln->names[currIndex] = kistrdup(name);
  aln->seqs[currIndex]  = NULL;
  aln->flags[currIndex] = KI_SEQ_NEW;
  
  char *seqSkip = " \t\r\n";
  int nKeep = 0;
  char *p, *q;
  for (p = seq; *p != '\0'; p++) {
    for (q = seqSkip; *q != '\0'; q++) {
      if (*p == *q)
        break;
    }
    if (*p != *q)
      nKeep++;
  }
  int nOld = (aln->seqs[currIndex] == NULL) ? 0 : strlen(aln->seqs[currIndex]);
  aln->seqs[currIndex] = (char*)kirealloc(aln->seqs[currIndex], nOld, nOld+nKeep+1, /*copy*/false);
  if (nOld + nKeep > aln->nPos)
    aln->nPos = nOld + nKeep;
  char *out = aln->seqs[currIndex] + nOld;
  for (p = seq; *p != '\0'; p++) {
    for (q = seqSkip; *q != '\0'; q++) {
      if (*p == *q)
        break;
    }
    if (*p != *q) {
      *out = *p;
      out++;
    }
  }
  assert(out - aln->seqs[currIndex] == nKeep + nOld);
  *out = '\0';
  
  seq_id_t id;
  id.id  = currIndex;
  id.cpu = ki_domain_rank;
  return id;
}
 
void kiAlignmentAdd(alignment_t* aln, char* name, char* seq, bool bCopy) {
  kiAlignmentGrow(aln);
  int currIndex = aln->nSeq - 1;
  aln->names[currIndex] = bCopy ? kiArenaStrdup(KI_ARENA_DEFAULT, name) : name;
  aln->seqs[currIndex]  = bCopy ? kiArenaStrdup(KI_ARENA_DEFAULT, seq)  : seq;
  aln->flags[currIndex] = KI_SEQ_NEW;
}



int kiAbsSeqId(int seqid) {
  return seqid >= 0 ? seqid : (-seqid-1);
}

bool kiAlignmentHasZeroOffset(alignment_t* aln) {
  bool hasZero = false;
  char *swap;
  int i;
  for (i = 0; i < aln->nSeq; ++i) {
    if (aln->flags[i] == 0) {
      hasZero = true;
      if (i > 0) {              /* move to top of the list */
        swap = aln->seqs[0];
        aln->seqs[0] = aln->seqs[i];
        aln->seqs[i] = swap;
        swap = aln->names[0];
        aln->names[0] = aln->names[i];
        aln->names[i] = swap;
        aln->flags[i] = aln->flags[0];
        aln->flags[0] = 0;
      }
      break;
    }
  }

  return hasZero;
}

alignment_t* kiReadFasta(/*IN*/FILE* fp) {
  int    nSeq  = 0;
  int    nPos  = 0;
  char** names = NULL;
  char** seqs  = NULL;
  int*   flags = NULL;
  char buf[KI_SEQ_BUF_SIZE] = "";
  if (fgets(buf, sizeof(buf), fp) == NULL) {
    fprintf(stderr, "Error reading header line\n");
    kiAbort(1);
  }
  int nSaved = 100;

  /* FASTA, truncate names at any of these */
  char *nameStop = "(),: \t\r\n";
  char *seqSkip  = " \t\r\n";

  seqs  = (char**)kimalloc(sizeof(char*) * nSaved);
  names = (char**)kimalloc(sizeof(char*) * nSaved);
  flags = kimalloc(sizeof(int) * nSaved);
  
  do {
    if (strlen(buf) >= sizeof(buf)) {
      fprintf(stderr, "Line too long in fasta file.\n");
      kiAbort(1);
    }
    /* loop over lines */
    if (buf[0] == '>') {
      /* truncate the name */
      char *p, *q;
      for (p = buf+1; *p != '\0'; p++) {
        for (q = nameStop; *q != '\0'; q++) {
          if (*p == *q) {
            *p = '\0';
            break;
          }
        }
        if (*p == '\0') break;
      }
      /* allocate space for another sequence */
      nSeq++;
      if (nSeq > nSaved) {
        int nNewSaved = nSaved * 2;
        seqs  = kirealloc(seqs,  sizeof(char*)*nSaved, sizeof(char*)*nNewSaved, /*copy*/false);
        names = kirealloc(names, sizeof(char*)*nSaved, sizeof(char*)*nNewSaved, /*copy*/false);
        flags = kirealloc(flags, sizeof(int)*nSaved,   sizeof(int)*nNewSaved,   /*copy*/false);
        nSaved = nNewSaved;
      }
      names[nSeq-1] = (char*)kimemdup(buf+1, strlen(buf));
      seqs[nSeq-1]  = NULL;
      flags[nSeq-1] = KI_SEQ_NEW;
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
      int nOld = (seqs[nSeq-1] == NULL) ? 0 : strlen(seqs[nSeq-1]);
      if (seqs[nSeq-1] == NULL) {
        seqs[nSeq-1] = (char*)kimalloc(nKeep+1);
      } else {
        seqs[nSeq-1] = (char*)kirealloc(seqs[nSeq-1], nOld+1, nOld+nKeep+1, /*copy*/false);
      }
      if (nOld+nKeep > nPos)
        nPos = nOld + nKeep;
      char *out = seqs[nSeq-1] + nOld;
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
      assert(out - seqs[nSeq-1] == nKeep + nOld);
      *out = '\0';
    }
  } while(fgets(buf,sizeof(buf),fp) != NULL);

  if (seqs[nSeq-1] == NULL) {
    fprintf(stderr, "No sequence data for last entry %s\n", names[nSeq-1]);
    kiAbort(1);
  }

  /* names = kirealloc(names, sizeof(char*)*nSaved, sizeof(char*)*nSeq, /\*copy*\/false); */
  /* seqs  = kirealloc(seqs,  sizeof(char*)*nSaved, sizeof(char*)*nSeq, /\*copy*\/false); */
  /* flags = kirealloc(flags, sizeof(int)*nSaved,   sizeof(int)*nSeq,   /\*copy*\/false); */
  /* nSaved = nSeq; */

  int i;
  for (i = 0; i < nSeq; i++) 
    kiUpcaseSeq(seqs[i]);

  if (ferror(fp)) {
    fprintf(stderr, "Error reading input file\n");
    kiAbort(1);
  }

  alignment_t* align = (alignment_t*)kimalloc(sizeof(alignment_t));
  align->nSeq   = nSeq;
  align->nPos   = nPos;
  align->names  = names;
  align->seqs   = seqs;
  align->flags  = flags;
  align->nSaved = nSaved;
  
  return align;
}


/* Origin structure to be used in combination with hashtable_t */
origin_list_t* kiAllocOriginList(int seqid, int offset) {
  origin_list_t* list = kiArenaMalloc(KI_ARENA_HASH, sizeof(origin_list_t));
  list->origin.seqid  = seqid;
  list->origin.offset = offset;
  list->next          = NULL;
  list->prev          = NULL;
  return list;
}

origin_list_t* kiFreeOriginList(origin_list_t* list) {
  if (list == NULL) return NULL;
  if (KI_USE_ARENA) return NULL;
  list->next = kiFreeOriginList(list->next);
  kifree(list, sizeof(origin_list_t));
  return NULL;
}

void kiPrintOriginList(origin_list_t* list) {
  origin_list_t* p = list;
  pmsg(0, "[ ");
  while (p) {
    pmsg(0, "(%d,%d) ", p->origin.seqid, p->origin.offset);
    p = p->next;
  }
  pmsg(0, "]\n");
}

