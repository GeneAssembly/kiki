#ifndef _SEQ_H_
#define _SEQ_H_

typedef struct {
  int id;
  int cpu;
} seq_id_t;

/* Sequence container  */
typedef struct {
  int nPos;
  int nSeq;
  char** names;
  char** seqs;
  int* flags;
  int nSaved; /* actual allocated size of names and seqs */
} alignment_t;


/* General sequence manipulation */
void kiUpcaseSeq(char* seq);
void kiReverseComplementSeq(char* seq, /*OUT*/char* comp);
void kiPartialReverseComplementSeq(char* seq, int beg, /*OUT*/char* comp);
void kiReverseComplementSeqN(char* seq, int n, /*OUT*/char* comp);
char kiComplementBase(char c);
int  kiBaseIndex(char c);
float kiSeqCmp(char* s1, char* s2); /* returns percent of discrepancy */
float kiSeqNCmp(char* s1, char* s2, int n);
float kiSeqNCmpX(char* s1, char* s2, int n, /*OUT*/int* iDiff);
bool kiIsProtein(char* s);
void kiDna2Protein(char* dna, char* aa, int frame);
float kiCalcGCContent(char* dna); /* in percentage */

/* Alignment routines */
float kiPairwiseAlignmentLocal(char* s1, char* s2);
float kiPairwiseAlignmentGlobal(char* s1, char* s2);
float kiPairwiseAlignmentSimpleGap(char* s1, char* s2);


/* Sequence container */
alignment_t* kiAllocAlignment();
alignment_t* kiReadFasta(/*IN*/FILE *fp);
void kiAlignmentGrow(alignment_t* aln);
seq_id_t kiAlignmentParseAndAdd(alignment_t* aln, char* name, char* seq);
void  kiAlignmentAdd(alignment_t* aln, char* name, char* seq, bool bCopy);
alignment_t* kiFreeAlignment(alignment_t*); /* returns NULL */
void kiClearAlignment(/*IN/OUT*/alignment_t*);
int kiAbsSeqId(int seqid);
bool kiAlignmentHasZeroOffset(alignment_t* aln);


/* Origin structure to be used in combination with hashtable_t */
typedef struct {
  int seqid;             /* index of sequence entry */
  int offset;            /* position of kmer relative to seq start */
} kmerorigin_t;

typedef struct origin_list_ {
  kmerorigin_t origin;
  struct origin_list_* next;
  struct origin_list_* prev;
} origin_list_t;
typedef origin_list_t* origin_iterator_t;

origin_list_t* kiAllocOriginList(int seqid, int offset);
origin_list_t* kiFreeOriginList(origin_list_t* list);
void kiPrintOriginList(origin_list_t* list);


#endif /* _SEQ_H_ */
