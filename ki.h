#ifndef _KI_H_
#define _KI_H_

#include <stdbool.h>
#include <stdarg.h>
 
#include "extern.h"


#define KI_UNDEF                    -1

#define KI_CMD_STOP                 -1
#define KI_CMD_MEM_INFO             1
#define KI_CMD_LOAD_SEQ             2
#define KI_CMD_LOAD_FASTA           3
#define KI_CMD_READ_FASTA           4
#define KI_CMD_READ_FASTA_OR_FASTQ  5
#define KI_CMD_LOAD_READ_SEQS       6
#define KI_CMD_GET_MIN_MAX_READLEN  7
#define KI_CMD_PROCESS_VARLEN_READS 8
#define KI_CMD_SET_K                9
#define KI_CMD_GET_KMER             10
#define KI_CMD_SET_HASH_MODE        11
#define KI_CMD_DUMMY_ASSEMBLE       12
#define KI_CMD_GET_SEED_SEQ         13
#define KI_CMD_GET_OVERLAPPING_SEQS 14
#define KI_CMD_GET_PREFIX_SEQS      15
#define KI_CMD_SEARCH_PROFILE       16
#define KI_CMD_RAIPHY_ORIGINAL      17
#define KI_CMD_STORE_STATE          18
#define KI_CMD_RESTORE_STATE        19
#define KI_CMD_TEST_ARGS            99


#define KI_OGV_NEW                  0
#define KI_OGV_TODO                 1
#define KI_OGV_EXTEND               2
#define KI_OGV_SOLID                3

#define KI_OG_MAX_SIZE              1000


/* Output macros */
#define KI_PM_FMT           "[%d/%d][%s %d/%d] "
#define KI_PM_ARG           ki_world_rank, ki_world_size, ki_domain_name, ki_domain_rank, ki_domain_size
#define pm0(...)            if (ki_domain_rank == 0) { pm(__VA_ARGS__); }
#define pmsg0(LVL, ...)     if (ki_domain_rank == 0) { pmsg(LVL, __VA_ARGS__); }
#define kipm(...)           { if (ki_world_size > 1) pm(KI_PM_FMT, KI_PM_ARG); pm(__VA_ARGS__); }
#define kipmsg(LVL, ...)    { if (ki_world_size > 1) pmsg(LVL, KI_PM_FMT, KI_PM_ARG); pmsg(LVL, __VA_ARGS__); }
#define kipm0(...)          if (ki_domain_rank == 0) { if (ki_world_size > 1) pm(KI_PM_FMT, KI_PM_ARG); pm(__VA_ARGS__); }
#define kipmsg0(LVL, ...)   if (ki_domain_rank == 0) { if (ki_world_size > 1) pmsg(LVL, KI_PM_FMT, KI_PM_ARG); pmsg(LVL, __VA_ARGS__); }


/* API */
void kiUserTestArgs(/*IN */int n, int* array, char* kmer, int kmer_id,
                    /*OUT*/int* new_n, int* new_array, char* new_kmer, int* new_kmer_id);
int  kiFarmerTestArgs();

void kiUserMemInfo(/*OUT*/long* used);
int  kiFarmerMemInfo();

void kiUserLoadSeq(char* seq, /*OUT*/seq_id_t* id);
int  kiFarmerLoadSeq();
 
void kiUserReadFasta(char* fileName, /*OUT*/long long* nSeqs);
int  kiFarmerReadFasta();

void kiUserReadFastaOrFastq(char* fileName, /*OUT*/long long* nSeqs);
int  kiFarmerReadFastaOrFastq();

void kiUserLoadFasta(char* fileName, /*OUT*/long long* nSeqs);
int  kiFarmerLoadFasta();

void kiUserLoadReadSeqs();
int  kiFarmerLoadReadSeqs();

void kiUserGetMinMaxReadLen(/*OUT*/int* minLen, int* maxLen);
int  kiFarmerGetMinMaxReadLen();

void kiUserProcessVarLenReads(int cutoff, int minLen, int maxLen, int minOverlap, /*OUT*/long long* nSeqsLeft);
int  kiFarmerProcessVarLenReads();

void kiUserSetK(int k);
int  kiFarmerSetK();

void kiUserGetKmer(char* kmer, /*OUT*/int* nSeqs);
int  kiFarmerGetKmer();

void kiUserSetHashMode(int hashMode);
int  kiFarmerSetHashMode();

void kiUserDummyAssemble(char* fileName);
int  kiFarmerDummyAssemble();

void kiUserGetSeedSeq(/*OUT*/char* seed);
int  kiFarmerGetSeedSeq();

void kiUserGetOverlappingSeqs(char* query, int minOverlap, float maxMismatch, bool bErase, /*OUT*/alignment_t* aln);
int  kiFarmerGetOverlappingSeqs();

void kiUserGetPrefixSeqs(char* query, int minOverlap, float maxMismatch, bool bErase, /*OUT*/alignment_t* aln);
int  kiFarmerGetPrefixSeqs();

void kiUserSearchProfile(char* profileName, float cutoff, /*OUT*/alignment_t* hits);
int  kiFarmerSearchProfile();

void kiUserRAIphyOriginal(char* dbName, char* binName);
int  kiFarmerRAIphyOriginal();

void kiUserStoreState(char* filePrefix, /*OUT*/long long* nProcessed);
int  kiFarmerStoreState();

void kiUserRestoreState(char* filePrefix, /*OUT*/long long* nProcessed);
int  kiFarmerRestoreState();

/* System functions */

void kiInit(int* argc, char** argv[]);
void kiFinalize();
void kiStart();
void kiStop();
void kiAbort(int error);
bool kiIsDomainRoot();
bool kiIsMyTurn();
void kiUserCall(int cmd);
void kiResetProgress();
void kiReportProgress(int name, char* info, long long count, long long total, int inteval);
void kiReportRootProgress(int name, char* info, int count, int total, int interval);


/* For internal use only */

typedef int (*ki_func_t)();

void kiProcessArgs(int* argc, char** argv[]);
bool kiIsParallel();
int  kiRunCommand(int cmd);

bool kiIsUserDomain();
bool kiIsFarmerDomain();

void kiDomainInit();
void kiDomainFinalize();

void kiDataInit();
void kiDataFinalize();

void kiUserStart();
void kiDealerStart();
void kiFarmerStart();

int  kiRandDealer();
int  kiRandFarmer();
int  kiNextFarmer();

void kiRegisterCommand(int cmd, ki_func_t func);
void kiRegisterCommands();
void kiClearCommands();

void kiPackArgs(void* buf, int bufSize, /*IN/OUT*/int* position, /*IN*/...);
void kiUnpackArgs(void* buf, int bufSize, /*IN/OUT*/int* position, /*OUT*/...);


/* Sequence packaging */
void kiPackMatchingSeq(void* buf, int bufSize, int* position, /*IN*/char* name, char* seq, int flag);
void kiUnpackAlignment(void* buf, int bufSize, int* position, /*OUT*/alignment_t* aln);

/* Collective hashtable routines */
void kiFreeHashtables();
void kiHashtableLoadSeq(int id);
void kiHashtableLoadAllSeqs();
void kiHashtableLoadAllSeqsMinLength();
void kiHashtableLoadAllSeqsFixedLength();

/* Utility functions */
bool kiIsPowerOfTwo(int x);
void kiPrintBuffer(char* buf, int len);
void kiTimeToString(double time, /*OUT*/char* str);


/* Implementation-specific data structures and methods */

/* Set implemented with simple int array */
typedef struct {
  int size;
  bool* array;
} int_array_set_t;
int_array_set_t* kiAllocIntArraySet();
int_array_set_t* kiFreeIntArraySet(int_array_set_t* set);
void kiIntArraySetAdd(int_array_set_t* set, int x);
bool kiIntArraySetHas(int_array_set_t* set, int x);
  

/* Overlap graph */
typedef struct {
  char* ext;
  float weight;
  struct ogv* v;
} ogedge_t;

typedef struct ogv {
  char* name;
  char* seq;
  ogedge_t* ies[4];             /* incoming edges */
  ogedge_t* oes[4];             /* outgoing edges */
  int occr;
  int offset;
  int flag;
  int id;
} ogvertex_t;

typedef struct int_list_ {
  int value;
  struct int_list_* next;
} int_list_t;

typedef struct {
  int minOverlap;
  int seqLen;                          /* for now, assume all reads are of equal length */
  int nVertices;
  int nSaved;
  ogvertex_t** vertices;
  int_list_t* tailSorted;                  /* list of vertex index tailSorted by offset */
  int_list_t** tailOffsetStarts;    /* array of pointers to elements in tailSorted */
  int_list_t** tailOffsetEnds;      /* array of pointers to elements in tailSorted */
  int nSavedTailOffsets;
  int maxTailOffset;
  int nTailOffsets;
} overlap_graph_t;

/* ogvertex_t* kiAllocOgVertex(overlap_graph_t* g); */
/* ogedge_t*   kiAllocOgEdge(overlap_graph_t* g); */

ogvertex_t* kiAllocOgVertex();
ogvertex_t* kiFreeOgVertex(ogvertex_t* v);

ogedge_t* kiAllocOgEdge();
ogedge_t* kiFreeOgEdge(ogedge_t* e);

int_list_t* kiAllocIntList(int val);
int_list_t* kiFreeIntList(int_list_t* list);
void kiPrintIntList(int_list_t* list);

overlap_graph_t* kiAllocOverlapGraph(int minOverlap);
overlap_graph_t* kiFreeOverlapGraph(overlap_graph_t* g);
void kiClearOverlapGraph(overlap_graph_t* g);

void kiOverlapGraphAppendVertex(overlap_graph_t* g, ogvertex_t* v);
void kiOverlapGraphAdd(overlap_graph_t* g, char* name, char* seq, int offset, int baseOffset, int minOverlap, bool bCopy);
void kiOverlapGraphExportDot(overlap_graph_t* g, char* filename, int step, char* postfix, bool big);
void kiOverlapGraphAddMatches(overlap_graph_t* g, alignment_t* matches, int baseOffset, int minOverlap);
void kiOverlapGraphReviseSeed(overlap_graph_t* g, char* query, /*OUT*/int* iNewSeed);
void kiOverlapGraphExplore(overlap_graph_t* g, int vi, /*OUT*/char* adventure);
bool kiOverlapGraphExtend(overlap_graph_t* g, /*IN/OUT*/int* vi);
void kiOverlapGraphGetContig(overlap_graph_t* g, int iSeed, /*OUT*/char* contig);
void kiOverlapGraphAppendContig(overlap_graph_t* g, int iSeed, /*OUT*/char* contig);
bool kiOgvUnextended(overlap_graph_t* g, int id);
bool kiOverlapGraphTooBig(overlap_graph_t* g);
void kiCombineElongation(char* seed, char* forward, char* backward, /*OUT*/char* contig);

long long kiStoreState(char* filePrefix);
long long kiRestoreState(char* filePrefix);

long long kiStoreStateToSeparateFiles(char* filePrefix);
long long kiRestoreStateFromSeparateFiles(char* filePrefix);


/* Parallel routines on the farm */
int  kiGetKmer(char* kmer);
int  kiReadFastaString(char* buf, alignment_t* seqs);
int  kiReadFastaShared(char* fileName, alignment_t* seqs);
long long kiReadFastaParallel(char* fileName); 
long long kiReadFastqParallel(char* fileName);
long long kiReadFastaSemiParallel(char* fileName); /* returns number of total seqs read */
bool kiIsFastq(char* fileName);


/* Obsolete routines */
void kiExtendSeq(/*IN/OUT*/char* seq, /*IN*/int sz, int beg, int root, /*OUT*/int* advance, int* leader); /* advance: advancement in number of chars; leader: rank of node that achieved greatest advancement */
void kiExactAssemble(char* outputName);


#endif /* _KI_H_ */
