
#include <locale.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "extern.h"


/* Packing and unpacking of function args */
#define KI_BUF_UNPACK       (kiIsParallel() ? ki_in_buf : ki_out_buf)

#define KI_PACK_ARGS(...)   kiPackArgs(ki_out_buf, KI_BUF_SIZE, &pos, __VA_ARGS__, KI_UNDEF)
#define KI_UNPACK_ARGS(...) kiUnpackArgs(KI_BUF_UNPACK, KI_BUF_SIZE, &pos, __VA_ARGS__, KI_UNDEF)

#define KI_FARMER_PACK(...) { int pos = 0;                           \
  kiPackArgs(ki_out_buf, KI_BUF_SIZE, &pos, __VA_ARGS__, KI_UNDEF);  \
  return pos; }

#define KI_FARMER_UNPACK(...) { int pos = 0; int cmd = 0;            \
  kiUnpackArgs(KI_BUF_UNPACK, KI_BUF_SIZE, &pos, KI_INT, &cmd, __VA_ARGS__, KI_UNDEF); }

#define KI_PACK_SEND_CMD(CMD, ...)   { int pos = 0; int cmd = CMD;                \
  kiPackArgs(ki_out_buf, KI_BUF_SIZE, &pos, KI_INT, &cmd, __VA_ARGS__, KI_UNDEF); \
  if (kiIsParallel())                                                             \
    KI_Send(ki_out_buf, pos, MPI_PACKED, kiRandDealer(), 0, ki_cmm_user_dealer);  \
  else                                                                            \
    kiRunCommand(cmd); }

#define KI_RECV_UNPACK_CMD(...) { int pos = 0;  MPI_Status status;  \
  KI_Recv(ki_in_buf, KI_BUF_SIZE, MPI_PACKED, MPI_ANY_SOURCE, 0, ki_cmm_user_dealer, &status); \
  kiUnpackArgs(KI_BUF_UNPACK, KI_BUF_SIZE, &pos, __VA_ARGS__, KI_UNDEF); }

#define KI_PACK_SEND_CMD_BARE(CMD)   { int pos = 0; int cmd = CMD;                \
  kiPackArgs(ki_out_buf, KI_BUF_SIZE, &pos, KI_INT, &cmd, KI_UNDEF);              \
  if (kiIsParallel())                                                             \
    KI_Send(ki_out_buf, pos, MPI_PACKED, kiRandDealer(), 0, ki_cmm_user_dealer);  \
  else                                                                            \
    kiRunCommand(cmd); }

#define KI_RECV_UNPACK_CMD_BARE() { MPI_Status status;  \
  KI_Recv(ki_in_buf, KI_BUF_SIZE, MPI_PACKED, MPI_ANY_SOURCE, 0, ki_cmm_user_dealer, &status);}
  

/* Global variables */

/* Communicators and groups */
MPI_Group ki_grp_world, ki_grp_domain, ki_grp_user, ki_grp_dealer, ki_grp_farmer;
MPI_Comm  ki_cmm_world, ki_cmm_domain, ki_cmm_user, ki_cmm_dealer, ki_cmm_farmer;

MPI_Group ki_grp_user_dealer, ki_grp_dealer_farmer;
MPI_Comm  ki_cmm_user_dealer, ki_cmm_dealer_farmer;

char ki_in_buf[KI_BUF_SIZE];
char ki_out_buf[KI_BUF_SIZE];
char ki_tmp_buf[KI_BUF_SIZE];

int  ki_world_rank  = KI_UNDEF;
int  ki_world_size  = KI_UNDEF;
int  ki_user_rank   = KI_UNDEF;
int  ki_user_size   = KI_UNDEF;
int  ki_dealer_rank = KI_UNDEF;
int  ki_dealer_size = KI_UNDEF;
int  ki_farmer_rank = KI_UNDEF;
int  ki_farmer_size = KI_UNDEF;
int  ki_top_dealer  = KI_UNDEF;
int  ki_top_farmer  = KI_UNDEF;
int  ki_domain_rank = KI_UNDEF;
int  ki_domain_size = KI_UNDEF;
int  ki_domain_id   = KI_UNDEF;
char ki_domain_name[20];

int  ki_round_robin = -1;

float ki_user_ratio = 0.02;     /* (0.02, 0.02, 0.96) */

ki_func_t ki_functions[KI_NUM_FUNCS];

/* Time */
double ki_time_beg;
double ki_time_end;

double ki_clock_start;
double ki_clock_stop;
double ki_clock_interval;
int    ki_clock_name    = KI_CLOCK_NONE;
long long ki_clock_count   = 0;
long long ki_clock_total   = 0;

int ki_kmer_len         = KI_DEF_KMER_LEN;
int ki_hash_mode        = KI_DEF_HASH_MODE;

int ki_exeution_time     = 0;
float ki_read_depletion  = -1.0;

/* Implementation-specific */
alignment_t*  ki_seqs   = NULL;
hashtable_t** ki_hash   = NULL;
int ki_seqs_index       = 0;
int ki_nseq_processed   = 0;


/* API */

void kiRegisterCommands() {
  kiClearCommands();
  kiRegisterCommand(KI_CMD_TEST_ARGS,            kiFarmerTestArgs);
  kiRegisterCommand(KI_CMD_MEM_INFO,             kiFarmerMemInfo);
  kiRegisterCommand(KI_CMD_LOAD_SEQ,             kiFarmerLoadSeq);
  kiRegisterCommand(KI_CMD_LOAD_FASTA,           kiFarmerLoadFasta);
  kiRegisterCommand(KI_CMD_READ_FASTA,           kiFarmerReadFasta);
  kiRegisterCommand(KI_CMD_READ_FASTA_OR_FASTQ,  kiFarmerReadFastaOrFastq);
  kiRegisterCommand(KI_CMD_LOAD_READ_SEQS,       kiFarmerLoadReadSeqs);
  kiRegisterCommand(KI_CMD_SET_K,                kiFarmerSetK);
  kiRegisterCommand(KI_CMD_GET_KMER,             kiFarmerGetKmer);
  kiRegisterCommand(KI_CMD_SET_HASH_MODE,        kiFarmerSetHashMode);
  kiRegisterCommand(KI_CMD_GET_MIN_MAX_READLEN,  kiFarmerGetMinMaxReadLen);
  kiRegisterCommand(KI_CMD_PROCESS_VARLEN_READS, kiFarmerProcessVarLenReads);
  kiRegisterCommand(KI_CMD_DUMMY_ASSEMBLE,       kiFarmerDummyAssemble);
  kiRegisterCommand(KI_CMD_GET_SEED_SEQ,         kiFarmerGetSeedSeq);
  kiRegisterCommand(KI_CMD_GET_OVERLAPPING_SEQS, kiFarmerGetOverlappingSeqs);
  kiRegisterCommand(KI_CMD_GET_PREFIX_SEQS,      kiFarmerGetPrefixSeqs);
  kiRegisterCommand(KI_CMD_SEARCH_PROFILE,       kiFarmerSearchProfile);
  kiRegisterCommand(KI_CMD_RAIPHY_ORIGINAL,      kiFarmerRAIphyOriginal);
  kiRegisterCommand(KI_CMD_STORE_STATE,          kiFarmerStoreState);
  kiRegisterCommand(KI_CMD_RESTORE_STATE,        kiFarmerRestoreState);
  kiRegisterCommand(KI_CMD_SET_TERMINATION,      kiFarmerSetTermination);
}

/* Demo function: packing and unpacking args */

void kiUserTestArgs(/*IN */ int n, int* array, char* kmer, int kmer_id,
                    /*OUT*/ int* new_n, int* new_array, char* new_kmer, int* new_kmer_id) {
  KI_PACK_SEND_CMD(KI_CMD_TEST_ARGS,
                   KI_V_INT,  n, array,
                   KI_STRING, kmer,
                   KI_INT,    &kmer_id);
  
  KI_RECV_UNPACK_CMD(KI_V_INT,  new_n, new_array,
                     KI_STRING, new_kmer,
                     KI_INT,    new_kmer_id);
}

int kiFarmerTestArgs() {
  int array[100], new_array[100];
  char kmer[100], new_kmer[100];
  int n, kmer_len, kmer_id, new_kmer_id;
  int i;
  
  KI_FARMER_UNPACK(KI_V_INT,  &n, array,
                   KI_STRING, kmer,
                   KI_INT,    &kmer_id);

  /* data manipulation */
  for (i = 0; i < n; ++i)
    new_array[i] = 10 + array[i];

  kmer_len = strlen(kmer);
  for (i = 0; i < kmer_len; ++i)
    new_kmer[i] = kmer[kmer_len - i - 1];
  new_kmer_id = -kmer_id;

  /* gather */

  KI_FARMER_PACK(KI_V_INT,  n, new_array,
                 KI_STRING, new_kmer,
                 KI_INT,    &new_kmer_id);
}

void kiUserCall(int cmd) {      /* demo only: unregistered functions */
  MPI_Status status;
  KI_Send(&cmd, 1, MPI_INT, kiRandDealer(), 0, ki_cmm_user_dealer);
  KI_Recv(ki_in_buf, KI_BUF_SIZE, MPI_PACKED, MPI_ANY_SOURCE, 0, ki_cmm_user_dealer, &status);
}

/* Show memory footprint */
void kiUserMemInfo(/*OUT*/long* used) {
  KI_PACK_SEND_CMD_BARE(KI_CMD_MEM_INFO);
  KI_RECV_UNPACK_CMD(KI_LONG, used);
}

int kiFarmerMemInfo() {
  long used = kiGetMallocUsed();
  kipm("Memory used = %ld KB\n", used >> 10);

  long totalUsed = 0;
  KI_Allreduce(&used, &totalUsed, 1, MPI_LONG, MPI_SUM, ki_cmm_domain);
  
  KI_FARMER_PACK(KI_LONG, &totalUsed);
}

/* Load a single sequence */
void kiUserLoadSeq(char* seq, /*OUT*/ seq_id_t* id) {
  KI_PACK_SEND_CMD(KI_CMD_LOAD_SEQ,
                   KI_STRING, seq);
  
  KI_RECV_UNPACK_CMD(KI_INT, &(id->id),
                     KI_INT, &(id->cpu));
}

int kiFarmerLoadSeq() {
  char seq[KI_SEQ_BUF_SIZE];
  KI_FARMER_UNPACK(KI_STRING, seq);

  seq_id_t id;
  kiNextFarmer();
  if (kiIsMyTurn()) {
    id = kiAlignmentParseAndAdd(ki_seqs, "", seq);
    kiHashtableLoadSeq(id.id);
  }
  
  KI_Bcast(&id, sizeof(id), MPI_BYTE, ki_round_robin, ki_cmm_domain);
  
  KI_FARMER_PACK(KI_INT, &(id.id),
                 KI_INT, &(id.cpu));
}

/* Load sequences from a FASTA file */
void kiUserLoadFasta(char* fileName, /*OUT*/long long* nSeqs) {
  KI_PACK_SEND_CMD(KI_CMD_LOAD_FASTA,
                   KI_STRING, fileName);
  KI_RECV_UNPACK_CMD(KI_LONG_LONG, nSeqs);
}

int kiFarmerLoadFasta() {
  char fileName[200];
  KI_FARMER_UNPACK(KI_STRING, fileName);
  
  long long nSeqs = kiReadFastaParallel(fileName);
  
  /* kiHashtableLoadAllSeqs(); */
  kiHashtableLoadAllSeqsFixedLength();

  KI_FARMER_PACK(KI_LONG_LONG, &nSeqs);
}

/* Load sequences from a FASTA file */
void kiUserReadFasta(char* fileName, /*OUT*/long long* nSeqs) {
  KI_PACK_SEND_CMD(KI_CMD_READ_FASTA,
                   KI_STRING, fileName);
  KI_RECV_UNPACK_CMD(KI_LONG_LONG, nSeqs);
}

int kiFarmerReadFasta() {
  char fileName[200];
  KI_FARMER_UNPACK(KI_STRING, fileName);
  
  long long nSeqs = kiReadFastaParallel(fileName);
  
  KI_FARMER_PACK(KI_LONG_LONG, &nSeqs);
}

/* Load sequences from a FASTA or FASTQ file */
void kiUserReadFastaOrFastq(char* fileName, /*OUT*/long long* nSeqs) {
  KI_PACK_SEND_CMD(KI_CMD_READ_FASTA_OR_FASTQ,
                   KI_STRING, fileName);
  KI_RECV_UNPACK_CMD(KI_LONG_LONG, nSeqs);
}

int kiFarmerReadFastaOrFastq() {
  char fileName[200];
  KI_FARMER_UNPACK(KI_STRING, fileName);
  
  long long nSeqs = 0;
  if (kiIsFastq(fileName)) 
    nSeqs = kiReadFastqParallel(fileName);
  else 
    nSeqs = kiReadFastaParallel(fileName);
    
  KI_FARMER_PACK(KI_LONG_LONG, &nSeqs);
}

/* Load sequences already read ki_seqs */
void kiUserLoadReadSeqs() {
  KI_PACK_SEND_CMD_BARE(KI_CMD_LOAD_READ_SEQS);
  KI_RECV_UNPACK_CMD_BARE();
}

int  kiFarmerLoadReadSeqs() {
  kiHashtableLoadAllSeqsFixedLength();
  return 0;
}

void kiUserGetMinMaxReadLen(/*OUT*/int* minLen, int* maxLen) {
  KI_PACK_SEND_CMD_BARE(KI_CMD_GET_MIN_MAX_READLEN);

  KI_RECV_UNPACK_CMD(KI_INT, minLen, KI_INT, maxLen);
}

int kiFarmerGetMinMaxReadLen() {
  int minLen = KI_SEQ_BUF_SIZE, maxLen = 0;
  int i, len;
  char** p = ki_seqs->seqs;
  for (i = 0; i < ki_seqs->nSeq; ++i, ++p) {
    len = strlen(*p);
    minLen = MIN(len, minLen);
    maxLen = MAX(len, maxLen);
  }

  int globalMin, globalMax;
  KI_Allreduce(&minLen, &globalMin, 1, MPI_INT, MPI_MIN, ki_cmm_domain);
  KI_Allreduce(&maxLen, &globalMax, 1, MPI_INT, MPI_MAX, ki_cmm_domain);
  /* kipmsg0(2, "Global read length: min = %d, max = %d\n", globalMin, globalMax); */
  
  KI_FARMER_PACK(KI_INT, &globalMin, KI_INT, &globalMax);
}

void kiUserProcessVarLenReads(int cutoff, int minLen, int maxLen, int minOverlap, /*OUT*/long long* nSeqsLeft) {
  KI_PACK_SEND_CMD(KI_CMD_PROCESS_VARLEN_READS,
                   KI_INT, &cutoff,
                   KI_INT, &minLen,
                   KI_INT, &maxLen,
                   KI_INT, &minOverlap);

  KI_RECV_UNPACK_CMD(KI_LONG_LONG, nSeqsLeft);
}

int kiFarmerProcessVarLenReads() {
  int cutoff, minLen, maxLen, minOverlap; /* cutoff = proposed fixed minReadLen */
  long long nSeqsLeft = 0;
  long long nFiltered = 0;

  KI_FARMER_UNPACK(KI_INT, &cutoff,
                   KI_INT, &minLen,
                   KI_INT, &maxLen,
                   KI_INT, &minOverlap);

  kipmsg0(2, "cutoff = %d, minLen = %d, maxLen = %d, minOverlap = %d\n", cutoff, minLen, maxLen, minOverlap);

  if (cutoff < minLen) cutoff = minLen;
  if (cutoff > maxLen) cutoff = maxLen;

  alignment_t* new_seqs = (alignment_t*)kimalloc(sizeof(alignment_t));

  new_seqs->nSeq   = 0;
  new_seqs->nPos   = 0;
  new_seqs->names  = NULL;
  new_seqs->seqs   = NULL;
  new_seqs->flags  = NULL;
  new_seqs->nSaved = 0;

  long long nSeqs = 0;
  char** p = ki_seqs->seqs;
  char** q = ki_seqs->names;
  char *name, *seq, *s;
  char postfix[100];
  int i, j, len, nRegions, stepSize;
  for (i = 0; i < ki_seqs->nSeq; ++i, ++p, ++q) {
    /* kiArenaStrdup(KI_ARENA_SEQS_EL,  */
    len = strlen(*p);
    if (len < cutoff) {
      nFiltered++;
      continue;
    }

    if (len == cutoff) {
      name = (char*)kiArenaStrdup(KI_ARENA_SEQS_EL, *q);
      seq  = (char*)kiArenaStrndup(KI_ARENA_SEQS_EL, *p, cutoff);
      kiAlignmentAdd(new_seqs, name, seq, /*copy*/false);
      nSeqs++;
    } else {
      nRegions = kiCeilingDevidedBy(len - minOverlap, cutoff - minOverlap);
      stepSize = kiCeilingDevidedBy(len - cutoff, nRegions - 1);

      /* kipm0("seq [all] = %s\n", *p); */
      for (j = 0; j < nRegions; ++j) {
        s = (j < nRegions-1) ? *p + j*stepSize : *p + len - cutoff;
        sprintf(postfix, "_%d/%d", j+1, nRegions);
        seq  = (char*)kiArenaStrndup(KI_ARENA_SEQS_EL, s, cutoff);
        name = (char*)kiArenaMalloc(KI_ARENA_SEQS_EL, strlen(*q) + strlen(postfix) + 1);
        strcpy(name, *q);
        strcpy(name+strlen(*q), postfix);
        kiAlignmentAdd(new_seqs, name, seq, /*copy*/false);
        nSeqs++;

        /* kipm0("%s\n", name); */
        /* int spc = s - *p; */
        /* kipm0("seq [%d/%d] = ", j+1, nRegions); */
        /* for (spc = 0; spc < s - *p; ++spc) kipm0(" "); */
        /* kipm0("%s\n", seq); */
      }
      /* kipm0("\n"); */
      /* kipmsg(2, "len = %d, nRegions = %d\n", len, nRegions); */
    }
  }

  kiFreeAlignment(ki_seqs);
  ki_seqs = new_seqs;

  long long totalFiltered;
  KI_Allreduce(&nFiltered, &totalFiltered, 1, MPI_LONG_LONG_INT, MPI_SUM,  ki_cmm_domain);
  kipmsg0(2, "Short reads filtered = %ld\n", totalFiltered);
  
  KI_Allreduce(&nSeqs, &nSeqsLeft, 1, MPI_LONG_LONG_INT, MPI_SUM,  ki_cmm_domain);
  KI_FARMER_PACK(KI_LONG_LONG, &nSeqsLeft);
}


/* Adjust kmer length */
void kiUserSetK(int k) {
  KI_PACK_SEND_CMD(KI_CMD_SET_K,
                   KI_INT, &k);
  KI_RECV_UNPACK_CMD_BARE();
}

int  kiFarmerSetK() {
  int k;
  KI_FARMER_UNPACK(KI_INT, &k);
  ki_kmer_len = k;
  return 0;
}

/* Adjust hash mode: ends only or complete seq */
void kiUserSetHashMode(int hashMode) {
  KI_PACK_SEND_CMD(KI_CMD_SET_HASH_MODE,
                   KI_INT, &hashMode);
  KI_RECV_UNPACK_CMD_BARE();
}

int kiFarmerSetHashMode() {
  int hashMode;
  KI_FARMER_UNPACK(KI_INT, &hashMode);
  ki_hash_mode = hashMode;
  return 0;
}

/* Get a list of sequences that own kmer */
void kiUserGetKmer(char* kmer, /*OUT*/int* nSeqs) {
  KI_PACK_SEND_CMD(KI_CMD_GET_KMER,
                   KI_STRING, kmer);
  KI_RECV_UNPACK_CMD(KI_INT, nSeqs);
}

int kiFarmerGetKmer() {
  char kmer[200];
  KI_FARMER_UNPACK(KI_STRING, kmer);

  int nSeqs = kiGetKmer(kmer);
  
  KI_FARMER_PACK(KI_INT, &nSeqs);
}

/* Dummy assembly */
void kiUserDummyAssemble(char* outputName) {
  KI_PACK_SEND_CMD(KI_CMD_DUMMY_ASSEMBLE,
                   KI_STRING, outputName);
  KI_RECV_UNPACK_CMD_BARE();
}

int kiFarmerDummyAssemble() {
  char outputName[200];
  KI_FARMER_UNPACK(KI_STRING, outputName);
                  
  kiExactAssemble(outputName);
  return 0;
}

/* Get seed sequence to extend on */
void kiUserGetSeedSeq(/*OUT*/char* seed, /*OUT*/int* termStatus) {
  KI_PACK_SEND_CMD_BARE(KI_CMD_GET_SEED_SEQ);
  KI_RECV_UNPACK_CMD(KI_STRING, seed,
                     KI_INT, termStatus);
}

int kiFarmerGetSeedSeq() {
  char seed[KI_SEQ_BUF_SIZE] = "";
  int termStatus = KI_TERMINATION_NONE;

  struct {
    int value;
    int index;
  } in, out;
  
  in.value = ki_seqs->nSeq - ki_nseq_processed;
  in.index = ki_domain_rank;
  
  KI_Allreduce(&in, &out, 1, MPI_2INT, MPI_MAXLOC, ki_cmm_domain);
  int cpu   = out.index;
  int nLeft = out.value;

  int interval = (ki_world_size > 8) ? 10 : 1;
  kiReportProgress(KI_CLOCK_MAIN, "Reads assembled", (long long)ki_nseq_processed, (long long)(ki_seqs->nSeq), interval);

  if (ki_domain_rank == cpu) {
    kipmsg(4, "Seqs processed = %d / %d\n", ki_nseq_processed, ki_seqs->nSeq);
    kipmsg(5, "nLeft = %d\n", nLeft);
  }
  
  if (ki_exeution_time > 0 && (int)(MPI_Wtime() - ki_time_beg) >= ki_exeution_time) {
    termStatus = KI_TERMINATION_TIME;
  }  
  if (ki_read_depletion > 0 && 1.0 * ki_clock_count / ki_clock_total >= ki_read_depletion) {
    kipmsg(3, "deplete = %.3f\n", 1.0 * ki_clock_count / ki_clock_total);
    termStatus = KI_TERMINATION_DEPLETION;
  }

  if (termStatus != KI_TERMINATION_NONE) { /* return */
    KI_Bcast_string(seed, cpu, ki_cmm_domain);
    KI_FARMER_PACK(KI_STRING, seed,
                   KI_INT, &termStatus);
  }

  if (nLeft > 0) {
    if (ki_domain_rank == cpu) {
      for (; ki_seqs_index < ki_seqs->nSeq; ++ki_seqs_index) {
        if (ki_seqs->flags[ki_seqs_index] == KI_SEQ_NEW) {      /* unprocessed */
          strcpy(seed, ki_seqs->seqs[ki_seqs_index]);
          ki_seqs->flags[ki_seqs_index] = KI_SEQ_BOOKED;        /* mark booked */
          ki_nseq_processed++;
          ki_seqs_index++;
          break;
        }
      }
    }
  }

  KI_Bcast_string(seed, cpu, ki_cmm_domain);
  KI_FARMER_PACK(KI_STRING, seed,
                 KI_INT, &termStatus);
}

void kiUserSetTermination(int executionTime, float readDepletionFraction) {
  KI_PACK_SEND_CMD(KI_CMD_SET_TERMINATION,
                   KI_INT,   &executionTime,
                   KI_FLOAT, &readDepletionFraction);
  KI_RECV_UNPACK_CMD_BARE();
}

int  kiFarmerSetTermination() {
  int executionTime;
  float readDepletionFraction;
  KI_FARMER_UNPACK(KI_INT,   &executionTime,
                   KI_FLOAT, &readDepletionFraction);

  ki_exeution_time  = executionTime;
  ki_read_depletion = readDepletionFraction;
  
  return 0;
}



/* Get overlapping sequences */
/* sequences that match query[t..$], for t in [0, strlen(query) - minOverlap] */
void kiUserGetOverlappingSeqs(char* query, int minOverlap, float maxMismatch, bool bErase,
                              /*OUT*/alignment_t* aln) {  
  KI_PACK_SEND_CMD(KI_CMD_GET_OVERLAPPING_SEQS,
                   KI_STRING, query,
                   KI_INT,    &minOverlap,
                   KI_FLOAT,  &maxMismatch,
                   KI_BOOL,   &bErase);

  int buflen, pos = 0;
  KI_RECV_UNPACK_CMD(KI_V_CHAR, &buflen, ki_tmp_buf);
  kipmsg(6, "kiUserGetOverlappingSeqs: buflen = %d\n", buflen);
  kiClearAlignment(aln);
  kiUnpackAlignment(ki_tmp_buf, buflen, &pos, aln);
}

int kiFarmerGetOverlappingSeqs() {
  char  query[KI_SEQ_BUF_SIZE];
  int   minOverlap;
  float maxMismatch;
  bool  bErase;
  
  KI_FARMER_UNPACK(KI_STRING, query,
                   KI_INT,    &minOverlap,
                   KI_FLOAT,  &maxMismatch,
                   KI_BOOL,   &bErase);

  kipmsg0(5, "getOverlappingSeqs: query = '%s', minOverlap = %d, maxMismatch = %f\n", query, minOverlap, maxMismatch);
  
  char* buf    = ki_tmp_buf;
  int   buflen = 0;

  hashtable_t* hash = ki_hash[ki_kmer_len];
  if (hash == NULL) {
    if (ki_seqs->nSeq > 0) {
      kipm("Hashtables not found for k = %d\n", ki_kmer_len);
      return 0;
    }
  } else {

    char  match[KI_SEQ_BUF_SIZE];
    char  kmer[KI_MAX_KMER_LEN+1];
    char* name  = NULL;
    int   qlen  = strlen(query);
    int   min_i = 0;
    int   max_i = qlen - MAX(ki_kmer_len, minOverlap);
    float diff;
    int   i;

    for (i = min_i; i <= max_i; ++i) {
      memcpy(kmer, query+i, ki_kmer_len);
      kmer[ki_kmer_len] = '\0';
      kipmsg0(5, "kmer = '%s'\n", kmer);
      hashiterator_t hi = kiFindMatch(hash, kmer);
      origin_iterator_t originIt = kiHashOrigins(hash, hi);
      while (originIt != NULL) {
        originIt = kiHashFetchUnusedMatchingSeq(hash, hi, originIt, bErase, /*OUT*/match, /*OUT*/&name);
        if (originIt == NULL) break;
        diff = kiSeqCmp(query+i, match);
        kipmsg(4, "Match [%2d] = '%s', diff = %f, name = %s\n", i, match, diff, name);
        if (diff <= maxMismatch) {
          kiPackMatchingSeq(buf, KI_BUF_SIZE, &buflen, name, match, i);
          originIt = kiHashEraseMatchingSeq(hash, hi, originIt, bErase);
        } else {
          originIt = originIt->next;
        }
      }
    }
  }
  
  int totLen;
  KI_Gatherv_buffer(buf, buflen, ki_in_buf, KI_BUF_SIZE, 0, ki_cmm_domain, &totLen);

  KI_FARMER_PACK(KI_V_CHAR, totLen, ki_in_buf);
}

/* Get overlapping sequences; (old: return raw buffer) */
void kiUserGetOverlappingSeqsRaw(char* query, int minOverlap, float maxMismatch,
                              /*OUT*/ char* buf, int* bufSize) {
  KI_PACK_SEND_CMD(KI_CMD_GET_OVERLAPPING_SEQS,
                   KI_STRING, query,
                   KI_INT,    &minOverlap,
                   KI_FLOAT,  &maxMismatch);

  KI_RECV_UNPACK_CMD(KI_V_CHAR, bufSize, buf);
}

int kiFarmerGetOverlappingSeqsRaw() {
  char  query[KI_SEQ_BUF_SIZE];
  int   minOverlap;
  float maxMismatch;
  
  KI_FARMER_UNPACK(KI_STRING, query,
                   KI_INT,    &minOverlap,
                   KI_FLOAT,  &maxMismatch);

  kipm0("getOverlappingSeqs: query = '%s', minOverlap = %d, maxMismatch = %f\n", query, minOverlap, maxMismatch);
  
  hashtable_t* hash = ki_hash[ki_kmer_len];
  if (hash == NULL) {
    kipm("Hashtables not found for k = %d\n", ki_kmer_len);
    return 0;
  }

  char* buf    = ki_tmp_buf;
  int   buflen = 0;

  char  match[KI_SEQ_BUF_SIZE];
  char* kmer  = kimalloc(sizeof(char) * (ki_kmer_len+1));
  int   qlen  = strlen(query);
  int   min_i = 0;
  int   max_i = qlen - MAX(ki_kmer_len, minOverlap);
  int   i, j;
  for (i = min_i; i <= max_i; ++i) {
    memcpy(kmer, query+i, ki_kmer_len);
    kmer[ki_kmer_len] = '\0';
    kipmsg0(5, "kmer = '%s'\n", kmer);
    hashiterator_t hi = kiFindMatch(hash, kmer);
    int nMatches = kiHashCount(hash, hi);
    origin_iterator_t originIt = kiHashOrigins(hash, hi);
    if (nMatches > 0) {
      for (j = 0; j < nMatches; ++j, originIt = originIt->next) {
        kiHashGetMatchingSeq(originIt, /*OUT*/match);
        float diff = kiSeqCmp(query+i, match);
        char* name = ki_seqs->names[kiAbsSeqId(originIt->origin.seqid)];
        kipmsg(5, "Match [i=%2d j=%2d] = '%s', diff = %f, name = %s\n", i, j, match, diff, name);
        if (diff <= maxMismatch) {
          strcpy(buf+buflen, match);
          buflen += (strlen(match)+1);
        }
      }
    }
  }
  kifree(kmer, sizeof(char) * (ki_kmer_len+1));

  int totLen;
  KI_Gatherv_buffer(buf, buflen, ki_in_buf, KI_BUF_SIZE, 0, ki_cmm_domain, &totLen);

  /* kipm("my buflen = %d\n", buflen); */
  /* kipm("my buf = "); */
  /* kiPrintBuffer(buf, buflen); */
  /* kipmsg0(1, "gathered buflen = %d\n", totLen); */
  /* kipmsg0(1, "gathered buf = "); */
  /* if (ki_domain_rank == 0)  */
    /* kiPrintBuffer(ki_in_buf, totLen); */

  KI_FARMER_PACK(KI_V_CHAR, totLen, ki_in_buf);
}

/* Get overlapping sequences */
/* sequences that match q_i (q_i = query[i..i+minOverlap-1]) for any i in [0..strlen(query)-minOverlap] */
void kiUserGetPrefixSeqs(char* query, int minOverlap, float maxMismatch, bool bErase,
                              /*OUT*/alignment_t* aln) {  
  KI_PACK_SEND_CMD(KI_CMD_GET_PREFIX_SEQS,
                   KI_STRING, query,
                   KI_INT,    &minOverlap,
                   KI_FLOAT,  &maxMismatch,
                   KI_BOOL,   &bErase);

  int buflen, pos = 0;
  KI_RECV_UNPACK_CMD(KI_V_CHAR, &buflen, ki_tmp_buf);
  kipmsg(6, "kiUserGetPrefixSeqs: buflen = %d\n", buflen);
  kiClearAlignment(aln);
  kiUnpackAlignment(ki_tmp_buf, buflen, &pos, aln);
}

int kiFarmerGetPrefixSeqs() {
  char  query[KI_SEQ_BUF_SIZE];
  int   minOverlap;
  float maxMismatch;
  bool  bErase;
  
  KI_FARMER_UNPACK(KI_STRING, query,
                   KI_INT,    &minOverlap,
                   KI_FLOAT,  &maxMismatch,
                   KI_BOOL,   &bErase);

  kipmsg0(5, "getPrefixSeqs: query = '%s', minOverlap = %d, maxMismatch = %f\n", query, minOverlap, maxMismatch);

  char* buf    = ki_tmp_buf;
  int   buflen = 0;
  
  hashtable_t* hash = ki_hash[ki_kmer_len];
  if (hash == NULL) {
    if (ki_seqs->nSeq > 0) {
      kipm("Hashtables not found for k = %d\n", ki_kmer_len);
      return 0;
    }
  } else {

    char  match[KI_SEQ_BUF_SIZE];
    char  kmer[KI_MAX_KMER_LEN+1];
    char* name  = NULL;
    int   qlen  = strlen(query);
    float diff;
    int   i;
  
    for (i = 0; i < qlen - minOverlap; ++i) {
      memcpy(kmer, query+i, ki_kmer_len);
      kmer[ki_kmer_len] = '\0';
      kipmsg0(5, "kmer = '%s'\n", kmer);
      hashiterator_t hi = kiFindMatch(hash, kmer);
      origin_iterator_t originIt = kiHashOrigins(hash, hi);
      while (originIt != NULL) {
        originIt = kiHashFetchUnusedMatchingSeq(hash, hi, originIt, bErase, /*OUT*/match, /*OUT*/&name);
        if (originIt == NULL) break;
        diff = kiSeqNCmp(query+i, match, minOverlap);
        kipmsg(4, "Match [%2d] = '%s', diff = %f, name = %s\n", i, match, diff, name);
        if (diff <= maxMismatch) {
          kiPackMatchingSeq(buf, KI_BUF_SIZE, &buflen, name, match, i);
          originIt = kiHashEraseMatchingSeq(hash, hi, originIt, bErase);
        } else {
          originIt = originIt->next;
        }
      }
    }
  }
  
  int totLen;
  KI_Gatherv_buffer(buf, buflen, ki_in_buf, KI_BUF_SIZE, 0, ki_cmm_domain, &totLen);

  KI_FARMER_PACK(KI_V_CHAR, totLen, ki_in_buf);
}

void kiUserSearchProfile(char* profileName, float cutoff, /*OUT*/alignment_t* hits) {
  KI_PACK_SEND_CMD(KI_CMD_SEARCH_PROFILE,
                   KI_STRING, profileName,
                   KI_FLOAT,  &cutoff);

  int buflen, pos = 0;
  KI_RECV_UNPACK_CMD(KI_V_CHAR, &buflen, ki_tmp_buf);
  kiClearAlignment(hits);
  kiUnpackAlignment(ki_tmp_buf, buflen, &pos, hits);
}

int  kiFarmerSearchProfile() {
  char profileName[200];
  float cutoff;
  
  KI_FARMER_UNPACK(KI_STRING, profileName,
                   KI_FLOAT,  &cutoff);

  ki_hash_mode = KI_HASH_COMPLETE;
  
  /* load profile */
  int startingIndex = ki_seqs->nSeq;
  int nProfileSeqs = kiReadFastaShared(profileName, ki_seqs);

  int i, j;
  for (i = startingIndex; i < startingIndex + nProfileSeqs; ++i) {
    kiHashtableLoadSeq(i);
  }

  char fileName1[] = "candidates.fa";
  char fileName2[] = "focus.fa";
  char fileName3[] = "match.txt";
  KI_File fh1, fh2, fh3;
  int rc;
  rc = KI_File_open(fileName1, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1);
  rc = KI_File_open(fileName2, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh2);
  rc = KI_File_open(fileName3, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh3);

  /* alignment_t* profile = kiAllocAlignment(); */
  /* kiReadFastaShared(profileName, profile); */
  hashtable_t* hash = ki_hash[ki_kmer_len];
  char kmer[KI_SEQ_BUF_SIZE];
  kmer[ki_kmer_len] = '\0';
  
  /* search in sequences */
  char aa[KI_SEQ_BUF_SIZE];
  int frame;
  char* match = kimalloc(sizeof(char) * KI_SEQ_BUF_SIZE);

  for (i = 0; i < ki_seqs->nSeq; ++i) {
    if (i >= startingIndex && i < startingIndex + nProfileSeqs) continue; /* skip prof seqs */
    /* compare seqs: 1. kmer lookup  -->  2. local alignment */
    /* kipm("%s\n", ki_seqs->seqs[i]); */
    for (frame = 1; frame <= 6; ++frame) {
      kiDna2Protein(ki_seqs->seqs[i], aa, frame);
      /* kipm("strlen(aa) = %d  aa = %s\n", strlen(aa), aa); */
      /* kipm("%3d %d  %s\n", i, frame, aa); */
      for (j = 0; j < (int)strlen(aa) - ki_kmer_len; ++j) {
        memcpy(kmer, aa+j, ki_kmer_len);
        hashiterator_t hi = kiFindMatch(hash, kmer);
        origin_iterator_t originIt = kiHashOrigins(hash, hi);
        if (originIt != NULL) {
          /* kipm("%s, i = %d, j = %d\n", ki_seqs->names[i], i, j); */
          /* kifprintf(fh, ">%s | frame %d | %s\n%s\n", ki_seqs->names[i], frame, ki_seqs->seqs[i], aa); */
          kiHashGetMatchingSeq(originIt, /*OUT*/match);
          float score = kiPairwiseAlignmentSimpleGap(aa, match);
          if (score > cutoff) {
            kmerorigin_t *p = &(originIt->origin);
            int id = kiAbsSeqId(p->seqid);
            kifprintf(fh1, ">%s | frame %d | kmer %s | score %.2f\n%s\n", ki_seqs->names[i], frame, kmer, score, aa);
            kifprintf(fh2, ">%s\n%s\n", ki_seqs->names[i], ki_seqs->seqs[i]);
            kifprintf(fh3, "%s\t%s\t%.2f\t%d\t%s\n", ki_seqs->names[i], ki_seqs->names[id], score, frame, kmer);
          }
          break;
        }
      }
    }

    int interval = (ki_world_size > 8) ? 10 : 1;
    kiReportRootProgress(KI_CLOCK_MAIN, "Reads searched", i+1, ki_seqs->nSeq - nProfileSeqs, interval);

  }
  KI_File_close(&fh1);
  KI_File_close(&fh2);
  KI_File_close(&fh3);

  kifree(match, sizeof(char) * KI_SEQ_BUF_SIZE);
  
  /* gather hits alignment in ki_in_buf */
  int totLen = 0;
  /* kiFreeAlignment(profile); */
  
  KI_FARMER_PACK(KI_V_CHAR, totLen, ki_in_buf);
}

void kiUserRAIphyOriginal(char* dbName, char* binName) {
  KI_PACK_SEND_CMD(KI_CMD_RAIPHY_ORIGINAL,
                   KI_STRING, dbName,
                   KI_STRING, binName);
  
  KI_RECV_UNPACK_CMD_BARE();
}

int  kiFarmerRAIphyOriginal() {
  char dbName[200];
  char binName[200];

  KI_FARMER_UNPACK(KI_STRING, dbName, KI_STRING, binName);

  KI_File fh;
  int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
  int rc = KI_File_open(binName, amode, MPI_INFO_NULL, &fh);
  if (rc != 0) {
    fprintf(stderr, "Could not open output file.\n");
    kiAbort(-1);
  }
  
  int bufSize = 1024*1024;
  char buf[bufSize];
  char* bufTop = buf;
  MPI_Status status;
  
  rai_db_t* db = NULL;
  db = loadDatabase(dbName);

  double margin;
  int i, class;

  for (i = 0; i < ki_seqs->nSeq; ++i) {
    margin = 0.;
    class = classifySequenceOriginal(ki_seqs->seqs[i], db, &margin); 
    sprintf(bufTop, ">%s\n%s\n", ki_seqs->names[i], db->names[class]);
    bufTop += strlen(bufTop);
    if (bufTop-buf > bufSize/2) {
      KI_File_write_shared(fh, buf, bufTop-buf, MPI_CHAR, &status);
      bufTop = buf;
    }
  }

  KI_File_write_shared(fh, buf, bufTop-buf, MPI_CHAR, &status);
  KI_File_close(&fh);

  kiFreeRaiDb(db);

  return 0;
}


/* Persistence: store and restore ki_seqs state */
void kiUserStoreState(char* filePrefix, /*OUT*/long long* nProcessed) {
  KI_PACK_SEND_CMD(KI_CMD_STORE_STATE, KI_STRING, filePrefix);
  KI_RECV_UNPACK_CMD(KI_LONG_LONG, nProcessed);
}

int  kiFarmerStoreState() {
  char filePrefix[KI_NAME_BUF_SIZE];
  long long nProcessed;
  
  KI_FARMER_UNPACK(KI_STRING, filePrefix);

  nProcessed = kiStoreState(filePrefix);

  KI_FARMER_PACK(KI_LONG_LONG, &nProcessed);
}

void kiUserRestoreState(char* filePrefix, /*OUT*/long long* nProcessed) {
  KI_PACK_SEND_CMD(KI_CMD_RESTORE_STATE, KI_STRING, filePrefix);
  KI_RECV_UNPACK_CMD(KI_LONG_LONG, nProcessed);
}

int  kiFarmerRestoreState() {
  char filePrefix[KI_NAME_BUF_SIZE];
  long long nProcessed;

  KI_FARMER_UNPACK(KI_STRING, filePrefix);
  
  nProcessed = kiRestoreState(filePrefix);

  KI_FARMER_PACK(KI_LONG_LONG, &nProcessed);
}


/* System functions */
void kiInit(int* argc, char** argv[]) {
  MPI_Init(argc, argv);
  
  ki_time_beg = MPI_Wtime();

  /* create a private world for KI */
  MPI_Comm_dup(MPI_COMM_WORLD, &ki_cmm_world);
  MPI_Comm_group(ki_cmm_world, &ki_grp_world);

  MPI_Comm_rank(ki_cmm_world, &ki_world_rank);
  MPI_Comm_size(ki_cmm_world, &ki_world_size);
  assert(ki_world_rank >= 0 && ki_world_size > 0);

  /* more debug info with higher level  */
  /* pmSetLevel(1); */
  /* pmSetLevel(3); */
  /* pmSetLevel(4); */

  kiProcessArgs(argc, argv);    /* may override pmLevel */

  kiRegisterCommands();
  
  pmsg(4, "Hello world. I am %d of %d\n", ki_world_rank, ki_world_size);
  
  if (kiIsParallel()) kiDomainInit();
  else {
    MPI_Comm_dup(MPI_COMM_WORLD, &ki_cmm_domain);
    ki_domain_size = 1;
    ki_domain_rank = 0;
  }
  
  kiDataInit();
}

void kiPrintLocalMemoryUsage() {
  long used = kiGetMallocUsed();
  
  char* locale = setlocale(LC_NUMERIC, "en_US.utf-8");
  
  if (locale) {
    kipmsg0(3, "Local memory used = %'ld KB\n", used >> 10);
  } else {
    kipmsg0(3, "Local memory used = %ld KB\n", used >> 10);
  }
}


void kiFinalize() {

  long used = kiGetMallocUsed();
  long totalUsed = 0;

  char* locale = setlocale(LC_NUMERIC, "en_US.utf-8");
  
  if (locale) {
    kipmsg0(3, "Local memory used = %'ld KB\n", used >> 10);
  } else {
    kipmsg0(3, "Local memory used = %ld KB\n", used >> 10);
  }
    
  KI_Allreduce(&used, &totalUsed, 1, MPI_LONG, MPI_SUM, ki_cmm_world);
  if (ki_world_rank == 0) {
    if (locale) {
      pmsg(1, "World memory used = %'ld KB\t\t\n", totalUsed >> 10);
    } else {
      pmsg(1, "World memory used = %ld KB\t\t\n", totalUsed >> 10);
    }
  }

#ifdef KI_TRACK_MEMORY
  long peak = kiGetMaxMallocHeap();
  long totalPeak = 0;
  long maxPeak = 0;
  if (locale) {
    kipmsg0(3, "Local memory peak = %'ld KB\n", peak >> 10);
  } else {
    kipmsg0(3, "Local memory peak = %ld KB\n", peak >> 10);
  }

  KI_Allreduce(&peak, &totalPeak, 1, MPI_LONG, MPI_SUM, ki_cmm_world);
  KI_Allreduce(&peak, &maxPeak,   1, MPI_LONG, MPI_MAX, ki_cmm_world);
  if (ki_world_rank == 0) {
    if (locale) {
      pmsg(1, "World memory peak = %'ld KB (max %'ld KB/proc)\n", totalPeak >> 10, maxPeak >> 10);
    } else {
      pmsg(1, "World memory peak = %ld KB (max %'ld KB/proc)\n", totalPeak >> 10, maxPeak >> 10);
    }
  }
#endif

  kiDataFinalize();
  
  if (kiIsParallel()) kiDomainFinalize();
  else MPI_Comm_free(&ki_cmm_domain);

  pmsg(4, "Bye.\n");

  used = kiGetMallocUsed();
  totalUsed = 0;
  if (used > 0) {
    if (locale) {
      kipmsg(3, "Local memory unfreed = %'ld B\n", used);
    } else {
      kipmsg(3, "Local memory unfreed = %'ld B\n", used);
    }
  }
  
  KI_Allreduce(&used, &totalUsed, 1, MPI_LONG, MPI_SUM, ki_cmm_world);
  if (ki_world_rank == 0) {
    if (locale) {
      pm("World memory leak = %'ld bytes\n", totalUsed);
    } else {
      pm("World memory leak = %ld bytes\n", totalUsed);
    }
  }
  
  ki_time_end = MPI_Wtime();
  char timeStr[100];
  kiTimeToString(ki_time_end - ki_time_beg, timeStr);
  
  if (ki_world_rank == 0) kipmsg(1, "Time elapsed = %s\n", timeStr);

  MPI_Group_free(&ki_grp_world);
  MPI_Comm_free(&ki_cmm_world);
  MPI_Finalize();
}

void kiStart() {
  switch (ki_domain_id) {
    case KI_DOMAIN_USER:   kiUserStart();   break;
    case KI_DOMAIN_DEALER: kiDealerStart(); break;
    case KI_DOMAIN_FARMER: kiFarmerStart(); break;
  }
}

void kiStop() {
  KI_Barrier(ki_cmm_domain);
  if (ki_domain_id == KI_DOMAIN_USER && ki_domain_rank == 0) {
    int cmd = KI_CMD_STOP;
    int i;
    for (i = 0; i < ki_dealer_size; ++i) {
      int dealer = ki_top_dealer + i;
      KI_Send(&cmd, 1, MPI_INT, dealer, 0, ki_cmm_user_dealer);
    }
  }
}

void kiAbort(int error) {
  MPI_Abort(ki_cmm_world, error);
}

bool kiIsDomainRoot() {
  return (ki_domain_rank == 0 || !kiIsParallel());
}


/* Internal routines */
void kiDomainInit() {
  /* float ratio[3]; */
  /* float ratio[3] = {0.2, 0.2, 0.6}; */
  /* float ratio[3] = {0.0, 0.0, 1.0};  /\* 1 user, 1 dealer, n-2 farmers *\/ */
  float ratio[3] = {0.02, 0.02, 0.96};  /* 2% user, %2 dealer, 96% farmers */
  /* float ratio[3] = {0.5, 0.0, 0.5};  /\* 1 user, 1 dealer, n-2 farmers *\/ */
  int nNodes[3];
  int sumNodes[4] = {0};
  int* ranks[3];
  int i, j;
  
  assert(ki_world_size >= 3);
  /* if (ki_world_size >= 8) { */
  /*   ratio[0] = 0.1; */
  /*   ratio[1] = 0.1; */
  /*   ratio[2] = 0.8; */
  /* } */

  if (ki_user_ratio < 0.5) {
    ratio[0] = ki_user_ratio;
    ratio[1] = ki_user_ratio;
    ratio[2] = 1 - ki_user_ratio * 2;
  }
  
  for (i = 0; i < 3; ++i) {
    nNodes[i] = (i < 2) ? MAX(1, (int)(ki_world_size * ratio[i]))
      : ki_world_size - sumNodes[i];
    sumNodes[i+1] = nNodes[i] + sumNodes[i];
    ranks[i] = malloc(sizeof(int) * nNodes[i]);
    for (j = 0; j < nNodes[i]; ++j)
      ranks[i][j] = sumNodes[i] + j;
    if (ki_world_rank < sumNodes[i+1])
      if (ki_domain_id == KI_UNDEF) ki_domain_id = i;
  }

  MPI_Comm_split(ki_cmm_world, ki_domain_id, ki_world_rank, &ki_cmm_domain);
  MPI_Comm_group(ki_cmm_domain, &ki_grp_domain);
  MPI_Comm_rank(ki_cmm_domain, &ki_domain_rank);
  MPI_Comm_size(ki_cmm_domain, &ki_domain_size);
  assert(ki_domain_rank >= 0 && ki_domain_size > 0);

  switch (ki_domain_id) {
    case KI_DOMAIN_USER:
      sprintf(ki_domain_name, "User");
      ki_cmm_user = ki_cmm_domain;
      break;
    case KI_DOMAIN_DEALER:
      sprintf(ki_domain_name, "Dealer");
      ki_cmm_dealer = ki_cmm_domain;
      break;
    case KI_DOMAIN_FARMER:
      sprintf(ki_domain_name, "Farmer");
      ki_cmm_farmer = ki_cmm_domain;
      break;
    default:
      sprintf(ki_domain_name, "Undef");
  }

  ki_user_size   = nNodes[0];
  ki_dealer_size = nNodes[1];
  ki_farmer_size = nNodes[2];

  ki_top_dealer  = nNodes[0];    /* rank in ki_cmm_user_dealer */
  ki_top_farmer  = nNodes[1];    /* rank in ki_cmm_dealer_farmer */

  if (kiIsParallel() && ki_world_rank == 0)
    pmsg(3, "Users = %d, Dealers = %d, Farmers = %d\n", ki_user_size, ki_dealer_size, ki_farmer_size);

  kipmsg(4, "top dealer = %d\n", ki_top_dealer);
  kipmsg(4, "top farmer = %d\n", ki_top_farmer);

  MPI_Group_incl(ki_grp_world, nNodes[0], ranks[0], &ki_grp_user);
  MPI_Group_incl(ki_grp_world, nNodes[1], ranks[1], &ki_grp_dealer);
  MPI_Group_incl(ki_grp_world, nNodes[2], ranks[2], &ki_grp_farmer);

  MPI_Group_union(ki_grp_user, ki_grp_dealer, &ki_grp_user_dealer);
  MPI_Comm_create(ki_cmm_world, ki_grp_user_dealer, &ki_cmm_user_dealer);

  MPI_Group_union(ki_grp_dealer, ki_grp_farmer, &ki_grp_dealer_farmer);
  MPI_Comm_create(ki_cmm_world, ki_grp_dealer_farmer, &ki_cmm_dealer_farmer);

  for (i = 0; i < 3; ++i) free(ranks[i]);

  kipmsg(4, "Three domains set up.\n");
}

void kiDomainFinalize() {

  switch (ki_domain_id) {
    case KI_DOMAIN_USER:
      break;
    case KI_DOMAIN_DEALER:
      break;
    case KI_DOMAIN_FARMER:
      kiGetHashCollisionRatio();
      break;
  }

  if (ki_cmm_user_dealer != MPI_COMM_NULL) MPI_Comm_free(&ki_cmm_user_dealer);
  if (ki_cmm_dealer_farmer != MPI_COMM_NULL) MPI_Comm_free(&ki_cmm_dealer_farmer);
  
  MPI_Group_free(&ki_grp_user_dealer);
  MPI_Group_free(&ki_grp_dealer_farmer);

  MPI_Group_free(&ki_grp_user);
  MPI_Group_free(&ki_grp_dealer);
  MPI_Group_free(&ki_grp_farmer);

  MPI_Group_free(&ki_grp_domain);
  MPI_Comm_free(&ki_cmm_domain);

  kipmsg(4, "Domain finalized.\n");
}

void kiDataInit() {
  kiArenaInit(KI_ARENA_MAX);
  if (kiIsFarmerDomain()) {
    ki_seqs         = (alignment_t*)kimalloc(sizeof(alignment_t));
    ki_seqs->nSeq   = 0;
    ki_seqs->nPos   = 0;
    ki_seqs->names  = NULL;
    ki_seqs->seqs   = NULL;
    ki_seqs->flags  = NULL;
    ki_seqs->nSaved = 0;

    ki_hash = (hashtable_t**)kimalloc(sizeof(hashtable_t) * (KI_MAX_KMER_LEN+1));
    memset(ki_hash, 0, sizeof(hashtable_t) * (KI_MAX_KMER_LEN+1));
  }
}

void kiDataFinalize() {
  if (kiIsFarmerDomain()) {
    kiFreeAlignment(ki_seqs);
    kiFreeHashtables();
  }
  kiArenaFinalize(KI_ARENA_MAX);
}

void kiUserStart() {
}

void kiDealerStart() {
  MPI_Status status;
  char timestr[50];
  double time, wait;
  
  while (true) {
    int pos = 0;
    int cmd = 0;
    int cnt = 0;
    int src = 0;
    
    getCurrentTime(timestr);
    kipmsg(3, "[%s] Expecting cmd from user\n", timestr);

    /* relays user command message to top farmer */
    KI_Recv(ki_in_buf, KI_BUF_SIZE, MPI_PACKED, MPI_ANY_SOURCE, 0, ki_cmm_user_dealer, &status);
    MPI_Get_count(&status, MPI_PACKED, &cnt);
    src = status.MPI_SOURCE;
    getCurrentTime(timestr);
    time = MPI_Wtime();
    
    MPI_Unpack(ki_in_buf, KI_BUF_SIZE, &pos, &cmd, 1, MPI_INT, ki_cmm_world);
    kipmsg(3, "[%s] Received cmd = %d from %d\n", timestr, cmd, src);
    
    /* only top dealer sends STOP signal to top farmer */
    if (cmd != KI_CMD_STOP || ki_domain_rank == 0)
      KI_Send(ki_in_buf, cnt, MPI_PACKED, ki_top_farmer, 0, ki_cmm_dealer_farmer);

    if (cmd == KI_CMD_STOP) break;

    getCurrentTime(timestr);
    wait = MPI_Wtime() - time;
    kipmsg(3, "[%s] Forwarded cmd = %d to farmer, waited for %.2f s\n", timestr, cmd, wait);
    time = MPI_Wtime();

    /* receives function output message even if its length is zero */
    KI_Recv(ki_in_buf, KI_BUF_SIZE, MPI_PACKED, MPI_ANY_SOURCE, 0, ki_cmm_dealer_farmer, &status);
    MPI_Get_count(&status, MPI_PACKED, &cnt);
    getCurrentTime(timestr);
    wait = MPI_Wtime() - time;
    kipmsg(3, "[%s] Received results (count = %d) to be sent to user %d, compute time  = %.2f s\n", timestr, cnt, src, wait);

    KI_Send(ki_in_buf, cnt, MPI_PACKED, src, 0, ki_cmm_user_dealer);
  }
}

void kiFarmerStart() {
  MPI_Status status;
  
  while (true) {
    int pos = 0;
    int cmd = 0;
    int cnt = 0;
    int src = 0;

    if (kiIsDomainRoot()) {
      KI_Recv(ki_in_buf, KI_BUF_SIZE, MPI_PACKED, MPI_ANY_SOURCE, 0, ki_cmm_dealer_farmer, &status);
      MPI_Get_count(&status, MPI_PACKED, &cnt);
      src = status.MPI_SOURCE;
    }
    
    KI_Bcast(&cnt, 1, MPI_INT, 0, ki_cmm_domain);
    KI_Bcast(ki_in_buf, cnt, MPI_PACKED, 0, ki_cmm_domain);

    MPI_Unpack(ki_in_buf, KI_BUF_SIZE, &pos, &cmd, 1, MPI_INT, ki_cmm_world);
  
    pos = kiRunCommand(cmd);
    if (pos < 0) break;         /* KI_CMD_STOP */
    
    if (kiIsDomainRoot())
      KI_Send(ki_out_buf, pos, MPI_PACKED, src, 0, ki_cmm_dealer_farmer);
  }
}

bool kiIsUserDomain() {
  return (ki_domain_id == KI_DOMAIN_USER || !kiIsParallel());
}

bool kiIsFarmerDomain() {
  return (ki_domain_id == KI_DOMAIN_FARMER || !kiIsParallel());
}

void kiProcessArgs(int* argc, char** argv[]) {
  int i;
  for (i = 0; i < *argc; ++i) {
    if (ki_world_rank == 0) pmsg(3, "argv[%d] = '%s'\n", i, (*argv)[i]);
    if (strcmp((*argv)[i], "-v") == 0) { /* set debug output level based on command line */
      int pmLevel = 3;                   
      ++i;
      if ((i < *argc) && strlen((*argv)[i]) == 1 && (*argv)[i][0] >= '0' && (*argv)[i][0] <= '9') {
        pmLevel = (*argv)[i][0] - '0';
      } else --i;
      pmsg(4, "commandline supplied pmLevel = %d\n", pmLevel);
      pmSetLevel(pmLevel);             /* works unless NDEBUG is set */
    }
    else if (strcmp((*argv)[i], "-u") == 0) {
      ++i;
      if ((i < *argc) && (*argv)[i][0] >= '0' && (*argv)[i][0] <= '9') {
        ki_user_ratio = atof((*argv)[i]);
        if (ki_user_ratio >= 1) ki_user_ratio /= ki_world_size;
        if (ki_world_rank == 0) pmsg(2, "ki_user_ratio = %.3f\n", ki_user_ratio);
      }
    }
  }
}

bool kiIsParallel() {
  return (ki_world_size > 1);
}

bool kiIsMyTurn() {
  return (ki_domain_rank == ki_round_robin || !kiIsParallel());
}

int kiRunCommand(int cmd) {
  /* check if cmd is valid */
  if (cmd == KI_CMD_STOP) return -1;

  int outputSize = 0;
  ki_func_t func = ki_functions[cmd];
  if (func == NULL) {
    if (kiIsDomainRoot()) 
      fprintf(stderr, "Error: command %d not registered.\n", cmd);
  } else {
    outputSize = (*func)();
  }
  return outputSize;
}

int kiRandDealer() {
  return ki_top_dealer + (rand() % ki_dealer_size);
}

int kiRandFarmer() {
  return ki_top_farmer + (rand() % ki_farmer_size);
}

int kiNextFarmer() {
  if (kiIsParallel()) {
    ++ki_round_robin;
    if (ki_round_robin >= ki_domain_size)
      ki_round_robin = 0;
  } else {
    ki_round_robin = 0;
  }
  return ki_round_robin;
}

void kiRegisterCommand(int cmd, ki_func_t func) {
  assert(cmd > 0 && cmd < KI_NUM_FUNCS);
  ki_functions[cmd] = func;
}

void kiClearCommands() {
  memset(ki_functions, 0, sizeof(ki_func_t) * KI_NUM_FUNCS);
}

void kiPrintBuffer(char* buf, int len) {
  int i;
  char *p = buf;
  for (i = 0; i < len; ++i, ++p) {
    if (*p >= 0x20 && *p <= 0x7e) /* printable ascii */
      fprintf(stderr, "%c", *p);
    else
      fprintf(stderr, "%c", *p > 0 ? '.' : '|');
  }
  fprintf(stderr, "\n");
  /* fprintf(stderr, "%d\n", ki_world_rank); */
}

void kiPackArgs(void* buf, int bufSize, /*IN/OUT*/int* position, /*IN*/...) {
  int type, count;
  void* arg;

  va_list vl;
  va_start(vl, position);
  while ((type = va_arg(vl, int)) != KI_UNDEF) {
    if (type == KI_STRING) {
      arg = va_arg(vl, char*);
      count = strlen(arg) + 1;
      MPI_Pack(&count, 1, MPI_INT, buf, bufSize, position, ki_cmm_world);
      kipmsg(6, "[STR_SIZE=%d] ", count);
      type = MPI_CHAR;
    } else {
      if ((type & KI_VECTOR) == KI_VECTOR) {
        count = va_arg(vl, int);
        MPI_Pack(&count, 1, MPI_INT, buf, bufSize, position, ki_cmm_world);
        kipmsg(6, "[VEC_SIZE=%d] ", count);
        type ^= KI_VECTOR;      /* convert to MPI_Datatype */
      } else {
        count = 1;
      }
      arg = va_arg(vl, void*);
    }
    MPI_Pack(arg, count, type, buf, bufSize, position, ki_cmm_world);
    {                           /* debugging output */
      int i;
      pmsg(6, "{");
      for (i = 0; i < count; ++i) {
        if (type == MPI_CHAR) pmsg(6, "%c", ((char*)arg)[i]);
        else if (type == MPI_INT) pmsg(6, "%d ", ((int*)arg)[i]);
        else if (type == MPI_FLOAT) pmsg(6, "%f ", ((int*)arg)[i]);
      }
      pmsg(6, "}\n");
    }
  }
  va_end(vl);
}

void kiUnpackArgs(void* buf, int bufSize, int* position, /*OUT*/...) {
  int type, count;
  int* argcount;
  void* arg;

  va_list vl;
  va_start(vl, position);
  while ((type = va_arg(vl, int)) != KI_UNDEF) {
    if (type == KI_STRING) {
      MPI_Unpack(buf, bufSize, position, &count, 1, MPI_INT, ki_cmm_world);      
      type = MPI_CHAR;
    } else if ((type & KI_VECTOR) == KI_VECTOR) {
      argcount = va_arg(vl, int*);
      MPI_Unpack(buf, bufSize, position, argcount, 1, MPI_INT, ki_cmm_world);      
      count = *argcount;
      type ^= KI_VECTOR;        /* convert to MPI_Datatype */
    } else {
      count = 1;
    }
    arg = va_arg(vl, void*);
    MPI_Unpack(buf, bufSize, position, arg, count, type, ki_cmm_world);
  }
  va_end(vl);
}


/* Small int set implemented with array */
int_array_set_t* kiAllocIntArraySet() {
  int_array_set_t* set = (int_array_set_t*)kiArenaMalloc(KI_ARENA_TEMP, sizeof(int_array_set_t));
  set->size = 0;
  set->array = NULL;
  return set;
}

int_array_set_t* kiFreeIntArraySet(int_array_set_t* set) {
  if (KI_USE_ARENA) {
    kiArenaClear(KI_ARENA_TEMP);
    return NULL;
  }
    
  if (set->size > 0) {
    assert(set->array);
    kifree(set->array, sizeof(bool) * set->size);
  }
  kifree(set, sizeof(int_array_set_t));
  return NULL;
}

void kiIntArraySetAdd(int_array_set_t* set, int x) {
  assert(x < 6553500); /* upper bound beyond which we should use hashtable */
  /* assert(x < 655350); /\* upper bound beyond which we should use hashtable *\/ */
  /* assert(x < KI_OG_MAX_SIZE); /\* upper bound beyond which we should use hashtable *\/ */
  if (x >= set->size) {
    int newSize = x * 2;
    if (set->size == 0) {
      newSize = MAX(newSize, 1000);
      set->array = (bool*)kiArenaMalloc(KI_ARENA_TEMP, sizeof(bool) * newSize);
    } else {
      set->array = (bool*)kiArenaRealloc(KI_ARENA_TEMP, set->array, sizeof(bool)*set->size, sizeof(bool)*newSize, /*copy*/false);
    }
    memset(set->array + set->size, 0, sizeof(bool) * (newSize - set->size));
    set->size = newSize;
  }
  set->array[x] = true;
}

bool kiIntArraySetHas(int_array_set_t* set, int x) {
  if (set->size <= x) return false;
  return (set->array[x]);
}


/* Sequence manipulation */
void kiPackMatchingSeq(void* buf, int bufSize, int* position, /*IN*/char* name, char* seq, int flag) {
  /* FIXME: temporary */
  if (*position >= (bufSize / ki_domain_size) - 200) {
    kipmsg(4, "skip kiPackMatchingSeq to prevent recv buf overflow\n");
    return;
  }
  
  kiPackArgs(buf, bufSize, position, KI_STRING, name, KI_STRING, seq, KI_INT, &flag, KI_UNDEF);
  kipmsg(7, "name = '%s', seq = '%s', flag = %d --> pos = %d\n", name, seq, flag, *position);
}

void kiUnpackAlignment(void* buf, int bufSize, int* position, /*OUT*/alignment_t* aln) {
  char name[KI_SEQ_BUF_SIZE];
  char seq[KI_SEQ_BUF_SIZE];
  int flag;
  int pos = *position;
  while (*position - pos < bufSize) {
    kiUnpackArgs(buf, bufSize, position, KI_STRING, name, KI_STRING, seq, KI_INT, &flag, KI_UNDEF);
    kipmsg(7, "pos = %d <-- name = '%s', seq = '%s', flag = %d \n", *position, name, seq, flag);
    kiAlignmentAdd(aln, name, seq, /*copy*/true);
    aln->flags[aln->nSeq-1] = flag;
  }
}


/* Collective hashtable routines */
void kiFreeHashtables() {
  int i;
  for (i = 0; i <= KI_MAX_KMER_LEN; ++i) {
    ki_hash[i] = kiFreeHashtable(ki_hash[i]); /* returns NULL */
  } 
  ki_hash = kifree(ki_hash, (sizeof(hashtable_t) * (KI_MAX_KMER_LEN+1)));
}

void kiHashtableLoadSeq(int id) {
  assert(strlen(ki_seqs->seqs[id]) >= ki_kmer_len);

  hashtable_t* hash = ki_hash[ki_kmer_len];
  if (hash == NULL) {
    bool bProtein = kiIsProtein(ki_seqs->seqs[id]);
    hash = kiMakeHashtable(NULL, 0, /*use default nBuckets*/0, /*copy*/true, /*protein*/bProtein);
    ki_hash[ki_kmer_len] = hash;
  }
  if (ki_hash_mode != KI_HASH_NONE) {
    char* seq = ki_seqs->seqs[id];
    int i,len = strlen(seq);
    char *p, *q; 
    if (ki_hash[ki_kmer_len]->bProtein) {
      if (ki_hash_mode == KI_HASH_COMPLETE) {
        for (i = 0, p = seq; i+ki_kmer_len <= len; ++i, p++) {
          kiHashtableAddCopy(hash, p, ki_kmer_len,  id,   i);
        }
      } else if (ki_hash_mode == KI_HASH_ENDS) {
        kiHashtableAddCopy(hash, seq,  ki_kmer_len,  id,   0);
      }
    } else {
      if (strchr(seq, 'N') != NULL) return;
      char* comp = (char*)kiArenaMalloc(KI_ARENA_FARMER, len + 1); 
      kiReverseComplementSeq(seq, comp);
      pmsg(5, "seq  = %s\ncomp = %s\n", seq, comp);

      if (ki_hash_mode == KI_HASH_COMPLETE) {
        for (i = 0, p = seq, q = comp; i+ki_kmer_len <= len; ++i, p++, q++) {
          kiHashtableAddCopy(hash, p, ki_kmer_len,  id,   i);
          kiHashtableAddCopy(hash, q, ki_kmer_len, -id-1, i);
        }
      } else if (ki_hash_mode == KI_HASH_ENDS) {
        kiHashtableAddCopy(hash, seq,  ki_kmer_len,  id,   0);
        kiHashtableAddCopy(hash, comp, ki_kmer_len, -id-1, 0);
      }
      kiArenaFree(KI_ARENA_FARMER, comp, len + 1);
    }
  }
}

void kiHashtableLoadAllSeqs() {
  int i;
  for (i = 0; i < ki_seqs->nSeq; ++i) {
    if (ki_seqs->flags[i] != KI_SEQ_NEW) continue;
    kiHashtableLoadSeq(i);
  }
  /* kiHashtableOptimize(); */
  kiArenaClear(KI_ARENA_FARMER);
}


void kiHashtableLoadAllSeqsMinLength() {
  float ratio = kiGetHashLoadFactor(ki_hash[ki_kmer_len]); /* should be after hash is initialized */
  kipmsg0(2, "Hashtable load factor = 1 : %.2f \n", 1.0/ratio);

  long used1 = kiGetMallocUsed();

  kipmsg0(3, "Before loading seqs to hash table:\n");
  kiPrintLocalMemoryUsage();

  double timeStart = MPI_Wtime();

  int i, len, minLen, maxLen = KI_SEQ_BUF_SIZE;
  char** p = ki_seqs->seqs;
  for (i = 0; i < ki_seqs->nSeq; ++i, ++p) {
    len = strlen(*p);
    minLen = MIN(len, minLen);
    maxLen = MAX(len, maxLen);
  }
  int globalMin, globalMax;
  KI_Allreduce(&minLen, &globalMin, 1, MPI_INT, MPI_MIN, ki_cmm_domain);
  KI_Allreduce(&maxLen, &globalMax, 1, MPI_INT, MPI_MAX, ki_cmm_domain);
  kipmsg0(2, "Global read length: min = %d, max = %\n", globalMin, globalMax);

  ki_seqs->nPos = globalMin;
  
  p = ki_seqs->seqs;
  for (i = 0; i < ki_seqs->nSeq; ++i, ++p) {
    if (ki_seqs->flags[i] != KI_SEQ_NEW) continue;
    /* FIXME: memory leak */
    (*p)[globalMin] = '\0';
    kiHashtableLoadSeq(i);
  }
  kiHashtableOptimize(ki_hash[ki_kmer_len]);

  kiArenaClear(KI_ARENA_FARMER);

  KI_Barrier(ki_cmm_domain);
  
  double time = MPI_Wtime() - timeStart;
  char timeStr[100];
  kiTimeToString(time, timeStr);

  kipmsg0(3, "After loading seqs to hash table:\n");
  kiPrintLocalMemoryUsage();

  long used2 = kiGetMallocUsed();
  long used0 = kiGetHashTableMemory();
  long used  = used2 - used1 - used0;  
  kipmsg0(3, "Extra local memory for loading seqs to hash = %ld MB\n", used >> 20);
  
  float hashSizePerSeq  = (float)used / ki_seqs->nSeq;
  float hashSizePerBase = (float)hashSizePerSeq / ki_seqs->nPos;
  kipmsg0(2, "Memory for sequences in hash = %.1f bytes / seq, %.1f bytes / base\n", hashSizePerSeq, hashSizePerBase);

  kipmsg0(2, "Loaded seqs to hashtable in %s.\n", timeStr);
}

/* ONIT */
void kiHashtableLoadAllSeqsFixedLength() {
  float ratio = kiGetHashLoadFactor(ki_hash[ki_kmer_len]); /* should be after hash is initialized */
  kipmsg0(2, "Hashtable load factor = 1 : %.2f \n", 1.0/ratio);

  long used1 = kiGetMallocUsed();

  kipmsg0(3, "Before loading seqs to hash table:\n");
  kiPrintLocalMemoryUsage();

  double timeStart = MPI_Wtime();

  int i, len;
  int minLen = KI_SEQ_BUF_SIZE, maxLen = 0;
  char** p = ki_seqs->seqs;
  for (i = 0; i < ki_seqs->nSeq; ++i, ++p) {
    len = strlen(*p);
    minLen = MIN(len, minLen);
    maxLen = MAX(len, maxLen);
  }
  int globalMin, globalMax;
  KI_Allreduce(&minLen, &globalMin, 1, MPI_INT, MPI_MIN, ki_cmm_domain);
  KI_Allreduce(&maxLen, &globalMax, 1, MPI_INT, MPI_MAX, ki_cmm_domain);
  /* kipmsg0(2, "Global read length: min = %d, max = %d\n", globalMin, globalMax); */

  ki_seqs->nPos = globalMin;
  
  p = ki_seqs->seqs;
  for (i = 0; i < ki_seqs->nSeq; ++i, ++p) {
    if (ki_seqs->flags[i] != KI_SEQ_NEW) continue;
    /* FIXME: memory leak */
    (*p)[globalMin] = '\0';
    kiHashtableLoadSeq(i);
  }
  kiHashtableOptimize(ki_hash[ki_kmer_len]);

  kiArenaClear(KI_ARENA_FARMER);

  KI_Barrier(ki_cmm_domain);
  
  double time = MPI_Wtime() - timeStart;
  char timeStr[100];
  kiTimeToString(time, timeStr);

  kipmsg0(3, "After loading seqs to hash table:\n");
  kiPrintLocalMemoryUsage();

  long used2 = kiGetMallocUsed();
  long used0 = kiGetHashTableMemory();
  long used  = used2 - used1 - used0;  
  kipmsg0(3, "Extra local memory for loading seqs to hash = %ld MB\n", used >> 20);
  
  float hashSizePerSeq  = (float)used / ki_seqs->nSeq;
  float hashSizePerBase = (float)hashSizePerSeq / ki_seqs->nPos;
  kipmsg0(2, "Memory for sequences in hash = %.1f bytes / seq, %.1f bytes / base\n", hashSizePerSeq, hashSizePerBase);

  kipmsg0(2, "Loaded seqs to hashtable in %s.\n", timeStr);
}


/* Utility functions */
bool kiIsPowerOfTwo(int x) {
  return (x & (x - 1)) == 0;
}

/* Utility collective subroutines */
void kiResetProgress() {
  ki_clock_start = -1.;
}

void kiReportProgress(int name, char* info, long long count, long long total, int interval) {
  if (ki_clock_total < 0) return;

  if (ki_clock_start < 0 || name != ki_clock_name) {
    ki_clock_start    = MPI_Wtime();
    ki_clock_stop     = ki_clock_start;
    ki_clock_name     = name;
    ki_clock_interval = 0.;
    ki_clock_total    = 0LL;
  }
  if (total > 0 && ki_clock_total <= 0) {
    KI_Allreduce(&total, &ki_clock_total, 1, MPI_LONG_LONG_INT, MPI_SUM,  ki_cmm_domain);
  }
  KI_Reduce(&count, &ki_clock_count, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, ki_cmm_domain);

  if (!kiIsDomainRoot()) return;
  /* kipm("ki_clock_start = %f, ki_clock_stop = %f, ki_clock_interval = %f\n", ki_clock_start, ki_clock_stop, ki_clock_interval); */
  /* kipm("ki_clock_count = %d, ki_clock_total = %d\n", ki_clock_count, ki_clock_total); */
  
  char eol = interval < 10 ? '\r' : '\n';
  double now = MPI_Wtime();
  char timeStr[100];
  ki_clock_interval += (now - ki_clock_stop);
  ki_clock_stop = now;
  if (ki_clock_interval > (double)interval) {
    ki_clock_interval = .0;
    if (ki_clock_total > ki_clock_count) {
      float eta = (ki_clock_total-ki_clock_count) / (ki_clock_count*1./(ki_clock_stop-ki_clock_start));
      kiTimeToString(eta, timeStr);
      pm(">> %s = %lld / %lld  (%.1f/s)  ETA: %s %c", info, ki_clock_count, ki_clock_total, ki_clock_count*1./(ki_clock_stop-ki_clock_start), timeStr, eol);
    }
  } else if (ki_clock_total == ki_clock_count) {
    kiTimeToString(ki_clock_stop - ki_clock_start, timeStr);
    pm(">> %s = %lld / %lld  (%.1f/s) Time: %s \n", info, ki_clock_count, ki_clock_total, ki_clock_count*1./(ki_clock_stop-ki_clock_start), timeStr);
    ki_clock_total = -1;
  }
}


void kiReportRootProgress(int name, char* info, int count, int total, int interval) {
  if (!kiIsDomainRoot()) return;
  
  if (ki_clock_total < 0) return;
  
  if (ki_clock_start < 0 || name != ki_clock_name) {
    ki_clock_start    = MPI_Wtime();
    ki_clock_stop     = ki_clock_start;
    ki_clock_name     = name;
    ki_clock_interval = 0.;
    ki_clock_total    = 0;
  }

  ki_clock_total = total;
  ki_clock_count = count;
  
  char eol = interval < 10 ? '\r' : '\n';
  double now = MPI_Wtime();
  char timeStr[100];
  ki_clock_interval += (now - ki_clock_stop);
  ki_clock_stop = now;
  if (ki_clock_interval > (double)interval) {
    ki_clock_interval = .0;
    if (ki_clock_total > ki_clock_count) {
      float eta = (ki_clock_total-ki_clock_count) / (ki_clock_count*1./(ki_clock_stop-ki_clock_start));
      kiTimeToString(eta, timeStr);
      pm(">> %s = %d / %d (~ %.1f/s) x %d  ETA: %s %c", info, ki_clock_count, ki_clock_total, ki_clock_count*1./(ki_clock_stop-ki_clock_start), ki_domain_size, timeStr, eol);
    }
  } else if (ki_clock_total == ki_clock_count) {
    kiTimeToString(ki_clock_stop - ki_clock_start, timeStr);
    pm(">> %s = %d / %d (~ %.1f/s) x %d  Time: %s \n", info, ki_clock_count, ki_clock_total, ki_clock_count*1./(ki_clock_stop-ki_clock_start), ki_domain_size, timeStr);
    ki_clock_total = -1;
  }
}



/* Parallel routines implemented on farmers */
int kiGetKmer(char* kmer) {
  int nMatches = 0;
  int len = strlen(kmer);
  hashtable_t* hash = ki_hash[len];
  if (hash == NULL) {
    kipm("Hashtables not found for k = %d\n", len);
  } else {
    hashiterator_t hi = kiFindMatch(hash, kmer);
    nMatches = kiHashCount(hash, hi);
    kipmsg(4, "%d local matches found.\n", nMatches);
  }
  int totalMatches = 0;
  KI_Allreduce(&nMatches, &totalMatches, 1, MPI_INT, MPI_SUM, ki_cmm_domain);
  return totalMatches;
}

void kiExactAssembleOld() {
  int totalSeq = 0;
  KI_Allreduce(&(ki_seqs->nSeq), &totalSeq, 1, MPI_INT, MPI_SUM, ki_cmm_domain);
  pm0("totalSeq = %d\n", totalSeq);
  
  bool* seqFlags = kimalloc(sizeof(bool) * ki_seqs->nSeq);
  memset(seqFlags, 0, sizeof(bool) * ki_seqs->nSeq);

  char* seq   = kimalloc(sizeof(char) * KI_SEQ_BUF_SIZE);
  char* match = kimalloc(sizeof(char) * KI_SEQ_BUF_SIZE);
  char* kmer  = kimalloc(sizeof(char) * (ki_kmer_len+1));
  memset(kmer, 0, sizeof(char) * (ki_kmer_len+1));

  hashtable_t* seenH = kiMakeHashtable(NULL, 0, 8000, /*copy*/true, /*protein*/false);
  
  int i, j;
  int totalConsidered = 0;
  int cpu  = 0;
  int len = 0;
  while (cpu < ki_domain_size) {
    int more = 1;
    while (more > 0) {
      int si = 0;
      if (cpu == ki_domain_rank) {
        while (seqFlags[si] && si < ki_seqs->nSeq) ++si;
        if (si >= ki_seqs->nSeq) more = 0;
        else strcpy(seq, ki_seqs->seqs[si]);
      }
      KI_Bcast(&more, 1, MPI_INT, cpu, ki_cmm_domain);
      if (more == 0) break;
      totalConsidered++;

      KI_Bcast(seq, KI_SEQ_BUF_SIZE, MPI_CHAR, cpu, ki_cmm_domain);
      len = strlen(seq);
      /* kipm0("seq='%s'\n", seq); */
      
      hashtable_t* hash = ki_hash[ki_kmer_len];
      if (hash == NULL) {
        kipm("Hashtables not found for k = %d\n", len);
        return;
      }
      
      kiClearHashtable(seenH);

      for (i = 1; i+ki_kmer_len <= len && len < KI_SEQ_BUF_SIZE-ki_kmer_len; ++i) {
        memcpy(kmer, seq+i, ki_kmer_len);
        
        hashiterator_t hi;
        hi = kiFindMatch(seenH, kmer);
        if (kiHashCount(seenH, hi) > 0) break; /* avoid loops */
        kiHashtableAdd(seenH, kmer, si);

        /* kipm0("kmer='%s'\n", kmer); */
        hi = kiFindMatch(hash, kmer);
        int nMatches = kiHashCount(hash, hi);
        if (nMatches > 0) {
          origin_iterator_t originIt = kiHashOrigins(hash, hi);
          /* kipmsg(2, "i=%d: %d matches found.\n", i, nMatches); */
          /* kipmsg(2, "query = %s\n", seq+i); */
          for (j = 0; j < nMatches && len < KI_SEQ_BUF_SIZE-ki_kmer_len; ++j, originIt = originIt->next) {
            kiHashGetMatchingSeq(originIt, /*OUT*/match);
            /* kipmsg(2, "match = %s\n", match); */
            float diff = kiSeqCmp(seq+i, match);
            if (diff < 1e-3) {
              int matchLen = strlen(match);
              if (matchLen > len-i) {
                if (len + i+matchLen-len < KI_SEQ_BUF_SIZE) {
                  strcpy(seq+len, match+len-i);
                  len += i+matchLen-len;
                  seqFlags[kiAbsSeqId(originIt->origin.seqid)] = true;
                }
              }
            }
          }
        }
      }

      /* if (len > 500) kipmsg(2, "local contig = %s\n", seq); */
      if (len > 1000) kipmsg(2, "local contig len = %d\n", strlen(seq));
      
      

      if (cpu == ki_domain_rank)
        seqFlags[si] = true;
    }
    ++cpu;
  }
  pm0("totalConsidered = %d\n", totalConsidered);

  kiFreeHashtable(seenH);
  
  kifree(kmer, sizeof(char) * (ki_kmer_len+1));
  kifree(match, sizeof(char) * KI_SEQ_BUF_SIZE);
  kifree(seq, sizeof(char) * KI_SEQ_BUF_SIZE);
  kifree(seqFlags, sizeof(bool) * ki_seqs->nSeq);

}

void kiExtendSeqOld(/*IN/OUT*/char* seq, /*IN*/int sz, int beg, int root, /*OUT*/int* advance, int* leader) {
/* advance: advancement in number of chars;
   leader:  rank of node that achieved greatest advancement */  
  int len = strlen(seq) + 1;
  KI_Bcast(&len, 1, MPI_INT, root, ki_cmm_domain);
  KI_Bcast(seq, len, MPI_CHAR, root, ki_cmm_domain);
  int oldLen = len = len - 1;
  /* kipm("len = %d\n", len); */
  
  hashtable_t* hash = ki_hash[ki_kmer_len];
  if (hash == NULL) {
    kipm("Hashtables not found for k = %d\n", len);
    return;
  }

  char* match = kimalloc(sizeof(char) * KI_SEQ_BUF_SIZE);
  char* kmer  = kimalloc(sizeof(char) * (ki_kmer_len+1));
  memset(kmer, 0, sizeof(char) * (ki_kmer_len+1));

  hashtable_t* seenH = kiMakeHashtable(NULL, 0, 8000, /*copy*/true, /*protein*/false);

  hashiterator_t hi;
  int i, j;
  char* p = seq;
  for (i = 0; i < beg && i+ki_kmer_len < len; ++i, ++p) {
    memcpy(kmer, p, ki_kmer_len);
    kiHashtableAdd(seenH, kmer, 0); /* any seqid >= 0  */
  }
  for (; i+ki_kmer_len <= len && len < sz-ki_kmer_len-1; ++i) { 
    memcpy(kmer, seq+i, ki_kmer_len);
        
    hi = kiFindMatch(seenH, kmer);
    if (kiHashCount(seenH, hi) > 0) break; /* avoid loops */
    kiHashtableAdd(seenH, kmer, 0);

    /* kipm0("kmer='%s'\n", kmer); */
    hi = kiFindMatch(hash, kmer);
    int nMatches = kiHashCount(hash, hi);
    if (nMatches > 0) {
      origin_iterator_t originIt = kiHashOrigins(hash, hi);
      /* kipmsg(2, "i=%d: %d matches found.\n", i, nMatches); */
      /* kipmsg(2, "query = %s\n", seq+i); */
      for (j = 0; j < nMatches && len < sz-ki_kmer_len; ++j, originIt = originIt->next) {
        kiHashGetMatchingSeq(originIt, /*OUT*/match);
        /* kipmsg(2, "match = %s\n", match); */
        float diff = kiSeqCmp(seq+i, match);
        if (diff < 1e-3) {
          int matchLen = strlen(match);
          if (matchLen > len - i) {
            if (len + i+matchLen-len < sz) {
              strcpy(seq+len, match+len-i);
              len += i+matchLen-len;
              ki_seqs->flags[kiAbsSeqId(originIt->origin.seqid)] = KI_SEQ_USED;
            }
          }
        }
      }
    }
  }

  kiFreeHashtable(seenH);
  
  kifree(kmer,  sizeof(char) * (ki_kmer_len+1));
  kifree(match, sizeof(char) * KI_SEQ_BUF_SIZE);

  /* kipm("after len = %d, oldLen = %d\n", len, oldLen); */
  /* if (len > 1000) kipm("after len = %d\n", len); */

  struct {
    int value;
    int index;
  } in, out;
  
  in.value = len - oldLen;
  in.index = ki_domain_rank;
  
  KI_Allreduce(&in, &out, 1, MPI_2INT, MPI_MAXLOC, ki_cmm_domain);
  *advance = out.value;
  *leader  = out.index;

}

void kiExtendSeq(/*IN/OUT*/char* seq, /*IN*/int sz, int beg, int root, /*OUT*/int* advance, int* leader) {
/* advance: advancement in number of chars;
   leader:  rank of node that achieved greatest advancement */  
  int len = strlen(seq) + 1;
  KI_Bcast(&len, 1, MPI_INT, root, ki_cmm_domain);
  KI_Bcast(seq, len, MPI_CHAR, root, ki_cmm_domain);
  int oldLen = len = len - 1;
  /* kipm("len = %d\n", len); */
  
  hashtable_t* hash = ki_hash[ki_kmer_len];
  if (hash == NULL) {
    kipm("Hashtables not found for k = %d\n", len);
    return;
  }

  char* match = kimalloc(sizeof(char) * KI_SEQ_BUF_SIZE);
  char* kmer  = kimalloc(sizeof(char) * (ki_kmer_len+1));
  memset(kmer, 0, sizeof(char) * (ki_kmer_len+1));

  /* hashtable_t* seenH = kiMakeHashtable(NULL, 0, 8000, /\*copy*\/true, /\*protein*\/false); */

  /* hashiterator_t hi; */
  /* int i, j; */
  /* char* p = seq; */
  /* for (i = 0; i < beg && i+ki_kmer_len < len; ++i, ++p) { */
  /*   memcpy(kmer, p, ki_kmer_len); */
  /*   kiHashtableAdd(seenH, kmer, 0); /\* any seqid >= 0  *\/ */
  /* } */
  int hashSize = 8000;
  hashtable_t* seenH = kiMakeHashtable(NULL, 0, hashSize, /*copy*/true, /*protein*/false);

  hashiterator_t hi;
  int i = MAX(0, beg - hashSize/8), j;
  char* p = seq + i;
  for (; i < beg && i+ki_kmer_len < len; ++i, ++p) {
    memcpy(kmer, p, ki_kmer_len);
    kiHashtableAdd(seenH, kmer, 0); /* any seqid >= 0  */
  }

  for (; i+ki_kmer_len <= len && len-oldLen < hashSize/2 && len < sz-ki_kmer_len-1; ++i) { 
    /* kipm("i=%d\n", i); */
    memcpy(kmer, seq+i, ki_kmer_len);
     
    hi = kiFindMatch(seenH, kmer);
    if (kiHashCount(seenH, hi) > 0) break; /* avoid loops */
    kiHashtableAdd(seenH, kmer, 0);
 
    /* kipm0("kmer='%s'\n", kmer); */
    hi = kiFindMatch(hash, kmer);
    int nMatches = kiHashCount(hash, hi);
    if (nMatches > 0) {
      origin_iterator_t originIt = kiHashOrigins(hash, hi);
      /* kipmsg(2, "i=%d: %d matches found.\n", i, nMatches); */
      /* kipmsg(2, "query = %s\n", seq+i); */
      for (j = 0; j < nMatches && len < sz-ki_kmer_len; ++j, originIt = originIt->next) {
        /* ki_seqs->flags[kiAbsSeqId(originIt->origin.seqid)] = 1; */
        kiHashGetMatchingSeq(originIt, /*OUT*/match);
        /* kipmsg(2, "match = %s\n", match); */
        float diff = kiSeqCmp(seq+i, match);
        if (diff < 1e-3) {
          int matchLen = strlen(match);
          if (matchLen > len - i) {
            if (len + i+matchLen-len < sz) {
              strcpy(seq+len, match+len-i);
              len += i+matchLen-len;
              ki_seqs->flags[kiAbsSeqId(originIt->origin.seqid)] = KI_SEQ_USED;
              /* kipm("len = %d\n", len); */
            }
          }
        }
      }
    }
  }

  kiFreeHashtable(seenH);
  
  kifree(kmer,  sizeof(char) * (ki_kmer_len+1));
  kifree(match, sizeof(char) * KI_SEQ_BUF_SIZE);

  /* kipm("after len = %d, oldLen = %d\n", len, oldLen); */
  /* if (len > 1000) kipm("after len = %d\n", len); */

  struct {
    int value;
    int index;
  } in, out;
  
  in.value = len - oldLen;
  in.index = ki_domain_rank;
  
  KI_Allreduce(&in, &out, 1, MPI_2INT, MPI_MAXLOC, ki_cmm_domain);
  *advance = out.value;
  *leader  = out.index;

}

void kiExactAssemble(char* outputName) {
  char contigFile1[200];
  char contigFile2[200];
  char contigFile3[200];

  sprintf(contigFile1, "%s.short.contig",  outputName);
  sprintf(contigFile2, "%s.medium.contig", outputName);
  sprintf(contigFile3, "%s.long.contig",   outputName);

  FILE *fp1, *fp2, *fp3;
  fp1 = fopen(contigFile1, "w");
  fp2 = fopen(contigFile2, "w");
  fp3 = fopen(contigFile3, "w");

  if (fp1 == NULL || fp2 == NULL || fp3 == NULL) 
    fprintf(stderr, "Error: cannot write to output contig file\n");

  char seqName[200];

  int totalSeq = 0;
  KI_Allreduce(&(ki_seqs->nSeq), &totalSeq, 1, MPI_INT, MPI_SUM, ki_cmm_domain);
  pm0("totalSeq = %d\n", totalSeq);

  double time0 = MPI_Wtime();
  double time;
  
  char* seq = kimalloc(sizeof(char) * KI_CONTIG_SIZE);

  int totalConsidered = 0;
  int cpu  = 0;
  while (cpu < ki_domain_size) {
    int more = 1;
    while (more > 0) {
      int si = 0;
      if (cpu == ki_domain_rank) {
        while (ki_seqs->flags[si] != KI_SEQ_NEW && si < ki_seqs->nSeq) ++si;
        if (si >= ki_seqs->nSeq) more = 0;
        else strcpy(seq, ki_seqs->seqs[si]);
      }
      KI_Bcast(&more, 1, MPI_INT, cpu, ki_cmm_domain);
      if (more == 0) break;
      totalConsidered++;

      int advance = -1;
      int leader  = cpu;
      int beg     = 1;
      /* int loopi = 0; */
      while (advance != 0) {
        kiExtendSeq(seq, KI_CONTIG_SIZE, beg, leader, &advance, &leader);
        beg += advance;
        /* if (advance > 0) kipm0("beg = %d, advance = %d, leader = %d\n", beg, advance, leader); */
        /* kipm0("loopi = %d, advance = %d\n", loopi++, advance); */
        /* kipm0("+ beg = %d, advance = %d, seqlen = %d\n", beg, advance, strlen(seq)); */
      }
      /* kipm0("+ beg = %d, advance = %d, seqlen = %d\n", beg, advance, strlen(seq)); */
      /* kipm0("seq+  = %s\n", seq); */
      
      /* kiAbort(1); */

      int len = strlen(seq);
      /* int oldLen = len - advance + 1; */
      char* copy = kimemdup(seq, len+1);
      kiReverseComplementSeq(copy, seq);
      kifree(copy, len+1);
      /* kipm0("seq-  = %s\n", seq); */
      
      advance = -1;
      leader = 0;
      beg = beg;
      while (advance != 0) {
        kiExtendSeq(seq, KI_CONTIG_SIZE, beg, leader, &advance, &leader);
        beg += advance;
        /* if (advance > 0) kipm0("beg = %d, advance = %d, leader = %d\n", beg, advance, leader); */
        /* kipm0("loopi = %d, advance = %d\n", loopi++, advance); */
        /* kipm0("- beg = %d, advance = %d, seqlen = %d\n", beg, advance, strlen(seq)); */
      }
      /* kipm0("- beg = %d, advance = %d, seqlen = %d\n", beg, advance, strlen(seq)); */
      /* kipm0("seq-- = %s\n\n", seq); */
      

      /* /\* if (len > 500) kipmsg(2, "local contig = %s\n", seq); *\/ */
      /* if (len > 1000) kipmsg(2, "local contig len = %d\n", strlen(seq)); */
      int seqlen = strlen(seq);
      /* if (seqlen > 1500) kipm0("seqlen = %d\n", seqlen); */
      if (seqlen >= 300) {
        FILE* fp = (seqlen < 600)  ? fp1 :
                   (seqlen < 2000) ? fp2 : fp3;
        strcpy(seqName, ki_seqs->names[si]);
        KI_Bcast(seqName, 200, MPI_CHAR, cpu, ki_cmm_domain);
        fprintf(fp, ">%s %d\n", seqName, cpu);
        fprintf(fp, "%s\n", seq);
      }

      if (totalConsidered % 1000 == 10) {
        time = MPI_Wtime() - time0;
        pm0("Expanded %d seqs in %.0f seconds (%.1f seq/s)\n", totalConsidered, time, totalConsidered/time);
      }
      
      if (cpu == ki_domain_rank)
        ki_seqs->flags[si] = KI_SEQ_USED;
    }
    ++cpu;
  }
  pm0("totalConsidered = %d\n", totalConsidered);

  kifree(seq, sizeof(char) * KI_CONTIG_SIZE);

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}

int kiReadFastaString(char* buf, alignment_t* seqs) {
  char readBuf[KI_SEQ_BUF_SIZE];
  /* char* nameStop = " ,\t\r\n"; */
  char* nameStop = "\t\r\n";
  char *name, *seq;
  char *pp, *p, *qq, *q;
  int nSeqs = 0;

  pp = buf;
  while (*pp != '>' && *pp != '\0') ++pp;

  while (*pp != '\0') {
    /* move pp to start of name */
    ++pp;

    /* move qq to start of seq */
    qq = pp+1;
    while (*qq != '\r' && *qq != '\n') ++qq;
    *qq = '\0';
    while (!((*qq >= 'a' && *qq <= 'z') || (*qq >= 'A' && *qq <= 'Z'))) ++qq;

    /* truncate the name */
    for (p = pp; *p != '\0'; ++p) {
      for (q = nameStop; *q != '\0'; ++q) {
        if (*p == *q) {
          *p = '\0';
          break;
        }
      }
      if (*p == '\0') break;
    }

    /* handles multiline fasta */
    for (q = readBuf, qq++; *qq != '\0' && *qq != '>'; ++qq) {
      if (*qq != '\r' && *qq != '\n' && *qq != '\t' && *qq != ' ' && *qq != '\0') {
        *(q++) = *qq;
      }
    }
    *q = '\0';
        
    name = (char*)kiArenaMemdup(KI_ARENA_DEFAULT, pp, p-pp+1); /* TODO: consider using a dedicated arena? */
    seq  = (char*)kiArenaMemdup(KI_ARENA_DEFAULT, readBuf, q-readBuf+1);

    kipmsg(7, "name = %s, len = %d, seq = %s\n", name, strlen(seq), seq);

    kiAlignmentAdd(seqs, name, seq, /*copy*/false);
    nSeqs++;

    /* position next pp */
    pp = qq;
    assert(*pp == '\0' || *pp == '>');
  }
  
  return nSeqs;
  
}

int kiReadFastaShared(char* fileName, alignment_t* seqs) {
  int nSeqs;
  int bufSize = KI_BUF_SIZE;
  char* buf = kimalloc(bufSize);

  if (kiIsDomainRoot()) {
    FILE* fp = fopen(fileName, "r");
    kipmsg(6, "fp = %x  fileName = %s\n", fp, fileName);
    int count = fread(buf, 1, bufSize, fp);
    if (count == bufSize) {
      fprintf(stderr, "Profile too big: %s\n", fileName);
      kiAbort(-1);
    }
    kipmsg(6, "read count = %d\n", count);
  }
  KI_Bcast(buf, bufSize, MPI_CHAR, 0, ki_cmm_domain);

  nSeqs = kiReadFastaString(buf, seqs);
  pmsg0(2, "Each farmer read %d seqs from '%s'\n", nSeqs, fileName);
  kifree(buf, bufSize);
  return nSeqs;
}


bool kiIsFastq(char* fileName) {
  int isFastq = 0;
  if (kiIsDomainRoot()) {
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) {
      fprintf(stderr, "Could not open trainning input file %s.\n", fileName);
      kiAbort(1);
    }
    char buf[KI_NAME_BUF_SIZE] = "";
    if (fgets(buf, sizeof(buf), fp) != NULL) {
      isFastq = (buf[0] == '@');
    }
    fclose(fp);
  }
  KI_Bcast(&isFastq, 1, MPI_INT, 0, ki_cmm_domain);
  kipmsg(7, "isFastq = %d\n", isFastq);
  return isFastq;
}


long long kiReadFastqParallel(char* fileName) {
  KI_File fh;
  MPI_Offset fileSize;
  MPI_Offset beg;
  MPI_Status status;
  
  int rc;
  char* buf;  
  /* char* nameStop   = "(),: \t\r\n"; */
  char* nameStop   = " ,\t\r\n";
  char *name, *seq;
  char *p, *q, *pp, *qq;
  /* char readBuf[KI_SEQ_BUF_SIZE]; */
  int chunk;
  double timeStart = MPI_Wtime();

  rc = KI_File_open(fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (rc != MPI_SUCCESS) {
    if (kiIsDomainRoot()) fprintf(stderr, "Error: unable to open file '%s'\n", fileName);
    kiAbort(1);
  }

  KI_File_get_size(fh, &fileSize);
  kipmsg0(5, "fileSize = %ld\n", fileSize);

  buf = (char*)kimalloc(KI_BUF_SIZE);

  chunk = fileSize / ki_domain_size;
  beg = ki_domain_rank * chunk;
  if (ki_domain_rank == ki_domain_size - 1)
    chunk = fileSize - beg;

  kipmsg(6, "Reading %ld: %ld\n", beg, chunk);
  
  int   maxEntrySize = 800; 
  int   bytesToRead  = MIN(chunk + maxEntrySize, fileSize - beg);
  int   bufSize      = MIN(bytesToRead, KI_BUF_SIZE-1);
  char* pos          = buf;
  bool  lastChunk    = (fileSize - beg <= chunk);

  int moveSize = 0;
  int count = 0;
  int totalCount = 0;
  long long nSeqsAdded = 0;

  if (chunk > 0 && beg < fileSize) {
    KI_File_seek(fh, beg, MPI_SEEK_SET);
    pos = buf;
    while (bytesToRead > 0) {
      count = MIN(bytesToRead, (bufSize - moveSize));
      assert(count > 0);
      int countCopy = count;
      int ret = KI_File_read(fh, pos, count, MPI_BYTE, &status);
      if (ret != MPI_SUCCESS) {
        fprintf(stderr, "MPI_read error\n");
      }
      if (count != status.count) {
        kipm("bytesToRead(l) = %ld, bytesToRead(d) = %d\n", bytesToRead, bytesToRead);
        kipm("bufSize = %d, moveSize = %d\n", bufSize, moveSize);
        kipm("countCopy = %d, count = %d, status.count = %d\n", countCopy, count, status.count);
      }
      assert(count == status.count);
      bytesToRead -= count;

      if (lastChunk && bytesToRead == 0)
        pos[count] = '@'; /* add fake fastq to mark the end of fastq */

      pp = buf;
      bool foundStart = 0;
      while (!foundStart) {
        while (*pp != '@') ++pp;
        for (qq = pp+1; *qq != '\r' && *qq != '\n'; ++qq);
        for (qq++;      *qq != '\r' && *qq != '\n'; ++qq);
        if (*(++qq) == '+')
          foundStart = 1;
        else
          for (pp++; *pp != '\r' && *pp != '\n'; ++pp);
      }

      while (true) {
        /* move pp to start of name */
        ++pp;

        /* move qq to start of seq */
        qq = pp+1;              
        while (*qq != '\r' && *qq != '\n') ++qq;
        *qq = '\0';
        while (!((*qq >= 'a' && *qq <= 'z') || (*qq >= 'A' && *qq <= 'Z'))) ++qq;

        /* truncate the name */
        for (p = pp; *p != '\0'; ++p) {
          for (q = nameStop; *q != '\0'; ++q) {
            if (*p == *q) {
              *p = '\0';
              break;
            }
          }
          if (*p == '\0') break; 
        }

        /* handles multiline fasta */
        for (q = qq; *q != '\r' && *q != '\n'; ++q);
        *q = '\0';
        
        name = (char*)kiArenaMemdup(KI_ARENA_SEQS, pp, p-pp+1);
        seq  = (char*)kiArenaMemdup(KI_ARENA_SEQS, qq, q-qq+1);

        kipmsg(8, "name = %s, len = %d, seq = %s\n", name, strlen(seq), seq);

        kiAlignmentAdd(ki_seqs, name, seq, /*copy*/false);
        nSeqsAdded++;

        /* position next pp */
        for (pp = q; *pp != '\r' && *pp != '\n'; ++pp); /* skip + line */
        for (pp++;   *pp != '\r' && *pp != '\n'; ++pp); /* skip quality line */

        pp++;
        assert(*pp == '@');

        if (totalCount + (pp - pos) > chunk) {
          bytesToRead = 0;
          break; /* done my share */
        }

        /* look ahead to see another read is possible */
        p = pp+1;
        int neol = 0;
        while ((p <= pos+count) && neol < 4) {
          if (*p == '\r' || *p == '\n') neol++;
          ++p;
        }
        
        if (p > pos+count) {
          moveSize = pos+count - pp;
          memcpy(buf, pp, moveSize);
          pos = buf + moveSize;
          break;
        }

        if (*p != '@') {
          kipm("lastChunk = %d, bytesToRead = %d\n", lastChunk, bytesToRead);
          kipm("wrong: p - (pos+count) = %d, count = %d, *p = '%c'\n", p - (pos+count), count, *p);
          kipm("*(p-1) = '%c'\n", *(p-1));
          
          /* char tmpStr[1002]; */
          /* char* beg = (p-buf) > 500 ? p-500 : buf; */
          /* char* end = (pos+count-p) > 500 ? p+500 : pos+count; */
          
          *(p+300) = '\0';// *(p-1) = 'B';
          kipm("Context='%s'\n", p-2);

          /* memset(tmpStr, 0, 1001); */
          /* memcpy(tmpStr, beg, end-beg); */
          /* kipm("Context='%s'\n", tmpStr); */
        }
        
        assert(*p == '@');

      }
      totalCount += count;
    }
  }
  kifree(buf, KI_BUF_SIZE);
  
  KI_File_close(&fh);

  int i;
  for (i = 0; i < ki_seqs->nSeq; i++)
    kiUpcaseSeq(ki_seqs->seqs[i]);

  double readTime = MPI_Wtime() - timeStart;
  char timeStr[100];
  kiTimeToString(readTime, timeStr);

  long long addedSum = 0;
  KI_Allreduce(&nSeqsAdded, &addedSum, 1, MPI_LONG_LONG, MPI_SUM, ki_cmm_domain);
  kipmsg0(2, "All nodes: read %d sequences from fastq (%s).\n", addedSum, timeStr);

  return addedSum;
}

long long kiReadFastaParallel(char* fileName) {
  KI_File fh;
  MPI_Offset fileSize;
  MPI_Offset beg;
  MPI_Status status;
  
  int rc;
  char* buf;  
  /* char* nameStop   = "(),: \t\r\n"; */
  char* nameStop   = " ,\t\r\n";
  char *name, *seq;
  char *p, *q, *pp, *qq;
  char readBuf[KI_SEQ_BUF_SIZE];
  long chunk;
  double timeStart = MPI_Wtime();

  rc = KI_File_open(fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (rc != MPI_SUCCESS) {
    if (kiIsDomainRoot()) fprintf(stderr, "Error: unable to open file '%s'\n", fileName);
    kiAbort(1);
  }

  KI_File_get_size(fh, &fileSize);
  kipmsg0(5, "fileSize = %ld\n", fileSize);

  buf = (char*)kimalloc(KI_BUF_SIZE);

  chunk = fileSize / ki_domain_size;
  beg = ki_domain_rank * chunk;
  if (ki_domain_rank == ki_domain_size - 1)
    chunk = fileSize - beg;

  kipmsg(6, "Reading %ld: %ld\n", beg, chunk);
  
  int   maxEntrySize = 800; 
  long  bytesToRead  = MIN(chunk + maxEntrySize, fileSize - beg);
  int   bufSize      = MIN(bytesToRead, KI_BUF_SIZE-1);
  char* pos          = buf;
  bool  lastChunk    = (fileSize - beg <= chunk);

  int moveSize = 0;
  int count = 0;
  int totalCount = 0;
  long long nSeqsAdded = 0;

  if (beg < fileSize) {
    KI_File_seek(fh, beg, MPI_SEEK_SET);
    pos = buf;
    while (bytesToRead > 0) {
      count = MIN(bytesToRead, bufSize - moveSize);
      KI_File_read(fh, pos, count, MPI_BYTE, &status);
      assert(count == status.count);
      bytesToRead -= count;

      if (lastChunk && bytesToRead == 0)
        pos[count] = '>'; /* add fake fasta to mark the end of fasta */

      pp = buf;
      while (*pp != '>') ++pp;

      while (true) {
        /* move pp to start of name */
        ++pp;

        /* move qq to start of seq */
        qq = pp+1;              
        while (*qq != '\r' && *qq != '\n') ++qq;
        *qq = '\0';
        while (!((*qq >= 'a' && *qq <= 'z') || (*qq >= 'A' && *qq <= 'Z'))) ++qq;

        /* truncate the name */
        for (p = pp; *p != '\0'; ++p) {
          for (q = nameStop; *q != '\0'; ++q) {
            if (*p == *q) {
              *p = '\0';
              break;
            }
          }
          if (*p == '\0') break; 
        }

        /* handles multiline fasta */
        int nChar = 0;
        for (q = readBuf; *qq != '>' && nChar < KI_SEQ_BUF_SIZE-1; ++qq, ++nChar) {
          if (*qq != '\r' && *qq != '\n' && *qq != '\t' && *qq != ' ' && *qq != '\0') {
            *(q++) = *qq;
          }
        }
        if (*qq != '>') {
          kipmsg(5, "Seq longer than limit (%d): %s\n", KI_SEQ_BUF_SIZE, pp);
          while (*qq != '>') ++qq;
        }
        *q = '\0';
        
        name = (char*)kiArenaMemdup(KI_ARENA_SEQS, pp, p-pp+1);
        seq  = (char*)kiArenaMemdup(KI_ARENA_SEQS, readBuf, q-readBuf+1);

        kipmsg(8, "name = %s, len = %d, seq = %s\n", name, strlen(seq), seq);

        kiAlignmentAdd(ki_seqs, name, seq, /*copy*/false);
        nSeqsAdded++;

        /* position next pp */
        pp = qq;
        assert(*pp == '>');

        if (totalCount + (pp - pos) > chunk) {
          bytesToRead = 0;
          break; /* done my share */
        }

        /* look ahead to see if another read is possible */
        p = pp+1;
        while ((p <= pos+count) && *p != '>') ++p;
        if (p > pos+count) {
          moveSize = pos+count - pp;
          memcpy(buf, pp, moveSize);
          pos = buf + moveSize;
          break;
        }

      }
      totalCount += count;
    }
  }
  kifree(buf, KI_BUF_SIZE);
  
  KI_File_close(&fh);

  int i;
  for (i = 0; i < ki_seqs->nSeq; i++)
    kiUpcaseSeq(ki_seqs->seqs[i]);

  double readTime = MPI_Wtime() - timeStart;
  char timeStr[100];
  kiTimeToString(readTime, timeStr);

  long long addedSum = 0;
  KI_Allreduce(&nSeqsAdded, &addedSum, 1, MPI_LONG_LONG, MPI_SUM, ki_cmm_domain);
  kipmsg0(2, "All nodes: read %ld sequences from fasta (%s).\n", addedSum, timeStr);

  return addedSum;
  
  /* http://kb.iu.edu/data/aqpe.html */
}

long long kiReadFastaSemiParallel(char* fileName) {
  char buf[KI_SEQ_BUF_SIZE*100] = "";
  int   chunk      = 100;
  int   nContSeqs  = 16;
  char* nameStop   = "(),: \t\r\n";
  char* seqSkip    = " \t\r\n";
  FILE* fp         = NULL;
  char* pp         = NULL;
  char* qq         = NULL;
  char* nextLine   = NULL;
  char* str        = NULL;
  int   bufSize    = 0;
  int   more       = 1;
  int   ith        = 0;

  long long nSeqsRead  = 0;
  long long nSeqsAdded = 0;

  int minLen = -1;
  int maxLen = -1;
  
  double timeStart = MPI_Wtime();
  
  if (kiIsDomainRoot()) {
    fp = fopen(fileName, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: cannot read %s\n", fileName);
      kiAbort(1);
    }
    if (fgets(buf, sizeof(buf), fp) == NULL) {
      fprintf(stderr, "Error: reading header line\n");
      kiAbort(1);
    }
    qq = buf + (strlen(buf) + 1); 
  }

  while (more != 0) {
    if (kiIsDomainRoot()) {
      while ((str = fgets(qq, sizeof(buf)-bufSize, fp)) != NULL) {
        if (*qq == '>') {
          nSeqsRead++;
          kipmsg(5, "nSeqsRead = %d\n", nSeqsRead);
          if (nSeqsRead % chunk == 0) break;
        }
        qq += (strlen(qq) + 1);
        bufSize = qq - buf;
      } 
      if (str == NULL) 
        more = 0;
    }
    KI_Bcast(&bufSize, 1, MPI_INT, 0, ki_cmm_domain);
    kipmsg0(5, "bufSize = %d\n", bufSize);
    KI_Bcast(buf, bufSize, MPI_CHAR, 0, ki_cmm_domain);
    kipmsg0(5, "buf = '%s'\n", buf);  
    KI_Bcast(&more, 1, MPI_INT, 0, ki_cmm_domain);

    pp = buf;
    while (pp - buf < bufSize) {
      /* loop over lines */
      nextLine = pp + (strlen(pp) + 1);
      if (*pp == '>') {
        /* paired-end reads better be handled by the same farmer */
        if (ith++ % nContSeqs == 0) kiNextFarmer();
        kipmsg(5, "ith = %d, round_robin = %d\n", ith, ki_round_robin);
        if (kiIsMyTurn()) {
          kipmsg(5, "owner = %d\n", ki_round_robin);
          /* truncate the name */
          char *p, *q;
          for (p = pp+1; *p != '\0'; p++) {
            for (q = nameStop; *q != '\0'; q++) {
              if (*p == *q) {
                *p = '\0';
                break;
              }
            }
            if (*p == '\0') break;
          }
          /* allocate space for another sequence */
          nSeqsAdded++;
          kiAlignmentGrow(ki_seqs);
          
          ki_seqs->names[ki_seqs->nSeq-1] = (char*)kimemdup(pp+1, p-pp); /* p-pp == strlen(pp) */
          ki_seqs->seqs[ki_seqs->nSeq-1]  = NULL;
          ki_seqs->flags[ki_seqs->nSeq-1] = KI_SEQ_NEW;
          kipmsg(5, "name = '%s'\n", ki_seqs->names[ki_seqs->nSeq-1]);
        }
      } else {
        if (kiIsMyTurn()) {
          /* count non-space characters and append to sequence */          
          int nKeep = 0;
          char *p, *q;
          for (p = pp; *p != '\0'; p++) {
            for (q = seqSkip; *q != '\0'; q++) {
              if (*p == *q)
                break;
            }
            if (*p != *q)
              nKeep++;
          }
          int nOld = (ki_seqs->seqs[ki_seqs->nSeq-1] == NULL) ? 0 : strlen(ki_seqs->seqs[ki_seqs->nSeq-1]);
          ki_seqs->seqs[ki_seqs->nSeq-1] = (char*)kirealloc(ki_seqs->seqs[ki_seqs->nSeq-1], nOld, nOld+nKeep+1, /*copy*/false);
          if (nOld+nKeep > ki_seqs->nPos)
            ki_seqs->nPos = nOld + nKeep;
          char *out = ki_seqs->seqs[ki_seqs->nSeq-1] + nOld;
          for (p = pp; *p != '\0'; p++) {
            for (q = seqSkip; *q != '\0'; q++) {
              if (*p == *q)
                break;
            }
            if (*p != *q) {
              *out = *p;
              out++;
            }
          }
          assert(out - ki_seqs->seqs[ki_seqs->nSeq-1] == nKeep + nOld);
          *out = '\0';
          kipmsg(5, "seq  = '%s'\n", ki_seqs->seqs[ki_seqs->nSeq-1]);

          int len = strlen(ki_seqs->seqs[ki_seqs->nSeq-1]);
          if (minLen == -1 || len < minLen) minLen = len;
          if (maxLen == -1 || len > maxLen) maxLen = len;
        }
      }
      pp = nextLine;
    } 
    if (kiIsDomainRoot()) {
      if (more != 0) {
        bufSize = strlen(qq) + 1;
        memcpy(buf, qq, bufSize);
        qq = buf + bufSize;
      }
    }
  }
  kipmsg(3, "Read %d sequences from fasta.\n", nSeqsAdded);
  kipmsg(3, "minLen = %d, maxLen = %d\n", minLen, maxLen);

  if (kiIsDomainRoot()) {
    fclose(fp);
  }

  int i;
  for (i = 0; i < ki_seqs->nSeq; i++) 
    kiUpcaseSeq(ki_seqs->seqs[i]);

  double readTime = MPI_Wtime() - timeStart;
  char timeStr[100];
  kiTimeToString(readTime, timeStr);

  long long addedSum = 0;
  KI_Allreduce(&nSeqsAdded, &addedSum, 1, MPI_LONG_LONG, MPI_SUM, ki_cmm_domain);
  kipmsg0(2, "All nodes: read %d sequences from fasta (%s).\n", addedSum, timeStr);
  return addedSum;
}


/* Overlap graph */
ogvertex_t* kiAllocOgVertex() {
  if (KI_USE_ARENA)
    return (ogvertex_t*)kiArenaMalloc(KI_ARENA_DEFAULT, sizeof(ogvertex_t));
    
  ogvertex_t* v = (ogvertex_t*)kimalloc(sizeof(ogvertex_t));
  v->occr   = 0;
  v->offset = 0;
  v->flag   = 0;
  v->id     = 0;
  v->name   = NULL;
  v->seq    = NULL;
  int i;
  for (i = 0; i < 4; ++i) {
    v->ies[i] = NULL;
    v->oes[i] = NULL;
  }
  return v;
}

ogedge_t* kiAllocOgEdge() {
  if (KI_USE_ARENA)
    return (ogedge_t*)kiArenaMalloc(KI_ARENA_DEFAULT, sizeof(ogedge_t));
  
  ogedge_t* e = (ogedge_t*)kimalloc(sizeof(ogedge_t));
  e->ext      = NULL;
  e->v        = NULL;
  e->weight   = .0;
  return e;
}

ogvertex_t* kiFreeOgVertex(ogvertex_t* v) {
  if (KI_USE_ARENA) return NULL;

  if (v == NULL) return NULL;

  /* FIXME: temp debug, v may be corrupted otherwise double free */
  /* if (v->seq == NULL && v->name != NULL) return NULL; */

  v->name = kifreestr(v->name);
  v->seq  = kifreestr(v->seq);
  int i;
  for (i = 0; i < 4; ++i) {
    v->ies[i] = kiFreeOgEdge(v->ies[i]);
    v->oes[i] = kiFreeOgEdge(v->oes[i]);
  }
  kifree(v, sizeof(ogvertex_t));
  return NULL;
}

ogedge_t* kiFreeOgEdge(ogedge_t* e) {
  if (KI_USE_ARENA) return NULL;

  if (e == NULL) return NULL;
  e->ext = kifreestr(e->ext);
  kifree(e, sizeof(ogedge_t));
  return NULL;
}

int_list_t* kiAllocIntList(int val) {
  int_list_t* list = kiArenaMalloc(KI_ARENA_DEFAULT, sizeof(int_list_t));
  list->value = val;
  list->next = NULL;
  return list;
}

int_list_t* kiFreeIntList(int_list_t* list) {
  if (list == NULL) return NULL;
  if (KI_USE_ARENA) return NULL;
  
  list->next = kiFreeIntList(list->next);
  kifree(list, sizeof(int_list_t));
  return NULL;
}

void kiPrintIntList(int_list_t* list) {
  int_list_t* p = list;
  pmsg(0, "[ ");
  while (p) {
    pmsg(0, "%d ", p->value);
    p = p->next;
  }
  pmsg(0, "]\n");
}

overlap_graph_t* kiAllocOverlapGraph(int minOverlap) {
  assert(minOverlap > 8);
  overlap_graph_t* g   = (overlap_graph_t*)kimalloc(sizeof(overlap_graph_t));
  g->minOverlap        = minOverlap;
  g->seqLen            = 0;
  g->nVertices         = 0;
  g->nSaved            = 0;
  g->vertices          = NULL;
  g->tailSorted        = NULL;
  g->tailOffsetStarts  = NULL;
  g->tailOffsetEnds    = NULL;
  g->maxTailOffset     = 0;
  g->nTailOffsets      = 0;
  g->nSavedTailOffsets = 0;
  return g;
}

overlap_graph_t* kiFreeOverlapGraph(overlap_graph_t* g) {
  if (g == NULL) return NULL;
  if (g->nSaved > 0) {
    /* if (!KI_USE_ARENA) { */
      ogvertex_t** p;
      int i;
      for (i = 0, p = g->vertices; i < g->nVertices; ++i, ++p) { 
        *p = kiFreeOgVertex(*p);
      }
    /* } */
    g->vertices = kifree(g->vertices, sizeof(ogvertex_t*) * g->nSaved);
    g->tailOffsetStarts = kifree(g->tailOffsetStarts, sizeof(int_list_t*) * g->nSavedTailOffsets);
    g->tailOffsetEnds   = kifree(g->tailOffsetEnds,   sizeof(int_list_t*) * g->nSavedTailOffsets);
    g->tailSorted = kiFreeIntList(g->tailSorted); /* frees tailOffsetStarts/Ends as well */
  }
  kifree(g, sizeof(overlap_graph_t));
  return NULL;
}

void kiClearOverlapGraph(overlap_graph_t* g) {
  assert(g);

  if (!KI_USE_ARENA) {
    ogvertex_t** p;
    int i;
    for (i = 0, p = g->vertices; i < g->nVertices; ++i, ++p) { 
      *p = kiFreeOgVertex(*p);
    }
  } else {
    kiArenaClear(KI_ARENA_DEFAULT);
  }
  
  memset(g->vertices, 0, sizeof(ogvertex_t*) * g->nSaved);
  memset(g->tailOffsetStarts, 0, sizeof(int_list_t*) * g->nSavedTailOffsets);
  memset(g->tailOffsetEnds,   0, sizeof(int_list_t*) * g->nSavedTailOffsets);
  
  g->tailSorted = kiFreeIntList(g->tailSorted); /* frees tailOffsetStarts/Ends as well */

  g->seqLen        = 0;
  g->nVertices     = 0;
  g->maxTailOffset = 0;
  g->nTailOffsets  = 0;
}

void kiOverlapGraphAppendVertex(overlap_graph_t* g, ogvertex_t* v) {
  kipmsg(7, "Adding new v: name = '%s'\n", v->name);
  if (v->id > 0) return;        /* v already in g */
  
  int i;
  if (g->nSaved == 0) {
    g->seqLen = strlen(v->seq);
    g->nSaved = 100;
    g->vertices = kimalloc(g->nSaved * sizeof(ogvertex_t*));
    for (i = 0; i < g->nSaved; ++i) {
      g->vertices[i] = NULL;
    }
    g->nSavedTailOffsets = 500;
    g->tailOffsetStarts  = kimalloc(g->nSavedTailOffsets * sizeof(int_list_t*));
    g->tailOffsetEnds    = kimalloc(g->nSavedTailOffsets * sizeof(int_list_t*));
    for (i = 0; i < g->nSavedTailOffsets; ++i) {
      g->tailOffsetStarts[i] = NULL;
      g->tailOffsetEnds[i]   = NULL;
    }
  } else if (g->nVertices == g->nSaved) {
    int nNewSaved = g->nSaved * 2;
    kipmsg(7, "new g->nSaved = %d\n", nNewSaved);
    g->vertices = kirealloc(g->vertices, sizeof(ogvertex_t*)*g->nSaved, sizeof(ogvertex_t*)*nNewSaved, /*copy*/false);
    for (i = g->nSaved; i < nNewSaved; ++i) {
      g->vertices[i] = NULL;
    }
    g->nSaved = nNewSaved;
  }
  /* if (g->nVertices == 8429) { */
    /* kipm("Watch out!\n"); */
  /* } */

  v->id = g->nVertices;
  g->vertices[g->nVertices] = v;
  g->nVertices++;
  /* kipm("v->id = %d\n", v->id); */
  /* if ((v->id == 12313 || v->id == 12314 || v->id == 12315) && ki_domain_rank == 0) { */
  /*   kipm("v->id = %d, v->name = '%s', v->seq = %s\n", v->id, v->name, v->seq); */
  /* } */

  if (v->seq == NULL) return;   /* intermediate nodes */

  if (v->offset >= g->nSavedTailOffsets) {
    int nNewSaved = MAX(g->nSavedTailOffsets * 2, v->offset + 1);
    kipmsg(7, "v->offset = %d\n", v->offset);
    kipmsg(7, "new g->nSavedTailOffsets = %d\n", nNewSaved);
    g->tailOffsetStarts = kirealloc(g->tailOffsetStarts, sizeof(int_list_t*)*g->nSavedTailOffsets, sizeof(int_list_t*)*nNewSaved, /*copy*/false);
    g->tailOffsetEnds = kirealloc(g->tailOffsetEnds,   sizeof(int_list_t*)*g->nSavedTailOffsets, sizeof(int_list_t*)*nNewSaved, /*copy*/false);
    for (i = g->nSavedTailOffsets; i < nNewSaved; ++i) {
      g->tailOffsetStarts[i] = NULL;
      g->tailOffsetEnds[i] = NULL;
    }
    g->nSavedTailOffsets = nNewSaved;
  }

  assert(v->offset < g->nSavedTailOffsets);
  if (v->id == 0) {
    g->tailSorted = kiAllocIntList(v->id);
    g->tailOffsetStarts[v->offset] = g->tailSorted;
    g->tailOffsetEnds[v->offset]   = g->tailSorted;
  } else {
    int_list_t *newNode, *endNode;
    i = v->offset;
    newNode = kiAllocIntList(v->id);
    if (g->tailOffsetEnds[i] != NULL) {
      endNode = g->tailOffsetEnds[i];
      newNode->next = endNode->next;
      endNode->next = newNode;
      g->tailOffsetEnds[i] = newNode;
    } else {
      g->tailOffsetStarts[i] = newNode;
      g->tailOffsetEnds[i] = newNode;
      for (--i; i > 0 && g->tailOffsetEnds[i] == NULL; --i);
      endNode = g->tailOffsetEnds[i];
      assert(endNode != NULL);
      newNode->next = endNode->next;
      endNode->next = newNode;
    }
  }
  if (v->offset > g->maxTailOffset)
    g->maxTailOffset = v->offset;
  g->nTailOffsets++;
}

void kiOverlapGraphExportDot(overlap_graph_t* g, char* filename, int step, char* postfix, bool big) {
  if (g->nVertices == 0) return;
  kipm("In kiOverlapGraphExportDot\n");
  

  /* if (g->nVertices < 15) return; */

  char sstep[20] = "";
  if (postfix) sprintf(sstep, "-%d", step);
  
  char dotname[200];
  sprintf(dotname, "%s%s%s.dot", filename, sstep, postfix ? postfix : "");
  
  /* FILE* fp = stderr; */
  FILE* fp = fopen(dotname, "w");
  assert(fp);

  char vbuf[100000];
  char ebuf[100000];
  bool seen[1000];
  memset(seen, 0, sizeof(bool) * 1000);

  int i, j;
  char* pv = vbuf;
  char* pe = ebuf;
  ogvertex_t** vs = g->vertices;
  ogvertex_t* v;
  ogedge_t* e;

  /* int min_i = 0; */
  int min_i = MAX(g->nVertices - 1000, 0);
  /* if (g->nVertices > 200) */
    kipm("Graph has %d vertices. Showing the last 1000. \n", g->nVertices);


  for (i = min_i; i < g->nVertices; ++i, ++vs) {
    v = *vs;
    if (!seen[v->id]) {
      seen[v->id] = true;

      /* vertices */
      char vAttr[100] = "";
      if (v->name == NULL) {    /* intermediate node */
        /* sprintf(vAttr, "shape=point"); */
        sprintf(vAttr, "label=\"\" margin=\"0,0\" width=0.2 height=0.2");
      } else {
        char label[100];
        sprintf(label, "%s", v->name);
        if (v->occr > 1) {
          sprintf(label, "%s(%d)", v->name, v->occr);
        }
        sprintf(vAttr, "label=\"%s\"", label);
        /* sprintf(vAttr, "label=\"%s/%d\"", label, v->flag); */
        /* if (v->occr > 1) { */
        /*   sprintf(vAttr, "%s peripheries=2", vAttr); */
        /* } */
      }
      if (v->flag & KI_OGV_EXTEND) {
        sprintf(vAttr, "%s color=blue penwidth=2.0", vAttr);
      }
      if (v->flag == KI_OGV_SOLID) {
        sprintf(vAttr, "%s style=filled fillcolor=lightskyblue", vAttr);
      }
      sprintf(pv, "  V%d\t[%s];\n", v->id, vAttr); pv += strlen(pv);

      /* edges */
      for (j = 0; j < 4; ++j) {
        if ((e = v->oes[j]) != NULL) {
          char eAttr[100] = "";
          if (big) {
            sprintf(eAttr, "label=\"%s\"", e->ext);
            sprintf(eAttr, "label=\"%s/%.0f\"", e->ext, e->weight);
            sprintf(eAttr, "%s \tminlen=%d", eAttr, (int)strlen(e->ext));
          } else {
            /* sprintf(eAttr, "label=\"%s/%.0f\"", e->ext, e->weight); */
            sprintf(eAttr, "label=\"%d %.0f.\"", (int)strlen(e->ext), e->weight);
          }
          sprintf(pe, "  V%d -> V%d\t[%s];\n", v->id, e->v->id, eAttr); pe += strlen(pe);
        }
      }
    }
  }
  
  /* rank */
  char rbuf[100000] = "";
  char* pr = rbuf;
  int_list_t* pl;
  if (big) {
    for (j = 0; j < g->maxTailOffset; ++j) {
      if ((pl = g->tailOffsetStarts[j]) != NULL) {
        v = g->vertices[pl->value];
        int offset = v->offset;
        sprintf(pr, "  { rank=same;"); pr += strlen(pr);
        for (; pl != NULL; pl = pl->next) {
          if (g->vertices[pl->value]->offset != offset)
            break;
          sprintf(pr, " V%d;", pl->value); pr += strlen(pr);
        }
        sprintf(pr, " }\n"); pr += strlen(pr);
      }
    }
  }
  
  char fbuf[100] = "fontname=Droid fontsize=16";

  fprintf(fp, "digraph G {\n");
  fprintf(fp, "  rankdir=LR;\n");
  fprintf(fp, "  ranksep=0.2;\n");
  fprintf(fp, "  node [%s margin=\"0,0\"];\n", fbuf);
  fprintf(fp, "  edge [%s arrowsize=0.75];\n\n", fbuf);
  fprintf(fp, "%s\n%s\n%s", vbuf, rbuf, ebuf);
  fprintf(fp, "}\n");
  
  fclose(fp);

  kipm("Out of kiOverlapGraphExportDot\n");

}

void kiOverlapGraphAdd(overlap_graph_t* g, char* name, char* seq, int offset, int baseOffset, int minOverlap, bool bCopy) {
  ogvertex_t* v = kiAllocOgVertex();
  v->offset = offset + baseOffset;
  v->occr   = 1;
  v->flag   = 0;
  if (bCopy) {
    v->name = kiArenaStrdup(KI_ARENA_DEFAULT, name);
    v->seq  = kiArenaStrdup(KI_ARENA_DEFAULT, seq);
  } else {
    v->name = name;
    v->seq  = seq;
  }
  /* if (strcmp(v->name, "427760") == 0 || */
  /*     strcmp(v->name, "16433") == 0 || */
  /*     strcmp(v->name, "10542") == 0 ||       */
  /*     strcmp(v->name, "12378") == 0 ||       */
  /*     strcmp(v->name, "7389") == 0 ||       */
  /*     strcmp(v->name, "1854") == 0 ||       */
  /*     strcmp(v->name, "13655") == 0) { */
  /*   kipm("adding v->name = %s, offset = %d\n", v->name, offset); */
  /* } */

  /* kipm("size = %d\n", sizeof(void*)); */
  if (g->nVertices == 0) {
    kiOverlapGraphAppendVertex(g, v);
    assert(offset == 0 && baseOffset == 0);
  } else {
    int seqLen = strlen(seq);
    if (baseOffset >= 0) {
      /* TODO */
      int maxOffset = baseOffset + offset;
      if (g->nSavedTailOffsets > 0 && g->nSavedTailOffsets-1 < maxOffset) maxOffset = g->nSavedTailOffsets-1;
      /* TODO */
      int minOffset = MAX(baseOffset + offset - seqLen + g->minOverlap, 0);
      int offset = minOffset;
      /* the following has the problem of ignoring the intermediate nodes that with minOffset */
      while (offset <= maxOffset && g->tailOffsetStarts[offset] == NULL) offset++; 
      /* while (offset <= maxOffset && g->tailOffsetStarts[offset] == NULL) offset--; */
      kipmsg(6, "offset = %d\n", offset);
      if (offset > maxOffset || offset >= g->nSavedTailOffsets) {
        kipm("seqLen = %d, minOffset = %d, maxOffset = %d\n", seqLen, minOffset, maxOffset);
        kipm("offset = %d, baseOffset = %d, g->nSavedTailOffsets = %d, g->maxTailOffset = %d\n", offset, baseOffset, g->nSavedTailOffsets, g->maxTailOffset);
      }
      
      assert(offset <= maxOffset);
      assert(offset < g->nSavedTailOffsets);
      assert(g->tailOffsetStarts[offset]);

      int_array_set_t* ancestors = kiAllocIntArraySet();
      int_list_t* it;
      int used = 0;
      for (it = g->tailOffsetStarts[offset]; it != NULL; it = it->next) {
        kipmsg(7, "it->value = %d\n", it->value);

        if (kiIntArraySetHas(ancestors, it->value))
          continue;             /* has traversed along this branch before */

        ogvertex_t* u = g->vertices[it->value];

        if (u == v) break;     
        if (u->offset > v->offset) break; /* past possible ancestors */

        kipmsg(6, "Consider potential ancestor of %s (%d) <-- %s (%d)\n", v->name, v->id, u->name, u->id);

        int d = v->offset - u->offset;
        int exlen = 0;
        int iDiff = 0;            
        
        float diff = kiSeqNCmpX(u->seq + d, v->seq, seqLen - d, &iDiff);
        kipmsg(7, "diff = %f\n", diff);
        float weight = 1. - diff;
        if (iDiff >= minOverlap) { /* iDiff set only if seqs differ */
          /* FIXME */
        } else if (diff < 1e-5) {
          kiIntArraySetAdd(ancestors, u->id);
          bool skip = false;
          char* s = v->seq + (seqLen - d);
          ogedge_t* e = u->oes[kiBaseIndex(s[0])];
          while (e != NULL && d > 0) {
            exlen = strlen(e->ext);
            iDiff = 0;
            if (kiSeqNCmpX(s, e->ext, exlen, &iDiff) > 1e-5) break;
            if (d < exlen) {
              iDiff = d; break;
            }
            e->weight += weight;
            d -= exlen;
            s += exlen;
            u = e->v;
            e = u->oes[kiBaseIndex(s[0])];
            weight = 1.;
            if (kiIntArraySetHas(ancestors, u->id)) {
              skip = true;
              break;
            }
            if (d > 0)
              kiIntArraySetAdd(ancestors, u->id);
          }
          if (skip) continue;

          used = 1;
          if (d == 0) {           /* duplicate */
            kipmsg(7, "duplicate\n");
            u->occr++;
            v = kiFreeOgVertex(v);
            break;                /* multipe precursors should have been handled already */
          } else if (e == NULL) { /* new read */
            kipmsg(7, "new branch\n");
            e = kiAllocOgEdge();
            e->ext    = kiArenaStrdup(KI_ARENA_DEFAULT, s);
            e->v      = v;
            e->weight = weight;
            u->oes[kiBaseIndex(s[0])] = e;
            kiOverlapGraphAppendVertex(g, v);
          } else {                /* split an existing edge */
            kipmsg(7, "split\n");

            d -= iDiff;
            char* t = e->ext;
            ogvertex_t* w = (d > 0) ? kiAllocOgVertex() : v;

            /* assert(w->oes[kiBaseIndex(t[iDiff])] == NULL || d > 0); */

            ogedge_t* e1;
            if (w->oes[kiBaseIndex(t[iDiff])] == NULL) {
              e1 = kiAllocOgEdge();
              w->oes[kiBaseIndex(t[iDiff])] = e1;
              e1->ext    = kiArenaStrdup(KI_ARENA_DEFAULT, t + iDiff);
              e1->v      = e->v;
              e1->weight = e->weight;
            } else {
              e1 = w->oes[kiBaseIndex(t[iDiff])];
              e1->weight += e->weight;
            }
            
            w->flag = e->v->flag;

            e->ext = kiArenaStrndup(KI_ARENA_DEFAULT, t, iDiff);
            e->v = w;
            e->weight += 1.;

            if (d > 0) {         /* add a new branch to the new intermediate vertex */
              kipmsg(7, "new node\n");
              ogedge_t* e2 = kiAllocOgEdge();
              w->oes[kiBaseIndex(s[iDiff])] = e2;
              e2->ext = kiArenaStrdup(KI_ARENA_DEFAULT, s + iDiff);
              e2->v = v;
              e2->weight = 1.;
              kiOverlapGraphAppendVertex(g, v);
            }
            kiOverlapGraphAppendVertex(g, w);
            kiArenaFreestr(KI_ARENA_DEFAULT, t);
          }

          /* kiOverlapGraphExportDot(g, "tmp"); */
          /* kiPrintIntList(g->tailSorted); */
        }
      }
      kiFreeIntArraySet(ancestors);
      if (used == 0) {
        /* kipm("used = 0\n"); */
        kiFreeOgVertex(v);
      }
    }
  }
}

void kiOverlapGraphAddMatches(overlap_graph_t* g, alignment_t* matches, int baseOffset, int minOverlap) {
  int i;

  /* temporary slow sorting */
  int swapFlag, j;
  char* swap;
  for (j = 0; j < matches->nSeq-1; ++j) {
    for (i = j+1; i < matches->nSeq; ++i) {
      if (matches->flags[i] < matches->flags[j]) {
          swap = matches->seqs[j];
          matches->seqs[j] = matches->seqs[i];
          matches->seqs[i] = swap;
          swap = matches->names[j];
          matches->names[j] = matches->names[i];
          matches->names[i] = swap;
          swapFlag = matches->flags[j];
          matches->flags[j] = matches->flags[i];
          matches->flags[i] = swapFlag;
      }
    }
  }

  for (i = 0; i < matches->nSeq; ++i) {
  /* for (i = 0; i < matches->nSeq && g->nVertices < KI_OG_MAX_SIZE; ++i) { */
    kipmsg(4, "Adding matches[%d] name = '%s', seq = '%s', flag = %d\n", i, matches->names[i], matches->seqs[i], matches->flags[i]);
    kiOverlapGraphAdd(g, matches->names[i], matches->seqs[i], matches->flags[i], baseOffset, minOverlap, /*copy*/false);
    matches->seqs[i] = matches->names[i] = NULL; /* if copy = false, strings taken over by graph */
  }
}

void kiOverlapGraphReviseSeed(overlap_graph_t* g, char* query, int* iNewSeed) {
  int nBest = 0;
  int iBest = 0;
  int best  = 0;
  int iSeed = 0;
  ogvertex_t* v;
  int_list_t* listIt = g->tailOffsetStarts[0];
  while (listIt != NULL) {
    v = g->vertices[listIt->value];
    if (strcmp(v->seq, query) == 0) {
      v->flag = KI_OGV_EXTEND;
      iSeed = v->id;
    }
    if (v->occr > best) {
      best  = v->occr;
      iBest = v->id;
      nBest = 1;
    } else if (v->occr == best) {
      nBest++;
    }
    if (listIt == g->tailOffsetEnds[0]) break;
    listIt = listIt->next;
  }
  v = g->vertices[iSeed];
  if (v->occr == best) {
    v->flag = KI_OGV_SOLID;
  } else {
    kipmsg(5, "Revised seed seq to extend.\n");
    v = g->vertices[iBest];
    v->flag = KI_OGV_TODO;
  }
  *iNewSeed = v->id;
}

bool kiOgvUnextended(overlap_graph_t* g, int id) {
  return (g->vertices[id]->flag & KI_OGV_EXTEND) == 0;
}

void kiOverlapGraphExplore(overlap_graph_t* g, int vi, /*OUT*/char* adventure) {
  adventure[0] = '\0';
  int i;
  float bestWeight;
  ogedge_t *bestEdge, *e;
  char* p = adventure;
  ogvertex_t* v = g->vertices[vi];
  while (v != NULL) {
    bestEdge = NULL;
    bestWeight = 0.0;
    for (i = 0; i < 4; ++i) {
      if ((e = v->oes[i]) != NULL) {
        if (e->weight > bestWeight) {
          bestWeight = e->weight;
          bestEdge = e;
        }
      }      
    }
    v->flag |= KI_OGV_EXTEND;
    if (bestEdge == NULL) break;
    strcpy(p, bestEdge->ext); p += strlen(p);
    kipmsg(5, "adventure len = %d\n", p - adventure);
    if (p - adventure > 200) break;
    v = bestEdge->v;
  }
}

bool kiOverlapGraphExtend(overlap_graph_t* g, /*IN/OUT*/int* vi) {
  bool progress = false;
  int i;
  float bestWeight;
  ogedge_t *bestEdge, *e;
  ogvertex_t* v = g->vertices[*vi];
  while (v != NULL && (v->flag & KI_OGV_EXTEND)) {
    bestEdge = NULL;
    bestWeight = 0.0;
    for (i = 0; i < 4; ++i) {
      if ((e = v->oes[i]) != NULL) {
        if (e->weight > bestWeight) {
          bestWeight = e->weight;
          bestEdge = e;
        }
      }      
    }
    v->flag = KI_OGV_SOLID;
    if (v->name != NULL) *vi = v->id;
    if (bestEdge == NULL) break;
    v = bestEdge->v;
    progress = true;
  }
  return progress;
}

void kiOverlapGraphGetContig(overlap_graph_t* g, int iSeed, /*OUT*/char* contig) {
  int i;
  char* p = contig;
  p[0] = '\0';
  ogvertex_t* v = g->vertices[iSeed];
  ogedge_t* e;
  while (v != NULL) {
    for (i = 0; i < 4; ++i) {
      if ((e = v->oes[i]) != NULL) {
        if (e->v->flag == KI_OGV_SOLID) {
          strcpy(p, e->ext); p += strlen(p);
          v = e->v;
          break;
        }
      }
    }
    if (e == NULL) break;
  }
}


void kiOverlapGraphAppendContig(overlap_graph_t* g, int iSeed, /*OUT*/char* contig) {
  int i;
  char* p = contig; p += strlen(p);
  ogvertex_t* v = g->vertices[iSeed];
  ogedge_t* e;
  while (v != NULL) {
    for (i = 0; i < 4; ++i) {
      if ((e = v->oes[i]) != NULL) {
        if (e->v->flag == KI_OGV_SOLID) {
          strcpy(p, e->ext); p += strlen(p);
          v = e->v;
          break;
        }
      }
    }
    if (e == NULL) break;
  }
}


int ki_max_tail_offset = 0;
int ki_og_max_size     = 0;
bool kiOverlapGraphTooBig(overlap_graph_t* g) {
  /* if (g->maxTailOffset > ki_max_tail_offset) { */
  /*   ki_max_tail_offset = g->maxTailOffset; */
  /*   kipmsg(2, "ki_max_tail_offset = %d\n", ki_max_tail_offset); */
  /* } */
  /* if (g->nVertices > ki_og_max_size) { */
  /*   ki_og_max_size = g->nVertices; */
  /*   kipmsg(2, "ki_og_max_size = %d\n", ki_og_max_size); */
  /* } */
  
  return (g->nVertices >= KI_OG_MAX_SIZE || g->maxTailOffset >= KI_CONTIG_SIZE / 2 - 1000);
}

void kiCombineElongation(char* seed, char* forward, char* backward, /*OUT*/char* contig) {
  char* p = contig;
  kiReverseComplementSeq(backward, p); p += strlen(p);
  strcpy(p, seed); p += strlen(p);
  strcpy(p, forward);
}

void stateToTxtFile(char* filePrefix) {
  char filename[KI_NAME_BUF_SIZE];
  sprintf(filename, "%s.%d", filePrefix, ki_domain_rank);
  FILE* fp = fopen(filename, "w");
  fprintf(fp, "%d\n", ki_seqs->nSeq);
  int i, bit;
  for (i = 0; i < ki_seqs->nSeq; ++i) {
    bit = (ki_seqs->flags[i] == KI_SEQ_NEW) ? 0 : 1;
    fprintf(fp, "%d ", bit);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

long long kiStoreState(char* filePrefix) {
  if (kiIsDomainRoot()) pm("\n");
  
  char filename[KI_NAME_BUF_SIZE];
  sprintf(filename, "%s.store", filePrefix);

  KI_File_backup(filename);

  long long nProcessed;
  long long nLocal = ki_nseq_processed;
  KI_Reduce(&nLocal, &nProcessed, 1, MPI_LONG_LONG, MPI_SUM, 0, ki_cmm_domain);
  
  int nSeq = ki_seqs->nSeq;
  int compactSize = kiCeilingDevidedBy(nSeq, 8); /* each seq flag is encoded in 1 bit, one additional int for nSeq */
  size_t alignedSize = kiAlignedSize(compactSize + sizeof(int), 4); 
  kipmsg(3, "nSeq = %d, compactSize = %d, alignedSize = %d\n", nSeq, compactSize, alignedSize);
  
  int *cnts = 0, *displs = 0;
  if (kiIsDomainRoot()) {
    cnts   = (int*)kiArenaMalloc(KI_ARENA_FARMER, sizeof(int) * ki_domain_size);
    displs = (int*)kiArenaMalloc(KI_ARENA_FARMER, sizeof(int) * ki_domain_size);
  }

  KI_Gather(&alignedSize, 1, MPI_INT, cnts, 1, MPI_INT, 0, ki_cmm_domain);

  KI_File fh;
  MPI_Status status;
  KI_File_open(filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  if (kiIsDomainRoot()) {
    int j, displsSum = (int)kiAlignedSize(sizeof(int)*(ki_domain_size+1), 4);    
    for (j = 0; j < ki_domain_size; ++j) {
      displs[j] = displsSum;
      displsSum += cnts[j];
      kipmsg(3, "Storing: cnts[%d] = %d, displs[%d] = %d\n", j, cnts[j], j, displs[j]);
    }
    KI_File_write(fh, &ki_domain_size, 1,     MPI_INT, &status);
    KI_File_write(fh, displs, ki_domain_size, MPI_INT, &status);
  }

  int displ;
  KI_Scatter(displs, 1, MPI_INT, &displ, 1, MPI_INT, 0, ki_cmm_domain);

  if (!KI_USE_ARENA) {
    if (kiIsDomainRoot()) {
      kifree(cnts,   (sizeof(int) * ki_domain_size));
      kifree(displs, (sizeof(int) * ki_domain_size));
    }
  } else {
    kiArenaClear(KI_ARENA_FARMER);
  }

  kipmsg(3, "displ = %d\n", displ);
  KI_File_seek(fh, displ, MPI_SEEK_SET);
  KI_File_write(fh, &nSeq, 1, MPI_INT, &status);

  int i;
  int bit;
  int nLeft = compactSize;      /* in bytes */
  char* p = ki_tmp_buf;
  memset(ki_tmp_buf, 0, KI_BUF_SIZE);
  for (i = 0; i < nSeq; ++i) {
    bit = (ki_seqs->flags[i] == KI_SEQ_NEW) ? 0 : 1;
    *p |= (bit << (7 - i % 8));
    if ((i+1) % 8 == 0) ++p;
    if ((i+1) % (KI_BUF_SIZE*8) == 0) {
      KI_File_write(fh, ki_tmp_buf, KI_BUF_SIZE, MPI_BYTE, &status);
      memset(ki_tmp_buf, 0, KI_BUF_SIZE);
      nLeft -= KI_BUF_SIZE;
      p = ki_tmp_buf;
      kipm("nLeft = %d KI_BUF_SIZE = %d\n", nLeft, KI_BUF_SIZE);
    }
  }
  KI_File_write(fh, ki_tmp_buf, nLeft, MPI_BYTE, &status);

  KI_File_close(&fh);

  /* stateToTxtFile("txt.store"); */
  return nProcessed;
  
}

long long kiRestoreState(char* filePrefix) {
  char filename[KI_NAME_BUF_SIZE];  
  sprintf(filename, "%s.store", filePrefix);

  if (!KI_File_exists(filename)) return 0;
  
  KI_File fh;
  MPI_Status status;
  int rc;
  
  rc = KI_File_open(filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (rc != MPI_SUCCESS) {
    if (kiIsDomainRoot()) fprintf(stderr, "Error: unable to open file '%s'\n", filename);
    kiAbort(1);
  }

  int *displs = 0;
  if (kiIsDomainRoot()) {
    int nParts;
    displs = (int*)kiArenaMalloc(KI_ARENA_FARMER, sizeof(int) * ki_domain_size);
    KI_File_read(fh, &nParts, 1, MPI_INT, &status);
    assert(status.count = 1);
    
    if (nParts != ki_domain_size) {
      fprintf(stderr, "Number of parts in .store file does not match domain size (old = %d vs now = %d).\nRestore failed.\n", nParts, ki_domain_size);
      KI_File_close(&fh);
      kiAbort(-1);
    }
    KI_File_read(fh, displs, ki_domain_size, MPI_INT, &status);
    assert(status.count = ki_domain_size);
  }

  int displ;
  KI_Scatter(displs, 1, MPI_INT, &displ, 1, MPI_INT, 0, ki_cmm_domain);
  kipmsg(3, "Restoring: displ = %d\n", displ);
  
  if (!KI_USE_ARENA) {
    if (kiIsDomainRoot()) {
      kifree(displs, (sizeof(int) * ki_domain_size));
    }
  } else {
    kiArenaClear(KI_ARENA_FARMER);
  }

  int nSeq;
  KI_File_seek(fh, displ, MPI_SEEK_SET);
  KI_File_read(fh, &nSeq, 1, MPI_INT, &status);
  assert(status.count = 1);
  
  if (nSeq != ki_seqs->nSeq) {
    fprintf(stderr, "Mismatch in number of sequences (old = %d vs now = %d).\nRestore failed.\n", nSeq, ki_seqs->nSeq);
    KI_File_close(&fh);
    kiAbort(-1);
  }
  

  int nToRead;
  int i = 0, nLeft = nSeq;      /* in bits */
  char* p;
  int bit;
  while (nLeft > 0) {
    nToRead = kiCeilingDevidedBy(nLeft, 8); /* in bytes */
    if (nToRead > KI_BUF_SIZE) nToRead = KI_BUF_SIZE;
    KI_File_read(fh, ki_tmp_buf, nToRead, MPI_BYTE, &status);
    assert(status.count = nToRead);
    nLeft -= nToRead * 8;
    for (p = ki_tmp_buf; i < nSeq; ++i) {
      bit = (*p >> (7 - i % 8)) & 1;
      ki_seqs->flags[i] = bit ? KI_SEQ_USED : KI_SEQ_NEW;
      ki_nseq_processed += bit;
      if ((i+1) % 8 == 0) ++p;
    }
  }

  long long nProcessed;
  long long nLocal = ki_nseq_processed;
  KI_Reduce(&nLocal, &nProcessed, 1, MPI_LONG_LONG, MPI_SUM, 0, ki_cmm_domain);
  
  /* stateToTxtFile("txt.restore"); */
  return nProcessed;
  
}

long long kiStoreStateToSeparateFiles(char* filePrefix) {
  char filename[KI_NAME_BUF_SIZE];  
  
  sprintf(filename, "store.%d", ki_domain_rank);
  FILE* fp = fopen(filename, "wb");
  /* KI_File* fp = NULL; */
  /* KI_File_open(filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, fp); */
  fwrite(&(ki_seqs->nSeq), sizeof(int), 1, fp);
  memset(ki_tmp_buf, 0, KI_BUF_SIZE);
  
  int i;
  int bit;
  char* p = ki_tmp_buf;
  for (i = 0; i < ki_seqs->nSeq; ++i) {
    bit = (ki_seqs->flags[i] == KI_SEQ_NEW) ? 0 : 1;
    *p |= (bit << (7 - i % 8));
    if ((i+1) % 8 == 0) ++p;
    if ((i+1) % (KI_BUF_SIZE*8) == 0) {
      fwrite(ki_tmp_buf, 1, KI_BUF_SIZE, fp);
      p = ki_tmp_buf;
    }
    /* fprintf(fp, "%d ", ); */
  }
  kipm("p-start = %d\n",(p-ki_tmp_buf));
  
  fwrite(ki_tmp_buf, 1, (p-ki_tmp_buf+1), fp);
  fclose(fp);
  
  /* KI_File_close(fp); */

  stateToTxtFile("txt.store");

  return 0;
}

long long kiRestoreStateFromSeparateFiles(char* filePrefix) {
  char filename[KI_NAME_BUF_SIZE];  
  FILE* fp;
  
  sprintf(filename, "store.%d", ki_domain_rank);
  fp = fopen(filename, "rb");
  if (fp == NULL) return 0;

  int nSeq;
  int nRead, nToRead;
  nRead = fread(&nSeq, sizeof(int), 1, fp);
  assert(nRead == 1);
  
  int i = 0, nLeft = nSeq;
  char* p;
  int bit;
  while (nLeft > 0) {
    nToRead = (nLeft-1) / 8 + 1;
    if (nToRead > KI_BUF_SIZE) nToRead = KI_BUF_SIZE;
    nRead = fread(ki_tmp_buf, 1, nToRead, fp);
    kipm("nRead = %d, nToRead = %d\n", nRead, nToRead);
    assert(nRead == nToRead);
    nLeft -= nRead * 8;
    
    for (p = ki_tmp_buf; i < nSeq; ++i) {
      bit = (*p >> (7 - i % 8)) & 1;
      ki_seqs->flags[i] = bit ? KI_SEQ_USED : KI_SEQ_NEW;
      if ((i+1) % 8 == 0) ++p;
    }
  }
  
  stateToTxtFile("txt.restore");
  return 0;
}
