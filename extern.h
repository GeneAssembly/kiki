#ifndef _EXTERN_H_
#define _EXTERN_H_

/* Common headers */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "mpi.h"
#include "debug.h"
#include "mem.h"
#include "io.h"
#include "comm.h"
#include "seq.h"
#include "hash.h"
#include "raiphy.h"
#include "util.h"

#include "ki.h"

#define KI_DOMAIN_USER   0
#define KI_DOMAIN_DEALER 1
#define KI_DOMAIN_FARMER 2

#define KI_STRING        0x00010000
#define KI_VECTOR        0x00100000
#define KI_CHAR          MPI_CHAR
#define KI_BYTE          MPI_BYTE
#define KI_BOOL          MPI_BYTE
#define KI_INT           MPI_INT
#define KI_LONG          MPI_LONG
#define KI_FLOAT         MPI_FLOAT
#define KI_DOUBLE        MPI_DOUBLE
#define KI_LONG_LONG     MPI_LONG_LONG
#define KI_V_CHAR        (MPI_CHAR      | KI_VECTOR)
#define KI_V_BYTE        (MPI_BYTE      | KI_VECTOR)
#define KI_V_INT         (MPI_INT       | KI_VECTOR)
#define KI_V_LONG        (MPI_LONG      | KI_VECTOR)
#define KI_V_FLOAT       (MPI_FLOAT     | KI_VECTOR)
#define KI_V_DOUBLE      (MPI_DOUBLE    | KI_VECTOR)
#define KI_V_LONG_LONG   (MPI_LONG_LONG | KI_VECTOR)

#define KI_NUM_FUNCS     100
/* #define KI_BUF_SIZE      0x08000000 /\* 128MB *\/ */
#define KI_BUF_SIZE      0x04000000 /* 64MB */
/* #define KI_BUF_SIZE      0x02000000 /\* 32MB *\/ */
/* #define KI_BUF_SIZE      0x020000 /\* 128KB *\/ */
#define KI_NAME_BUF_SIZE 1000
#define KI_SEQ_BUF_SIZE  5000
#define KI_CONTIG_SIZE   10000000   /* 10M bps */

#define KI_ARENA_DEFAULT 0
#define KI_ARENA_TEMP    1
#define KI_ARENA_USER    2
#define KI_ARENA_FARMER  3
#define KI_ARENA_SEQS    4
#define KI_ARENA_SEQS_EL 5      /* equal length */
#define KI_ARENA_HASH    6
#define KI_ARENA_MAX     10

#define KI_SEQ_NEW       0
#define KI_SEQ_BOOKED    1
#define KI_SEQ_USED      2

#define KI_DEF_KMER_LEN  15
#define KI_MAX_KMER_LEN  63

#define KI_HASH_NONE     0
#define KI_HASH_ENDS     1
#define KI_HASH_COMPLETE 2

#define KI_DEF_HASH_MODE KI_HASH_ENDS

#define KI_CLOCK_NONE    0
#define KI_CLOCK_LOAD    1
#define KI_CLOCK_MAIN    2

#define KI_FASTA_ID_ONLY 0
#define KI_FASTA_ID_DESC 1

#define KI_TERMINATION_NONE      0
#define KI_TERMINATION_TIME      1
#define KI_TERMINATION_DEPLETION 2

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


/* Global variables */
extern MPI_Group ki_grp_world, ki_grp_domain, ki_grp_user, ki_grp_dealer, ki_grp_farmer;
extern MPI_Comm  ki_cmm_world, ki_cmm_domain, ki_cmm_user, ki_cmm_dealer, ki_cmm_farmer;

extern MPI_Group ki_grp_user_dealer, ki_grp_dealer_farmer;
extern MPI_Comm  ki_cmm_user_dealer, ki_cmm_dealer_farmer;

extern int  ki_world_rank;
extern int  ki_world_size;
extern int  ki_user_rank;
extern int  ki_user_size;
extern int  ki_dealer_rank;
extern int  ki_dealer_size;
extern int  ki_farmer_rank;
extern int  ki_farmer_size;
extern int  ki_top_dealer;
extern int  ki_top_farmer;
extern int  ki_domain_rank;
extern int  ki_domain_size;
extern int  ki_domain_id;
extern char ki_domain_name[20];

extern int ki_kmer_len;
extern int ki_hash_mode;

/* Implementation-specific */
extern alignment_t*  ki_seqs;
extern hashtable_t** ki_hash;

extern int ki_seqs_index;
extern int ki_nseq_processed;

/* Tables */
extern int base2complement[117];
extern int base2int[117];
extern int rai_base2int[117];

#endif /* _EXTERN_H_ */
