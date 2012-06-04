#if defined(SIMPLE_HASH) 
#else

#include "extern.h"
#include "debug.h"
#include "mem.h"
#include "hash_chain.h"

#define KI_HASH_PREFIX_LEN  12
#define KI_HASH_NUM_CHAINS  (1 << (2*KI_HASH_PREFIX_LEN))

/* Local functions */
inline bool hashFound(hashiterator_t hi) {
  return (hi.bucketNo >= 0);
}

/* double local_time_count = 0.; */
inline int dna2int(char* seq, int prefixLen) {
  /* double t1 = MPI_Wtime(); */
  int i, val = 0;
  for (i = 0; i < prefixLen; i++) {
    val |= ((base2int[(int)seq[i]]) << (i << 1));
  }
  /* double t2 = MPI_Wtime(); */
  /* local_time_count += (t2-t1); */
  return val;
}

inline long long dna2ll(char* seq, int prefixLen) {
  char* start = seq+prefixLen;
  int i, len = MIN(strlen(seq), ki_kmer_len) - prefixLen;
  assert(len > 0);
  long long val = 0L;
  for (i = 0; i < len; i++) {
    val |= ((base2int[(int)start[i]]) << (i << 1));
  }
  return val;
}

inline int aa2hashint(char c) {
  if (c == '*') return 31;
  if (c < 'A' || c > 'Z') return 'X';
  return (int)(c - 'A');
}


inline int pep2int(char* seq, int prefixLen) {
  int i, val = 0;
  for (i = 0; i < prefixLen; i++) {
    val |= (aa2hashint(seq[i]) << (i * 5));
  }
  return val;
}

inline long long pep2ll(char* seq, int prefixLen) {
  char* start = seq+prefixLen;
  int i, len = MIN(strlen(seq), ki_kmer_len) - prefixLen;
  assert(len > 0);
  long long val = 0L;
  for (i = 0; i < len; i++) {
    val |= (aa2hashint(start[i]) << (i * 5));
  }
  return val;
}

inline int seq2int(char* seq, hashtable_t* hash) {
  return hash->bProtein ? pep2int(seq, hash->prefixLen) : dna2int(seq, hash->prefixLen);
}

inline long long seq2ll(char* seq, hashtable_t* hash) {
  return hash->bProtein ? pep2ll(seq, hash->prefixLen) : dna2ll(seq, hash->prefixLen);
}

/* Hashtable functions */

static inline int calc_hash(char *str)
{
  int i, ret;

  ret = 0;
  for (i = 0; i < strlen(str); i++) {
    ret += (unsigned int) (str[i]);
    ret += ((1 << (i % 4)) * (unsigned int) (str[i]));
  }

  return ret;
}

static inline void create_string_hash(hashtable_t* hash, int hi)
{
  /* hash->bucket_string_hashes[hi] = calc_hash(hash->buckets[hi].string); */
}

hashchain_t* kiAllocHashChain() {
  hashchain_t* chain = kiArenaMalloc(KI_ARENA_HASH, sizeof(hashchain_t));
  memset(chain, 0, sizeof(hashchain_t));
  return chain;
}

int kiHashChainGrow(hashchain_t* chain) {
  chain->nBuckets++;
  
  if (chain->nSaved == 0) {
    chain->nSaved = 1;
    chain->buckets = (hashbucket_t*)kiArenaMalloc(KI_ARENA_HASH, sizeof(hashbucket_t) * chain->nSaved);
    memset(chain->buckets, 0, sizeof(hashbucket_t) * chain->nSaved);
  }
  if (chain->nBuckets > chain->nSaved) {
    int nNewSaved = chain->nSaved * 2;
    chain->buckets = kiArenaRealloc(KI_ARENA_HASH, chain->buckets, sizeof(hashbucket_t)*(chain->nSaved), sizeof(hashbucket_t)*nNewSaved, /*copy*/false);
    memset(chain->buckets + chain->nSaved, 0, sizeof(hashbucket_t) * (nNewSaved - chain->nSaved));
    chain->nSaved = nNewSaved;
  }
  return chain->nBuckets - 1;
}

inline hashbucket_t* kiGetHashBucket(hashtable_t* hash, hashiterator_t hi) {
  if (!hash->chains[hi.chainNo] || !hashFound(hi) || hash->chains[hi.chainNo]->nBuckets <= hi.bucketNo) 
    return NULL;
  else
    return (hash->chains[hi.chainNo]->buckets + hi.bucketNo);
}

/* Make hashtable from list of strings and make no copies of them */
hashtable_t* kiMakeHashtable(char** strings, int nStrings, int nBuckets, bool bCopy, bool bProtein) {
  assert(nStrings <= nBuckets);

  long used1 = kiGetMallocUsed();

  hashtable_t* hash = (hashtable_t*)kimalloc(sizeof(hashtable_t));
  hash->bCopy = bCopy;
  hash->bProtein = bProtein;
  kipmsg0(3, "nBuckets = %d\n", nBuckets);
  
  if (nBuckets > 0) {
    hash->prefixLen = 1;
    int temp = nBuckets - 1;
    while ((temp >>= 2) > 0) {
      hash->prefixLen++;
    } 
  } else {
    hash->prefixLen = KI_HASH_PREFIX_LEN;
  }
  kipmsg0(4, "prefixLen = %d\n", hash->prefixLen);

  hash->prefixLen = MIN(hash->prefixLen, KI_HASH_PREFIX_LEN);
  hash->nChains   = 1 << (2 * hash->prefixLen);
  if (bProtein) {
    hash->prefixLen = MAX(hash->prefixLen/3, 1);
    hash->nChains   = 1 << (5 * hash->prefixLen);
  } 
  kipmsg0(3, "nChains = %d, prefixLen = %d\n", hash->nChains, hash->prefixLen);
  
  hash->chains    = (hashchain_t**)kimalloc(sizeof(hashchain_t*) * hash->nChains);
  memset(hash->chains, 0, sizeof(hashchain_t*) * hash->nChains);

  kipmsg0(8, "prefixLen = %d\n", hash->prefixLen);
  
  int i;
  for (i = 0; i < nStrings; i++) {
    kiHashtableAddCopy(hash, strings[i], strlen(strings[i]), i, /*offset*/0);
  }

  long used2 = kiGetMallocUsed();
  kipmsg0(3, "Hash table size = %ld MB\n", (used2 - used1) >> 20);

  return hash;
}

/* Add a string to hashtable and make no copy */
void kiHashtableAdd(hashtable_t* hash, char* str, int seqid) {
}

/* Add a replicated string to hashtable and record its origin */
void kiHashtableAddCopy(hashtable_t* hash, char* str, size_t len, int seqid, int offset) {
  /* reverse seq if seqid < 0 */
  char ch = str[len];
  str[len] = '\0';
  hash->nCount++;
  
  hashiterator_t hi = kiFindMatch(hash, str);
  hashbucket_t* bucket;
  if (!hashFound(hi)) {
    /* save a unique entry */
    if (hash->chains[hi.chainNo] == NULL) 
      hash->chains[hi.chainNo] = kiAllocHashChain();
    hi.bucketNo = kiHashChainGrow(hash->chains[hi.chainNo]);
    bucket = kiGetHashBucket(hash, hi);
    bucket->stringHash = hi.ll;
    bucket->nCount = 1;
    bucket->first = kiAbsSeqId(seqid);
  } else {
    /* record a duplicate entry */
    bucket = kiGetHashBucket(hash, hi);
    bucket->nCount++;
  }
  /* insert to the top of the list */
  origin_list_t* old = bucket->origins;
  origin_list_t* new = kiAllocOriginList(seqid, offset);
  new->next = old;
  if (old) old->prev = new;
  bucket->origins = new;

  /* restore input string */
  str[len] = ch;
}

int bucketCompare(const void *p1, const void *p2) {
  long long l1 = ((hashbucket_t *)p1)->stringHash;
  long long l2 = ((hashbucket_t *)p2)->stringHash;
  
  if (l1 > l2) return +1;
  if (l1 < l2) return -1;
  return 0;
}

void kiHashtableOptimize(hashtable_t* hash) {
  /* return; */
  hashchain_t* chain;  
  int i;
  for (i = 0; i < hash->nChains; ++i) {
    chain = hash->chains[i];
    if (chain != NULL && chain->nBuckets > 10) {
      /* kipm("i = %d, nBuckets = %d\n", i, hash->chains[i]->nBuckets); */
      /* kipm("before: "); */
      /* int j; */
      /* for (j = 0; j < chain->nBuckets; ++j) { */
      /*   pm("%lld ", chain->buckets[j].stringHash); */
      /* } */
      /* pm("\n"); */
      qsort(chain->buckets, chain->nBuckets, sizeof(hashbucket_t), bucketCompare);
      chain->sorted = true;
      /* kipm("after: "); */
      /* for (j = 0; j < chain->nBuckets; ++j) { */
      /*   pm("%lld ", chain->buckets[j].stringHash); */
      /* } */
      /* pm("\n\n"); */
    }
  }
}

void kiClearHashtable(hashtable_t* hash) {
  if (hash == NULL) return;  
  hashbucket_t* p = NULL;
  hashchain_t** q = hash->chains;
  int i, j;
  for (i = 0; i < hash->nChains; ++i, ++q) {
    if (*q != NULL) {
      p = (*q)->buckets;
      for (j = 0; j < (*q)->nBuckets; ++j, ++p) {
        if (hash->bCopy && !KI_USE_ARENA)
          kifreestr(p->string);
        p->string  = NULL;
        p->nCount  = 0;
        p->first   = -1;
        p->stringHash = -1;
        p->origins = kiFreeOriginList(p->origins);
      }
    }
  }
  hash->nCount = 0;
}

hashtable_t* kiFreeHashtable(hashtable_t* hash) {
  if (hash == NULL) return NULL;  
  hashbucket_t* p = NULL;
  hashchain_t** q = hash->chains;
  int i, j;
  for (i = 0; i < hash->nChains; ++i, ++q) {
    if (*q != NULL) {
      p = (*q)->buckets;
      for (j = 0; j < (*q)->nBuckets; ++j, ++p) {
        if (p->string != NULL) {
          if (hash->bCopy && !KI_USE_ARENA)
            kifreestr(p->string);
        }
        p->origins = kiFreeOriginList(p->origins);
      }
      kiArenaFree(KI_ARENA_HASH, (*q)->buckets, sizeof(hashbucket_t) * (*q)->nSaved);
      kiArenaFree(KI_ARENA_HASH, *q, sizeof(hashchain_t));
    }
  }
  kifree(hash->chains, sizeof(hashchain_t*) * hash->nChains);
  kifree(hash, sizeof(hashtable_t));
  return NULL;
}


long long hash_collision_count = 0;
long long hash_query_count = 0;

hashiterator_t kiFindMatch(hashtable_t* hash, char *string) {
  hashiterator_t hi;
  hashchain_t* chain;
  hi.chainNo = seq2int(string, hash);
  hi.ll = seq2ll(string, hash);
  hi.bucketNo = -1;
  hash_query_count++;
  chain = hash->chains[hi.chainNo];
  long long target = hi.ll;
  if (chain != NULL) {
    int i;
    hashbucket_t* p = chain->buckets;
    if (chain->sorted) {        /* binary search */
      int beg = 0, end = chain->nBuckets - 1, mid = end/2;
      while (beg < end) {
        if (target < p[mid].stringHash) {
          end = mid - 1;
        } else if (target > p[mid].stringHash) {
          beg = mid + 1;
        } else break;
        mid = (beg + end) / 2;
        hash_collision_count++;
      }
      if (p[mid].stringHash == target)
        hi.bucketNo = mid;
    } else {
      for (i = 0; i < chain->nBuckets; ++i, ++p) {
        if (p->stringHash == target) {
          hi.bucketNo = i;
          break;
        }
        hash_collision_count++;
      }
    }
  }
  kipmsg(7, "string = '%s', chainNo = %d, ll = %d, bucketNo = %d\n", string, hi.chainNo, hi.ll, hi.bucketNo);
  return hi;
}

float kiGetHashCollisionRatio() {
  /* kipm("time spent in seq2int: %.3f s\n", local_time_count); */
  
  if (hash_query_count > 0)
    kipmsg0(3, "Hash collision ratio = %.3f\n", 1.0*hash_collision_count / hash_query_count);

  return 0.;
}

long kiGetHashTableMemory() {
  return KI_HASH_NUM_CHAINS * 8;
}

float kiGetHashLoadFactor(hashtable_t* hash) {
  /* float ratio = 1. * hash->nCount / hash->nChains; */
  /* return ratio; */
  return 1. * ki_seqs->nSeq / KI_HASH_NUM_CHAINS;
}

char *kiGetHashString(hashtable_t* hash, hashiterator_t hi) {
  hashbucket_t* bucket = kiGetHashBucket(hash, hi);
  return bucket->string;
}

int kiHashCount(hashtable_t* hash, hashiterator_t hi) {
  hashbucket_t* bucket = kiGetHashBucket(hash, hi);
  return bucket->nCount;
}

int kiHashFirst(hashtable_t* hash, hashiterator_t hi) {
  hashbucket_t* bucket = kiGetHashBucket(hash, hi);
  return bucket->first;
}

origin_iterator_t kiHashOrigins(hashtable_t* hash, hashiterator_t hi) {
  hashbucket_t* bucket = kiGetHashBucket(hash, hi);
  return bucket ? bucket->origins : NULL;
}

void kiHashGetMatchingSeq(origin_iterator_t it, /*OUT*/char* match) {  /* assumes global var ki_seqs */
  assert(it != NULL);
  kmerorigin_t *p = &(it->origin);
  if (p->seqid >= 0) {
    strcpy(match, ki_seqs->seqs[p->seqid] + p->offset);
  } else {                      /* reverse sequence */
    kiPartialReverseComplementSeq(ki_seqs->seqs[kiAbsSeqId(p->seqid)], p->offset, match);
  }
}

origin_iterator_t kiHashFetchUnusedMatchingSeq(hashtable_t* hash, hashiterator_t hi, origin_iterator_t it, bool bErase, /*OUT*/char* match, /*OUT*/char** name) {  /* removes used seqs from ki_seqs */

  assert(it != NULL);
  *match = '\0';
  *name = NULL;
  
  /* find the first unused seq and remove the used ones */
  kmerorigin_t *p = &(it->origin);
  int id = kiAbsSeqId(p->seqid);
  while (ki_seqs->flags[id] == KI_SEQ_USED) {
    it = kiHashEraseMatchingSeq(hash, hi, it, bErase);
    if (it == NULL) return NULL;
    p = &(it->origin);
    id = kiAbsSeqId(p->seqid);
  }

  /* get matching seq and name */
  if (p->seqid >= 0) {
    strcpy(match, ki_seqs->seqs[p->seqid] + p->offset);
  } else {                      /* reverse sequence */
    kiPartialReverseComplementSeq(ki_seqs->seqs[id], p->offset, match);
  }
  *name = ki_seqs->names[id];

  return it;
}

origin_iterator_t kiHashEraseMatchingSeq(hashtable_t* hash, hashiterator_t hi, origin_iterator_t it, bool bErase) {  /* assumes global var ki_seqs */
  
  assert(it != NULL);
  int id = kiAbsSeqId(it->origin.seqid);
  if (ki_seqs->flags[id] != KI_SEQ_USED) {
    if (ki_seqs->flags[id] != KI_SEQ_BOOKED) {
      ki_nseq_processed++;
      /* kipm("erased, ki_nseq_processed = %d\n", ki_nseq_processed); */
    }
    
    ki_seqs->flags[id] = KI_SEQ_USED; 
  } 

  if (bErase) {
    hashbucket_t* bucket = kiGetHashBucket(hash, hi);
    bucket->flag = -1;
  }
  
  return it->next;
}

#endif
