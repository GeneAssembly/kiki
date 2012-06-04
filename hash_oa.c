#if defined(SIMPLE_HASH) 

#include "extern.h"
#include "debug.h"
#include "mem.h"
#include "hash_oa.h"

/* #define KI_NUM_BUCKETS   8000000 */
/* #define KI_NUM_BUCKETS   16000000 */
#define KI_NUM_BUCKETS   24000000
/* #define KI_NUM_BUCKETS   30000000 */
/* #define KI_NUM_BUCKETS   40000000 */
/* #define KI_NUM_BUCKETS   80000000 */


/* Hashtable functions */
/* Open addressing with linear probing (n = 1) */

/* pavan */
/* static inline int calc_hash(char *str) */
/* { */
/*     int i, ret; */

/*     ret = 0; */
/*     for (i = 0; i < strlen(str); i++) { */
/*         ret += (unsigned int) (str[i]); */
/*         ret += ((1 << (i % 4)) * (unsigned int) (str[i])); */
/*     } */

/*     return ret; */
/* } */

inline calc_hash(char* str) {
  int i, val = 0;
  
  int len = MIN(strlen(str), 10);
  for (i = 0; i < len; i++) {
    val |= ((base2int[(int)str[i]]) << (i << 1));
  }
  return val;
}

static inline void create_string_hash(hashtable_t* hash, int hi)
{
    hash->bucket_string_hashes[hi] = calc_hash(hash->buckets[hi].string);
}

/* Make hashtable from list of strings and make no copies of them */
hashtable_t* kiMakeHashtable(char **strings, int nStrings, int nBuckets, bool bCopy, bool bProtein) {
  hashtable_t* hash = (hashtable_t*)kimalloc(sizeof(hashtable_t));
  hash->bCopy = bCopy;
  hash->nBuckets = nBuckets > 0 ? nBuckets : KI_NUM_BUCKETS;  /* 8*nStrings; */
  hash->buckets = (hashbucket_t*)kimalloc(sizeof(hashbucket_t) * hash->nBuckets);
  hash->bucket_string_hashes = (int *)kimalloc(sizeof(int) * hash->nBuckets);

  int i;

  for (i = 0; i < hash->nBuckets; i++)
      hash->bucket_string_hashes[i] = 0;

  hashbucket_t* p = hash->buckets;
  for (i = 0; i < hash->nBuckets; ++i, ++p) {
    p->string  = NULL;
    p->nCount  = 0;
    p->first   = -1;
    p->origins = NULL;
  }
  for (i = 0; i < nStrings; i++) {
    hashiterator_t hi = kiFindMatch(hash, strings[i]);
    if (hash->buckets[hi].string == NULL) {
      /* save a unique entry */
      assert(hash->buckets[hi].nCount == 0);
      hash->buckets[hi].string = strings[i];
      create_string_hash(hash, hi);
      hash->buckets[hi].nCount = 1;
      hash->buckets[hi].first = i;
    } else {
      /* record a duplicate entry */
      assert(hash->buckets[hi].string != NULL);
      assert(strcmp(hash->buckets[hi].string, strings[i]) == 0);
      assert(hash->buckets[hi].first >= 0);
      hash->buckets[hi].nCount++;
    }
  }
  return(hash);
}

/* Add a string to hashtable and make no copy */
void kiHashtableAdd(hashtable_t* hash, char* str, int seqid) {
  hashiterator_t hi = kiFindMatch(hash, str);
  if (hash->buckets[hi].string == NULL) {
    /* save a unique entry */
    assert(hash->buckets[hi].nCount == 0);
    hash->buckets[hi].string  = hash->bCopy ? (char*)kimemdup(str, strlen(str)+1) : str;
    create_string_hash(hash, hi);
    hash->buckets[hi].nCount  = 1;
    hash->buckets[hi].first   = seqid;
    hash->buckets[hi].origins = NULL;
  } else {
    /* record a duplicate entry */
    assert(hash->buckets[hi].string != NULL);
    assert(strcmp(hash->buckets[hi].string, str) == 0);
    assert(hash->buckets[hi].first >= 0);
    hash->buckets[hi].nCount++;
  }
}

/* Add a replicated string to hashtable and record its origin */
void kiHashtableAddCopy(hashtable_t* hash, char* str, size_t len, int seqid, int offset) {
  /* reverse seq if seqid < 0 */
  char ch = str[len];
  str[len] = '\0';
  hashiterator_t hi = kiFindMatch(hash, str);
  
  if (hash->buckets[hi].string == NULL) {
    /* save a unique entry */
    assert(hash->buckets[hi].nCount == 0);
    hash->buckets[hi].string = (char*)kiArenaMemdup(KI_ARENA_HASH, str, len+1);
    create_string_hash(hash, hi);
    kipmsg(5, "str = %s\n", hash->buckets[hi].string);
    hash->buckets[hi].nCount = 1;
    hash->buckets[hi].first = kiAbsSeqId(seqid);
  } else {
    /* record a duplicate entry */
    assert(hash->buckets[hi].string != NULL);
    assert(strcmp(hash->buckets[hi].string, str) == 0);
    assert(hash->buckets[hi].first >= 0);
    hash->buckets[hi].nCount++;
  }
  /* insert to the top of the list */
  origin_list_t* old = hash->buckets[hi].origins;
  origin_list_t* new = kiAllocOriginList(seqid, offset);
  new->next = old;
  if (old) old->prev = new;
  hash->buckets[hi].origins = new;
  
  /* restore input string */
  str[len] = ch;
}

void kiClearHashtable(hashtable_t* hash) {
  if (hash == NULL) return;  
  hashbucket_t* p = hash->buckets;
  int i;
  for (i = 0; i < hash->nBuckets; ++i, ++p) {
    if (p->string != NULL) {
      if (hash->bCopy && !KI_USE_ARENA)
        kifreestr(p->string);
      p->string  = NULL;
      p->nCount  = 0;
      p->first   = -1;
      p->origins = kiFreeOriginList(p->origins);
    }
    hash->bucket_string_hashes[i] = 0;
  }
}

hashtable_t* kiFreeHashtable(hashtable_t* hash) {
  if (hash != NULL) {
    hashbucket_t* p = hash->buckets;
    int i;
    for (i = 0; i < hash->nBuckets; ++i, ++p) {
      if (p->string != NULL) {
        if (hash->bCopy && !KI_USE_ARENA)
          kifreestr(p->string);
        p->origins = kiFreeOriginList(p->origins);
      }
    }
    kifree(hash->buckets, sizeof(hashbucket_t) * hash->nBuckets);
    kifree(hash->bucket_string_hashes, sizeof(int) * hash->nBuckets);
    kifree(hash, sizeof(hashtable_t));
  }
  return NULL;
}

void kiHashtableOptimize(hashtable_t* hash) {
}

#define CRC_HASH_MASK    0x0000000000ffffffL

static int crc_table[256] = {
  0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419,
  0x706af48f, 0xe963a535, 0x9e6495a3, 0x0edb8832, 0x79dcb8a4,
  0xe0d5e91e, 0x97d2d988, 0x09b64c2b, 0x7eb17cbd, 0xe7b82d07,
  0x90bf1d91, 0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de,
  0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7, 0x136c9856,
  0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9,
  0xfa0f3d63, 0x8d080df5, 0x3b6e20c8, 0x4c69105e, 0xd56041e4,
  0xa2677172, 0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
  0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940, 0x32d86ce3,
  0x45df5c75, 0xdcd60dcf, 0xabd13d59, 0x26d930ac, 0x51de003a,
  0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423, 0xcfba9599,
  0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
  0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d, 0x76dc4190,
  0x01db7106, 0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f,
  0x9fbfe4a5, 0xe8b8d433, 0x7807c9a2, 0x0f00f934, 0x9609a88e,
  0xe10e9818, 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
  0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e, 0x6c0695ed,
  0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950,
  0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3,
  0xfbd44c65, 0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2,
  0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a,
  0x346ed9fc, 0xad678846, 0xda60b8d0, 0x44042d73, 0x33031de5,
  0xaa0a4c5f, 0xdd0d7cc9, 0x5005713c, 0x270241aa, 0xbe0b1010,
  0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
  0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17,
  0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad, 0xedb88320, 0x9abfb3b6,
  0x03b6e20c, 0x74b1d29a, 0xead54739, 0x9dd277af, 0x04db2615,
  0x73dc1683, 0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8,
  0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1, 0xf00f9344,
  0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb,
  0x196c3671, 0x6e6b06e7, 0xfed41b76, 0x89d32be0, 0x10da7a5a,
  0x67dd4acc, 0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
  0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252, 0xd1bb67f1,
  0xa6bc5767, 0x3fb506dd, 0x48b2364b, 0xd80d2bda, 0xaf0a1b4c,
  0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55, 0x316e8eef,
  0x4669be79, 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
  0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe,
  0xb2bd0b28, 0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31,
  0x2cd99e8b, 0x5bdeae1d, 0x9b64c2b0, 0xec63f226, 0x756aa39c,
  0x026d930a, 0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713,
  0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38, 0x92d28e9b,
  0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21, 0x86d3d2d4, 0xf1d4e242,
  0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1,
  0x18b74777, 0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c,
  0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45, 0xa00ae278,
  0xd70dd2ee, 0x4e048354, 0x3903b3c2, 0xa7672661, 0xd06016f7,
  0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc, 0x40df0b66,
  0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
  0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605,
  0xcdd70693, 0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8,
  0x5d681b02, 0x2a6f2b94, 0xb40bbe37, 0xc30c8ea1, 0x5a05df1b,
  0x2d02ef8d
};

int crc32_v(const char *buf) {
  int crc;
  
  if (buf == NULL)
    return 0;

  crc = 0xffffffff;
  while (*buf) {
    crc = crc_table[((int) crc ^ (*buf++)) & 0xff] ^ (crc >> 8);
  }
  crc ^= 0xffffffff;
  
  return crc & CRC_HASH_MASK;
}

int crc32_dna(const char *buf) {
  int crc;
  
  if (buf == NULL)
    return 0;

  crc = 0xffffffff;

  int i = 0, a = 0;
  
  while (*buf) {
    a = (a << 2) | (int)*buf++;
    if (++i == 4) {
      crc = crc_table[((int) crc ^ a) & 0xff] ^ (crc >> 8);
      i = 0;
    }
  }
  crc = crc_table[((int) crc ^ a) & 0xff] ^ (crc >> 8);
  crc ^= 0xffffffff;
  
  return crc & CRC_HASH_MASK;
}

long long hash_collision_count = 0;
long long hash_query_count = 0;

#define MAXADLER 65521
hashiterator_t kiFindMatch(hashtable_t* hash, char *string) {
  /* Adler-32 checksum */
  unsigned int hashA = 1;
  unsigned int hashB = 0;

  int local_hash, bucket_count = hash->nBuckets;
  hashbucket_t *buckets = hash->buckets;
  int *bucket_string_hashes = hash->bucket_string_hashes;

  /* char *p; */
  /* for (p = string; *p != '\0'; p++) { */
  /*   hashA += (unsigned int)*p; */
  /*   hashB += hashA; */
  /* } */
  /* hashA %= MAXADLER; */
  /* hashB %= MAXADLER; */

  /* hashiterator_t hi = (hashB*65536+hashA) % bucket_count; */

  hashiterator_t hi = crc32_v(string) % bucket_count;
  /* hashiterator_t hi = crc32_dna(string) % bucket_count; */

  local_hash = calc_hash(string);
  hash_query_count++;
  while (bucket_string_hashes[hi] &&
         ((bucket_string_hashes[hi] != local_hash) || strcmp(buckets[hi].string, string))) {
    hi++;
    hash_collision_count++;
    if (hi >= bucket_count)
      hi = 0;
  }
  int first = kiHashFirst(hash, hi);
  kipmsg(5, "string = %s, first = %d\n", string, first);
  return hi;
}

float kiGetHashLoadFactor(hashtable_t* hash) {
  return 1. * ki_seqs->nSeq / KI_NUM_BUCKETS;
}

float kiGetHashCollisionRatio() {
  kipm("Hash (open addressing) collision ratio = %.3f\n", 1.0*hash_collision_count / hash_query_count);
  return 0.;
}

long kiGetHashTableMemory() {
  /* TODO */
  return 0;
}

char *kiGetHashString(hashtable_t* hash, hashiterator_t hi) {
  return hash->buckets[hi].string;
}

int kiHashCount(hashtable_t* hash, hashiterator_t hi) {
  return hash->buckets[hi].nCount;
}

int kiHashFirst(hashtable_t* hash, hashiterator_t hi) {
  return hash->buckets[hi].first;
}

origin_iterator_t kiHashOrigins(hashtable_t* hash, hashiterator_t hi) {
  return hash->buckets[hi].origins;
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
    if (it->prev != NULL) {       /* not root of the list */
      it->prev->next = it->next;
      if (it->next != NULL) {
        it->next->prev = it->prev;
      }
    } else {    /* root */
      hash->buckets[hi].origins = it->next;
      if (it->next != NULL) {
        it->next->prev = NULL;
      }
    }
    hash->buckets[hi].nCount--;
    origin_iterator_t rv = it->next;
    kifree(it, sizeof(origin_list_t));
    return rv;
  }
  
  return it->next;
}

#else
#endif
