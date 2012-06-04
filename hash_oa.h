#ifndef _HASH_H_
#define _HASH_H_

#include "ki.h"

/* Hashtable functions */
/* Open addressing with linear probing (n = 1) */

typedef struct {
  char* string;
  int nCount;			/* number of times this entry was seen */
  int first;			/* index of first entry with this value */
  origin_list_t* origins;
  /* kmerorigin_t* origins; */
} hashbucket_t;

typedef struct { 
  int nBuckets;
  /* hashvalue -> bucket. Or look in bucket + 1, +2, etc., till you hit a NULL string */
  hashbucket_t *buckets;
  int *bucket_string_hashes;
  bool bCopy;           /* string in bucket by ref or copy */
  bool bProtein;
} hashtable_t;
typedef int hashiterator_t;


hashtable_t *kiMakeHashtable(char** strings, int nStrings, int nBuckets, bool bCopy, bool bProtein);
hashtable_t *kiFreeHashtable(hashtable_t* hash); /*returns NULL*/
hashiterator_t kiFindMatch(hashtable_t* hash, char* str);
void kiClearHashtable(hashtable_t* hash);
void kiHashtableAdd(hashtable_t* hash, char* str, int seqid);
void kiHashtableAddCopy(hashtable_t* hash, char* str, size_t sz, int seqid, int offset); /* reverse seq if seqid < 0 */
void kiHashtableOptimize(hashtable_t* hash);

char* kiGetHashString(hashtable_t* hash, hashiterator_t hi); /* returns NULL if we have run out of values */
int   kiHashCount(hashtable_t* hash, hashiterator_t hi);
int   kiHashFirst(hashtable_t* hash, hashiterator_t hi);


void  kiHashGetMatchingSeq(origin_iterator_t it, /*OUT*/char* match);
origin_iterator_t kiHashFetchUnusedMatchingSeq(hashtable_t* hash, hashiterator_t hi, origin_iterator_t it, bool bErase, /*OUT*/char* match, /*OUT*/char** name);
origin_iterator_t kiHashEraseMatchingSeq(hashtable_t* hash, hashiterator_t hi, origin_iterator_t it, bool bErase); /* returns the iterator pointing to the next element */
origin_iterator_t kiHashOrigins(hashtable_t* hash, hashiterator_t hi);
float kiGetHashLoadFactor(hashtable_t* hash);
float kiGetHashCollisionRatio();
long kiGetHashTableMemory();


#endif /* _HASH_H_ */
