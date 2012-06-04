#ifndef _MEM_H_
#define _MEM_H_

#ifndef __APPLE__
#define KI_TRACK_MEMORY  1
#endif

#ifdef KI_TRACK_MEMORY
/* malloc.h apparently doesn't exist on MacOS */
#include <malloc.h>
#endif

#define KI_USE_ARENA     1

#include <stdbool.h>

/* Memory */
void* kimalloc(size_t sz);       /* prints "Out of memory" and exits on failure */
void* kifree(void *, size_t sz); /* always returns NULL */
char* kifreestr(char* str);
void* kimemdup(void *data, size_t sz);
char* kistrdup(char* str);
char* kistrndup(char* str, int n);
void* kirealloc(void *data, size_t szOld, size_t szNew, bool bCopy);

long kiGetMallocUsed();
long kiGetMallocAll();
long kiGetMaxMallocHeap();

/* Customized memory allocator */
typedef struct {
  void** pages;
  int pageSize;
  int pageCount;
  int nSaved;                   /* number of allocated pages */
  int pageNo;
  int top;                      /* (pageNo:top): pointer to the unallocated memory */
} arena_t;

void  kiArenaInit(int maxArena);
void  kiArenaFinalize(int maxArena);
void* kiArenaMalloc(int id, size_t sz);
void* kiArenaFree(int id, void* data, size_t sz);
void  kiArenaClear(int id);

void* kiArenaRealloc(int id, void *data, size_t szOld, size_t szNew, bool bCopy);
void* kiArenaMemdup(int id, void *data, size_t sz);
char* kiArenaStrdup(int id, char* str);
char* kiArenaStrndup(int id, char* str, int n);
char* kiArenaFreestr(int id, char* str);

size_t kiAlignedSize(size_t sz, int bits);
int kiCeilingDevidedBy(int x, int d);

#endif /* _MEM_H_ */
