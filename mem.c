
#include "extern.h"
#include "debug.h"
#include "ki.h"
#include "mem.h"

/* Memory */
long ki_malloc_all      = 0;
long ki_malloc_used     = 0;	/* useful allocations by kimalloc */
long ki_max_malloc_heap = 0;    /* maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */


/* Memory operations */
void *kimalloc(size_t sz) {
  if (sz == 0) return NULL;
  void* new = malloc(sz);
  if (new == NULL) {
    fprintf(stderr, "Out of memory\n");
    kiAbort(1);
  }
  ki_malloc_all  += sz;
  ki_malloc_used += sz;
  kipmsg(6, "kimalloc = +%d\n", sz);
#ifdef KI_TRACK_MEMORY
  struct mallinfo mi = mallinfo();
  if (mi.arena+mi.hblkhd > ki_max_malloc_heap)
    ki_max_malloc_heap = mi.arena+mi.hblkhd;
#endif
  /* gcc malloc should always return 16-byte-aligned values... */
  /* assert(IS_ALIGNED(new)); */
  return new;
}

void *kimemdup(void *data, size_t sz) {
  if (data == NULL) return NULL;
  void* new = kimalloc(sz);
  memcpy(/*to*/new, /*from*/data, sz);
  return new;
}

char* kistrdup(char* str) {
  if (str == NULL) return NULL;
  return kimemdup(str, strlen(str)+1);
}

char* kistrndup(char* str, int n) {
  if (str == NULL) return NULL;
  assert(n <= strlen(str));
  char* s = kimemdup(str, n+1);
  s[n] = '\0';
  return s;
}

void *kirealloc(void *data, size_t szOld, size_t szNew, bool bCopy) {
  if (data == NULL && szOld == 0)
    return(kimalloc(szNew));
  if (data == NULL || szOld == 0 || szNew == 0) {
    fprintf(stderr,"Empty kirealloc\n");
    kiAbort(1);
  }
  if (szOld == szNew)
    return(data);
  void *new = NULL;
  if (bCopy) {
    /* Try to reduce memory fragmentation by allocating anew and copying
       Seems to help in practice */
    new = kimemdup(data, szNew);
    kifree(data, szOld);
  } else {
    new = realloc(data,szNew);
    if (new == NULL) {
      fprintf(stderr, "Out of memory\n");
      kiAbort(1);
    }
    /* assert(IS_ALIGNED(new)); */
    ki_malloc_all  += (szNew-szOld);
    ki_malloc_used += (szNew-szOld);
    kipmsg(6, "kirealloc = +%d\n", szNew-szOld);

#ifdef KI_TRACK_MEMORY
    struct mallinfo mi = mallinfo();
    if (mi.arena+mi.hblkhd > ki_max_malloc_heap)
      ki_max_malloc_heap = mi.arena + mi.hblkhd;
#endif
  }
  return(new);
}

void *kifree(void *p, size_t sz) {
  if (p==NULL) return(NULL);
  free(p);
  ki_malloc_used -= sz;
  kipmsg(6, "kifree = -%d\n", sz);
  
  return(NULL);
}

char* kifreestr(char* str) {
  if (str == NULL) return NULL;
  kifree(str, strlen(str)+1);
  return NULL;
}


long kiGetMallocUsed() {
  return ki_malloc_used;
}

long kiGetMallocAll() {
  return ki_malloc_all;
}

long kiGetMaxMallocHeap() {
  return ki_max_malloc_heap;
}


/* Customized memory allocator */
#define KI_ARENA_PAGE_SIZE 0x01000000   /* 16 MB */
/* #define KI_ARENA_NUM_PAGES 64 */
/* #define KI_ARENA_NUM_PAGES 256 */
#define KI_ARENA_NUM_PAGES 6400

arena_t ki_arena[KI_ARENA_MAX];

void kiArenaInit(int n) {
  if (!KI_USE_ARENA) return;

  int id;
  for (id = 0; id < n; ++id) {
    arena_t* a   = &(ki_arena[id]);
    a->pageSize  = KI_ARENA_PAGE_SIZE;
    a->pageCount = KI_ARENA_NUM_PAGES;
    a->nSaved    = 0;
    a->pageNo    = 0;
    a->top       = 0;
    a->pages     = (void**)kimalloc(sizeof(void*) * a->pageCount);
    memset(a->pages, 0, sizeof(void*) * a->pageCount);
  }
}

void kiArenaFinalize(int n) {
  if (!KI_USE_ARENA) return;

  int id, i;
  for (id = 0; id < n; ++id) {
    arena_t* a = &(ki_arena[id]);
    for (i = 0; i < a->nSaved; ++i) {
      a->pages[i] = kifree(a->pages[i], a->pageSize);
    }
    a->pages  = kifree(a->pages, sizeof(void*) * a->pageCount);
    a->nSaved = 0;
    a->pageNo = 0;
    a->top    = 0;
  }
}

void* kiArenaMalloc(int id, size_t sz) {
  if (!KI_USE_ARENA) return kimalloc(sz);

  arena_t* a = &(ki_arena[id]);
  assert(sz < a->pageSize);
  if (a->pageNo >= a->pageCount - 1 && a->top + sz >= a->pageSize) {
    fprintf(stderr, "Out of arena memory!\n");
    kiAbort(1);
  }
  if (a->nSaved == 0) {
    a->pages[a->nSaved++] = kimalloc(a->pageSize);
  } else if (a->top + sz >= a->pageSize) {
    if (a->pageNo == a->nSaved - 1) {
      void* p = kimalloc(a->pageSize);
      memset(p, 0, a->pageSize);
      a->pages[a->nSaved++] = p;
    }
    a->pageNo++;
    a->top = 0;
  }

  /* always allocate from the top of the stack */
  void* rv = a->pages[a->pageNo] + a->top;
  
  size_t alignedSize = sz;         /* 16-byte-aligned values */
  if (sz > 0) alignedSize = (((sz - 1) >> 4) + 1) << 4;
  
  a->top += alignedSize;
  
  /* kipm("pageNo = %d\t top = %d\t sz = %d\t alignedSize = %d\n", a->pageNo, a->top, sz, alignedSize); */
  return rv;
}

void* kiArenaFree(int id, void* data, size_t sz) {
  if (!KI_USE_ARENA) return kifree(data, sz);

  return NULL;
}

void kiArenaClear(int id) {
  if (!KI_USE_ARENA) return;

  arena_t* a = &(ki_arena[id]);
  int i;
  for (i = 0; i < a->pageNo; ++i) {
    memset(a->pages[i], 0, a->pageSize);
  }
  memset(a->pages[a->pageNo], 0, a->top);
  a->pageNo = 0;
  a->top = 0;
}

void* kiArenaRealloc(int id, void *data, size_t szOld, size_t szNew, bool bCopy) {
  if (!KI_USE_ARENA) return kirealloc(data, szOld, szNew, bCopy);
  
  void* new = kiArenaMalloc(id, szNew);
  memcpy(new, data, szOld);
  return new;
}

void* kiArenaMemdup(int id, void *data, size_t sz) {
  if (data == NULL) return NULL;
  void* new = kiArenaMalloc(id, sz);
  memcpy(/*to*/new, /*from*/data, sz);
  return new;
}

char* kiArenaStrdup(int id, char* str) {
  if (str == NULL) return NULL;
  return kiArenaMemdup(id, str, strlen(str)+1);
}

char* kiArenaStrndup(int id, char* str, int n) {
  if (str == NULL) return NULL;
  assert(n <= strlen(str));
  char* s = kiArenaMemdup(id, str, n+1);
  s[n] = '\0';
  return s;
}

char* kiArenaFreestr(int id, char* str) {
  if (!KI_USE_ARENA) return kifreestr(str);
  return NULL;
}

size_t kiAlignedSize(size_t sz, int bits) {
  assert(bits > 0);
  if (sz > 0) return (((sz - 1) >> bits) + 1) << bits;
  return 0;
}

int kiCeilingDevidedBy(int x, int d) {
  assert(d > 0);
  if (x > 0) return (((x - 1) / d) + 1);
  return 0;
}

