#include "debug.h"
#include <stdio.h>
#include <stdarg.h>


#if defined(NDEBUG) && defined(__GNUC__)
/* Nothing. pmesg has been "defined away" in debug.h already. */
#else
/* the higher, the more messages... */
int debug_pmsg_level     = 2;
int debug_pmsg_level_def = 2;

void pmsg(int level, char* format, ...) {
#ifdef NDEBUG
  /* Empty body, so a good compiler will optimise calls
     to pmesg away */
#else
  va_list args;

  if (level > debug_pmsg_level)
    return;
  
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
#endif /* NDEBUG */
}
#endif /* NDEBUG && __GNUC__ */


#if defined(NDEBUG) && defined(__GNUC__)
#else
void pm(char* format, ...) {
#ifdef NDEBUG
#else
  va_list args;

  if (debug_pmsg_level < debug_pmsg_level_def)
    return;
  
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
#endif /* NDEBUG */
}
#endif /* NDEBUG && __GNUC__ */


#if defined(NDEBUG) && defined(__GNUC__)
#else
void pmSetLevel(int level) {
#ifdef NDEBUG
#else
  debug_pmsg_level = level;
#endif /* NDEBUG */
}
#endif /* NDEBUG && __GNUC__ */


