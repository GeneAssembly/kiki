#ifndef _DEBUG_H_
#define _DEBUG_H_

/* Usage:

   #include "debug.h"
   pmSetLevel(2);
   ...
   pmsg(3, "var = %d", var);
   pm("str = $s", str);
*/


#if defined(NDEBUG) && defined(__GNUC__)
/* gcc's cpp has extensions; it allows for macros with a variable number of
   arguments. We use this extension here to preprocess pmesg away. */
#define pmsg(level, format, args...) ((void)0)
#define pm(format, args...)          ((void)0)
#define pmSetLevel(level)            ((void)0)

#else
void pmsg(int level, char* format, ...);
void pm(char* format, ...);
void pmSetLevel(int level);

/* print a message, if it is considered significant enough.
   Adapted from [K&R2], p. 174 */
#endif

#endif /* _DEBUG_H_ */
