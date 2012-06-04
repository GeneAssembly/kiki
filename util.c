#include "util.h"

#include <sys/time.h>
#include <time.h>

void getCurrentTime(char* timestr) {
  struct timeval tv;
  time_t curtime;
  gettimeofday(&tv, NULL);
  curtime=tv.tv_sec;
  strftime(timestr, 30, "%m-%d-%Y %T", localtime(&curtime));
}
