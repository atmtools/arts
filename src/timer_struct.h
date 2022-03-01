#ifndef timer_struct_h
#define timer_struct_h

#ifdef TIME_SUPPORT
#include <sys/times.h>
#endif

struct Timer {
  bool running{false};
  bool finished{false};
#ifdef TIME_SUPPORT
  struct tms cputime_start;
  clock_t realtime_start;
  struct tms cputime_end;
  clock_t realtime_end;
#endif
};

#endif
