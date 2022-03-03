#include "stdarg.h"
#include "stdio.h"

#include "log.h"

void _log(PrintVerbosity level, const char* func, const char *file, int line, const char *fmt, ...) {
  if (level <= VERB_LEVEL) {
    va_list args;
    fprintf(stderr, "%s:%i (%s): ", file, line, func);
    va_start(args, fmt);

    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");

    va_end(args);
  }
}
