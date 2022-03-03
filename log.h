#ifndef LOGGING_H
#define LOGGING_H

typedef enum {
  VERB_ERR = 0,
  VERB_WRN,
  VERB_DBG,
} PrintVerbosity;

#define VERB_LEVEL VERB_ERR

#define log_msg(...) _log(VERB_DBG, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define log_wrn(...) _log(VERB_WRN, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define log_err(...) _log(VERB_ERR, __func__, __FILE__, __LINE__, __VA_ARGS__)
void _log(PrintVerbosity level, const char* func, const char *file, int line, const char *fmt, ...);

#endif // LOGGING_H
