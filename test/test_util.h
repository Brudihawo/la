#include "stdio.h"
#include "stdlib.h"

#define RED   "\033[1;31m"
#define RED_N "\033[0;31m"
#define GRN   "\033[1;32m"
#define GRN_N "\033[0;32m"
#define RST   "\033[0m"

float rand_float() { return (float)rand() / (float)RAND_MAX; }

#define TEST_PASS(name) \
  fprintf(stdout, GRN "PASSED:"RST" %s in %s\n  %s:%i\n",\
      name, __func__, __FILE__, __LINE__)
#define TEST_FAIL(name) \
  fprintf(stdout, RED "FAILED:"RST" %s in %s\n  %s:%i\n",\
      name, __func__, __FILE__, __LINE__)

#define TEST_FAIL_MSG(name, msg, ...) \
  fprintf(stdout, RED "FAILED:"RST" %s in %s\n  %s:%i\n  " msg "\n",\
      name, __func__, __FILE__, __LINE__, __VA_ARGS__)

#define TEST_PASS_MSG(name, msg, ...) \
  fprintf(stdout, GRN "PASSED:"RST" %s in %s\n  %s:%i\n  " msg "\n",\
      name, __func__, __FILE__, __LINE__, __VA_ARGS__)

#define TEST(name, cond) cond; cond ? TEST_PASS(name) : TEST_FAIL(name);
