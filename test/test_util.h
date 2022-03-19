#include "stdio.h"

#define RED   "\033[1;31m"
#define RED_N "\033[0;31m"
#define GRN   "\033[1;32m"
#define GRN_N "\033[0;32m"
#define RST   "\033[0m"

#define TEST_PASS(name) \
  fprintf(stdout, GRN "PASSED: %s in %s\n"GRN_N"  %s:%i"RST"\n",\
      name, __func__, __FILE__, __LINE__)
#define TEST_FAIL(name) \
  fprintf(stdout, RED "FAILED: %s in %s\n"RED_N"  %s:%i"RST"\n",\
      name, __func__, __FILE__, __LINE__)

#define TEST_FAIL_MSG(name, msg, ...) \
  fprintf(stdout, RED "FAILED: %s in %s\n"RED_N"  %s:%i"RST"\n  " msg "\n",\
      name, __func__, __FILE__, __LINE__, __VA_ARGS__)

#define TEST(name, cond) cond; cond ? TEST_PASS(name) : TEST_FAIL(name);
