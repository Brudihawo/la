#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"

#define RED   "\033[1;31m"
#define RED_N "\033[0;31m"
#define GRN   "\033[1;32m"
#define GRN_N "\033[0;32m"
#define RST   "\033[0m"


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

typedef bool (*test_func)(void);

void _run_tests(test_func* fn, long n_tests, const char* file) {
  fprintf(stderr, "========================================================\n");
  fprintf(stderr, "Running %ld tests in %s\n", n_tests, file);
  for (long i = 0; i < n_tests; ++i) {
    if (!fn[i]()) {
      fprintf(stderr, "Failed at Test %ld / %ld\n", i, n_tests);
      exit(EXIT_FAILURE);
    }
  }
  fprintf(stderr, "Passed all (%ld) tests in %s\n", n_tests, file);
}

#define run_tests(fns) _run_tests(fns, sizeof(fns) / sizeof(fns[0]), __FILE__);

