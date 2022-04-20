#include "la.h"
#include "limits.h"
#include "test_util.h"

#define N_ITER 9
#define MAX_SZ 100
#define MAX_VAL 10000

void print_arr(long *arr, long sz) {
  printf("%4ld | ", (long)0);
  for (long i = 0; i < sz; ++i) {
    if ((i % 10) == 0 && i != 0)
      printf("\n%4ld | ", i / 10);
    printf("%3ld ", arr[i]);
  }
  printf("\n");
}

long random_arr(long **arr, long size) {
  *arr = malloc(size * sizeof(long));

  for (long i = 0; i < size; ++i) {
    (*arr)[i] = rand() % MAX_VAL;
  }

  radix_sort(*arr, size, NULL);

  long real_sz = size;
  for (long i = 1; i < size; ++i) {
    if ((*arr)[i] == (*arr)[i - 1]) {
      --real_sz;
      (*arr[i - 1]) = MAX_VAL + 1;
    }
  }

  radix_sort(*arr, size, NULL);
  return real_sz;
}

bool test_bin_search_random() {
  bool fail = false;
  for (int i = 0; i < N_ITER && !fail; ++i) {
    const long size = rand() % MAX_SZ;

    long *arr;
    const long real_sz = random_arr(&arr, size);

    for (long i = 0; i < N_ITER && !fail; ++i) {
      const long idx = rand() % real_sz;
      const long val = arr[idx];

      const long out = binary_search(arr, real_sz, val);
      if (out != idx) {
        TEST_FAIL_MSG("binary_search",
                      "could not find value %ld at position %ld", val, idx);
        fail = true;
        print_arr(arr, real_sz);
      }
    }

    free(arr);
  }

  if (!fail)
    TEST_PASS("binary_search");

  return !fail;
}

bool test_bin_search_edge() {
  long *arr;
  const long real_sz = random_arr(&arr, 10);

  if (binary_search(arr, real_sz, arr[0]) < 0) {
    TEST_FAIL("binary_search first element");
    return false;
  }

  if (binary_search(arr, real_sz, arr[real_sz - 1]) < 0) {
    TEST_FAIL("binary_search first element");
    return false;
  }

  if (binary_search(arr, real_sz, arr[real_sz / 2]) < 0) {
    TEST_FAIL("binary_search middle element");
    return false;
  }
  TEST_PASS("binary_search edge cases");
  return true;
}

int main() {
  srand(69);

  test_func funcs[] = {
      test_bin_search_random,
      test_bin_search_edge,
  };

  run_tests(funcs);
}
