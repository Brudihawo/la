#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"

#include "sparse.h"
#include "test_util.h"

#define N_DIAGS 3
#define SIZE 10
#define N_VALS 30

float rand_float() { return (float)rand() / (float)RAND_MAX; }

void gen_randoms(long *rows, long *cols, float *vals, long n_vals) {
  for (long i = 0; i < n_vals; ++i) {
    bool repeat = true;
    while (repeat) {
      vals[i] = rand_float();
      cols[i] = rand() % SIZE;
      rows[i] = rand() % SIZE;
      repeat = false;
      for (long j = 0; j < i; j++) {
        if (cols[j] == cols[i] && rows[j] == rows[i]) {
          repeat = true;
        }
      }
    }
  }
}

void test_random_pos_prod() {
  float *vals = malloc(N_VALS * sizeof(float));
  long *pos_row = malloc(N_VALS * sizeof(long));
  long *pos_col = malloc(N_VALS * sizeof(long));

  gen_randoms(pos_row, pos_col, vals, N_VALS);

  SM_from_pos_with(SIZE, SIZE, N_VALS, pos_row, pos_col, vals);

  free(vals);
  free(pos_row);
  free(pos_row);
}

void test_tridiag_vec_prod() {
  SMatF tridiag =
      SM_diag_regular((long[N_DIAGS]){-1, 0, 1},
                      (float[N_DIAGS]){1.0f, 2.0f, 3.0f}, N_DIAGS, SIZE);

  SMatF vec = SM_vec_empty(SIZE);

  for (int i = 0; i < SIZE; ++i) {
    if (i == 4)
      SM_set_or_panic(vec, i, 0, 1.0f);
    else
      SM_set_or_panic(vec, i, 0, 0.0f);
  }

  SMatF res = SM_prod_prepare(tridiag, vec);

  SM_prod(tridiag, vec, res);

  bool fail = false;
  long failed_at = -1;
  float expected;
  for (long i = 0; i < SIZE; ++i) {
    if (i > 2 && i < 6) {
      if (SM_at(res, i, 0) != (float)(6 - i)) {
        fail = true;
        expected = (float)(6 - i);
        TEST_FAIL_MSG(
            "Matrix vector multiplication",
            "Unexpected value in output at position %ld: %.2f. Expected %.2f.",
            failed_at, SM_at(vec, i, 0), expected);
        break;
      }
    } else {
      if (SM_at(res, i, 0) != 0.0f) {
        break;
        failed_at = i;
        expected = 0.0f;
        TEST_FAIL_MSG(
            "Matrix vector multiplication",
            "Unexpected value in output at position %ld: %.2f. Expected %.2f.",
            failed_at, SM_at(vec, i, 0), expected);
        fail = true;
      }
    }
  }

  if (!fail)
    TEST_PASS("Matrix vector multiplication");

  SM_free(vec);
  SM_free(tridiag);
  SM_free(res);
}

int main(void) {
  test_tridiag_vec_prod();
  return EXIT_SUCCESS;
}
