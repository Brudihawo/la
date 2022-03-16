#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"

#include "sparse.h"

#define N_DIAGS 3
#define SIZE 10
#define N_VALS 30

float rand_float() { return (float)rand() / (float)RAND_MAX; }

void gen_randoms(long* cols, long* rows, float* vals, long n_vals) {
  for (long i = 0; i < n_vals; ++i) {
    bool repeat = true;
    while (repeat) {
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
}

void test_tridiag_vec_prod() {
  SMatF tridiag = SM_diag_regular((long[N_DIAGS]){-1, 0, 1},
                                  (float[N_DIAGS]){1.0f, 2.0f, 3.0f}, N_DIAGS, SIZE);

  SMatF vec = SM_vec_empty(SIZE); 

  for (int i = 0; i < SIZE; ++i) {
    if (i == 4) SM_set_or_panic(vec, i, 0, 1.0f);
    else SM_set_or_panic(vec, i, 0, 0.0f);
  }

  SMatF res = SM_prod_prepare(tridiag, vec);

  SM_prod(tridiag, vec, res);

  // SM_print_nonzero(tridiag);
  SM_print(tridiag);
  SM_print(vec);
  SM_print(res);
}

int main(void) {
  test_tridiag_vec_prod();
  return EXIT_SUCCESS;
}
