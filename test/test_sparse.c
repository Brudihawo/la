#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"

#include "sparse.h"
#include "test_util.h"

#define N_DIAGS 3
#define SIZE 10
#define N_VALS 30

float rand_float() { return (float)rand() / (float)RAND_MAX; }

int comp_long(const void* a, const void* b) {
  long va = *(long*)a;
  long vb = *(long*)b;
  return va > vb;
}

void gen_randoms(long *rows, long *cols, float *vals, long n_vals) {
  long* idcs = malloc(n_vals * sizeof(long));

  for (long i = 0; i < n_vals; ++i) {
    bool repeat = true;
    while (repeat) {
      vals[i] = rand_float();
      idcs[i] = rand() % (SIZE * SIZE);

      repeat = false;
      for (long j = 0; j < i; j++) {
        if (idcs[j] == idcs[i]) {
          repeat = true;
        }
      }
    }
  }
  qsort(idcs, n_vals, sizeof(long), comp_long);

  for (long i = 0; i < n_vals; ++i) {
    cols[i] = idcs[i] % SIZE;
    rows[i] = idcs[i] / SIZE;
  }
  free(idcs);
}

SMatF gen_random_square(long size, long n_vals) {
  float *vals = malloc(n_vals * sizeof(float));
  long *pos_row = malloc(n_vals * sizeof(long));
  long *pos_col = malloc(n_vals * sizeof(long));

  gen_randoms(pos_row, pos_col, vals, n_vals);

  SMatF ret = SM_from_pos_with(size, size, n_vals, pos_row, pos_col, vals);

  free(vals);
  free(pos_row);
  free(pos_col);

  return ret;
}

// Test matrix product of random matrix with identity. If result and input matrix are
// identical, the product should be working correctly.
void test_random_pos_prod() {
  SMatF identity = SM_diag_regular((long[1]){ 0 }, (float[1]){ 1.0f }, 1, SIZE);

  SMatF A = gen_random_square(SIZE, N_VALS);
  SMatF target = SM_prod_prepare(A, identity);

  SM_prod(A, identity, target);

  if (SM_eq(A, target))
    TEST_PASS("Random matrix product with identity");
  else TEST_FAIL("Random matrix product with identity");

  SM_free(A);
  SM_free(identity);
  SM_free(target);
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

void test_equality_self() {
  const SMatF A = gen_random_square(SIZE, N_VALS);
  
  if (SM_eq(A, A))
    TEST_PASS("Random matrix self equality");
  else TEST_FAIL("Random matrix self equality");

  SM_free(A);
}

void test_structure_equality_self() {
  const SMatF A = gen_random_square(SIZE, N_VALS);
  
  if (SM_structure_eq(A, A))
    TEST_PASS("Random structural self equality");
  else TEST_FAIL("Random structureal self equality");

  SM_free(A);
}

int main(void) {
  srand(69);
  test_structure_equality_self();
  test_equality_self();
  test_tridiag_vec_prod();
  test_random_pos_prod();
  return EXIT_SUCCESS;
}
