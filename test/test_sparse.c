#include "stdbool.h"
#include "stdio.h"

#include "matf.h"
#include "sparse.h"

#include "test_util.h"

#define N_DIAGS 3
#define SIZE 10
#define N_VALS 30

int comp_long(const void *a, const void *b) {
  const long va = *(long *)a;
  const long vb = *(long *)b;
  return va > vb;
}

void gen_randoms(long *rows, long *cols, float *vals, long n_vals, long n_rows,
                 long n_cols) {
  long *idcs = malloc(n_vals * sizeof(long));

  for (long i = 0; i < n_vals; ++i) {
    bool repeat = true;
    while (repeat) {
      vals[i] = rand_float();
      idcs[i] = rand() % (n_rows * n_cols);

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
    cols[i] = idcs[i] % n_cols;
    rows[i] = idcs[i] / n_cols;
  }
  free(idcs);
}

SMatF gen_random(long n_vals, long n_rows, long n_cols) {
  float *vals = malloc(n_vals * sizeof(float));
  long *pos_row = malloc(n_vals * sizeof(long));
  long *pos_col = malloc(n_vals * sizeof(long));

  gen_randoms(pos_row, pos_col, vals, n_vals, n_rows, n_cols);

  SMatF ret = SM_from_pos_with(n_rows, n_cols, n_vals, pos_row, pos_col, vals);

  free(vals);
  free(pos_row);
  free(pos_col);

  return ret;
}

// Test matrix product of random matrix with identity. If result and input
// matrix are identical, the product should be working correctly.
void test_random_pos_prod() {
  SMatF identity = SM_diag_regular((long[1]){0}, (float[1]){1.0f}, 1, SIZE);

  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  SMatF target = SM_prod_prepare(A, identity);

  SM_prod(A, identity, target);

  if (SM_eq(A, target))
    TEST_PASS("Random matrix product with identity");
  else
    TEST_FAIL("Random matrix product with identity");

  SM_free(A);
  SM_free(identity);
  SM_free(target);
}

void test_random_pos_dif_sizes_prod_ver_matf() {
  float *values_A = malloc(N_VALS * sizeof(float));
  long *row_pos_A = malloc(N_VALS * sizeof(long));
  long *col_pos_A = malloc(N_VALS * sizeof(long));

  float *values_B = malloc(N_VALS * sizeof(float));
  long *row_pos_B = malloc(N_VALS * sizeof(long));
  long *col_pos_B = malloc(N_VALS * sizeof(long));

  long a_rows = SIZE;
  long a_cols = SIZE + 5;
  long b_rows = SIZE + 5;
  long b_cols = SIZE;

  gen_randoms(row_pos_A, col_pos_A, values_A, N_VALS, a_rows, a_cols);
  gen_randoms(row_pos_B, col_pos_B, values_B, N_VALS, b_rows, b_cols);

  SMatF s_A =
      SM_from_pos_with(a_rows, a_cols, N_VALS, row_pos_A, col_pos_A, values_A);
  SMatF s_B =
      SM_from_pos_with(b_rows, b_cols, N_VALS, row_pos_B, col_pos_B, values_B);

  SMatF s_T = SM_prod_prepare(s_A, s_B);
  SM_prod(s_A, s_B, s_T);

  MatF m_A = MF_with(a_rows, a_cols, 0.0f);
  MatF m_B = MF_with(b_rows, b_cols, 0.0f);
  MatF m_T = MF_empty(a_rows, b_cols);

  for (long i = 0; i < N_VALS; ++i) {
    *MF_PTR(m_A, row_pos_A[i], col_pos_A[i]) = values_A[i];
    *MF_PTR(m_B, row_pos_B[i], col_pos_B[i]) = values_B[i];
  }

  MF_prod(m_A, m_B, m_T);

  bool fail = false;
  for (long row = 0; row < m_T.rows && !fail; ++row) {
    for (long col = 0; col < m_T.cols && !fail; ++col) {
      if (SM_at(s_T, row, col) != MF_AT(m_T, row, col)) {
        TEST_FAIL_MSG("Matrix multiplication verified with MatF",
                      "Mismatch at (%ld, %ld)", row, col);
        fail = true;
      }
    }
  }

  if (!fail)
    TEST_PASS("Matrix multiplication verified with MatF");

  SM_free(s_A);
  SM_free(s_B);
  SM_free(s_T);

  MF_FREE(m_A);
  MF_FREE(m_B);
  MF_FREE(m_T);

  free(values_A);
  free(values_B);
  free(row_pos_A);
  free(row_pos_B);
  free(col_pos_A);
  free(col_pos_B);
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
  const SMatF A = gen_random(N_VALS, SIZE, SIZE);

  if (SM_eq(A, A))
    TEST_PASS("Random matrix self equality");
  else
    TEST_FAIL("Random matrix self equality");

  SM_free(A);
}

void test_structure_equality_self() {
  const SMatF A = gen_random(N_VALS, SIZE, SIZE);

  if (SM_structure_eq(A, A))
    TEST_PASS("Random structural self equality");
  else
    TEST_FAIL("Random structureal self equality");

  SM_free(A);
}

void test_add_sub_random() {
  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  SMatF B = gen_random(N_VALS, SIZE, SIZE);

  SMatF target = SM_addsub_prepare(A, B);
  SM_add(A, B, target);

  bool fail = false;
  for (long row = 0; row < SIZE && !fail; ++row) {
    for (long col = 0; col < SIZE && !fail; ++col) {
      float t = SM_at(target, row, col);
      float a = SM_at(A, row, col);
      float b = SM_at(B, row, col);

      if (t != a + b) {
        TEST_FAIL_MSG("Matrix Addition", "Failed at (%ld, %ld) (%f + %f != %f)",
                      row, col, a, b, t);
        fail = true;
      }
    }
  }

  if (!fail)
    TEST_PASS("Matrix Addition");

  SM_sub(A, B, target);

  fail = false;
  for (long row = 0; row < SIZE && !fail; ++row) {
    for (long col = 0; col < SIZE && !fail; ++col) {
      float t = SM_at(target, row, col);
      float a = SM_at(A, row, col);
      float b = SM_at(B, row, col);

      if (t != a - b) {
        TEST_FAIL_MSG("Matrix Subtraction",
                      "Failed at (%ld, %ld) (%f - %f != %f)", row, col, a, b,
                      t);
        fail = true;
      }
    }
  }

  if (!fail)
    TEST_PASS("Matrix Subtraction");

  float scale_fac = rand_float() * 10.0f;
  SMatF target2 = SM_empty_like(A);
  SM_scl(A, scale_fac, target2);

  fail = false;
  for (long row = 0; row < SIZE && !fail; ++row) {
    for (long col = 0; col < SIZE && !fail; ++col) {
      float t = SM_at(target2, row, col);
      float a = SM_at(A, row, col);

      if (t != a * scale_fac) {
        TEST_FAIL_MSG("Matrix Scaling", "Failed at (%ld, %ld) (%f * %f != %f)",
                      row, col, a, scale_fac, t);
        fail = true;
      }
    }
  }
  if (!fail)
    TEST_PASS("Matrix Scaling");

  SM_free(A);
  SM_free(B);
  SM_free(target);
  SM_free(target2);
}

void test_transpose() {
  float *values = malloc(N_VALS * sizeof(float));
  long *row_pos = malloc(N_VALS * sizeof(long));
  long *col_pos = malloc(N_VALS * sizeof(long));
  long a_rows = SIZE;
  long a_cols = SIZE + 5;

  gen_randoms(row_pos, col_pos, values, N_VALS, a_rows, a_cols);
  SMatF A = SM_from_pos_with(a_rows, a_cols, N_VALS, row_pos, col_pos, values);

  SMatF T = SM_transpose(A);

  bool fail = false;
  for (long i = 0; i < A.nrows && !fail; ++i) {
    for (long j = 0; j < A.ncols && !fail; ++j) {
      const float a = SM_at(A, i, j);
      const float t = SM_at(T, j, i);
      if (a != t) {
        fail = true;
        TEST_FAIL_MSG(
            "Matrix Transposition",
            "Value at (%ld, %ld) in A does not equal value at (%ld, %ld) "
            "in T (%.2f != %.2f)",
            i, j, j, i, a, t);
      }
    }
  }

  if (!fail)
    TEST_PASS_MSG("Matrix Transposition", "Value Equality passed%s", "");

  SMatF TT = SM_transpose(T);

  fail = !SM_eq(A, TT) || fail;

  if (!fail) {
    TEST_PASS("Matrix Transposition");
  } else {
    TEST_FAIL_MSG("Matrix Transposition",
                  "Double does not equal original Matrix %s", "");
  }

  SM_free(A);
  SM_free(T);
  SM_free(TT);
  free(values);
  free(row_pos);
  free(col_pos);
}

void test_scalar_prod_self() {
  SMatF vec = SM_vec_empty(SIZE);
  for (long i = 0; i < vec.nvals; ++i) {
    vec.vals[i] = rand_float();
  }

  float abs2 = SM_scalar(vec, vec);

  float tst = 0;
  for (long i = 0; i < vec.nvals; ++i) {
    tst += vec.vals[i] * vec.vals[i];
  }

  if (tst != abs2)
    TEST_FAIL("Scalar Product");
  else
    TEST_PASS("Scalar Product");

  SM_free(vec);
}

int main(void) {
  srand(69);
  test_structure_equality_self();
  test_equality_self();
  test_tridiag_vec_prod();
  test_random_pos_dif_sizes_prod_ver_matf();
  test_random_pos_prod();
  test_add_sub_random();
  test_transpose();
  test_scalar_prod_self();
  return EXIT_SUCCESS;
}
