#include "stdbool.h"
#include "stdio.h"
#include "assert.h"

#include "matf.h"
#include "sparse.h"

#include "test_util.h"

#define N_DIAGS 3
#define SIZE 96
#define N_VALS 5

void gen_randoms(long *rows, long *cols, float *vals, long n_vals, long n_rows,
                 long n_cols) {
  assert(n_vals < n_rows * n_cols);
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
  qsort(idcs, n_vals, sizeof(long), long_gt);

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

bool test_init_from_pos_with(void) {
  float *vals = malloc(N_VALS * sizeof(float));
  long *pos_row = malloc(N_VALS * sizeof(long));
  long *pos_col = malloc(N_VALS * sizeof(long));

  const long nrows = SIZE;
  const long ncols = SIZE;
  gen_randoms(pos_row, pos_col, vals, N_VALS, nrows, ncols);

  SMatF A = SM_from_pos_with(nrows, ncols, N_VALS, pos_row, pos_col, vals);

  bool pass = true;
  long idx = 0;
  for (long row = 0; row < A.nrows; ++row) {
    for (long col = 0; col < A.ncols; ++col) {
      if (SM_has_loc(A, row, col)) {
        if (SM_at(A, row, col) != vals[idx]) {
          TEST_FAIL_MSG("Initialisation from positions",
                        "Wrong value at index %ld", idx);
          pass = false;
        }
        ++idx;
      } else {
        if (SM_at(A, row, col) != 0.0f) {
          TEST_FAIL_MSG("Initialisation from positions",
                        "Nonzero Value at (%ld, %ld)", row, col);
          pass = false;
        }
      }
    }
  }

  if (idx != A.nvals) {
    TEST_FAIL_MSG("Initialisation from positions",
                  "Number of values is wrong%s", "");
  } else {
    if (pass)
      TEST_PASS("Initialisation from positions");
  }

  free(vals);
  free(pos_row);
  free(pos_col);
  SM_free(A);
  return pass;
}

// Test matrix product of random matrix with identity. If result and input
// matrix are identical, the product should be working correctly.
bool test_random_pos_prod(void) {
  SMatF identity = SM_diag_regular((long[1]){0}, (float[1]){1.0f}, 1, SIZE);

  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  SMatF target = SM_prod_prepare(A, identity);

  SM_prod(A, identity, target);

  bool pass = true;
  if (SM_structure_eq(A, target))
    TEST_PASS("Random matrix product with identity - Structural equality");
  else {
    TEST_FAIL("Random matrix product with identity - Structural equality");
    pass = false;
  }

  if (SM_eq(A, target))
    TEST_PASS("Random matrix product with identity - Value equality");
  else {
    TEST_FAIL("Random matrix product with identity - Value equality");
    pass = false;
  }

  SM_free(A);
  SM_free(identity);
  SM_free(target);
  return pass;
}

bool test_random_pos_dif_sizes_prod_ver_matf(void) {
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

  bool pass = true;
  for (long row = 0; row < m_T.rows && pass; ++row) {
    for (long col = 0; col < m_T.cols && pass; ++col) {
      if (SM_at(s_T, row, col) != MF_AT(m_T, row, col)) {
        TEST_FAIL_MSG("Matrix multiplication verified with MatF",
                      "Mismatch at (%ld, %ld)", row, col);
        pass = false;
      }
    }
  }

  if (pass)
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
  return pass;
}

bool test_tridiag_vec_prod(void) {
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

  bool pass = true;
  long failed_at = -1;
  float expected;
  for (long i = 0; i < SIZE; ++i) {
    if (i > 2 && i < 6) {
      if (SM_at(res, i, 0) != (float)(6 - i)) {
        expected = (float)(6 - i);
        TEST_FAIL_MSG(
            "Matrix vector multiplication",
            "Unexpected value in output at position %ld: %.2f. Expected %.2f.",
            failed_at, SM_at(vec, i, 0), expected);
        pass = false;
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
        pass = false;
      }
    }
  }

  if (pass)
    TEST_PASS("Matrix vector multiplication");

  SM_free(vec);
  SM_free(tridiag);
  SM_free(res);
  return pass;
}

bool test_equality_self(void) {
  const SMatF A = gen_random(N_VALS, SIZE, SIZE);

  bool pass = true;
  if (SM_eq(A, A))
    TEST_PASS("Random matrix self equality");
  else {
    TEST_FAIL("Random matrix self equality");
    pass = false;
  }

  SM_free(A);
  return pass;
}

bool test_structure_equality_self(void) {
  const SMatF A = gen_random(N_VALS, SIZE, SIZE);

  bool pass = true;
  if (SM_structure_eq(A, A))
    TEST_PASS("Random structural self equality");
  else {
    TEST_FAIL("Random structureal self equality");
    pass = false;
  }

  SM_free(A);
  return pass;
}

bool test_add_sub_random(void) {
  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  SMatF B = gen_random(N_VALS, SIZE, SIZE);

  SMatF target = SM_addsub_prepare(A, B);
  SM_add(A, B, target);

  bool pass = true;
  for (long row = 0; row < SIZE && pass; ++row) {
    for (long col = 0; col < SIZE && pass; ++col) {
      float t = SM_at(target, row, col);
      float a = SM_at(A, row, col);
      float b = SM_at(B, row, col);

      if (t != a + b) {
        TEST_FAIL_MSG("Matrix Addition", "Failed at (%ld, %ld) (%f + %f != %f)",
                      row, col, a, b, t);
        pass = false;
      }
    }
  }

  if (pass)
    TEST_PASS("Matrix Addition");

  SM_sub(A, B, target);

  for (long row = 0; row < SIZE && pass; ++row) {
    for (long col = 0; col < SIZE && pass; ++col) {
      float t = SM_at(target, row, col);
      float a = SM_at(A, row, col);
      float b = SM_at(B, row, col);

      if (t != a - b) {
        TEST_FAIL_MSG("Matrix Subtraction",
                      "Failed at (%ld, %ld) (%f - %f != %f)", row, col, a, b,
                      t);
        pass = false;
      }
    }
  }

  if (pass)
    TEST_PASS("Matrix Subtraction");

  float scale_fac = rand_float() * 10.0f;
  SMatF target2 = SM_empty_like(A);
  SM_scl(A, scale_fac, target2);

  for (long row = 0; row < SIZE && pass; ++row) {
    for (long col = 0; col < SIZE && pass; ++col) {
      float t = SM_at(target2, row, col);
      float a = SM_at(A, row, col);

      if (t != a * scale_fac) {
        TEST_FAIL_MSG("Matrix Scaling", "Failed at (%ld, %ld) (%f * %f != %f)",
                      row, col, a, scale_fac, t);
        pass = false;
      }
    }
  }
  if (pass)
    TEST_PASS("Matrix Scaling");

  SM_free(A);
  SM_free(B);
  SM_free(target);
  SM_free(target2);
  return pass;
}

bool test_transpose(void) {
  float *values = malloc(N_VALS * sizeof(float));
  long *row_pos = malloc(N_VALS * sizeof(long));
  long *col_pos = malloc(N_VALS * sizeof(long));
  long a_rows = SIZE;
  long a_cols = SIZE + 5;

  gen_randoms(row_pos, col_pos, values, N_VALS, a_rows, a_cols);
  SMatF A = SM_from_pos_with(a_rows, a_cols, N_VALS, row_pos, col_pos, values);

  SMatF T = SM_transpose(A);

  bool pass = true;
  for (long i = 0; i < A.nrows && pass; ++i) {
    for (long j = 0; j < A.ncols && pass; ++j) {
      const float a = SM_at(A, i, j);
      const float t = SM_at(T, j, i);
      if (a != t) {
        pass = false;
        TEST_FAIL_MSG(
            "Matrix Transposition",
            "Value at (%ld, %ld) in A does not equal value at (%ld, %ld) "
            "in T (%.2f != %.2f)",
            i, j, j, i, a, t);
      }
    }
  }

  if (pass)
    TEST_PASS_MSG("Matrix Transposition", "Value Equality passed%s", "");

  SMatF TT = SM_transpose(T);

  if (SM_eq(A, TT)) {
    TEST_PASS("Matrix Transposition");
  } else {
    TEST_FAIL_MSG("Matrix Transposition",
                  "Double does not equal original Matrix %s", "");
    pass = false;
  }

  SM_free(A);
  SM_free(T);
  SM_free(TT);
  free(values);
  free(row_pos);
  free(col_pos);
  return pass;
}

bool test_scalar_prod_self(void) {
  SMatF vec = SM_vec_empty(SIZE);
  for (long i = 0; i < vec.nvals; ++i) {
    vec.vals[i] = rand_float();
  }

  float abs2 = SM_scalar(vec, vec);

  float tst = 0;
  for (long i = 0; i < vec.nvals; ++i) {
    tst += vec.vals[i] * vec.vals[i];
  }

  bool pass = true;
  if (tst != abs2) {
    TEST_FAIL("Scalar Product");
    pass = false;
  } else
    TEST_PASS("Scalar Product");

  SM_free(vec);
  return pass;
}

bool test_energy_norm_identity(void) {
  SMatF vec = SM_vec_empty(SIZE);
  for (long i = 0; i < vec.nvals; ++i) {
    vec.vals[i] = rand_float();
  }

  float tst = SM_scalar(vec, vec);

  SMatF identity = SM_diag_regular((long[1]){0}, (float[1]){1.0f}, 1, SIZE);
  float ret = SM_energy_norm(identity, vec, NULL);
  bool pass = true;
  if (almost_eq(ret, tst)) {
    TEST_PASS("Energy Norm test with identity");
  } else {
    TEST_FAIL("Energy Norm test with identity");
    pass = false;
  }
  SM_free(identity);
  SM_free(vec);
  return pass;
}

bool test_is_triup(void) {
  SMatF is_triup = SM_diag_regular((long[3]){0, 1, 2},
                                   (float[3]){1.0f, -3.0f, 2.0f}, 3, SIZE);

  SMatF no_triup = SM_diag_regular((long[3]){-1, 0, 2},
                                   (float[3]){1.0f, -3.0f, 2.0f}, 3, SIZE);

  bool pass = true;
  if (!SM_is_triup(is_triup, true)) {
    TEST_FAIL("Triangular upper matrix test (shouldn't be triup)");
    SM_print_nonzero(is_triup);
    pass = false;
  } else {
    TEST_PASS("Triangular upper matrix test (shouldn't bee triup)");
  }

  if (SM_is_triup(no_triup, true)) {
    TEST_FAIL("Triangular upper matrix test (should be triup)");
    SM_print_nonzero(no_triup);
    pass = false;
  } else {
    TEST_PASS("Triangular upper matrix test (should be triup)");
  }

  SM_free(is_triup);
  SM_free(no_triup);
  return pass;
}

bool test_back_sub(void) {
  SMatF b = SM_vec_empty(SIZE);
  for (long i = 0; i < b.nvals; ++i) {
    b.vals[i] = rand_float();
  }
  SMatF target = SM_empty_like(b);
  SMatF A = SM_diag_regular((long[3]){0, 1, 2}, (float[3]){1.0f, -3.0f, 2.0f},
                            3, SIZE);

  SM_back_sub(A, b, target);

  SMatF b1 = SM_prod_prepare(A, target);
  SM_prod(A, target, b1);
  SMatF err = SM_addsub_prepare(b, b1);
  SM_sub(b, b1, err);

  float err_abs = SM_abs(err);
  bool pass = true;
  if (err_abs > 0.0001f) {
    TEST_FAIL_MSG("Backward substitution", "Error too large: %f", err_abs);
    pass = false;
  } else {
    TEST_PASS("Backward substitution");
  }

  SM_free(A);
  SM_free(b);
  SM_free(b1);
  SM_free(err);
  SM_free(target);
  return pass;
}

bool test_forw_sub(void) {
  SMatF b = SM_vec_empty(SIZE);
  for (long i = 0; i < b.nvals; ++i) {
    b.vals[i] = rand_float();
  }
  SMatF target = SM_empty_like(b);
  SMatF A = SM_diag_regular((long[3]){-2, -1, 0}, (float[3]){1.0f, -3.0f, 2.0f},
                            3, SIZE);

  SM_forw_sub(A, b, target);

  SMatF b1 = SM_prod_prepare(A, target);
  SM_prod(A, target, b1);
  SMatF err = SM_addsub_prepare(b, b1);
  SM_sub(b, b1, err);

  float err_abs = SM_abs(err);
  bool pass = true;
  if (err_abs > 0.0001f) {
    TEST_FAIL_MSG("Forward substitution", "Error too large: %f", err_abs);
    pass = false;
  } else {
    TEST_PASS("Forward substitution");
  }

  SM_free(A);
  SM_free(b);
  SM_free(b1);
  SM_free(err);
  SM_free(target);
  return pass;
}

bool test_subset_diag(void) {
  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  SMatF diag = SM_subset_diag(A);

  bool pass = true;
  for (long i = 0; i < A.ncols && pass; ++i) {
    float correct = SM_at(A, i, i);
    float found = SM_at(diag, i, i);
    if (correct != found) {
      pass = false;
      TEST_FAIL_MSG("Subset diagonal - diagonal values",
                    "Unexpected value on diagonal (%ld: %f), expected %f", i,
                    found, correct);
    }

    if (SM_has_loc(A, i, i) && !SM_has_loc(diag, i, i)) {
      pass = false;
      TEST_FAIL_MSG("Subset diagonal - diagonal values",
                    "Missing value at (%ld, %ld)", i, i);
    }

    if (!SM_has_loc(A, i, i) && SM_has_loc(diag, i, i)) {
      pass = false;
      TEST_FAIL_MSG("Subset diagonal - diagonal values",
                    "Additional value at (%ld, %ld)", i, i);
    }
  }

  if (pass) {
    TEST_PASS("Subset diagonal - diagonal values");
  }

  for (long row = 0; row < A.nrows && pass; ++row) {
    for (long col = 0; col < A.ncols && pass; ++col) {
      if (row == col)
        continue;

      if (SM_has_loc(diag, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset diagonal - off-diagonal values",
                      "Additional value at (%ld, %ld)", row, col);
      }
    }
  }
  if (pass)
    TEST_PASS("Subset diagonal - off-diagonal values");

  SM_free(A);
  SM_free(diag);
  return pass;
}

bool test_subset_trilo(void) {
  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  bool use_diag = true;
  SMatF trilo = SM_subset_trilo(A, use_diag);

  bool pass = true;
  for (long row = 0; row < trilo.nrows && pass; ++row) {
    for (long col = 0; col < trilo.ncols && pass; ++col) {
      if (use_diag && col > row)
        continue;
      if (!use_diag && col >= row)
        continue;

      float correct = SM_at(A, row, col);
      float found = SM_at(trilo, row, col);
      if (correct != found) {
        pass = false;
        TEST_FAIL_MSG("Subset trilo - triangular lower values",
                      "Unexpected value on at (%ld, %ld: %f), expected %f", row,
                      col, found, correct);
      }

      if (SM_has_loc(A, row, col) && !SM_has_loc(trilo, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset trilo - triangular lower values",
                      "Missing value at (%ld, %ld)", row, col);
      }

      if (!SM_has_loc(A, row, col) && SM_has_loc(trilo, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset trilo - triangular lower values",
                      "Additional value at (%ld, %ld)", row, col);
      }
    }
  }

  if (pass) {
    TEST_PASS("Subset trilo - triangular lower values");
  }

  for (long row = 0; row < A.nrows && pass; ++row) {
    for (long col = 0; col < A.ncols && pass; ++col) {
      if (use_diag && col <= row)
        continue;
      if (!use_diag && col < row)
        continue;

      if (SM_has_loc(trilo, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset trilo - upper values",
                      "Additional value at (%ld, %ld)", row, col);
      }
    }
  }
  if (pass)
    TEST_PASS("Subset trilo - upper values");

  SM_free(A);
  SM_free(trilo);
  return pass;
}

bool test_subset_triup(void) {
  SMatF A = gen_random(N_VALS, SIZE, SIZE);
  bool use_diag = true;
  SMatF triup = SM_subset_triup(A, use_diag);

  bool pass = true;
  for (long row = 0; row < triup.nrows && pass; ++row) {
    for (long col = 0; col < triup.ncols && pass; ++col) {
      if (use_diag && col < row)
        continue;
      if (!use_diag && col <= row)
        continue;

      float correct = SM_at(A, row, col);
      float found = SM_at(triup, row, col);
      if (correct != found) {
        pass = false;
        TEST_FAIL_MSG("Subset triup - triangular upper values",
                      "Unexpected value on at (%ld, %ld: %f), expected %f", row,
                      col, found, correct);
      }

      if (SM_has_loc(A, row, col) && !SM_has_loc(triup, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset triup - triangular upper values",
                      "Missing value at (%ld, %ld)", row, col);
      }

      if (!SM_has_loc(A, row, col) && SM_has_loc(triup, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset triup - triangular upper values",
                      "Additional value at (%ld, %ld)", row, col);
      }
    }
  }

  if (pass) {
    TEST_PASS("Subset triup - triangular upper values");
  }

  for (long row = 0; row < A.nrows && pass; ++row) {
    for (long col = 0; col < A.ncols && pass; ++col) {
      if (use_diag && col >= row)
        continue;
      if (!use_diag && col > row)
        continue;

      if (SM_has_loc(triup, row, col)) {
        pass = false;
        TEST_FAIL_MSG("Subset triup - triangular lower values",
                      "Additional value at (%ld, %ld)", row, col);
      }
    }
  }
  if (pass)
    TEST_PASS("Subset triup - triangular lower values");

  SM_free(A);
  SM_free(triup);
  return pass;
}

int main(void) {
  srand(69);
  test_func funcs[] = {
    test_init_from_pos_with,
    test_structure_equality_self,
    test_equality_self,
    test_random_pos_prod,
    test_random_pos_dif_sizes_prod_ver_matf,
    test_tridiag_vec_prod,
    test_add_sub_random,
    test_transpose,
    test_scalar_prod_self,
    test_energy_norm_identity,
    test_is_triup,
    test_back_sub,
    test_forw_sub,
    test_subset_diag,
    test_subset_trilo,
    test_subset_triup,
  };

  run_tests(funcs);

  return EXIT_SUCCESS;
}
