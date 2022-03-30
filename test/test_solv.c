#include "la.h"

#include "matf.h"
#include "stdbool.h"
#include "stdio.h"

#include "test_util.h"

#define TEST_SIZE 10

void create_test_mfs(MatF *A, MatF *A2, MatF *b, MatF *b2, MatF *u, MatF *c) {
  float *vals = malloc(TEST_SIZE * TEST_SIZE * sizeof(float));
  float *b_vals = malloc(TEST_SIZE * sizeof(float));
  for (long i = 0; i < TEST_SIZE * TEST_SIZE; i++) {
    vals[i] = rand_float() * 10.0f;
  }

  for (long i = 0; i < TEST_SIZE; i++) {
    b_vals[i] = rand_float();
  }

  *A = MF_from(TEST_SIZE, TEST_SIZE, vals);
  *b = MF_from(TEST_SIZE, 1, b_vals);

  *u = MF_EMPTY_LIKE(*b);
  *c = MF_EMPTY_LIKE(*b); // correction value

  *A2 = MF_clone(*A);
  *b2 = MF_clone(*b);
}

void cleanup_test(MatF *A, MatF *A2, MatF *b, MatF *b2, MatF *u, MatF *c) {
  MF_FREE(*A);
  MF_FREE(*A2);
  MF_FREE(*b);
  MF_FREE(*b2);
  MF_FREE(*u);
  MF_FREE(*c);
}

#define pass_fail(c, b2, name, pass)                                           \
  if (MF_eq(c, b2)) {                                                          \
    TEST_PASS(name);                                                           \
  } else {                                                                     \
    MatF tmp = MF_EMPTY_LIKE(c);                                               \
    MF_sub(b2, c, tmp);                                                        \
                                                                               \
    TEST_FAIL_MSG(name, "Relative Error: %f (OK: %f)\n",                       \
                  VEC_abs(tmp) / VEC_abs(b2), REL_ERROR);                      \
    pass = false;                                                              \
                                                                               \
    printf("Expected A * u:\n");                                               \
    MF_print(b2);                                                              \
    printf("And got:\n");                                                      \
    MF_print(c);                                                               \
                                                                               \
    MF_FREE(tmp);                                                              \
  }

bool test_gaussian_normal(void) {
  MatF A, b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);

  MF_gauss_elim(A, b, u); // A, b is changed
  MF_prod(A2, u, c);

  bool pass = true;
  pass_fail(c, b2, "gaussian normal", pass);

  cleanup_test(&A, &A2, &b, &b2, &u, &c);
  return pass;
}

bool test_gaussian_perm(void) {
  MatF A, b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);
  // modify A, A2
  *MF_PTR(A, 0, 1) = 0.0f;
  *MF_PTR(A2, 0, 1) = 0.0f;

  MF_gauss_elim(A, b, u); // A, b is changed
  MF_prod(A2, u, c);

  bool pass = true;
  pass_fail(c, b2, "gaussian permute", pass);

  cleanup_test(&A, &A2, &b, &b2, &u, &c);
  return pass;
}

bool test_lu_normal(void) {
  MatF A, b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);

  RowPerm *perms = MF_lu_decomp(A); // A is changed
  MF_lu_solv(A, b, perms, u);
  MF_prod(A2, u, c);

  bool pass = true;
  pass_fail(c, b2, "lu normal", pass);

  if (perms != NULL)
    RP_free(perms);
  cleanup_test(&A, &A2, &b, &b2, &u, &c);

  return pass;
}

bool test_lu_perm(void) {
  MatF A, b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);
  // modify A, A2
  *MF_PTR(A, 0, 0) = 0.0f;
  *MF_PTR(A2, 0, 0) = 0.0f;

  RowPerm *perms = MF_lu_decomp(A); // A is changed
  MF_lu_solv(A, b, perms, u);
  MF_prod(A2, u, c);

  bool pass = true;
  pass_fail(c, b2, "lu normal", pass);

  if (perms != NULL) {
    RP_free(perms);
    free(perms);
  }
  cleanup_test(&A, &A2, &b, &b2, &u, &c);
  return pass;
}

int main(void) {
  srand(69420);

  // TODO: Write proper tests for MatF functions
  test_func test_funcs[] = {
    test_gaussian_normal,
    test_gaussian_perm,
    test_lu_normal,
    test_lu_perm,
  };

  run_tests(test_funcs);
}
