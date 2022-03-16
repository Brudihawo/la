#include "la.h"

#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "matf.h"

// float between 0 and 1
float rand_float() { return (float)rand() / (float)RAND_MAX; }

#define TEST_SIZE 10
#define RED "\033[31;1m"
#define GRN "\033[32;1m"
#define RST "\033[0m"

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

void pass_fail(MatF c, MatF b2, const char* name) {

  if (MF_eq(c, b2)) {
    printf(GRN"Passed:"RST" %s\n", name);
  } else {
    printf(RED"Failed:"RST" %s\n", name);
    printf("Expected A * u:\n");
    MF_print(b2);
    printf("And got:\n");
    MF_print(c);

    MatF tmp = MF_EMPTY_LIKE(c);

    MF_sub(b2, c, tmp);
    printf("Relative Error: %f (OK: %f)\n", VEC_abs(tmp) / VEC_abs(b2), REL_ERROR);

    MF_FREE(tmp);
  }

}

void test_gaussian_normal() {
  MatF A , b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);

  MF_gauss_elim(A, b, u); // A, b is changed
  MF_prod(A2, u, c);

  pass_fail(c, b2, "gaussian normal");

  cleanup_test(&A, &A2, &b, &b2, &u, &c);
}

void test_gaussian_perm() {
  MatF A , b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);
  // modify A, A2
  *MF_PTR(A, 0, 1) = 0.0f;
  *MF_PTR(A2, 0, 1) = 0.0f;

  MF_gauss_elim(A, b, u); // A, b is changed
  MF_prod(A2, u, c);

  pass_fail(c, b2, "gaussian permute");

  cleanup_test(&A, &A2, &b, &b2, &u, &c);
}

void test_lu_normal() {
  MatF A , b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);
  
  RowPerm *perms = MF_lu_decomp(A); // A is changed
  MF_lu_solv(A, b, perms, u);
  MF_prod(A2, u, c);

  pass_fail(c, b2, "lu normal");
}

void test_lu_perm() {
  MatF A , b, u, c, A2, b2;
  create_test_mfs(&A, &A2, &b, &b2, &u, &c);
  // modify A, A2
  *MF_PTR(A, 0, 0) = 0.0f;
  *MF_PTR(A2, 0, 0) = 0.0f;
  
  RowPerm *perms = MF_lu_decomp(A); // A is changed
  MF_lu_solv(A, b, perms, u);
  MF_prod(A2, u, c);

  pass_fail(c, b2, "lu normal");
}

int main(void) {
  srand(69420);

  // TODO: Write proper tests
  test_gaussian_normal();
  test_gaussian_perm();
  test_lu_normal();
  test_lu_perm();
}
