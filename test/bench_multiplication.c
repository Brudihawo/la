#include "matf.h"
#include "math.h"
#include "sparse.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

float rand_float() { return (float)rand() / (float)RAND_MAX; }

int comp_long(const void *a, const void *b) {
  const long va = *(long *)a;
  const long vb = *(long *)b;
  return va > vb;
}

void gen_randoms(long *rows, long *cols, float *vals, long n_vals, long size) {
  long *idcs = malloc(n_vals * sizeof(long));

  for (long i = 0; i < n_vals; ++i) {
    bool repeat = true;
    while (repeat) {
      vals[i] = rand_float();
      idcs[i] = rand() % (size * size);

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
    cols[i] = idcs[i] % size;
    rows[i] = idcs[i] / size;
  }
  free(idcs);
}

void matrix_product_time(long size, long nvals) {
  float *values_A = malloc(nvals * sizeof(float));
  long *row_pos_A = malloc(nvals * sizeof(long));
  long *col_pos_A = malloc(nvals * sizeof(long));

  float *values_B = malloc(nvals * sizeof(float));
  long *row_pos_B = malloc(nvals * sizeof(long));
  long *col_pos_B = malloc(nvals * sizeof(long));

  gen_randoms(row_pos_A, col_pos_A, values_A, nvals, size);
  gen_randoms(row_pos_B, col_pos_B, values_B, nvals, size);

  MatF m_A = MF_with(size, size, 0.0f);
  MatF m_B = MF_with(size, size, 0.0f);

  for (long i = 0; i < nvals; ++i) {
    *MF_PTR(m_A, row_pos_A[i], col_pos_A[i]) = values_A[i];
    *MF_PTR(m_B, row_pos_B[i], col_pos_B[i]) = values_B[i];
  }

  clock_t m_start = clock();

  MatF m_T = MF_EMPTY_LIKE(m_A);
  MF_prod(m_A, m_B, m_T);

  clock_t m_end = clock();

  double m_time = (double)(m_end - m_start) / CLOCKS_PER_SEC;

  SMatF s_A =
      SM_from_pos_with(size, size, nvals, row_pos_A, row_pos_A, values_A);
  SMatF s_B =
      SM_from_pos_with(size, size, nvals, row_pos_B, row_pos_B, values_B);

  clock_t s_start = clock();

  SMatF s_T = SM_prod_prepare(s_A, s_B);
  SM_prod(s_A, s_B, s_T);

  clock_t s_end = clock();
  double s_time = (double)(s_end - s_start) / CLOCKS_PER_SEC;

  for (long row = 0; row < m_T.rows; ++row) {
    for (long col = 0; col < m_T.cols; ++col) {
      if (MF_AT(m_T, row, col) != SM_at(s_T, row, col)) {
        printf("Incorrect product\n");
        break;
      }
    }
  }

  fprintf(stderr, "%5ld %13.8f %13.8f %13.8f\n", size, m_time, s_time, s_time / m_time);

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

int main(void) {
  srand(69);
  fprintf(stderr, "# SIZE   MatF/s        SMatF/s       SMatF / MatF\n");
  for (int order = 4; order < 13; ++order) {
    const long size = (long)pow(2, order);
    const long nvals = size * 5;

    matrix_product_time(size, nvals);
  }
}
