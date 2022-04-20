#include "float.h"
#include "la.h"
#include "limits.h"
#include "math.h"
#include "sparse.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "stdbool.h"
#include "assert.h"

#define FPRC "%11.6f"


void gen_randoms(long *rows, long *cols, float *vals, long n_vals, long n_rows,
                 long n_cols) {
  long *idcs = malloc(n_vals * sizeof(long));

  for (long i = 0; i < n_vals; ++i) {
    vals[i] = rand_float();
    idcs[i] = rand() % (n_rows * n_cols);
  }
  radix_sort(idcs, n_vals, NULL);

  long invalid_indices = 0;
  for (long i = 1; i < n_vals; ++i) {
    if (idcs[i] == idcs[i - 1]) {
      idcs[i - 1] = LONG_MAX;
      ++invalid_indices;
    }
  }

  qsort(idcs, n_vals, sizeof(long), long_gt);

  long *tmp_arr = malloc(invalid_indices * sizeof(long));

  for (long i = n_vals - invalid_indices; i < n_vals; ++i) {
    const long cur_corr_idx = i + invalid_indices - n_vals;

    if (cur_corr_idx > 0)
      radix_sort(&idcs[n_vals - invalid_indices], cur_corr_idx, tmp_arr);

    bool repeat = true;
    while (repeat) {
      vals[i] = rand_float();
      idcs[i] = rand() % (n_rows * n_cols);

      repeat = false;

      long pos =
          binary_search(&idcs[n_vals - invalid_indices - 1], cur_corr_idx, idcs[i]);
      if (pos > 0) {
        repeat = true;
      }

      if (binary_search(idcs, n_vals - invalid_indices, idcs[i]) > 0) {
        repeat = true;
      }
    }
  }
  free(tmp_arr);

  radix_sort(idcs, n_vals, NULL);

  for (long i = 1; i < n_vals; ++i) {
    assert(idcs[i] != idcs[i - 1]);
  }

  for (long i = 0; i < n_vals; ++i) {
    cols[i] = idcs[i] % n_cols;
    rows[i] = idcs[i] / n_cols;
  }

  free(idcs);
}


void matrix_product_time(long size, long nvals, long n_runs) {
  float *values_A = malloc(nvals * sizeof(float));
  long *row_pos_A = malloc(nvals * sizeof(long));
  long *col_pos_A = malloc(nvals * sizeof(long));

  float *values_B = malloc(nvals * sizeof(float));
  long *row_pos_B = malloc(nvals * sizeof(long));
  long *col_pos_B = malloc(nvals * sizeof(long));

  gen_randoms(row_pos_A, col_pos_A, values_A, nvals, size, size);
  gen_randoms(row_pos_B, col_pos_B, values_B, nvals, size, size);

  SMatF s_A =
      SM_from_pos_with(size, size, nvals, row_pos_A, col_pos_A, values_A);
  SMatF s_B =
      SM_from_pos_with(size, size, nvals, row_pos_B, col_pos_B, values_B);

  double *times = malloc(n_runs * sizeof(double));

  for (long n = 0; n < n_runs; ++n) {
    clock_t s_start = clock();

    SMatF s_T = SM_prod_prepare(s_A, s_B);
    SM_prod(s_A, s_B, s_T);

    clock_t s_end = clock();
    times[n] = (double)(s_end - s_start) / CLOCKS_PER_SEC;

    SM_free(s_T);
  }

  double mean = 0.0f;
  double min = FLT_MAX;
  double max = 0.0f;
  for (long n = 0; n < n_runs; ++n) {
    mean += times[n];
    if (times[n] < min) {
      min = times[n];
    }

    if (times[n] > max) {
      max = times[n];
    }
  }

  mean /= n_runs;
  fprintf(stderr, "%6ld " FPRC " " FPRC " " FPRC "\n", size, min, max, mean);

  SM_free(s_A);
  SM_free(s_B);

  free(values_A);
  free(values_B);
  free(row_pos_A);
  free(row_pos_B);
  free(col_pos_A);
  free(col_pos_B);
}

int main(void) {
  srand(69);
  const long n_runs = 5;
  fprintf(stderr, "# SIZE %11s %11s %11s [%ld iterations per size]\n", "Min",
          "Max", "Mean", n_runs);
  for (int order = 6; order < 22; ++order) {
    const long size = (long)pow(2, order);
    const long nvals = size * 5;

    matrix_product_time(size, nvals, n_runs);
  }
}
