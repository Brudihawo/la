#include "sparse.h"

#include "assert.h"
#include "log.h"
#include "memory.h"
#include "stdio.h"
#include "stdlib.h"

#include "la.h"
#include <math.h>

#define check_size(rows, cols, n_vals)                                         \
  if (n_vals > rows * cols) {                                                  \
    log_err("%s:%i Number of values in SMatF (nvals = %ld) cannot exceed "     \
            "rows * columns (%ld * %ld = %ld).",                               \
            __FILE__, __LINE__, n_vals, rows, cols, rows *cols);               \
    exit(EXIT_FAILURE);                                                        \
  }

SMatF SM_empty(long rows, long cols, long n_vals) {
  check_size(rows, cols, n_vals);
  return (SMatF){
      .nrows = rows,
      .ncols = cols,
      .nvals = n_vals,

      .col_sizes = calloc(cols, sizeof(long)),
      .col_starts = malloc(cols * sizeof(long)),
      .col_idcs = malloc(n_vals * sizeof(long)),
      .col_pos = malloc(n_vals * sizeof(long)),
      .vals = malloc(n_vals * sizeof(float)),
      .row_starts = malloc(rows * sizeof(long)),
      .row_sizes = calloc(rows, sizeof(long)),
  };
}

SMatF SM_vec_empty(long rows) {
  SMatF ret = SM_empty(rows, 1, rows);

  ret.col_sizes[0] = rows;
  ret.col_starts[0] = 0;
  for (long i = 0; i < rows; ++i) {
    ret.col_idcs[i] = i;
    ret.col_pos[i] = 0;

    ret.row_sizes[i] = 1;
    ret.row_starts[i] = i;
  }
  return ret;
}

float SM_abs(SMatF A) {
  float acc = 0.0f;
  for (long i = 0; i < A.nvals; ++i) {
    acc += A.vals[i] * A.vals[i];
  }
  return sqrtf(acc);
}

SMatF SM_empty_like(SMatF A) {
  SMatF ret = {
      .nrows = A.nrows,
      .ncols = A.ncols,
      .nvals = A.nvals,

      .col_sizes = malloc(A.ncols * sizeof(long)),
      .col_starts = malloc(A.ncols * sizeof(long)),
      .col_idcs = malloc(A.nvals * sizeof(long)),
      .col_pos = malloc(A.nvals * sizeof(long)),
      .vals = malloc(A.nvals * sizeof(float)),
      .row_starts = malloc(A.nrows * sizeof(long)),
      .row_sizes = malloc(A.nrows * sizeof(long)),
  };

  // copy non-zero structure of A
  memcpy(ret.col_sizes, A.col_sizes, A.ncols * sizeof(long));
  memcpy(ret.col_starts, A.col_starts, A.ncols * sizeof(long));
  memcpy(ret.col_idcs, A.col_idcs, A.nvals * sizeof(long));
  memcpy(ret.col_pos, A.col_pos, A.nvals * sizeof(long));
  memcpy(ret.row_starts, A.row_starts, A.nrows * sizeof(long));
  memcpy(ret.row_sizes, A.row_sizes, A.nrows * sizeof(long));

  return ret;
}

SMatF SM_clone(SMatF A) {
  SMatF ret = SM_empty_like(A);

  // copy values from A
  memcpy(ret.vals, A.vals, A.nvals * sizeof(float));

  return ret;
}

/* @brief initialize row / column start arrays after row / column sizes have
 *        been initialized. For internal use only. This is not in the header
 */
void SM_init_start_arrs(SMatF A) {
  A.row_starts[0] = 0;
  for (long row = 1; row < A.nrows; ++row) {
    A.row_starts[row] = A.row_starts[row - 1] + A.row_sizes[row - 1];
  }

  A.col_starts[0] = 0;
  for (long col = 1; col < A.ncols; ++col) {
    A.col_starts[col] = A.col_starts[col - 1] + A.col_sizes[col - 1];
  }
}

/* @brief initialize col_idcs after col_sizes, row_sizes, row_starts, col_starts have
 * been initialized
 */
void SM_init_col_idcs(SMatF A) {
  long cur_val_idx = 0;
  for (long col = 0; col < A.ncols; ++col) { // iterate over columns first
    for (long row = 0; row < A.nrows; ++row) { // then rows
      if (SM_has_loc(A, row, col)) {
        A.col_idcs[cur_val_idx] = SM_idx(A, row, col);
        ++cur_val_idx;
      }
    }
  }
}

SMatF SM_empty_from_pos(long n_rows, long n_cols, long n_vals, long *row_pos,
                        long *col_pos) {
  check_size(n_rows, n_cols, n_vals);
  SMatF ret = {
      .nrows = n_rows,
      .ncols = n_cols,
      .nvals = n_vals,

      .col_sizes = calloc(n_cols, sizeof(long)),
      .col_starts = malloc(n_cols * sizeof(long)),
      .col_idcs = malloc(n_vals * sizeof(long)),
      .col_pos = malloc(n_vals * sizeof(long)),
      .vals = malloc(n_vals * sizeof(float)),
      .row_starts = malloc(n_rows * sizeof(long)),
      .row_sizes = calloc(n_rows, sizeof(long)),
  };

  memcpy(ret.col_pos, col_pos, n_vals * sizeof(long));

  for (long i = 0; i < n_vals; ++i) {
    assert((row_pos[i] > 0 && row_pos[i] < n_rows) &&
           "Row Position out of bounds");
    assert((col_pos[i] > 0 && col_pos[i] < n_cols) &&
           "Column Position out of bounds");
    // increment corresponding row and column sizes
    ++ret.col_sizes[col_pos[i]];
    ++ret.row_sizes[row_pos[i]];
  }

  SM_init_start_arrs(ret);

  // validate nonzero position order
  long last_idx = -1;
  for (long i = 0; i < n_vals; ++i) {
    long cur_idx = SM_idx(ret, row_pos[i], col_pos[i]);
    assert(last_idx < cur_idx &&
           "Positions need to be passed in row-major order");
    last_idx = cur_idx;
  }

  ret.col_pos = col_pos;
  for (long i = 0; i < n_vals; ++i) {
    ret.col_idcs[i] = SM_idx(ret, row_pos[i], col_pos[i]);
  }

  return ret;
}

SMatF SM_from_pos_with(long n_rows, long n_cols, long n_vals, long *row_pos,
                       long *col_pos, float *vals) {
  check_size(n_rows, n_cols, n_vals);

  SMatF ret = {
      .nrows = n_rows,
      .ncols = n_cols,
      .nvals = n_vals,

      .col_sizes = calloc(n_cols, sizeof(long)),
      .col_starts = malloc(n_cols * sizeof(long)),
      .col_idcs = malloc(n_vals * sizeof(long)),
      .col_pos = malloc(n_vals * sizeof(long)),
      .vals = malloc(n_vals * sizeof(float)),
      .row_starts = malloc(n_rows * sizeof(long)),
      .row_sizes = calloc(n_rows, sizeof(long)),
  };

  memcpy(ret.vals, vals, n_vals * sizeof(float));
  memcpy(ret.col_pos, col_pos, n_vals * sizeof(long));

  for (long i = 0; i < n_vals; ++i) {
    if ((row_pos[i] < 0 || row_pos[i] >= n_rows) ||
        (col_pos[i] < 0 || col_pos[i] >= n_cols)) {
      log_err(
          "Position (%ld, %ld) is out of bounds for matrix of size (%ld, %ld)");
      exit(EXIT_FAILURE);
    }

    // increment corresponding row and column sizes
    ++ret.col_sizes[col_pos[i]];
    ++ret.row_sizes[row_pos[i]];
  }

  SM_init_start_arrs(ret);

  // validate nonzero position order
  long last_idx = -1;
  for (long i = 0; i < n_vals; ++i) {
    long cur_idx = SM_idx(ret, row_pos[i], col_pos[i]);
    if (last_idx > cur_idx) {
      assert(i != 0);
      log_err("Positions need to be passed in row-major order "
              "(not satisfied by last and current positions %ld / %ld at (%ld, "
              "%ld) "
              "and (%ld, %ld) with row-major indices %ld and %ld: %ld > %ld)",
              i - 1, i, row_pos[i - 1], col_pos[i - 1], row_pos[i], col_pos[i],
              last_idx, cur_idx, last_idx, cur_idx);
      exit(EXIT_FAILURE);
    }

    last_idx = cur_idx;
  }

  SM_init_col_idcs(ret);

  return ret;
}

SMatF SM_empty_diag(long *diags, long n_diags, long size) {
  long n_vals = 0;
  for (long d = 0; d < n_diags; ++d) {
    const long cur_diag = diags[d];
    if (labs(cur_diag) >= size) {
      log_err("Failed trying to set diagonal %ld in matrix of size %ld. "
              "Cannot set values outside of matrix.",
              cur_diag, size);
      exit(EXIT_FAILURE);
    }

    n_vals += (size - labs(cur_diag));
  }

  SMatF ret = SM_empty(size, size, n_vals);

  // compute row / col sizes
  long count_row = 0;
  for (long rc_idx = 0; rc_idx < size; ++rc_idx) {
    for (long d = 0; d < n_diags; ++d) {
      const long cur_diag = diags[d];
      // check value present in row and diagonal
      // rc_idx represents a row here
      if (((cur_diag < 0) && (rc_idx >= -cur_diag)) || // diagonal below main
          ((cur_diag >= 0) && (rc_idx + cur_diag < size))) { // diagonal above
        // row sizes
        ++ret.row_sizes[rc_idx];

        // set column of current value position
        ret.col_pos[count_row] = rc_idx + cur_diag;
        ++count_row;
      }

      // check value present in column and diagonal
      // rc_idx represents a column here
      if (((cur_diag >= 0) && (rc_idx >= cur_diag)) || // diagonal above main
          ((cur_diag < 0) && (rc_idx - cur_diag < size))) { // diagonal below
        // column sizes
        ++ret.col_sizes[rc_idx];
      }
    }
  }

  // row / column starts
  SM_init_start_arrs(ret);

  SM_init_col_idcs(ret);

  return ret;
}

SMatF SM_diag_regular(long *diags, float *diag_vals, long n_diags, long size) {
  SMatF ret = SM_empty_diag(diags, n_diags, size);

  for (long d = 0; d < n_diags; ++d) {
    const long cur_diag = diags[d];
    long row = cur_diag > 0 ? 0 : -cur_diag;
    long col = cur_diag > 0 ? cur_diag : 0;
    for (long idx = 0; idx < size - labs(cur_diag); ++idx) {
      SM_set_or_panic(ret, row, col, diag_vals[d]);

      ++row;
      ++col;
    }
  }
  return ret;
}

// TODO: check correctness
bool SM_has_loc(SMatF A, long row, long col) {
  // Check if row is present in matrix
  if (SM_ROW_EMPTY(A, row))
    return false;

  const long row_start = A.row_starts[row];
  for (long i = 0; i < A.row_sizes[row]; i++) {
    if (A.col_pos[row_start + i] == col)
      return true;
  }
  return false;
}

long SM_idx(SMatF A, long row, long col) {
  if ((row >= A.nrows || col >= A.ncols) || (row < 0 || col < 0)) {
    log_err("Position (%ld, %ld) out of bounds in Matrix of size (%ld, %ld).",
            row, col, A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  if (SM_ROW_EMPTY(A, row))
    return SM_NOT_PRESENT;

  const long row_start = A.row_starts[row];
  for (long i = 0; i < A.row_sizes[row]; i++) {
    if (A.col_pos[row_start + i] == col)
      return row_start + i;
  }
  return SM_NOT_PRESENT;
}

bool SM_structure_eq(SMatF A, SMatF B) {
  // TODO: I think this can be optimized. I dont think we have to do all these
  // checks.
  if (A.ncols != B.ncols || A.nrows != B.nrows || A.nvals != B.nvals)
    return false;

  if (memcmp(A.row_sizes, B.row_sizes, A.nrows * sizeof(long)))
    return false;
  if (memcmp(A.row_starts, B.row_starts, A.nrows * sizeof(long)))
    return false;

  if (memcmp(A.col_sizes, B.col_sizes, A.ncols * sizeof(long)))
    return false;
  if (memcmp(A.col_starts, B.col_starts, A.ncols * sizeof(long)))
    return false;

  if (memcmp(A.col_pos, B.col_pos, A.nvals * sizeof(long)))
    return false;
  if (memcmp(A.col_idcs, B.col_idcs, A.nvals * sizeof(long)))
    return false;

  return true;
}

bool SM_eq(SMatF A, SMatF B) {
  // compare structure / short circuit
  if (!SM_structure_eq(A, B))
    return false;

  // compare values
  if (memcmp(A.vals, B.vals, A.nvals * sizeof(float)) != 0)
    return false;

  // if non-zero structure and values are equal, A and B are equal
  return true;
}

// TODO: do i want a non-panicking setter?
void SM_set_or_panic(SMatF A, long row, long col, float val) {
  long idx = SM_idx(A, row, col);
  assert(idx != SM_NOT_PRESENT &&
         "Can only assign to non-zero type position in SMatF");

  A.vals[idx] = val;
}

float SM_at(SMatF A, long row, long col) {
  assert(row < A.nrows && col < A.ncols && "Position out of bounds");

  long idx = SM_idx(A, row, col);
  return idx == SM_NOT_PRESENT ? 0.0f : A.vals[idx];
}

long SM_col(SMatF A, long row, long col_idx) {
  if (row >= A.nrows || col_idx >= A.ncols) {
    log_err("Row and column index (%ld, %ld) out of bounds "
            "for matrix of size %ld x %ld",
            row, col_idx, A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  if (SM_ROW_EMPTY(A, row) || col_idx >= A.row_sizes[row])
    return SM_NOT_PRESENT;

  return A.col_pos[A.row_starts[row] + col_idx];
}

long SM_col_or_panic(SMatF A, long row, long col_idx) {
  long col = SM_col(A, row, col_idx); // column in target
  assert(col != SM_NOT_PRESENT && col > -1 && "Bug in row->column-iteration");
  return col;
}

float *SM_ptr(SMatF A, long row, long col) {
  assert(row < A.nrows && col < A.ncols && "Position out of bounds");

  long idx = SM_idx(A, row, col);
  return idx == SM_NOT_PRESENT ? NULL : &A.vals[idx];
}

float *SM_ptr_or_panic(SMatF A, long row, long col) {
  assert(row < A.nrows && col < A.ncols && "Position out of bounds");
  long idx = SM_idx(A, row, col);
  assert(idx != SM_NOT_PRESENT &&
         "Can only get pointer to non-zero type element in SMatF");
  return &A.vals[idx];
}

SMatF SM_addsub_prepare(SMatF A, SMatF B) {
  assert((A.nrows == B.nrows && A.ncols == B.ncols) &&
         "Size mismatch. A and B need to have the same size!");

  // short circuit if a and b have the same structure
  if (SM_structure_eq(A, B))
    return SM_empty_like(A);
  // compute number of values in target
  long nvals = 0;
  for (long col = 0; col < B.ncols; ++col) {   // columns in target
    for (long row = 0; row < A.nrows; ++row) { // rows in target
      if (SM_has_loc(A, row, col) || SM_has_loc(B, row, col)) {
        ++nvals;
      }
    }
  }

  SMatF ret = SM_empty(A.nrows, A.ncols, nvals);
  long cur_val_idx = 0;
  for (long row = 0; row < A.nrows; ++row) {   // rows in target
    for (long col = 0; col < B.ncols; ++col) { // columns in target
      if (SM_has_loc(A, row, col) || SM_has_loc(B, row, col)) {
        ret.row_sizes[row] += 1;
        ret.col_sizes[col] += 1;
        ret.col_pos[cur_val_idx] = col;
        ++cur_val_idx;
      }
    }
  }

  SM_init_start_arrs(ret);

  SM_init_col_idcs(ret);

  return ret;
}

void SM_add(SMatF A, SMatF B, SMatF target) {
  assert((A.nrows == B.nrows && A.ncols == B.ncols) &&
         (target.ncols == A.ncols && target.nrows == A.nrows) &&
         "Size mismatch. A, B, and target need to have the same size!");

  for (long row = 0; row < target.nrows; ++row) {
    for (long col = 0; col < target.nrows; ++col) {
      if (SM_has_loc(target, row, col)) {
        // TODO: think about checking querying nonzero in A or B first, but i
        // think this is more efficient right now.
        SM_set_or_panic(target, row, col,
                        SM_at(A, row, col) + SM_at(B, row, col));
      }
    }
  }
}

void SM_sub(SMatF A, SMatF B, SMatF target) {
  assert((A.nrows == B.nrows && A.ncols == B.ncols) &&
         (target.ncols == A.ncols && target.nrows == A.nrows) &&
         "Size mismatch. A, B, and target need to have the same size!");

  for (long row = 0; row < target.nrows; ++row) {
    for (long col = 0; col < target.nrows; ++col) {
      if (SM_has_loc(target, row, col)) {
        // TODO: think about checking querying nonzero in A or B first, but i
        // think this is more efficient right now.
        SM_set_or_panic(target, row, col,
                        SM_at(A, row, col) - SM_at(B, row, col));
      }
    }
  }
}

void SM_scl(SMatF A, float s, SMatF target) {
  assert(SM_structure_eq(A, target) &&
         "A and target must have the same shape!");

  for (long i = 0; i < A.nvals; ++i) {
    target.vals[i] = A.vals[i] * s;
  }
}

SMatF SM_prod_prepare(SMatF A, SMatF B) {
  assert((A.ncols == B.nrows) && "Size mismatch, needs a.cols == b.rows.");

  SMatF ret = {
      .nrows = A.nrows,
      .ncols = B.ncols,
      .nvals = 0, // tbd

      .col_sizes = calloc(B.ncols, sizeof(long)),
      .col_starts = malloc(B.ncols * sizeof(long)),
      .col_idcs = NULL,
      .col_pos = NULL,
      .vals = NULL,
      .row_starts = malloc(A.nrows * sizeof(long)),
      .row_sizes = calloc(A.nrows, sizeof(long)),
  };

  // get memory requirements
  for (long t_row = 0; t_row < ret.nrows; t_row++) {
    for (long t_col = 0; t_col < ret.ncols;
         t_col++) { // iterate on row, column in target
      for (long test_idx = 0; test_idx < A.row_sizes[t_row]; test_idx++) {
        // Test idx in A and B for column / row respectively
        long idx = SM_col(A, t_row, test_idx);
        if (idx == SM_NOT_PRESENT)
          continue;

        if (SM_has_loc(A, t_row, SM_col_or_panic(A, t_row, test_idx)) &&
            SM_has_loc(B, SM_col_or_panic(A, t_row, test_idx), t_col)) {
          ret.nvals++;
          ret.col_sizes[t_col]++;
          ret.row_sizes[t_row]++;

          break;
        }
      }
    }
  }

  // Set row / column starts
  SM_init_start_arrs(ret);

  // allocate memory
  ret.vals = malloc(ret.nvals * sizeof(float));
  ret.col_idcs = malloc(ret.nvals * sizeof(long));
  ret.col_pos = malloc(ret.nvals * sizeof(long));

  // set positions for row-wise iteration, general manipulation
  long cur_idx = 0;
  for (long t_row = 0; t_row < ret.nrows; t_row++) {
    for (long t_col = 0; t_col < ret.ncols;
         t_col++) { // iterate on row, column in target
      for (long test_idx = 0; test_idx < A.row_sizes[t_row]; test_idx++) {
        // Test idx in A and B for column / row respectively
        const long test_pos = SM_col_or_panic(A, t_row, test_idx);
        if (SM_has_loc(A, t_row, test_pos) && SM_has_loc(B, test_pos, t_col)) {
          ret.col_pos[cur_idx] = t_col;
          cur_idx++;
          break;
        }
      }
    }
  }

  // TODO: Do i need column-wise iteration in SMatF?
  // Assignments for column-wise iteration
  SM_init_col_idcs(ret);

  return ret;
}

void SM_prod(SMatF A, SMatF B, SMatF target) {
  // size assertions
  assert(A.ncols == B.nrows && "Size mismatch between A and B");
  assert(target.nrows == A.nrows && "target.nrows == A.nrows");
  assert(target.ncols == B.ncols && "target.ncols == B.ncols");

  for (long t_row = 0; t_row < target.nrows; t_row++) { // rows in target
    for (long t_col_i = 0; t_col_i < target.row_sizes[t_row];
         t_col_i++) { // values in row of target (present columns' indices)
      const long t_col = SM_col_or_panic(target, t_row, t_col_i); // column in target
      SM_set_or_panic(target, t_row, t_col, 0.0f);

      for (long a_col_i = 0; a_col_i < A.row_sizes[t_row]; a_col_i++) {
        const long idx = SM_col_or_panic(A, t_row, a_col_i); // column in a

        if (SM_has_loc(B, idx, t_col)) { // check if value is non-zero in b
          *SM_ptr_or_panic(target, t_row, t_col) +=
              SM_at(A, t_row, idx) * SM_at(B, idx, t_col);
        }
      }
    }
  }
}

SMatF SM_jacobi(SMatF A, SMatF b, float rel_err_max, long n_iter,
                float stale_bound) {
  if (A.nrows != A.ncols) {
    log_err("Can only perform the jacobi-method for matrix inversion on "
            "invertible matrices. Here, A.nrows != A.ncols (%ld != %ld).",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  SMatF u_k = SM_empty_like(b);
  SMatF u_k1 = SM_empty_like(b);

  // random initialisation of target
  for (long i = 0; i < u_k.nvals; ++i) {
    u_k.vals[i] = (float)rand() / (float)RAND_MAX;
  }

  SMatF cur_result = SM_prod_prepare(A, b);
  SMatF err = SM_empty_like(b);
  float last_err = -1;

  for (long iter = 0; iter < n_iter; ++iter) {
    SM_prod(A, u_k, cur_result);
    SM_sub(b, cur_result, err);
    float rel_err = SM_abs(err) / SM_abs(b);

    if (rel_err < rel_err_max || rel_err / last_err > stale_bound) {
      break;
    } else {
      if (iter % 10 == 0)
        log_msg("Iteration %ld: %f %f", iter, rel_err, rel_err / last_err);

      last_err = rel_err;
      for (long i = 0; i < u_k.nvals; ++i) {
        u_k1.vals[i] = 0;
        for (long j = 0; j < u_k.nvals; ++j) {
          if (i == j)
            continue;

          if (SM_has_loc(A, i, j)) {
            u_k1.vals[i] -= SM_at(A, i, j) * u_k.vals[j];
          }
        }
        u_k1.vals[i] += SM_at(b, i, 0);
        u_k1.vals[i] /= SM_at(A, i, i);
      }

      float *tmp = u_k.vals;
      u_k.vals = u_k1.vals;
      u_k1.vals = tmp;
    }
  }

  // TODO: [SM_jacobi] reduce allocations or move outside for repeated calling
  SM_free(cur_result);
  SM_free(err);
  SM_free(u_k1);
  return u_k;
}

SMatF SM_transpose(SMatF A) {
  SMatF ret = SM_empty_like(A);

  long* tmp_sizes;
  tmp_sizes = ret.col_sizes;
  ret.col_sizes = ret.row_sizes;
  ret.row_sizes = tmp_sizes;

  for (long row = 0; row < A.nrows; ++row) {
    for (long col_idx = 0; col_idx < A.col_sizes[row]; ++col_idx) {
      long col = SM_col_or_panic(A, row, col_idx);

      SM_set_or_panic(ret, col, row, SM_at(A, row, col));
    }
  }

  SM_init_start_arrs(ret);

  for (long row = 0; row < ret.nrows; ++row) {
    for (long col_idx = 0; col_idx < ret.col_sizes[row]; ++col_idx) {
      long col = SM_col_or_panic(ret, row, col_idx);
    }
  }

  return ret;
}

void SM_gauss_seidel_or(SMatF A, SMatF b, SMatF target, float or, float rel_err,
                        long n_iter);

// SMatF SM_conjugate_gradients(SMatF A, SMatF b, float rel_err) {
// }

void SM_print(SMatF A) {
  printf("[");
  for (long row = 0; row < A.nrows; row++) {
    if (row != 0) {
      printf(" ");
    }
    printf(" ");
    for (long col = 0; col < A.ncols; col++) {
      printf(LA_PRINT_FMT, SM_at(A, row, col));
    }
    if (row != A.nrows - 1) {
      printf("\n");
    }
  }
  printf(" ]\n");
}

void SM_print_nonzero(SMatF A) {
  for (long row = 0; row < A.nrows; row++) {
    printf("| ");
    for (long col = 0; col < A.ncols; col++) {
      printf(SM_has_loc(A, row, col) ? "* " : "  ");
    }
    printf("|\n");
  }
}

void SM_print_shape(SMatF A) { printf("(%ld x %ld)\n", A.nrows, A.ncols); }

void SM_print_meta(SMatF A) {
  printf("SMatF (%ld x %ld), %ld values set.\n", A.nrows, A.ncols, A.nvals);
}

void SM_free(SMatF A) {
  free(A.vals);
  free(A.col_idcs);
  free(A.row_sizes);
  free(A.col_sizes);
  free(A.row_starts);
  free(A.col_starts);
  free(A.col_pos);
}
