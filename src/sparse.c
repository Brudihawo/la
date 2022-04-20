#include "sparse.h"

#include "assert.h"
#include "limits.h"
#include "log.h"
#include "memory.h"
#include "stdio.h"
#include "stdlib.h"

#include "float.h"
#include "la.h"
#include "math.h"

typedef struct {
  long row, col;
} SMPos;

typedef struct {
  SMPos pos;
  float val;
} SMEntry;

typedef struct {
  SMEntry *vals;
  long size;
  long cap;
} PosBuf; // resizable buffer of SMEntry

typedef struct {
  long *vals;
  long size;
  long cap;
} RSBufL; // resizable buffer of longs

/* @brief Comparison for row major sorting of SMPos
 */
int SMPos_comp_rm(const void *lhs, const void *rhs) {
  const SMPos l = *(SMPos *)lhs;
  const SMPos r = *(SMPos *)rhs;

  if (l.row == r.row) {
    return l.col > r.col;
  } else {
    return l.row > r.row;
  }
}

/* @brief Comparison for row major sorting of SMEntry
 */
int SMEntry_comp_rm(const void *lhs, const void *rhs) {
  const SMEntry *l = (SMEntry *)lhs;
  const SMEntry *r = (SMEntry *)rhs;

  return SMPos_comp_rm(l, r);
}

RSBufL RSBL_with_size(long cap) {
  return (RSBufL){
      .vals = malloc(cap * sizeof(long)),
      .cap = cap,
      .size = 0,
  };
}

PosBuf PB_with_size(long cap) {
  return (PosBuf){
      .vals = malloc(cap * sizeof(SMEntry)),
      .cap = cap,
      .size = 0,
  };
}

#define DYNBUF_CAP_SCL 1.8f
static inline void RSBL_push_val(RSBufL *buf, long val) {
  if (buf->size >= buf->cap) {
    buf->cap *= DYNBUF_CAP_SCL;
    buf->vals = realloc(buf->vals, buf->cap * sizeof(long));
  }

  buf->vals[buf->size] = val;
  ++buf->size;
}

static inline void PB_push_val(PosBuf *buf, long row, long col, float val) {
  if (buf->size >= buf->cap) {
    buf->cap *= DYNBUF_CAP_SCL;
    buf->vals = realloc(buf->vals, buf->cap * sizeof(SMEntry));
  }

  buf->vals[buf->size].pos.row = row;
  buf->vals[buf->size].pos.col = col;
  buf->vals[buf->size].val = val;
  ++buf->size;
}

void RSBL_free(RSBufL buf) { free(buf.vals); }

void PB_free(PosBuf buf) { free(buf.vals); }

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
    ret.col_pos[i] = 0;

    ret.row_sizes[i] = 1;
    ret.row_starts[i] = i;
  }
  return ret;
}

float SM_sse(SMatF A) {
  float acc = 0.0f;
  for (long i = 0; i < A.nvals; ++i) {
    acc += A.vals[i] * A.vals[i];
  }
  return acc;
}

float SM_abs(SMatF A) { return sqrtf(SM_sse(A)); }

SMatF SM_empty_like(SMatF A) {
  SMatF ret = SM_empty(A.nrows, A.ncols, A.nvals);

  // copy non-zero structure of A
  memcpy(ret.col_sizes, A.col_sizes, A.ncols * sizeof(long));
  memcpy(ret.col_starts, A.col_starts, A.ncols * sizeof(long));
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

SMatF SM_empty_from_pos(long n_rows, long n_cols, long n_vals, long *row_pos,
                        long *col_pos) {
  check_size(n_rows, n_cols, n_vals);
  SMatF ret = SM_empty(n_rows, n_cols, n_vals);
  PosBuf buf = PB_with_size(n_vals);

  for (long i = 0; i < n_vals; ++i) {
    if (row_pos[i] < 0 || row_pos[i] >= n_rows) {
      log_err("Row position %ld out of bounds for matrix with %ld rows.",
              row_pos[i], n_rows);
      exit(EXIT_FAILURE);
    }

    if (col_pos[i] < 0 || col_pos[i] >= n_cols) {
      log_err("Column position %ld out of bounds for matrix with %ld columns.",
              col_pos[i], n_cols);
      exit(EXIT_FAILURE);
    }

    // increment corresponding row and column sizes
    ++ret.col_sizes[col_pos[i]];
    ++ret.row_sizes[row_pos[i]];
    PB_push_val(&buf, row_pos[i], col_pos[i], 0.0f);
  }

  SM_init_start_arrs(ret);

  qsort(buf.vals, n_vals, sizeof(SMEntry), SMEntry_comp_rm);
  for (long i = 0; i < n_vals; ++i) {
    ret.col_pos[i] = buf.vals[i].pos.col;
  }

  PB_free(buf);

  return ret;
}

SMatF SM_from_pos_with(long n_rows, long n_cols, long n_vals, long *row_pos,
                       long *col_pos, float *vals) {
  check_size(n_rows, n_cols, n_vals);
  SMatF ret = SM_empty(n_rows, n_cols, n_vals);
  PosBuf buf = PB_with_size(n_vals);

  for (long i = 0; i < n_vals; ++i) {
    if (row_pos[i] < 0 || row_pos[i] >= n_rows) {
      log_err("Row position %ld out of bounds for matrix with %ld rows.",
              row_pos[i], n_rows);
      exit(EXIT_FAILURE);
    }

    if (col_pos[i] < 0 || col_pos[i] >= n_cols) {
      log_err("Column position %ld out of bounds for matrix with %ld columns.",
              col_pos[i], n_cols);
      exit(EXIT_FAILURE);
    }

    // increment corresponding row and column sizes
    ++ret.col_sizes[col_pos[i]];
    ++ret.row_sizes[row_pos[i]];
    PB_push_val(&buf, row_pos[i], col_pos[i], vals[i]);
  }

  SM_init_start_arrs(ret);

  qsort(buf.vals, n_vals, sizeof(SMEntry), SMEntry_comp_rm);
  for (long i = 0; i < n_vals; ++i) {
    ret.col_pos[i] = buf.vals[i].pos.col;
    ret.vals[i] = buf.vals[i].val;
  }

  PB_free(buf);

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
        const long row_pos = rc_idx + cur_diag;

        // set column of current value position
        ret.col_pos[count_row] = row_pos;
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

bool SM_has_loc(SMatF A, long row, long col) {
  if ((row >= A.nrows || col >= A.ncols) || (row < 0 || col < 0)) {
    log_err("Position (%ld, %ld) out of bounds in Matrix of size (%ld, %ld).",
            row, col, A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  // Check if row is present in matrix
  if (SM_COL_EMPTY(A, col))
    return false;

  const long row_start = A.row_starts[row];
  const long pos = binary_search(&A.col_pos[row_start], A.row_sizes[row], col);
  return pos >= 0;
}

long SM_idx(SMatF A, long row, long col) {
  if ((row >= A.nrows || col >= A.ncols) || (row < 0 || col < 0)) {
    log_err("Position (%ld, %ld) out of bounds in Matrix of size (%ld, %ld).",
            row, col, A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  // Check if row is present in matrix
  if (SM_COL_EMPTY(A, col))
    return SM_NOT_PRESENT;

  const long row_start = A.row_starts[row];
  const long pos = binary_search(&A.col_pos[row_start], A.row_sizes[row], col);
  return pos < 0 ? SM_NOT_PRESENT : row_start + pos;
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

void SM_swap_vals(SMatF *A, SMatF *B) {
  assert(A->nvals == B->nvals);
  assert(A->nrows == B->nrows);
  assert(A->ncols == B->ncols);

  // swap pointers
  float *tmp = A->vals;
  A->vals = B->vals;
  B->vals = tmp;

  long *tmp_l = A->row_sizes;
  A->row_sizes = B->row_sizes;
  B->row_sizes = tmp_l;
}

// TODO: do i want a non-panicking setter?
void SM_set_or_panic(SMatF A, long row, long col, float val) {
  const long idx = SM_idx(A, row, col);
  assert(idx != SM_NOT_PRESENT &&
         "Can only assign to non-zero type position in SMatF");

  A.vals[idx] = val;
}

float SM_at(const SMatF A, long row, long col) {
  if (row >= A.nrows || col >= A.ncols) {
    log_err(
        "Position (%ld, %ld) is out of bounds for matrix of size (%ld, %ld)",
        row, col, A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

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

  if (A.row_sizes[row] <= col_idx)
    return SM_NOT_PRESENT;

  return A.col_pos[A.row_starts[row] + col_idx];
}

long SM_col_or_panic(SMatF A, long row, long col_idx) {
  const long col = SM_col(A, row, col_idx); // column in target
  if (col < 0) {
    log_err("Bug in row->column-iteration");
    exit(EXIT_FAILURE);
  }
  return col;
}

float *SM_ptr(SMatF A, long row, long col) {
  assert(row < A.nrows && col < A.ncols && "Position out of bounds");

  const long idx = SM_idx(A, row, col);
  return idx == SM_NOT_PRESENT ? NULL : &A.vals[idx];
}

float *SM_ptr_or_panic(SMatF A, long row, long col) {
  assert(row < A.nrows && col < A.ncols && "Position out of bounds");
  const long idx = SM_idx(A, row, col);
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

  // TODO: Optimize this
  SMatF ret = SM_empty(A.nrows, A.ncols, nvals);
  long cur_val_idx = 0;
  for (long row = 0; row < A.nrows; ++row) {   // rows in target
    for (long col = 0; col < B.ncols; ++col) { // columns in target
      if (SM_has_loc(A, row, col) || SM_has_loc(B, row, col)) {
        ret.col_pos[cur_val_idx] = col;
        ++ret.row_sizes[row];
        ++ret.col_sizes[col];
        ++cur_val_idx;
      }
    }
  }

  SM_init_start_arrs(ret);

  return ret;
}

void SM_add(SMatF A, SMatF B, SMatF target) {
  assert((A.nrows == B.nrows && A.ncols == B.ncols) &&
         (target.ncols == A.ncols && target.nrows == A.nrows) &&
         "Size mismatch. A, B, and target need to have the same size!");

  for (long row = 0; row < target.nrows; ++row) {
    for (long col_idx = 0; col_idx < target.row_sizes[row]; ++col_idx) {
      const long col = SM_col_or_panic(target, row, col_idx);
      SM_set_or_panic(target, row, col,
                      SM_at(A, row, col) + SM_at(B, row, col));
    }
  }
}

void SM_add_scl(SMatF A, SMatF B, float s_a, float s_b, SMatF target) {
  assert((A.nrows == B.nrows && A.ncols == B.ncols) &&
         (target.ncols == A.ncols && target.nrows == A.nrows) &&
         "Size mismatch. A, B, and target need to have the same size!");

  for (long row = 0; row < target.nrows; ++row) {
    for (long col_idx = 0; col_idx < target.row_sizes[row]; ++col_idx) {
      const long col = SM_col_or_panic(target, row, col_idx);
      SM_set_or_panic(target, row, col,
                      s_a * SM_at(A, row, col) + s_b * SM_at(B, row, col));
    }
  }
}

void SM_sub(SMatF A, SMatF B, SMatF target) {
  assert((A.nrows == B.nrows && A.ncols == B.ncols) &&
         (target.ncols == A.ncols && target.nrows == A.nrows) &&
         "Size mismatch. A, B, and target need to have the same size!");

  for (long row = 0; row < target.nrows; ++row) {
    for (long col_idx = 0; col_idx < target.row_sizes[row]; ++col_idx) {
      const long col = SM_col_or_panic(target, row, col_idx);
      SM_set_or_panic(target, row, col,
                      SM_at(A, row, col) - SM_at(B, row, col));
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

void SM_scl_inplace(SMatF A, float s) {
  for (long i = 0; i < A.nvals; ++i) {
    A.vals[i] *= s;
  }
}

SMatF SM_prod_prepare(SMatF A, SMatF B) {
  assert((A.ncols == B.nrows) && "Size mismatch, needs a.cols == b.rows.");

  long tmp_buf_init_cap = (A.nvals + B.nvals) * 2;
  RSBufL col_pos = RSBL_with_size(tmp_buf_init_cap);

  SMatF ret = {
      .nrows = A.nrows,
      .ncols = B.ncols,
      .nvals = 0, // tbd

      .col_sizes = calloc(B.ncols, sizeof(long)),
      .col_starts = calloc(B.ncols, sizeof(long)),
      .col_pos = NULL,
      .vals = NULL,
      .row_starts = calloc(A.nrows, sizeof(long)),
      .row_sizes = calloc(A.nrows, sizeof(long)),
  };

  // get memory requirements
  for (long t_row = 0; t_row < ret.nrows; t_row++) {
    for (long test_idx = 0; test_idx < A.row_sizes[t_row]; test_idx++) {
      const long a_col = SM_col_or_panic(A, t_row, test_idx);
      // we get a column index in a and iterate over the corresponding column in
      // B then, we set all increment all appropriate sizes this will cause more
      // values per row than we actually have
      for (long t_col_i = 0; t_col_i < B.row_sizes[a_col]; t_col_i++) {
        const long t_col = SM_col_or_panic(B, a_col, t_col_i);

        ++ret.nvals;
        ++ret.col_sizes[t_col];
        ++ret.row_sizes[t_row];

        RSBL_push_val(&col_pos, t_col);
      }
    }
  }

  // Set row / colum starts with wrong row sizes so we can iterate more easily
  long max_row_size = ret.row_sizes[0];
  for (long row = 1; row < ret.nrows; ++row) {
    ret.row_starts[row] = ret.row_starts[row - 1] + ret.row_sizes[row - 1];
    if (max_row_size < ret.row_sizes[row])
      max_row_size = ret.row_sizes[row];
  }

#define RADIX
#ifdef RADIX
  long *tmp_arr = malloc(max_row_size * sizeof(long));
#endif

  // remove duplicate columns that we got from position collection earlier
  // sort row-subsets of col_pos
  for (long row = 0; row < ret.nrows; ++row) {
    if (ret.row_sizes[row] > 1) {
#ifdef RADIX
      radix_sort(&col_pos.vals[ret.row_starts[row]], ret.row_sizes[row],
                 tmp_arr);
#else
      qsort(&col_pos.vals[ret.row_starts[row]], ret.row_sizes[row],
            sizeof(long), long_lt);
#endif
    }
  }

#ifdef RADIX
  free(tmp_arr);
#endif

  // find repeated values in col_pos and remove them
  for (long row = 0; row < ret.nrows; ++row) {
    if (ret.row_sizes[row] > 1) {
      const long cur_row_size = ret.row_sizes[row];
      long to_remove = 0;
      for (long i = ret.row_starts[row] + 1;
           i < ret.row_starts[row] + cur_row_size; ++i) {
        if (col_pos.vals[i] == col_pos.vals[i - 1]) {
          ++to_remove;
          col_pos.vals[i - 1] = SM_NOT_PRESENT;
          --ret.col_sizes[col_pos.vals[i]];
        }
      }

      ret.row_sizes[row] -= to_remove;
      ret.nvals -= to_remove;
    }
  }

  ret.col_pos = malloc(ret.nvals * sizeof(long));
  long count = 0;
  for (long idx = 0; idx < col_pos.size; ++idx) {
    if (col_pos.vals[idx] != SM_NOT_PRESENT) {
      ret.col_pos[count] = col_pos.vals[idx];
      ++count;
    }
  }

  SM_init_start_arrs(ret);

  // allocate value array memory
  ret.vals = malloc(ret.nvals * sizeof(float));

  // free col_pos
  RSBL_free(col_pos);

  return ret;
}

void SM_prod(SMatF A, SMatF B, SMatF target) {
  // size assertions
  assert(A.ncols == B.nrows && "Size mismatch between A and B");
  assert(target.nrows == A.nrows && "target.nrows == A.nrows");
  assert(target.ncols == B.ncols && "target.ncols == B.ncols");

  memset(target.vals, (int)0x0, target.nvals * sizeof(float));

  for (long t_row = 0; t_row < target.nrows; t_row++) { // rows in target
    for (long a_col_i = 0; a_col_i < A.row_sizes[t_row]; a_col_i++) {
      // iterate over columns present in row of A
      // idx is column in a / row in b
      const long idx = SM_col_or_panic(A, t_row, a_col_i);
      const long a_idx = A.row_starts[t_row] + a_col_i;

      for (long t_col_i = 0; t_col_i < B.row_sizes[idx]; t_col_i++) {
        const long t_col = SM_col_or_panic(B, idx, t_col_i); // column in target
        const long b_idx = B.row_starts[idx] + t_col_i;

        *SM_ptr_or_panic(target, t_row, t_col) +=
            A.vals[a_idx] * B.vals[b_idx];
      }
    }
  }
}

void SM_prod_scl(SMatF A, SMatF B, float s, SMatF target) {
  // size assertions
  assert(A.ncols == B.nrows && "Size mismatch between A and B");
  assert(target.nrows == A.nrows && "target.nrows == A.nrows");
  assert(target.ncols == B.ncols && "target.ncols == B.ncols");

  for (long t_row = 0; t_row < target.nrows; t_row++) { // rows in target
    for (long t_col_i = 0; t_col_i < target.row_sizes[t_row];
         t_col_i++) { // values in row of target (present columns' indices)
      const long t_col =
          SM_col_or_panic(target, t_row, t_col_i); // column in target
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

  for (long i = 0; i < target.nvals; ++i) {
    target.vals[i] *= s;
  }
}

float SM_scalar(SMatF A, SMatF B) {
  if (A.nrows != B.nrows || A.ncols != B.ncols) {
    log_err("For scalar product, A and B need to have the same dimensions. "
            "A: (%ld x %ld), B: (%ld x %ld)",
            A.nrows, A.ncols, B.nrows, B.ncols);
    exit(EXIT_FAILURE);
  }

  float ret = 0.0f;
  for (long row = 0; row < A.nrows; row++) {
    for (long col_idx = 0; col_idx < A.row_sizes[row];
         col_idx++) { // iterate on row, column in A
      // Test position existance in B
      const long col = SM_col_or_panic(A, row, col_idx);
      if (SM_has_loc(B, row, col)) {
        ret += SM_at(A, row, col) * SM_at(B, row, col);
      }
    }
  }

  return ret;
}

SMatF SM_jacobi(SMatF A, SMatF b, float rel_tol, long n_iter) {
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
  float abs_tol = SM_abs(b) * rel_tol;

  for (long iter = 0; iter < n_iter; ++iter) {
    SM_prod(A, u_k, cur_result);
    SM_sub(b, cur_result, err);
    float abs_err = SM_abs(err);

    if (abs_err < abs_tol) {
      break;
    } else {
      if (iter % 10 == 0)
        log_msg("Iteration %ld: %f > %f", iter, abs_err, abs_tol);

      for (long i = 0; i < u_k.nvals; ++i) {
        u_k1.vals[i] = 0;
      }

      for (long row = 0; row < A.nrows; ++row) {
        for (long col_idx = 0; col_idx < A.row_sizes[row]; ++col_idx) {
          const long col = SM_col_or_panic(A, row, col_idx);
          if (row == col)
            continue;
          u_k1.vals[row] -= SM_at(A, row, col) * u_k.vals[col];
        }
      }

      for (long i = 0; i < u_k.nvals; ++i) {
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
  SMatF ret = SM_empty(A.ncols, A.nrows, A.nvals);

  memcpy(ret.col_sizes, A.row_sizes, A.nrows * sizeof(long));
  memcpy(ret.row_sizes, A.col_sizes, A.ncols * sizeof(long));

  memcpy(ret.col_starts, A.row_starts, A.nrows * sizeof(long));
  memcpy(ret.row_starts, A.col_starts, A.ncols * sizeof(long));

  // Create and sort position array
  SMPos *r_pos = malloc(A.nvals * sizeof(SMPos));
  long idx = 0;
  for (long a_row = 0; a_row < A.nrows; ++a_row) {
    for (long a_col_idx = 0; a_col_idx < A.row_sizes[a_row]; ++a_col_idx) {
      const long a_col = SM_col(A, a_row, a_col_idx);
      r_pos[idx] = (SMPos){a_col, a_row};
      ++idx;
    }
  }
  qsort(r_pos, A.nvals, sizeof(SMPos), SMPos_comp_rm);

  // set values in ret
  for (long i = 0; i < A.nvals; ++i) {
    ret.col_pos[i] = r_pos[i].col;
    SM_set_or_panic(ret, r_pos[i].row, r_pos[i].col,
                    SM_at(A, r_pos[i].col, r_pos[i].row));
  }

  free(r_pos);

  return ret;
}

float SM_energy_norm(SMatF A, SMatF x, SMatF *tmp_p) {
  if (A.nrows != A.ncols) {
    log_err("Can only compute energy norm on A if A is spd. "
            "A.nrows != A.ncols (%ld !=%ld)",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  if (x.nrows != A.nrows || x.ncols != 1) {
    log_err("Size mismatch, x.nrows != A.ncols (%ld != %ld) or x.ncols != 1 "
            "(%ld != 1).",
            x.nrows, A.nrows, x.ncols);
    exit(EXIT_FAILURE);
  }

  // initialisation of temporary storage if neccessary
  SMatF tmp;
  if (tmp_p == NULL) {
    tmp = SM_empty_like(x);
  } else {
    tmp = *tmp_p;
    assert(tmp.nrows == x.nrows && tmp.ncols == x.ncols);
  }

  for (long i = 0; i < tmp.nvals; ++i) {
    tmp.vals[i] = 0.0f;
  }

  for (long row = 0; row < A.nrows; ++row) {
    for (long col_idx = 0; col_idx < A.row_sizes[row]; ++col_idx) {
      const long col = SM_col_or_panic(A, row, col_idx);
      tmp.vals[col] += SM_at(A, row, col) * SM_at(x, row, 0);
    }
  }

  float ret = SM_scalar(tmp, x);

  if (tmp_p == NULL)
    SM_free(tmp);

  return ret;
}

SMatF SM_cg(SMatF A, SMatF b, float rel_tol, long n_iter) {
  assert(A.nrows == A.ncols);
  SMatF u_k = SM_empty_like(b);
  SMatF u_k1 = SM_empty_like(b);
  SMatF r_k = SM_empty_like(b);

  const float abs_tol = SM_abs(b) * rel_tol;

  // random initialisation of target
  for (long i = 0; i < u_k.nvals; ++i) {
    u_k.vals[i] = (float)rand() / (float)RAND_MAX;
  }

  SM_prod(A, u_k, r_k);
  SM_sub(b, r_k, r_k);

  SMatF d_k = SM_clone(r_k);
  SMatF tmp = SM_clone(r_k);

  // TODO: this diverges. Find out why
  // TODO: [verify] conjugate gradient takes at most N iterations to compute
  // solution
  for (long i = 0; i < n_iter; ++i) {
    // <r_k, r_k>
    const float rk_sse = SM_sse(r_k);

    if (sqrt(rk_sse) < abs_tol)
      break;
    if (i % 10 == 0)
      log_msg("Iteration %ld / %ld: Error %f > %f", i, n_iter, sqrt(rk_sse),
              abs_tol);

    // <r_k, r_k> / <d_k, d_k>_A
    const float a_k = rk_sse / SM_energy_norm(A, d_k, &u_k1);

    // done calculating u_k+1
    // u_k1 -> u_k + tmp = u_k + a_k d_k
    SM_add_scl(u_k, d_k, 1.0f, a_k, u_k);

    // tmp -> a_k A d_k
    SM_prod_scl(A, d_k, a_k, tmp);

    // r_k+1 -> r_k - tmp = r_k - a_k A d_k
    SM_sub(r_k, tmp, r_k);

    // <r_k+1, r_k+1> / <r_k, r_k>
    const float rk1_sse = SM_sse(r_k);
    const float b_k = rk1_sse / rk_sse;

    SM_add_scl(d_k, r_k, b_k, 1.0f, d_k);
  }

  SM_free(u_k1);
  SM_free(r_k);
  SM_free(tmp);
  SM_free(d_k);
  return u_k;
}

SMatF SM_subset_diag(SMatF A) {
  if (A.nrows != A.ncols) {
    log_err("Can only subset diagonal from square SMatF. A has dimensions "
            "(%ld, %ld)",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  RSBufL elements = RSBL_with_size(A.nrows);

  for (long i = 0; i < A.nrows; ++i) {
    if (SM_has_loc(A, i, i)) {
      RSBL_push_val(&elements, i);
    }
  }
  SMatF ret = SM_empty(A.nrows, A.ncols, elements.size);

  memcpy(ret.col_pos, elements.vals, elements.size * sizeof(long));

  for (long i = 0; i < elements.size; ++i) {
    long pos = ret.col_pos[i];
    ret.col_sizes[pos] = 1;
    ret.row_sizes[pos] = 1;
  }
  SM_init_start_arrs(ret);

  for (long i = 0; i < elements.size; ++i) {
    long pos = ret.col_pos[i];
    SM_set_or_panic(ret, pos, pos, SM_at(A, pos, pos));
  }

  RSBL_free(elements);
  return ret;
}

SMatF SM_subset_trilo(SMatF A, bool diag) {
  if (A.nrows != A.ncols) {
    log_err("Can only subset diagonal from square SMatF. A has dimensions "
            "(%ld, %ld)",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  PosBuf elements = PB_with_size(A.nvals);

  for (long row = 0; row < A.nrows; ++row) {
    for (long col_i = 0; col_i < A.row_sizes[row]; ++col_i) {
      const long col = SM_col(A, row, col_i);
      if (diag && row < col)
        continue;
      if (!diag && row <= col)
        continue;

      PB_push_val(&elements, row, col, SM_at(A, row, col));
    }
  }
  SMatF ret = SM_empty(A.nrows, A.ncols, elements.size);

  for (long i = 0; i < elements.size; ++i) {
    const long col = elements.vals[i].pos.col;
    const long row = elements.vals[i].pos.row;

    ++ret.col_sizes[col];
    ++ret.row_sizes[row];
    ret.col_pos[i] = col;
  }
  SM_init_start_arrs(ret);

  for (long i = 0; i < elements.size; ++i) {
    const long col = elements.vals[i].pos.col;
    const long row = elements.vals[i].pos.row;
    const float val = elements.vals[i].val;

    SM_set_or_panic(ret, row, col, val);
  }

  PB_free(elements);
  return ret;
}

SMatF SM_subset_triup(SMatF A, bool diag) {
  if (A.nrows != A.ncols) {
    log_err("Can only subset diagonal from square SMatF. A has dimensions "
            "(%ld, %ld)",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  PosBuf elements = PB_with_size(A.nvals);

  for (long row = 0; row < A.nrows; ++row) {
    for (long col_i = 0; col_i < A.row_sizes[row]; ++col_i) {
      const long col = SM_col(A, row, col_i);
      if (diag && row > col)
        continue;
      if (!diag && row >= col)
        continue;

      PB_push_val(&elements, row, col, SM_at(A, row, col));
    }
  }
  SMatF ret = SM_empty(A.nrows, A.ncols, elements.size);

  for (long i = 0; i < elements.size; ++i) {
    const long col = elements.vals[i].pos.col;
    const long row = elements.vals[i].pos.row;

    ++ret.col_sizes[col];
    ++ret.row_sizes[row];
    ret.col_pos[i] = col;
  }
  SM_init_start_arrs(ret);

  for (long i = 0; i < elements.size; ++i) {
    const long col = elements.vals[i].pos.col;
    const long row = elements.vals[i].pos.row;
    const float val = elements.vals[i].val;

    SM_set_or_panic(ret, row, col, val);
  }

  PB_free(elements);
  return ret;
}

bool SM_is_triup(SMatF A, bool diagonal) {
  if (A.ncols != A.nrows)
    return false;

  for (long row = 0; row < A.nrows; ++row) {
    for (long col_i = 0; col_i < A.row_sizes[row]; ++col_i) {
      const long col = SM_col(A, row, col_i);
      if (diagonal) {
        if (col < row)
          return false;
      } else {
        if (col <= row)
          return false;
      }
    }
  }
  return true;
}

bool SM_is_trilo(SMatF A, bool diagonal) {
  if (A.ncols != A.nrows)
    return false;

  for (long row = 0; row < A.nrows; ++row) {
    for (long col_i = 0; col_i < A.row_sizes[row]; ++col_i) {
      const long col = SM_col(A, row, col_i);
      if (diagonal) {
        if (col > row)
          return false;
      } else {
        if (col >= row)
          return false;
      }
    }
  }
  return true;
}

void SM_back_sub(SMatF tri_up, SMatF b, SMatF target) {
  // size checks
  if (!SM_is_triup(tri_up, true)) {
    log_err("Structure mismatch: tri_up needs to be a triangular upper matrix");
    exit(EXIT_FAILURE);
  }

  if (b.ncols != 1 || b.nrows != tri_up.nrows) {
    log_err("Size mismatch: Expected b to have dimensions (%ld x 1), found "
            "(%ld x %ld).",
            tri_up.nrows, b.nrows, b.ncols);
    exit(EXIT_FAILURE);
  }

  if (target.ncols != 1 || target.nrows != tri_up.ncols) {
    log_err("Size mismatch: Expected b to have dimensions (%ld x 1), found "
            "(%ld x %ld).",
            tri_up.ncols, target.nrows, target.ncols);
    exit(EXIT_FAILURE);
  }

  SM_set_or_panic(target, target.nrows - 1, 0,
                  SM_at(b, b.nrows - 1, 0) /
                      SM_at(tri_up, tri_up.nrows - 1, tri_up.ncols - 1));

  for (long i = 0; i < tri_up.nrows; ++i) {
    if (fabs(SM_at(tri_up, i, i)) < 0.0000001f) {
      log_err("Zero element on main diagonal in backward substitution. "
              "(Division by Zero)");
      exit(EXIT_FAILURE);
    }
  }
  for (long row = b.nrows - 2; row >= 0; --row) {
    SM_set_or_panic(target, row, 0, SM_at(b, row, 0));

    for (long col_i = tri_up.row_sizes[row] - 1; col_i >= 0; --col_i) {
      const long col = SM_col(tri_up, row, col_i);
      if (col == row)
        continue;

      *SM_ptr_or_panic(target, row, 0) -=
          SM_at(target, col, 0) * SM_at(tri_up, row, col);
    }
    float ann = SM_at(tri_up, row, row);
    if (fabs(ann) < 0.0000001f) {
      log_wrn("Division by very small number in backward substitution. This "
              "might lead to unexpected results");
    }
    *SM_ptr_or_panic(target, row, 0) /= ann;
  }
}

void SM_forw_sub(SMatF tri_lo, SMatF b, SMatF target) {
  // size checks
  if (!SM_is_trilo(tri_lo, true)) {
    log_err("Structure mismatch: tri_lo needs to be a triangular lower matrix");
    exit(EXIT_FAILURE);
  }

  if (b.ncols != 1 || b.nrows != tri_lo.nrows) {
    log_err("Size mismatch: Expected b to have dimensions (%ld x 1), found "
            "(%ld x %ld).",
            tri_lo.nrows, b.nrows, b.ncols);
    exit(EXIT_FAILURE);
  }

  if (target.ncols != 1 || target.nrows != tri_lo.ncols) {
    log_err("Size mismatch: Expected b to have dimensions (%ld x 1), found "
            "(%ld x %ld).",
            tri_lo.ncols, target.nrows, target.ncols);
    exit(EXIT_FAILURE);
  }

  SM_set_or_panic(target, 0, 0, b.vals[0] / SM_at(tri_lo, 0, 0));
  for (long row = 1; row < b.nrows; row++) {
    SM_set_or_panic(target, row, 0, SM_at(b, row, 0));

    for (long col_i = 0; col_i < tri_lo.row_sizes[row]; ++col_i) {
      const long col = SM_col(tri_lo, row, col_i);
      if (row == col)
        continue;

      *SM_ptr_or_panic(target, row, 0) -=
          SM_at(target, col, 0) * SM_at(tri_lo, row, col);
    }

    float ann = SM_at(tri_lo, row, row);
    if (fabs(ann) < 0.0000001f) {
      log_wrn("Division by very small number in backward substitution. This "
              "might lead to unexpected results");
    }
    *SM_ptr_or_panic(target, row, 0) /= ann;
  }
}

SMatF SM_sor(SMatF A, SMatF b, float or, float rel_err, long n_iter) {
  if (A.nrows != A.ncols) {
    log_err("Can only perform the jacobi-method for matrix inversion on "
            "invertible matrices. Here, A.nrows != A.ncols (%ld != %ld).",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  if (b.nrows != A.nrows || b.ncols != 1) {
    log_err("Size mismatch. Expected b to be of size (%ld, %ld), but got (%ld, "
            "%ld)",
            A.nrows, 1, b.nrows, b.ncols);
    exit(EXIT_FAILURE);
  }

  SMatF B = SM_subset_trilo(A, true);
  SMatF c = SM_empty_like(b);
  SMatF r = SM_empty_like(b);
  SMatF r1 = SM_empty_like(b);
  SMatF t = SM_empty_like(b);
  SMatF tmp = SM_empty_like(b);

  // randomly initialise target
  for (long i = 0; i < t.nvals; ++i) {
    t.vals[i] = rand_float();
  }

  float max_err = rel_err * SM_abs(b);
  float err_abs = FLT_MAX;

  SM_prod(A, t, r);
  SM_sub(b, r, tmp);
  SM_swap_vals(&r, &tmp);

  for (long i = 0; i < n_iter; ++i) {
    err_abs = SM_abs(r);
    if (err_abs < max_err) {
      break;
    }

    if (i % 10 == 0)
      log_msg("Iteration %ld: %f > %f", i, err_abs, max_err);

    SM_forw_sub(B, r, c);
    SM_add_scl(t, c, 1.0, or, tmp);
    SM_swap_vals(&t, &tmp);

    SM_prod(A, c, tmp);
    SM_add_scl(r, tmp, 1.0f, - or, r1);

    SM_swap_vals(&r, &r1);
  }

  SM_free(B);
  SM_free(c);
  SM_free(r);
  SM_free(r1);
  SM_free(tmp);
  return t;
}

SMatF SM_incom_lr(SMatF A, SMatF b, float rel_err, long n_iter) {
  if (A.nrows != A.ncols) {
    log_err("Can only perform the jacobi-method for matrix inversion on "
            "invertible matrices. Here, A.nrows != A.ncols (%ld != %ld).",
            A.nrows, A.ncols);
    exit(EXIT_FAILURE);
  }

  if (b.nrows != A.nrows || b.ncols != 1) {
    log_err("Size mismatch. Expected b to be of size (%ld, %ld), but got (%ld, "
            "%ld)",
            A.nrows, 1, b.nrows, b.ncols);
    exit(EXIT_FAILURE);
  }

  // TODO: This is wrong
  SMatF L = SM_subset_trilo(A, true);
  SMatF R = SM_subset_triup(A, true);

  SMatF c = SM_empty_like(b);
  SMatF r = SM_empty_like(b);
  SMatF r1 = SM_empty_like(b);
  SMatF t = SM_empty_like(b);
  SMatF tmp = SM_empty_like(b);

  // randomly initialise target
  for (long i = 0; i < t.nvals; ++i) {
    t.vals[i] = rand_float();
  }

  float max_err = rel_err * SM_abs(b);
  float err_abs = FLT_MAX;

  SM_prod(A, t, r);
  SM_sub(b, r, tmp);
  SM_swap_vals(&r, &tmp);

  for (long i = 0; i < n_iter; ++i) {
    err_abs = SM_abs(r);
    if (err_abs < max_err) {
      break;
    }

    if (i % 10 == 0)
      log_msg("Iteration %ld: %f > %f", i, err_abs, max_err);

    SM_forw_sub(L, r, tmp);
    SM_back_sub(R, tmp, c);

    SM_add(t, c, tmp);
    SM_swap_vals(&t, &tmp);

    SM_prod(A, c, tmp);
    SM_sub(r, tmp, r1);

    SM_swap_vals(&r, &r1);
  }

  SM_free(L);
  SM_free(R);
  SM_free(c);
  SM_free(r);
  SM_free(r1);
  SM_free(tmp);
  return t;
}

void SM_vec_iteration(SMatF A, SMatF target, float rel_tol) {
  assert(A.nrows == A.ncols);

  for (long i = 0; i < target.nvals; ++i) {
    target.vals[i] = rand_float();
  }

  SMatF y = SM_empty_like(target);

  float *t_vals = target.vals;

  float err = FLT_MAX;
  long iter = 0;
  while (err > rel_tol) {
    SM_prod(A, target, y);

    const float ev = SM_at(target, 0, 0) / SM_at(y, 0, 0);
    err = 0.0f;
    for (long i = 0; i < y.nvals; ++i) {
      err += powf(ev * SM_at(y, i, 0) - SM_at(target, i, 0), 2);
    }
    const float y_len = SM_abs(y);
    err = sqrtf(err) / y_len;
    if (iter % 10 == 0)
      log_msg("Iteration %ld, EV: %f, ERR: %f", iter, ev, err);

    SM_scl_inplace(y, 1.0f / y_len);
    SM_swap_vals(&y, &target);
    ++iter;
  }

  if (t_vals != target.vals) {
    SM_swap_vals(&target, &y);
  }

  SM_free(y);
}

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
  free(A.row_sizes);
  free(A.col_sizes);
  free(A.row_starts);
  free(A.col_starts);
  free(A.col_pos);
}
