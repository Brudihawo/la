#include "la.h"

#include "assert.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "log.h"
#include <string.h>

bool almost_eq(float a, float b) {
  if (a == b)
    return true;
  float dif = a < b ? a - b : b - a;
  if (a == 0.0f)
    return dif / b < REL_ERROR;
  return dif / a < REL_ERROR;
}

RowPerm RP_new(long rows) {
  long *vals = malloc(rows * sizeof(long));
  for (long i = 0; i < rows; i++) {
    vals[i] = i;
  }
  return (RowPerm){
      .order = vals,
      .size = rows,
      .n_swaps = 0,
  };
}

void RP_add(RowPerm *perms, long i, long j) {
  long tmp = perms->order[j];
  perms->order[j] = perms->order[i];
  perms->order[i] = tmp;
  perms->n_swaps++;
}

MatF MF_empty(long rows, long cols) {
  return (MatF){
      .vals = malloc(rows * cols * sizeof(float)),
      .rows = rows,
      .cols = cols,
  };
}

MatF MF_with(long rows, long cols, float init_val) {
  MatF ret = (MatF){
      .vals = malloc(rows * cols * sizeof(float)),
      .rows = rows,
      .cols = cols,
  };

  for (long row = 0; row < rows; row++) {
    for (long col = 0; col < cols; col++) {
      ret.vals[MF_IDX(ret, row, col)] = init_val;
    }
  }

  return ret;
}

MatF MF_from(long rows, long cols, float *vals) {
  // TODO: can i somehow check that rows * cols = sizeof vals?
  return (MatF){.rows = rows, .cols = cols, .vals = vals};
}

MatF MF_clone(MatF A) {
  MatF ret = {.rows = A.rows,
              .cols = A.cols,
              .vals = malloc(A.rows * A.cols * sizeof(float))};
  memcpy(ret.vals, A.vals, A.rows * A.cols * sizeof(float));
  return ret;
}

MatF MF_transpose(MatF A) {
  MatF ret = MF_EMPTY_LIKE_T(A);
  for (long row = 0; row < A.rows; row++) {
    for (long col = 0; col < A.cols; col++) {
      *MF_PTR(ret, col, row) = MF_AT(A, row, col);
    }
  }
  return ret;
}

void VEC_r_perm(MatF v, RowPerm perms) {
  if (perms.n_swaps == 0)
    return;

  long performed_perms = 0;
  for (long row = 0; row < v.rows; row++) {
    if (performed_perms == perms.n_swaps)
      break;
    if (perms.order[row] != row) {
      float tmp = VEC_AT(v, row);
      *VEC_PTR(v, row) = VEC_AT(v, perms.order[row]);
      *VEC_PTR(v, perms.order[row]) = tmp;
      performed_perms++;
    }
  }
}

bool MF_eq(MatF A, MatF B) {
  if (A.cols != B.cols || A.rows != B.rows)
    return false;

  // pointer to value arrays is the same
  if (A.vals == B.vals)
    return true;

  for (long i = 0; i < A.cols; i++) {
    for (long j = 0; j < A.rows; j++) {
      if (!almost_eq(MF_AT(A, i, j), MF_AT(B, i, j))) {
        return false;
      }
    }
  }
  return true;
}

float MF_at(MatF A, long row, long col) { return A.vals[MF_IDX(A, row, col)]; }

float *MF_ptr(MatF A, long row, long col) {
  return &A.vals[MF_IDX(A, row, col)];
}

MatF MF_subset_copy(MatF A, long start_row, long end_row, long start_col,
                    long end_col) {
  // TODO: Handle size checks nicely
  assert(start_row < end_row && start_col < end_col && start_row > 0 &&
         end_row > 0);

  assert(end_row <= A.rows && end_col <= A.cols);
  long rows = end_row - start_row;
  long cols = end_col - start_col;
  MatF ret = {
      .rows = rows,
      .cols = cols,
      .vals = malloc(rows * cols * sizeof(float)),
  };

  for (long row = start_row; row < end_row; row++) { // copy line by line
    memcpy(MF_PTR(ret, row - start_row, 0), MF_PTR(A, row, start_col),
           cols * sizeof(float));
  }

  return ret;
}

// Scalar and Elementwise Operations
float MF_det(MatF A) {
  // Gaussian elimination for creating triangular upper matrix
  // Then, determinant is the product of main diagonal items
  MatF tmp = MF_clone(A);
  MF_make_tri_up(tmp);

  float ret = 0.0f;
  for (long i = 0; i < tmp.rows; i++) {
    ret *= MF_AT(tmp, i, i);
  }

  return ret;
}

void MF_add(MatF a, MatF b, MatF target) {
  // TODO: Proper error reporting for size checks
  assert(((a.rows == b.rows && a.cols == b.cols) &&
          (a.rows == target.rows && a.cols == target.cols)) &&
         "a, b and target need to be of same size.");

  for (long row = 0; row < a.rows; row++) {
    for (long col = 0; col < a.cols; col++) {
      long idx = MF_IDX(a, row, col);
      target.vals[idx] = a.vals[idx] + b.vals[idx];
    }
  }
}

void MF_sub(MatF a, MatF b, MatF target) {
  assert(((a.rows == b.rows && a.cols == b.cols) &&
          (a.rows == target.rows && a.cols == target.cols)) &&
         "a, b and target need to be of same size.");

  for (long row = 0; row < a.rows; row++) {
    for (long col = 0; col < a.cols; col++) {
      long idx = MF_IDX(a, row, col);
      target.vals[idx] = a.vals[idx] - b.vals[idx];
    }
  }
}

void MF_scl(MatF a, float s, MatF target) {
  assert((a.rows == target.rows && a.cols == target.cols) &&
         "a and target need to be of same size.");

  for (long row = 0; row < a.rows; row++) {
    for (long col = 0; col < a.cols; col++) {
      long idx = MF_IDX(a, row, col);
      target.vals[idx] = s * a.vals[idx];
    }
  }
}

// TODO: parallelisation via openmp
float MF_scalar(MatF a, MatF b) {
  assert((a.rows == b.rows && a.cols == b.cols) &&
         "a and b need to be of same size.");

  float ret = 0;
  for (long row = 0; row < a.rows; row++) {
    for (long col = 0; col < a.cols; col++) {
      long idx = MF_IDX(a, row, col);
      ret += a.vals[idx] * b.vals[idx];
    }
  }

  return ret;
}

void MF_element(MatF a, MatF b, MatF target) {
  assert(((a.rows == b.rows && a.cols == b.cols) &&
          (a.rows == target.rows && a.cols == target.cols)) &&
         "a, b and target need to be of same size.");

  for (long row = 0; row < a.rows; row++) {
    for (long col = 0; col < a.cols; col++) {
      long idx = MF_IDX(a, row, col);
      target.vals[idx] = a.vals[idx] * b.vals[idx];
    }
  }
}

void MF_prod(MatF a, MatF b, MatF target) {
  assert(
      (a.cols == b.rows && (target.rows == a.rows && target.cols == b.cols)) &&
      "Size mismatch, needs a.cols == b.rows, target.rows == a.rows, "
      "target.cols == b.cols.");

  assert((a.vals != target.vals && b.vals != target.vals) &&
         "Writing into an operand matrix will cause unexpected results.");

  for (long col = 0; col < target.cols; col++) {
    for (long row = 0; row < target.rows; row++) {
      long target_idx = MF_IDX(target, row, col);
      target.vals[target_idx] = 0.0f;
      for (int idx = 0; idx < b.rows; idx++) {
        target.vals[target_idx] += MF_AT(a, row, idx) * MF_AT(b, idx, col);
      }
    }
  }
}

float VEC_abs(MatF v) {
  float ret = 0.0f;
  for (long i = 0; i < v.rows; i++) {
    ret += v.vals[i] * v.vals[i];
  }
  return sqrt(ret);
}

// Solving Equations
void MF_forw_sub(MatF tri_lo, MatF b, MatF target) {
  assert((target.rows == b.rows && target.cols == 1) &&
         (tri_lo.rows == tri_lo.cols && tri_lo.rows == b.rows));

  target.vals[0] = b.vals[0] / MF_AT(tri_lo, 0, 0);
  for (long i = 1; i < b.rows; i++) {
    target.vals[i] = b.vals[i];
    for (long j = 0; j < i; j++) {
      target.vals[i] -= target.vals[j] * MF_AT(tri_lo, i, j);
    }
    target.vals[i] /= MF_AT(tri_lo, i, i);
  }
}

void MF_forw_sub_lu(MatF tri_lo, MatF b, MatF target) {
  if (target.rows != b.rows || target.cols != 1) {
    log_err("Unexpected target size: Got (%ld x %ld), expected (%ld x %ld).\n",
            target.rows, target.cols, b.rows, 1);
    exit(1);
  }

  if (tri_lo.rows != tri_lo.cols || tri_lo.rows != b.rows) {
    log_err("Unexpected Operand sizes. Cannot compute for tri_lo: (%ld x %ld), "
            "b: (%ld x %ld).\n",
            tri_lo.rows, tri_lo.cols, b.rows, b.cols);
    exit(1);
  }

  target.vals[0] = b.vals[0];
  for (long i = 1; i < b.rows; i++) {
    target.vals[i] = b.vals[i];
    for (long j = 0; j < i; j++) {
      target.vals[i] += target.vals[j] * MF_AT(tri_lo, i, j);
    }
  }
}

void MF_back_sub(const MatF tri_up, MatF b, MatF target) {
  assert((target.rows == b.rows && target.cols == 1) &&
         (tri_up.rows == tri_up.cols && tri_up.rows == b.rows));

  target.vals[b.rows - 1] =
      b.vals[b.rows - 1] / MF_AT(tri_up, tri_up.rows - 1, tri_up.cols - 1);
  for (long i = b.rows - 2; i >= 0; i--) {
    target.vals[i] = b.vals[i];
    for (long j = i + 1; j < b.rows; j++) {
      target.vals[i] -= target.vals[j] * MF_AT(tri_up, i, j);
    }
    target.vals[i] /= MF_AT(tri_up, i, i);
  }
}

static void swap_cols(MatF A, long i, long j, float *cpybuf) {
  // k -> buf
  // n -> k
  // buf -> n
  memcpy(cpybuf, MF_PTR(A, i, 0), A.cols * sizeof(float));
  memcpy(MF_PTR(A, i, 0), MF_PTR(A, j, 0), A.cols * sizeof(float));
  memcpy(MF_PTR(A, j, 0), cpybuf, A.cols * sizeof(float));
}

RowPerm *MF_lu_decomp(MatF A) {
  assert(A.rows == A.cols && "A needs to be a square Matrix");

  RowPerm *rperms = malloc(sizeof(RowPerm));
  *rperms = RP_new(A.rows);
  // TODO: implement stability measures in gauss eliminiation
  // Transform to triangular upper matrix
  // Save Operations into lower triangular half
  float *cpybuf = malloc(A.cols * sizeof(float)); // copy buffer
  for (long i = 0; i < A.cols - 1; i++) {         // zeroing of i-th column
    for (long n = i + 1; n < A.rows; n++) {       // row n
      for (long j = i + 1; j < A.cols; j++) {     // column j
        if (MF_AT(A, i, i) == 0.0f) {             // pivot
          long k = 0;
          for (; MF_AT(A, k, i) == 0.0f; k++)
            ;
          // swap rows of a
          swap_cols(A, i, k, cpybuf);

          // increase permutation count and record permutation
          RP_add(rperms, i, k);
        }
        *MF_PTR(A, n, j) -= MF_AT(A, n, i) / MF_AT(A, i, i) * MF_AT(A, i, j);
      }
      *MF_PTR(A, n, i) = -MF_AT(A, n, i) / MF_AT(A, i, i);
    }
  }

  if (rperms->n_swaps == 0)
    return NULL;
  return rperms;
}

void MF_lu_solv(const MatF LU, MatF b, RowPerm *perms, MatF target) {
  MatF y = MF_EMPTY_LIKE(b);

  if (perms != NULL)
    VEC_r_perm(b, *perms);

  MF_forw_sub_lu(LU, b, y);
  MF_back_sub(LU, y, target);

  if (perms != NULL) {
    VEC_r_perm(b, *perms);
  }

  MF_FREE(y);
}

void MF_gauss_elim(MatF A, MatF b, MatF target) {
  assert(A.rows == A.cols && "A needs to be a square Matrix");
  assert((target.cols == 1 && b.cols == 1) &&
         "target and b need to be column vectors");
  assert((A.rows == b.rows && target.rows == b.rows) &&
         "Size mismatch: Needed A.rows == b.rows, target.rows == b.rows");

  RowPerm perms = RP_new(A.rows);

  // TODO: implement pivoting and stability measures in gauss eliminiation
  // TODO: do i want to zero the right lower half of A?
  // Transform to triangular upper matrix
  float *cpybuf = malloc(A.cols * sizeof(float)); // copy buffer
  for (long i = 0; i < A.cols - 1; i++) {         // zeroing of i-th column
    for (long n = i + 1; n < A.rows; n++) {       // row n
      if (MF_AT(A, i, i) == 0.0f) {               // pivot
        long k = 0;
        for (; MF_AT(A, k, i) == 0.0f; k++)
          ;
        // swap rows of a
        swap_cols(A, i, k, cpybuf);
        // swap values in b
        float tmp = VEC_AT(b, k);
        *VEC_PTR(b, k) = VEC_AT(b, i);
        *VEC_PTR(b, i) = tmp;

        // increase permutation count and record permutation
        RP_add(&perms, i, k);
      }
      for (long j = i + 1; j < A.cols; j++) { // column j
        *MF_PTR(A, n, j) -= MF_AT(A, n, i) / MF_AT(A, i, i) * MF_AT(A, i, j);
      }
      *VEC_PTR(b, n) -= MF_AT(A, n, i) / MF_AT(A, i, i) * VEC_AT(b, i);
      *MF_PTR(A, n, i) = 0.0f;
    }
  }
  free(cpybuf);

  // invertibility check
  // TODO: check for invertibility of A
  for (long i = 0; i < A.cols; i++) {
    if (MF_at(A, i, i) == 0.0f) {
      log_err("Diagonal of triangular upper version contains zero-element at "
              "(%ld, %ld)!",
              i, i);
      exit(1);
    }
  }

  MF_back_sub(A, b, target);

  RP_FREE(perms);
}

void MF_make_tri_up(MatF A) {
  assert(A.rows == A.cols && "A needs to be a square Matrix");

  // TODO: implement pivoting and stability measures in gauss eliminiation
  // TODO: do i want to zero the right lower half of A?
  // Transform to triangular upper matrix
  for (long i = 0; i < A.cols - 1; i++) {     // zeroing of i-th column
    for (long n = i + 1; n < A.rows; n++) {   // row n
      for (long j = i + 1; j < A.cols; j++) { // column j
        *MF_PTR(A, n, j) -= MF_AT(A, n, i) / MF_AT(A, i, i) * MF_AT(A, i, j);
      }
    }
  }

  for (long i = 1; i < A.rows; i++) {
    for (long j = 0; j < i; j++) {
      *MF_PTR(A, i, j) = 0.0f;
    }
  }
}

bool MF_invertible(MatF A) {
  MatF tmp = MF_clone(A);
  MF_make_tri_up(A);

  for (long i = 0; i < tmp.rows; i++) {
    if (MF_AT(tmp, i, i) == 0.0f)
      return false;
  }

  return true;
}

void MF_print(MatF A) {
  printf("[");
  for (long row = 0; row < A.rows; row++) {
    if (row != 0) {
      printf(" ");
    }
    printf(" ");
    for (long col = 0; col < A.cols; col++) {
      printf("%9.2f", MF_AT(A, row, col));
    }
    if (row != A.rows - 1) {
      printf("\n");
    }
  }
  printf(" ]\n");
}

void MF_print_shape(MatF A) { printf("(%ld x %ld)", A.rows, A.cols); }

// typedef struct {
//   long* col_sizes;  //< Number of elements in columns
//   long* col_starts; //< Inidices to starts of cols (in col_idcs array)
//   long* col_idcs;   //< Indices into values array (in column-major order)
//   long* col_pos;    //< Column positions of values in vals array
//   float* vals;      //< Values present in matrix
//   long* row_starts; //< Inidices to starts of rows (in vals/ cols array)
//   long* row_sizes;  //< Number of elements in rows
//
//   long nrows, ncols, nvals; //< Number of rows, cols and values in matrix
// } SMatF; // Sparse Matrix

// #define SM_col_iter(mat, col_idx_var_name, col_var_name)
// #define SM_row_iter(mat, row_idx_var_name, row_var_name)

SMatF SM_empty(long rows, long cols, long n_vals) {
  return (SMatF){
      .nrows = rows,
      .ncols = cols,
      .nvals = n_vals,

      .col_sizes = malloc(cols * sizeof(long)),
      .col_starts = malloc(cols * sizeof(long)),
      .col_idcs = malloc(n_vals * sizeof(long)),
      .col_pos = malloc(n_vals * sizeof(long)),
      .vals = malloc(n_vals * sizeof(float)),
      .row_starts = malloc(rows * sizeof(long)),
      .row_sizes = malloc(rows * sizeof(long)),
  };
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

// TODO: check correctness
bool SM_has_loc(SMatF A, long row, long col) {
  // Check if row is present in matrix
  if (SM_ROW_EMPTY(A, row))
    return false;

  const long row_start = A.row_starts[row];
  for (long i = 0; i < A.row_sizes[row] && i < A.col_pos[row_start + i]; i++) {
    if (A.col_pos[row_start + i] == col)
      return true;
  }
  return false;
}

long SM_idx(SMatF A, long row, long col) {
  assert(row < A.nrows && col < A.ncols && "Position out of bounds");

  if (SM_ROW_EMPTY(A, row))
    return SM_NOT_PRESENT;

  const long row_start = A.row_starts[row];
  for (long i = 0; i < A.row_sizes[row] && i < A.col_pos[row_start + i]; i++) {
    if (A.col_pos[row_start + i] == col)
      return row_start + i;
  }
  return SM_NOT_PRESENT;
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
  assert(row < A.nrows && col_idx < A.ncols && "Position out of bounds");

  if (SM_ROW_EMPTY(A, row) || col_idx < A.row_sizes[row])
    return SM_NOT_PRESENT;

  return A.col_pos[A.row_starts[row] + col_idx];
}

long SM_col_or_panic(SMatF A, long row, long col_idx) {
  long col = SM_col(A, row, col_idx); // column in target
  assert(col != SM_NOT_PRESENT && "Bug in row->column-iteration");
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

/* @brief perform matrix multiplication with two sparse matrices
 *        assumes target matches required size and positions
 *        (created with SM_prod_prepare ...)
 *
 * saves output into target
 */
void SM_prod(SMatF A, SMatF B, SMatF target) {
  // size assertions
  assert(A.ncols == B.nrows && "Size mismatch between A and B");
  assert(target.nrows == A.nrows && "target.nrows == A.nrows");
  assert(target.ncols == B.ncols && "target.ncols == B.ncols");

  for (long t_row = 0; t_row < target.nrows; t_row++) { // rows in target
    for (long t_col_i = 0; t_col_i < target.row_sizes[t_row];
         t_col_i++) { // values in row of target (present columns' indices)
      long t_col = SM_col_or_panic(target, t_row, t_col_i); // column in target
      SM_set_or_panic(target, t_row, t_col, 0.0f);

      for (long a_col_i = 0; a_col_i < A.row_sizes[t_row]; a_col_i++) {
        long idx = SM_col_or_panic(A, t_row, a_col_i); // column in a

        if (SM_has_loc(B, idx, t_col)) { // check if value is non-zero in b
          *SM_ptr_or_panic(target, t_row, t_col) +=
              SM_at(A, t_row, idx) * SM_at(B, idx, t_col);
        }
      }
    }
  }
}

/* @brief prepare target for matrix multiplication
 *        allocates memory in target to accomodate neccesary values
 */
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
        if (SM_has_loc(A, t_row, test_idx) && SM_has_loc(B, test_idx, t_col)) {
          ret.nvals++;
          ret.col_sizes[t_col]++;
          ret.row_sizes[t_row]++;

          break;
        }
      }
    }
  }

  // Set row starts
  ret.row_starts[0] = 0;
  for (long t_row = 1; t_row < ret.nrows; t_row++) {
    ret.row_starts[t_row] =
        ret.row_starts[t_row - 1] + ret.row_sizes[t_row - 1];
  }

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

        if (SM_has_loc(A, t_row, test_idx) && SM_has_loc(B, test_idx, t_col)) {
          ret.col_pos[cur_idx] = t_col;
          cur_idx++;
          break;
        }
      }
    }
  }

  // TODO: Do i need column-wise iteration in SMatF?
  // Assignments for column-wise iteration

  // set column starts
  ret.col_starts[0] = 0;
  for (long t_col = 1; t_col < ret.ncols; t_col++) {
    ret.col_starts[t_col] =
        ret.col_starts[t_col - 1] + ret.col_sizes[t_col - 1];
  }

  // set col_idcs
  cur_idx = 0;
  for (long t_col = 0; t_col < ret.ncols; t_col++) {
    for (long t_row = 0; t_row < ret.nrows;
         t_row++) { // iterate on row, column in target
      for (long test_idx = 0; test_idx < A.row_sizes[t_row]; test_idx++) {
        // Test idx in A and B for column / row respectively
        if (SM_has_loc(ret, t_row, t_col)) {
          ret.col_idcs[cur_idx] = SM_idx(ret, t_row, t_col);
          cur_idx++;
          break;
        }
      }
    }
  }

  return ret;
}

void SM_print(SMatF A) {
  printf("[");
  for (long row = 0; row < A.nrows; row++) {
    if (row != 0) {
      printf(" ");
    }
    printf(" ");
    for (long col = 0; col < A.ncols; col++) {
      printf("%9.2f", SM_at(A, row, col));
    }
    if (row != A.nrows - 1) {
      printf("\n");
    }
  }
  printf(" ]\n");
}
