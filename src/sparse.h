#ifndef SPARSE_H
#define SPARSE_H
#include "stdbool.h"

// TODO: Update Sparse Matrix structure documentation
/* Sparse Matrix. All rows are present, but can be empty.
 *
 * vals: holds all values in matrix
 *
 * row_starts: An array of size nrows and holds indices into vals array
 * indicating the first value in a new column.
 *
 * row_sizes: An array containing the number of elements in each row
 *
 * col_pos: array matching the size of vals, holding column positions of values
 *
 * nrows, ncols, nvals: number of rows, cols and values in sparse matrix
 *
 */
typedef struct {
  long *col_sizes;  //< Number of elements in columns
  long *col_starts; //< Inidices to starts of cols (in col_idcs array)
  long *col_idcs;   //< Indices into values array (in column-major order)
  long *col_pos;    //< Column positions of values in vals array
  float *vals;      //< Values present in matrix
  long *row_starts; //< Inidices to starts of rows (in vals/ cols array)
  long *row_sizes;  //< Number of elements in rows

  long nrows, ncols, nvals; //< Number of rows, cols and values in matrix
} SMatF;                    // Sparse Matrix

#define SM_NOT_PRESENT -1

#define SM_ROW_EMPTY(mat, row) ((mat).row_sizes[row] == 0)

SMatF SM_empty(long rows, long cols, long n_vals);
SMatF SM_empty_like(SMatF A);

// TODO: Implement SM_empty_from_pos
SMatF SM_empty_from_pos(long n_rows, long n_cols, long n_vals, long *row_pos,
                        long *col_pos);

// TODO: Implement SM_empty_from_pos_vals
SMatF SM_from_pos_with(long n_rows, long n_cols, long n_vals,
                       long *row_pos, long *col_pos, float *vals);

/* @brief create empty sparse matrix with diagonal non-zero elements.
 *
 * @param diags diagonals which can be non-zero.
 *        Main diagonal is denoted by zero.
 *        Diagonals below are denoted by negative sign, ones above by positive sign.
 * @param n_diags number of diagonals to populate
 * @param size size of matrix (square)
 *
 * @return Diagonal Matrix with no values set (uninitialized memory)
 */
SMatF SM_empty_diag(long* diags, long n_diags, long size);

/* @brief Create regular sparse matrix with diagonal non-zero elements.
 *
 * @param diags diagonals which can be non-zero.
 *        Main diagonal is denoted by zero.
 *        Diagonals below are denoted by negative sign, ones above by positive sign.
 * @param n_diags number of diagonals to populate
 * @param size size of matrix (square)
 *
 * @return Diagonal Matrix with values set to values specified in diag_vals
 */
SMatF SM_diag_regular(long *diags, float* diag_vals, long n_diags, long size);

SMatF SM_vec_empty(long rows);

bool SM_has_loc(SMatF A, long row, long col);
long SM_idx(SMatF A, long row, long col);

void SM_set_or_panic(SMatF A, long row, long col, float val);
float SM_at(SMatF A, long row, long col);

long SM_col(SMatF A, long row, long col_idx);
long SM_col_or_panic(SMatF A, long row, long col_idx);

SMatF SM_prod_prepare(SMatF A, SMatF B);
void SM_prod(SMatF A, SMatF B, SMatF target);

void SM_print(SMatF A);
void SM_print_nonzero(SMatF A);
void SM_print_meta(SMatF A);
void SM_print_shape(SMatF A);
#endif // SPARSE_H
