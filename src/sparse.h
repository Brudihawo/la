#ifndef SPARSE_H
#define SPARSE_H
#include "stdbool.h"

// TODO: Update Sparse Matrix structure documentation
/* Sparse Matrix. All rows are present, but can be empty.
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

/* @brief Create empty sparse matrix
 *        Only allocates memory. No information on position of values is present
 *
 * @param rows number of rows
 * @param cols number of columns
 * @param n_vals number of non-zero values
 *
 * @return Matrix with malloced (but uninitialized) memory
 */
SMatF SM_empty(long rows, long cols, long n_vals);

/* @brief Create empty sparse matrix by copying size and non-zero structure from
 * other. Does not set actual values.
 *
 * @param A SMatF to copy structure from
 *
 * @return Matrix with malloced (but uninitialized) memory
 */
SMatF SM_empty_like(SMatF A);

/* @brief Create a clone of A
 *
 * @param A SMatF to clone
 */
SMatF SM_clone(SMatF A);

/* @brief Create empty SMatF from ordered list of non-zero positions
 *
 * @param n_rows number of rows
 * @param n_cols number of columns
 * @param n_vals number of values
 * @param row_pos row positions that can be non-zero (in row-major order)
 * @param col_pos column positions that can be non-zero (in row-major order)
 *
 * @return SMatF with corresponding non-zero structure but empty values
 */
SMatF SM_empty_from_pos(long n_rows, long n_cols, long n_vals, long *row_pos,
                        long *col_pos);

/* @brief Create SMatF from ordered list of non-zero positions
 *
 * @param n_rows number of rows
 * @param n_cols number of columns
 * @param n_vals number of values
 * @param row_pos row positions that can be non-zero (in row-major order)
 * @param col_pos column positions that can be non-zero (in row-major order)
 * @param vals values of positions (value array is set to this at the moment)
 *
 * @return SMatF with corresponding non-zero structure and set values
 */
SMatF SM_from_pos_with(long n_rows, long n_cols, long n_vals, long *row_pos,
                       long *col_pos, float *vals);

/* @brief create empty sparse matrix with diagonal non-zero elements.
 *
 * @param diags diagonals which can be non-zero.
 *        Main diagonal is denoted by zero.
 *        Diagonals below are denoted by negative sign, ones above by positive
 * sign.
 * @param n_diags number of diagonals to populate
 * @param size size of matrix (square)
 *
 * @return Diagonal Matrix with no values set (uninitialized memory)
 */
SMatF SM_empty_diag(long *diags, long n_diags, long size);

/* @brief Create regular sparse matrix with diagonal non-zero elements
 *        and set each diagonal to a value.
 *
 * @param diags diagonals which can be non-zero.
 *        Main diagonal is denoted by zero.
 *        Diagonals below are denoted by negative sign, ones above by positive
 * sign.
 * @param n_diags number of diagonals to populate
 * @param size size of matrix (square)
 *
 * @return Diagonal Matrix with values set to values specified in diag_vals
 */
SMatF SM_diag_regular(long *diags, float *diag_vals, long n_diags, long size);

/* @brief Create an empty SMatF vector
 *
 * @param rows rows of vector
 *
 * @return SMatF with all positions present in (rows x 1)
 */
SMatF SM_vec_empty(long rows);

/* @brief compute square root of sum of squared elements in A
 *        This is not a standard Matrix norm!
 *        Use this to compute the absolute value of a vector as SMatF
 */
float SM_abs(SMatF A);

/* @brief determine whether a position in SMatF can be non-zero
 *
 * @param A SMatF to test
 * @param row row to test
 * @param col column to test
 *
 * @return value can be non-zero
 */
bool SM_has_loc(SMatF A, long row, long col);

/* @brief Compute index into value array for position in SMatF
 *
 * @param A SMatF to query
 * @param row row
 * @param col column
 *
 * @return index into value array or SM_NOT_PRESENT if not present.
 */
long SM_idx(SMatF A, long row, long col);

/* @brief Check for equal non-zero structure of A and B
 *
 * @param A
 * @param B
 *
 * @return whether A and B have the same non-zero structure.
 */
bool SM_structure_eq(SMatF A, SMatF B);

/* @brief Check for equality of A and B
 *
 * @param A
 * @param B
 *
 * @return whether A and B are equal (element-wise comparison)
 */
bool SM_eq(SMatF A, SMatF B);

/* @brief Get column from column index and row
 *
 * @param A SMatF to query
 * @param row row
 * @param col_idx column index
 *
 * @return column index or SM_NOT_PRESENT if not present
 */
long SM_col(SMatF A, long row, long col_idx);

/* @brief Get column from column index and row. Crashes if not present
 *
 * @param A SMatF to query
 * @param row row
 * @param col_idx column index
 *
 * @return column
 */
long SM_col_or_panic(SMatF A, long row, long col_idx);

/* @brief Set value in SMatF. Crash on trying to set non-zero or out of bounds
 * value
 *
 * @param A SMatF to modify
 * @param row row to set
 * @param col_idx column to set
 * @param val value to set
 */
void SM_set_or_panic(SMatF A, long row, long col, float val);

/* @brief return value at position in SMatF
 *
 * @param A SMatF to query
 * @param row row
 * @param col column
 *
 * @return value at position in SMatF
 */
float SM_at(SMatF A, long row, long col);

/* @brief prepare target for element wise sum of A and B
 *
 * @param A summand A
 * @param B summand B
 *
 * @return SMatF with expected structure to hold A+B
 */
SMatF SM_addsub_prepare(SMatF A, SMatF B);

/* @brief Elementwise Addition of A and B
 *
 * @param A summand A
 * @param B summand B
 * @param target location to save to (can only be A or B if A and B have the
 * same non-zero structure.)
 */
void SM_add(SMatF A, SMatF B, SMatF target);

/* @brief Elementwise subtraction of B from A
 *
 * @param A SMatF to subtract from
 * @param B SMatF to subtract
 * @param target location to save to (can only be A or B if A and B have the
 * same non-zero structure.)
 */
void SM_sub(SMatF A, SMatF B, SMatF target);

/* @brief Elementwise scalar multiplication of A by s
 *
 * @param A SMatF to Scale
 * @param s Scaling factor
 * @param target where to save results. Can be A
 */
void SM_scl(SMatF A, float s, SMatF target);

/* @brief prepare target for matrix product of A and B
 *
 * @param A left factor of matrix product
 * @param B right factor of matrix product
 *
 * @return SMatF with expected structure to hold A*B
 */
SMatF SM_prod_prepare(SMatF A, SMatF B);

/* @brief SMatF matrix multiplication
 *
 * @param A left factor of matrix product
 * @param B right factor of matrix product
 * @param target SMatF to write results into (create by SM_prod_prepare)
 *
 * @return SMatF with expected structure to hold A*B
 */
void SM_prod(SMatF A, SMatF B, SMatF target);

/* @brief Print contents of A
 *
 * @param A SMatF to print
 */
void SM_print(SMatF A);

/* @brief Print non-zero structure of A
 *
 * @param A SMatF to print
 */
void SM_print_nonzero(SMatF A);

/* @brief Print meta information of A (for debug purposes)
 *
 * @param A SMatF to print
 */
void SM_print_meta(SMatF A);

/* @brief Print shape of A
 *
 * @param A SMatF to print
 */
void SM_print_shape(SMatF A);

/* @brief Free SMatF memory
 *
 * @param A SMatF to free
 */
void SM_free(SMatF A);
#endif // SPARSE_H
