#include "stdbool.h"
#include "stdlib.h"
#include "string.h"

#define REL_ERROR 0.00001
#define LA_PRINT_FMT "%5.2f"

typedef struct {
  float *vals;
  long rows, cols;
} MatF;

typedef struct {
  long *order;
  long size;
  long n_swaps;
} RowPerm;

RowPerm RP_new(long rows);

void RP_add(RowPerm *perms, long i, long j);

#define RP_FREE(perms) free(perms.order)

#define MF_IDX(mat, row, col) (mat.cols * (row) + (col))

/* @brief Create an empty MatF
 *
 * @param rows number of rows
 * @param cols number of cols
 *
 * @return MatF initalised to 0
 */
MatF MF_empty(long rows, long cols);

/* @brief Create an empty MatF with dimensions of mat
 *
 * @param mat matrix with dimensions to copy
 */
#define MF_EMPTY_LIKE(mat) MF_empty((mat).rows, (mat).cols)

/* @brief Create an empty MatF with dimensions of mat Transposed
 *
 * @param mat matrix with dimensions to copy
 */
#define MF_EMPTY_LIKE_T(mat) MF_empty((mat).cols, (mat).rows)

/* @brief clone A into new matrix
 *
 * @param A: matrix to clone
 *
 * @return Clone of A
 */
MatF MF_clone(MatF A);

/* @brief Create an initialised MatF
 *
 * @param rows number of rows
 * @param cols number of cols
 *
 * @return MatF initialised with init_val
 */
MatF MF_with(long rows, long cols, float init_val);

/* @brief Shorthand for empty vector (Nx1 MatF)
 *
 * @param size length of vector
 *
 * @return Nx1 MatF initilised to 0
 */
#define VEC_EMPTY(size) MF_empty(size, 1);

/* @brief Shorthand for initialised vector (Nx1 MatF)
 *
 * @param size length of vector
 * @param init_val value to initialise vector with
 *
 * @return Nx1 MatF initilised to init_val elementwise
 */
#define VEC_WITH(size, init_val) MF_with(size, 1, init_val);

/* @brief create MatF from float array
 *
 * @param rows number of rows
 * @param cols number of columns
 * @param vals array of values
 */
MatF MF_from(long rows, long cols, float *vals);

/* @brief Transpose MatF
 *
 * @param A: Matrix to transpose
 * @return Transposed Matrix
 */
MatF MF_transpose(MatF A);

/* @brief Apply Row permutations to A
 *
 * @param v: vector (Nx1 MatF) to permute
 * @param perms: RowPermutations to apply
 */
void VEC_r_perm(MatF v, RowPerm perms);

// /* @brief Apply Row permutations to A
//  * @param A: MatF to permute
//  * @param perms: RowPermutations to apply
//  */
// void MF_r_perm(MatF A, RowPerm perms);

/* @brief: free value array of MatF
 *
 * @param mat MatF to free
 */
#define MF_FREE(mat) free((mat).vals);

bool MF_eq(MatF A, MatF B);

/* @brief Get value in MatF
 *
 * @param A matrix to query
 * @param row row
 * @param col column
 *
 * @return Value at A[row, col]
 */
#define MF_AT(A, row, col) A.vals[MF_IDX(A, row, col)]

/* @brief Get value in vector (Nx1 MatF)
 *
 * @param vec vector to query
 * @param idx row
 *
 * @return Value at vec[idx]
 */
#define VEC_AT(vec, idx) MF_AT(vec, idx, 0);

/* @brief Get pointer to location in MatF
 *
 * @param A matrix to query
 * @param row row
 * @param col column
 *
 * @return pointer to A[row, col]
 */
#define MF_PTR(A, row, col) (&A.vals[MF_IDX(A, row, col)])

/* @brief Get pointer to value in vector (Nx1 MatF)
 *
 * @param vec vector to query
 * @param idx row
 *
 * @return pointer to vec[idx]
 */
#define VEC_PTR(vec, idx) MF_PTR(vec, idx, 0)

/* @brief Create a subset of a matrix
 *
 * @param A Matrix to subset
 * @param start_row inclusive row to start at
 * @param end_row exclusive row to end at
 * @param start_col inclusive column to start at
 * @param end_col exclusive column to end at
 *
 * @return Subset MatF (copying!)
 */
MatF MF_subset_copy(MatF A, long start_row, long end_row, long start_col,
                    long end_col);

// Scalar and Elementwise Operations

/* @brief compute determinant of A
 *
 * @param A matrix to compute determinant for
 *
 * @return determinant
 */
float MF_det(MatF A);

/* @brief Elementwise addition of a and b
 *
 * @param a first summand
 * @param b second summand
 * @param target where to put results
 */
void MF_add(MatF a, MatF b, MatF target);

/* @brief Elementwise subtraction of a and b (a - b)
 *
 * @param a matrix
 * @param b matrix to subtract
 * @param target where to put results
 */
void MF_sub(MatF a, MatF b, MatF target);

/* @brief Elementwise scaling of a by s
 *
 * @param a MatF to scale
 * @param s scaling factor
 * @param target where to put results
 */
void MF_scl(MatF a, float s, MatF target);

/* @brief Scalar Product of a and b (sum of products of elements)
 *
 * @param a first operand
 * @param b second operand
 *
 * @return scalar product of a and b
 */
float MF_scalar(MatF a, MatF b);

/* @brief element wise product of a and b
 *
 * @param a first operand
 * @param b second operand
 * @param target where to put results
 */
void MF_element(MatF a, MatF b, MatF target);

/* @brief Matrix product of a and b
 *
 * @param a first operand
 * @param b second operand
 * @param target where to put results
 */
void MF_prod(MatF a, MatF b, MatF target);

/* @brief compute length of vector
 *
 * @param v: vector to compute length for
 */
float VEC_abs(MatF v);

// Solving Equations

/* @brief perform forward substitution for right hand side b (tri_lo * target =
 * b)
 *
 * @param tri_lo triangular lower matrix for forward substitution
 * @param b right hand side of equation
 * @param target where to put results
 */
void MF_forw_sub(MatF tri_lo, MatF b, MatF target);

/* @brief perform forward substitution for right hand side b (tri_up * target =
 * b)
 *
 * @param tri_up triangular upper matrix for backwards substitution
 * @param b right hand side of equation
 * @param target where to put results
 */
void MF_back_sub(const MatF tri_up, MatF b, MatF target);

/* @brief Perform LU-Decomposition for A using Gaussian elimination (inplace)
 *
 * @param A matrix to decompose
 * @return pointer to applied permutations
 */
RowPerm *MF_lu_decomp(MatF A);

/* @brief Calculate solution to LU * target = b
 *
 * @param LU: matrix containing LU decomposition entries
 * @param b: right hand side
 * @param perms: permutations to apply to b
 * @param target: where to put results
 */
void MF_lu_solv(const MatF LU, MatF b, RowPerm *perms, MatF target);

/* @brief Calculate solution to A target = b using Gauss-Elimination
 *        This converts A to a triangular upper matrix and changes b to match
 * the new system of equations. Invalidates A and b
 *
 * @param A: matrix
 * @param b: right hand side
 * @param target where to put results
 */
void MF_gauss_elim(MatF A, MatF b, MatF target);

/* @brief Turn A into a triangular upper matrix using gaussian elimination
 *
 * @param A: Matrix to transform
 */
void MF_make_tri_up(MatF A);

/* @brief Check if A is invertible
 *
 * @param A Matrix to check
 *
 * @return bool indicating whether a is invertible
 */
bool MF_invertible(MatF A);

/* @brief print string representation of A to stdout
 */
void MF_print(MatF A);

/* @brief print string representation of shape of A to stdout
 */
void MF_print_shape(MatF A);

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
SMatF SM_empty_from_pos_vals(long n_rows, long n_cols, long n_vals,
                             long *row_pos, long *col_pos, float *vals);

bool SM_has_loc(SMatF A, long row, long col);
long SM_idx(SMatF A, long row, long col);

void SM_set_or_panic(SMatF A, long row, long col, float val);
float SM_at(SMatF A, long row, long col);

long SM_col(SMatF A, long row, long col_idx);
long SM_col_or_panic(SMatF A, long row, long col_idx);

SMatF SM_prod_prepare(SMatF A, SMatF B);
void SM_prod(SMatF A, SMatF B, SMatF target);

void SM_print(SMatF A);
void SM_print_shape(SMatF A);
