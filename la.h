#include "stdbool.h"

typedef struct {
  float* vals;
  long rows, cols;
} MatF;

/* @brief Create an empty MatF
 *
 * @param rows number of rows
 * @param cols number of cols
 *
 * @return MatF initalised to 0
 */
MatF MF_empty(long rows, long cols);

/* @brief Create an initialised MatF
 *
 * @param rows number of rows
 * @param cols number of cols
 *
 * @return MatF initialised with init_val
 */
MatF MF_with(long nx, long ny, float init_val);

/* @brief Shorthand for empty vector (Nx1 MatF)
 *
 * @param size length of vector
 *
 * @return Nx1 MatF initilised to 0
 */
#define VEC_EMPTY(size) MF_empty(nx, 1, init_val);

/* @brief Shorthand for initialised vector (Nx1 MatF)
 *
 * @param size length of vector
 * @param init_val value to initialise vector with
 *
 * @return Nx1 MatF initilised to init_val elementwise
 */
#define VEC_NEW(size, init_val) MF_with(nx, 1, init_val);

/* @brief Get value in MatF
 *
 * @param A matrix to query
 * @param row row
 * @param col column
 *
 * @return Value at A[row, col]
 */
float MF_at(MatF A, long row, long col);

/* @brief Get value in vector (Nx1 MatF)
 *
 * @param vec vector to query
 * @param idx row
 *
 * @return Value at vec[idx]
 */
#define VEC_AT(vec, idx) MF_at(vec, idx, 0);

/* @brief Get pointer to location in MatF
 *
 * @param A matrix to query
 * @param row row
 * @param col column
 *
 * @return pointer to A[row, col]
 */
float MF_ptr(MatF A, long row, long col);

/* @brief Get pointer to value in vector (Nx1 MatF)
 *
 * @param vec vector to query
 * @param idx row
 *
 * @return pointer to vec[idx]
 */
#define VEC_PTR(vec, idx) MF_ptr(vec, idx, 0);

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
MatF MF_subset_copy(MatF A, long start_row, long end_row, long start_col, long end_col);

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

// Solving Equations

/* @brief perform forward substitution for right hand side b (tri_lo * target = b)
 *
 * @param tri_lo triangular lower matrix for forward substitution
 * @param b right hand side of equation
 * @param target where to put results
 */
void MF_forw_sub(MatF tri_lo, MatF b, MatF target);

/* @brief perform forward substitution for right hand side b (tri_up * target = b)
 *
 * @param tri_up triangular upper matrix for backwards substitution
 * @param b right hand side of equation
 * @param target where to put results
 */
void MF_back_sub(MatF tri_up, MatF b, MatF target);

/* @brief Perform LU-Decomposition for A
 *
 * @param A matrix to decompose
 * @param target where to put results
 */
void MF_lu_decomp(MatF A, MatF target);

/* @brief Calculate solution to LU * target = b
 *
 * @param LU matrix containing LU decomposition entries
 * @param b right hand side
 * @param target where to put results
 */
void MF_lu_calc(MatF LU, MatF b, MatF target);

/* @brief Calculate solution to A target = b using Gauss-Elimination
 *
 * @param A matrix
 * @param b right hand side
 * @param target where to put results
 */
void MF_gauss_elim(MatF A, MatF b, MatF target);

/* @brief Check if A is invertible
 *
 * @param A Matrix to check
 * 
 * @return bool indicating whether a is invertible
 */
bool MF_invertible(MatF A);

