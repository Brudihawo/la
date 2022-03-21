#include "stdio.h"

#include "math.h"
#include "sparse.h"
#include "stdlib.h"

#define FIELD_SIZE 100
#define V_X 2.0f
#define V_Y -2.0f
#define ALPHA 10.0f
#define DELTA_X 1.0f
#define DELTA_Y 1.0f

#define PI 3.14159265359

#define FREQ 2.0f

#define IDX(row, col) (row) * FIELD_SIZE + (col)

/*          y=0
 *  +-----------------+
 *  |                 |
 *  |                 |
 * x|                 |x
 * =|                 |=
 * 0|                 |1
 *  |                 |
 *  |                 |
 *  +-----------------+
 *          y=1
 */

int main(void) {
  // Values to insert into matrix
  const float a =
      -2.0f * ALPHA / (DELTA_X * DELTA_X) - 2.0f * ALPHA / (DELTA_Y * DELTA_Y);
  const float b_1 = ALPHA / (DELTA_X * DELTA_X) + V_X / (2.0f * DELTA_X);
  const float b_2 = ALPHA / (DELTA_Y * DELTA_Y) + V_Y / (2.0f * DELTA_Y);
  const float c_1 = ALPHA / (DELTA_X * DELTA_X) - V_X / (2.0f * DELTA_X);
  const float c_2 = ALPHA / (DELTA_Y * DELTA_Y) - V_Y / (2.0f * DELTA_Y);

  SMatF mat = SM_diag_empty((long[5]){-FIELD_SIZE, -1, 0, 1, FIELD_SIZE}, 5,
                            FIELD_SIZE * FIELD_SIZE);

  SMatF field = SM_vec_empty(FIELD_SIZE * FIELD_SIZE);

  // Boundary Conditions in system of equations
  for (long i = 0; i < FIELD_SIZE * FIELD_SIZE; i++) {
    SM_set_or_panic(mat, i, i, 1.0f);
  }

  long *bound_x = malloc(4 * (FIELD_SIZE - 1) * sizeof(long));
  long *bound_y = malloc(4 * (FIELD_SIZE - 1) * sizeof(long));

  // boundary coordinates
  // start at top left
  // clockwise
  for (long i = 0; i < FIELD_SIZE - 1; i++) {
    // top row
    bound_x[i] = i;
    bound_y[i] = 0;

    // right column
    bound_x[i + FIELD_SIZE - 1] = FIELD_SIZE - 1;
    bound_y[i + FIELD_SIZE - 1] = i;

    // bottom row right to left
    bound_x[i + 2 * (FIELD_SIZE - 1)] = FIELD_SIZE - i - 1;
    bound_y[i + 2 * (FIELD_SIZE - 1)] = FIELD_SIZE - 1;

    bound_x[i + 3 * (FIELD_SIZE - 1)] = 0;
    bound_y[i + 3 * (FIELD_SIZE - 1)] = FIELD_SIZE - i - 1;
  }

  // Boundary conditions in right hand side
  SMatF rhs = SM_empty_like(field);
  // set sides
  for (long i = 0; i < 4 * (FIELD_SIZE - 1); i++) {
    SM_set_or_panic(
        rhs, IDX(bound_x[i], bound_y[i]), 0,
        (1.0f + sinf(FREQ * 2.0f * PI / (4.0f * (FIELD_SIZE - 1)) * (float)i)) / 2.0);
  }

  for (long i = 1; i < FIELD_SIZE - 1; i++) {   // row
    for (long j = 1; j < FIELD_SIZE - 1; j++) { // column
      SM_set_or_panic(mat, IDX(i, j), IDX(i - 1, j), b_1);
      SM_set_or_panic(mat, IDX(i, j), IDX(i, j - 1), b_2);
      SM_set_or_panic(mat, IDX(i, j), IDX(i, j), a);
      SM_set_or_panic(mat, IDX(i, j), IDX(i, j + 1), c_2);
      SM_set_or_panic(mat, IDX(i, j), IDX(i + 1, j), c_1);
    }
  }

  SMatF result = SM_jacobi(mat, rhs, 0.00001f, 10000, 0.9999f);
  for (long i = 0; i < FIELD_SIZE; i++) {
    for (long j = 0; j < FIELD_SIZE; j++) {
      printf("%ld %ld %7.4f\n", i, j, SM_at(result, IDX(i, j), 0));
    }
  }

  free(bound_x);
  free(bound_y);
  SM_free(mat);
  SM_free(rhs);
  SM_free(field);
}
