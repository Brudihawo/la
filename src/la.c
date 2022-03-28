#include "la.h"

#include "stdbool.h"
#include "stdlib.h"

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

void RP_free(RowPerm *perms) {
  free(perms->order);
}

float rand_float() { return (float)rand() / (float)RAND_MAX; }
