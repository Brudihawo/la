#ifndef LA_H
#define LA_H

#include "stdbool.h"

#define REL_ERROR 0.00001
#define LA_PRINT_FMT "%6.2f"

typedef struct {
  long *order;
  long size;
  long n_swaps;
} RowPerm;

RowPerm RP_new(long rows);

void RP_add(RowPerm *perms, long i, long j);

void RP_free(RowPerm *perms);

bool almost_eq(float a, float b);

float rand_float();

int long_gt(const void *va, const void *vb);

int long_lt(const void *va, const void *vb);

#endif // LA_H
