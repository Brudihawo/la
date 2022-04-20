#include "la.h"

#include "stdio.h"
#include "assert.h"
#include "stdbool.h"
#include "limits.h"
#include "memory.h"
#include "stdlib.h"

long binary_search(long *arr, long size, long val) {
  if (size < 1) return -1;
  if (size == 1) {
    if (val != arr[0])
      return -1;
    else
      return 0;
  }

  long end = size;
  long begin = 0;
  long cur_sz = end - begin;

  while (cur_sz > 0) {
    cur_sz = end - begin;
    const long mid_idx = begin + cur_sz / 2;
    const long mid_val = arr[mid_idx];

    if (val == mid_val)
      return mid_idx;

    const long ne = mid_idx - 1;
    const long nb = mid_idx + 1;
    const int cmp = val < mid_val;

    // hopefully this is faster?
    end = cmp * ne + (1 - cmp) * end;
    begin = cmp * begin + (1 - cmp) * nb;
    // if (val < mid_val) {
    //   end = mid_idx - 1;
    // } else {
    //   begin = mid_idx + 1;
    // }
  }

  return -1;
}

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

int long_gt(const void *va, const void *vb) {
  return *(long *)va > *(long *)vb;
}

int long_lt(const void *va, const void *vb) {
  return *(long *)va < *(long *)vb;
}

#define RADIX_MAX 100
void count_sort(long *arr, long *tmp, long size, long exp) {
  long counts[RADIX_MAX] = {0};
  for (long i = 0; i < size; ++i) {
    ++counts[(arr[i] / exp) % RADIX_MAX];
  }

  for (long i = 1; i < RADIX_MAX; ++i) {
    counts[i] += counts[i - 1];
  }

  for (long i = size - 1; i >= 0; --i) {
    tmp[counts[(arr[i] / exp) % RADIX_MAX] - 1] = arr[i];
    --counts[(arr[i] / exp) % RADIX_MAX];
  }

  memcpy(arr, tmp, size * sizeof(long));
}

long get_max(long *arr, long size) {
  long max = LONG_MIN;
  for (long i = 0; i < size; ++i) {
    if (arr[i] > max)
      max = arr[i];
  }
  return max;
}

void radix_sort(long *arr, long size, long *tmp_arr) {
  bool was_null = false;
  if (tmp_arr == NULL) {
    was_null = true;
    tmp_arr = malloc(size * sizeof(long));
  }

  long max = get_max(arr, size);
  for (long exp = 1; max / exp > 0; exp *= RADIX_MAX) {
    count_sort(arr, tmp_arr, size, exp);
    assert(exp < LONG_MAX / RADIX_MAX && "long overflow");
  }
  if (was_null)
    free(tmp_arr);
}

