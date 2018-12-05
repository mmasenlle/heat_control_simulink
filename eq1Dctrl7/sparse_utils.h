
#include "matrix.h"

struct sparse_t;

struct sparse_t *new_sparse(int m, int n);
double sparse_get(const struct sparse_t *sp, int i, int j);
void sparse_set(struct sparse_t *sp, int i, int j, double v);
void sparse_add(struct sparse_t *sp, int i, int j, double v);
void sparse2mxArray(const struct sparse_t *sp, mxArray *ma);
void sparse_free(struct sparse_t *sp);
