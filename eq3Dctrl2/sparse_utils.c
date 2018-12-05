
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparse_utils.h"

struct _sparse_nonzero
{
	int row;
	double value;
	struct _sparse_nonzero *next;
};

struct sparse_t
{
	int m,n;
	struct _sparse_nonzero *columns[0];
};

static struct _sparse_nonzero *new_nonzero(int j, double v)
{
    struct _sparse_nonzero *nn = malloc(sizeof(struct _sparse_nonzero));
    nn->row = j;
    nn->value = v;
    nn->next = 0;
    return nn;
}

struct sparse_t *new_sparse(int m, int n)
{
    struct sparse_t *sp = malloc(sizeof(struct sparse_t)+(n*sizeof(struct _sparse_nonzero *)));
    if (!sp) return sp;
    sp->m=m; sp->n=n;
    memset(sp->columns,0,(n*sizeof(struct _sparse_nonzero *)));
    return sp;
}

double sparse_get(const struct sparse_t *sp, int i, int j)
{
    if (i < sp->m && j < sp->n) {
        struct _sparse_nonzero *p = sp->columns[j];
        while (p && p->row <= i) {
            if (i == p->row) return p->value;
            p = p->next;
        }
    }
    return 0.0;
}

void sparse_set(struct sparse_t *sp, int i, int j, double v)
{
    if (i < sp->m && j < sp->n) {
        struct _sparse_nonzero **p = &sp->columns[j];
        while (*p && (*p)->row < i) p = &((*p)->next);
        if (!*p) {
            if (v) *p = new_nonzero(i, v);
        }
        else if ((*p)->row == i) {
            if (v) {
                (*p)->value = v;
            } else {
                free(*p);
                *p = (*p)->next;
            }
        }
        else if (v) { /* row > i */
            struct _sparse_nonzero *q = *p;
            *p = new_nonzero(i, v);
            (*p)->next = q;
        }
    }
}

void sparse_add(struct sparse_t *sp, int i, int j, double v)
{
    if (i < sp->m && j < sp->n) {
        struct _sparse_nonzero *p = sp->columns[j];
        while (p && p->row <= i) {
            if (i == p->row) {
                p->value += v;
                return;
            }
            p = p->next;
        }
        sparse_set(sp, i, j, v);
    }
}

void sparse_free(struct sparse_t *sp)
{
    int j;
    for(j=0; j < sp->n; j++) {
        struct _sparse_nonzero *p = sp->columns[j];
        while (p) {
            struct _sparse_nonzero *q = p;
            p = p->next;
            free(q);
        }
    }
    free(sp);
}

void sparse2mxArray(const struct sparse_t *sp, mxArray *ma)
{
    int j, k = 0;
    double *pr = mxGetPr(ma);
    mwIndex *jc = mxGetJc(ma);
    mwIndex *ir = mxGetIr(ma);
    mxSetM(ma, sp->m);
    mxSetN(ma, sp->n);
    for(j=0; j < sp->n; j++) {
        struct _sparse_nonzero *p = sp->columns[j];
        jc[j] = k;
        while (p) {
            ir[k] = p->row;
            pr[k] = p->value;
            k++;
            p = p->next;
        }
    }
    jc[j] = k;
}

/*
void sparse_print(const struct sparse_t *sp)
{
    int j, k = 0;
	printf("sparse_print -> (%d,%d)\n", sp->m, sp->n);
    for(j=0; j < sp->n; j++) {
        struct _sparse_nonzero *p = sp->columns[j];
        //jc[j] = k;
        while (p) {
          //  ir[k] = p->row;
            //pr[k] = p->value;
printf("(%d,%d) -> %f\n", p->row, j, p->value);
            k++;
            p = p->next;
        }
    }
   // jc[j] = k;
printf("sparse_print -> total: %d\n", k);
}
*/
