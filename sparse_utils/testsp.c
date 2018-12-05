
#include "mex.h"
#include "sparse_utils.h"


/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    mxArray *m = mxCreateSparse(5, 5, 15, mxREAL);
    struct sparse_t *sp = new_sparse(5, 5);
    mexPrintf("(3,3)-> %f\n", sparse_get(sp,3,3));
    sparse_set(sp, 1, 1, 1.1);
    sparse_set(sp, 2, 2, 2.2);
    sparse_set(sp, 3, 3, 3.3);
    sparse_set(sp, 0, 0, 0.05);
    sparse_set(sp, 4, 3, 4.3);
    sparse_set(sp, 3, 1, 3.1);
    sparse_set(sp, 7, 7, 7.17);
    sparse2mxArray(sp, m);
    plhs[0] = m;
mexPrintf("(3,3)-> %f\n", sparse_get(sp,3,3));
sparse_set(sp, 3, 3, 0.0);
mexPrintf("(3,3)-> %f\n", sparse_get(sp,3,3));
    sparse_free(sp);
    mexPrintf("Agur\n");
    return;
}
