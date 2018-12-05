
//#include "mex.h"
#include "sparse_utils.h"


/* gateway function 
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
*/
int main(int argc, char *argv[])
{
/*    mxArray *m = mxCreateSparse(5, 5, 15, mxREAL);*/
    struct sparse_t *sp = new_sparse(5, 5);
    sparse_set(sp, 1, 1, 1.1);
    sparse_set(sp, 2, 2, 2.2);
    sparse_set(sp, 3, 3, 3.3);
    sparse_set(sp, 0, 0, 5.5);
    sparse_set(sp, 2, 0, 5.5);
    sparse_set(sp, 2, 0, 5.5);
    sparse_set(sp, 1, 0, 5.5);
    sparse_set(sp, 4, 3, 4.3);
    sparse_set(sp, 3, 1, 3.1879787);
    sparse_set(sp, 3, 1, 3.1879787);
    sparse_set(sp, 3, 1, 3.1879787);
    sparse_set(sp, 4, 4, 0.4);
	sparse_print(sp);
    sparse_set(sp, 2, 2, 0.0);
    sparse_set(sp, 4, 3, 0.0);
    sparse_set(sp, 0, 0, 0.0);
    sparse_set(sp, 0, 0, 0.0);
    sparse_set(sp, 0, 0, 0.0);
/*    sparse2mxArray(sp, m);
    plhs[0] = m;
  */ 
	sparse_print(sp);
printf("testsp1-> (3,1)-> %f\n", sparse_get(sp, 3, 1));
printf("testsp1-> (3,3)-> %f\n", sparse_get(sp, 3, 3));
printf("testsp1-> (4,3)-> %f\n", sparse_get(sp, 4, 3));
printf("testsp1-> (2,0)-> %f\n", sparse_get(sp, 2, 0));
    sparse_free(sp);
printf("Agur\n");
 return;
}
