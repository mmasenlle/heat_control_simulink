
#include "mex.h"
#include "sparse_utils.h"

struct init_data_t
{
    int nel,*n;
    struct sparse_t *spC,*spK;
    double *Kek,*Cel,*Kev;
    double *vT,*vk,*f,*cpe1,*cpe2;
    double mk,v;
};
static void init_data(struct init_data_t *data, mxArray *plhs[], const mxArray *prhs[])
{
    double Tm,Teps;
    int nn = (int)mxGetScalar(mxGetField(prhs[0], 0, "n"));
    data->nel = nn-1;
    plhs[0] = mxCreateSparse(nn, nn, 15*nn, mxREAL);
    plhs[1] = mxCreateSparse(nn, nn, 15*nn, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nn, 1, mxREAL);
    data->f = mxGetPr(plhs[2]);
    data->vT = mxGetPr(prhs[2]);
    data->cpe1 = mxGetPr(prhs[4]);
	data->cpe2 = mxGetPr(prhs[5]);
    data->vk = mxGetPr(prhs[6]);
    data->v = mxGetScalar(prhs[7]);
    data->spC = new_sparse(nn, nn);
    data->spK = new_sparse(nn, nn);
    data->n = mxGetData(mxGetField(prhs[0], 0, "cnodes"));
    data->Kek = mxGetPr(mxGetField(prhs[0], 0, "Kek"));
    data->Kev = mxGetPr(mxGetField(prhs[0], 0, "Kev"));
    data->Cel = mxGetPr(mxGetField(prhs[0], 0, "Cel"));
}

static void update(struct init_data_t *data)
{
	data->n += 2;
    data->Kek += 4;
    data->Cel += 1;
    data->Kev += 4;
}
static int *calc(struct init_data_t *data)
{
    data->mk = (data->vk[data->n[0]]+data->vk[data->n[1]])*0.5;
    return data->n;
}

/*
 * [C K f] = fem1d_runc(data,params,T,t,cp,k,v)
 */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    int j,k,p,*n;
    struct init_data_t d;
    init_data(&d, plhs, prhs);
    for (p=0;p<d.nel;p++,update(&d)) {
        n = calc(&d);
/* matriz K*/
        for (j=0;j<2;j++) {
            for (k=0;k<2;k++) {
                int idx = j*2+k;
                sparse_set(d.spK,n[k],n[j], sparse_get(d.spK,n[k],n[j]) + d.Kek[idx]*d.mk
                       + d.Kev[idx]*d.cpe2[n[j]]*d.v);
            }
/* matriz C*/sparse_set(d.spC,n[j],n[j], sparse_get(d.spC,n[j],n[j]) + *d.Cel*d.cpe1[n[j]]);
        }
    }
    sparse2mxArray(d.spC, plhs[0]);
    sparse2mxArray(d.spK, plhs[1]);
    sparse_free(d.spC);
    sparse_free(d.spK);
    return;
}
