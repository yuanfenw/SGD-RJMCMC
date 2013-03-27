#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int r, j, i, num_species, num_rxns = 8;
    double *coeffs, *s;
    double *F;
 
    num_species = mxGetN(prhs[0]);
    s = mxGetPr(prhs[0]);
    coeffs = mxGetPr(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(1,num_rxns, mxREAL);
    F=mxGetPr(plhs[0]);
    
    
    for (r=0; r<num_rxns; r++){
        F[r] = 1;
        for (j=0; j<num_species; j++){
            if (*(coeffs+r+j*num_rxns)> 0) {
                for (i=0; i< *(coeffs+r+j*num_rxns); i++) 
                 F[r] *= s[j]-i;
            }
            F[r] *= (s[j]<coeffs[r+j*num_rxns]? 0: 1);
        }
    }
//    Fr = &F;
//    
    return;
}