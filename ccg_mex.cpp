/*
 *This is a mex file to compute cross correlations as in Kohn et al. 2005
 * HAS A BUG WHEN USED FOR MULTIPLE TRIALS AND IS NOT FASTER THAN MATLAB CORR!!!!
 */
#include "mex.h"
#include "matrix.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, mxArray *prhs[])
{
    if (nrhs < 2)
        mexErrMsgTxt("Wrong number of input arguments.");

    
#define SP1 prhs[0]
#define SP2 prhs[1]
    
    double* sp1;
    double* sp2;
    double *geo_mean_rates;
    double *ccg;
    double u1, u2;
    int shift, N, M, maxdelay,i,t,w;

    // Pointer to input arguments
    sp1 = mxGetPr(SP1);
    sp2 = mxGetPr(SP2);
    M = mxGetM(SP1);
    N = mxGetN(SP1);
    if (M != mxGetM(SP2))
        mexErrMsgTxt("Number of trials must be the same in the 2 binary spiketrains.");
    if (N != mxGetN(SP2))
        mexErrMsgTxt("Binary spiketrains must have the same length.");

    maxdelay = floor(N/2.0);
    if (nrhs > 2)
        maxdelay = mxGetScalar(prhs[2]);

    //Compute mean firing rates
    u1 = 0.0;
    u2 = 0.0;
    for (i = 0; i<M; i++)
        for (t = 0; t<N; t++) {
            u1 += sp1[t+i*N];
            u2 += sp2[t+i*N];
        }
    #ifdef DEBUG
        printf("Total number of spikes is: %lf and %lf\n",u1,u2);
    #endif
    u1 /= N*M;
    u2 /= N*M;
    #ifdef DEBUG
        mexPrintf("Working with a %d x %d matrix.\n",M,N);
        mexPrintf("Max delay is %d, %d\n",maxdelay,maxdelay*2+1);
    #endif
    // Compute the correlation
    plhs[0] = mxCreateDoubleMatrix(1, maxdelay*2+1, mxREAL);
    ccg = mxGetPr(plhs[0]);
    for (w = -maxdelay; w <= maxdelay; w++) {
        ccg[w + maxdelay] = 0;
        for (i = 0; i<M; i++) {
            for (t = 0; t<N; t++) {
                shift = t + w;
                if (shift >= 0 && shift <= N)  
                    // There is a bug with the shift!!!
                    ccg[w + maxdelay] += sp1[t + i*N] * sp2[shift + i*N];          
            }
        }
        ccg[w + maxdelay] /= M*(N - abs(w))*sqrt(u1*u2);
    }
    mxSetPr(plhs[0],ccg);
    //Output rates
    plhs[1] = mxCreateDoubleMatrix(1, 2, mxREAL);
    geo_mean_rates = mxGetPr(plhs[1]);
    geo_mean_rates[0] = u1;
    geo_mean_rates[1] = u2;
    mxSetPr(plhs[1],geo_mean_rates);
    
}