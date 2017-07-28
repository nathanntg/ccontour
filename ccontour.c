/**
 * 
 */

#include "mex.h"
#include "consensus_contour.h"

static double getScalar(const mxArray *in, const char *err_id, const char *err_str) {
    /* check scalar */
    if (!mxIsDouble(in) || mxIsComplex(in) || mxGetN(in) * mxGetM(in) != 1) {
        mexErrMsgIdAndTxt(err_id, err_str);
    }
    
    /* get the scalar input */
    return mxGetScalar(in);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *s, *t;
    size_t sn, sm, sl;
    double fs;
    
    /* ARGUMENT 1: vector, signal */
    
    /*  check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidNumInputs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:ccc:invalidNumOutputs", "One output required.");
    }
    
    /*  get the dimensions of the matrix input s */
    sn = mxGetN(prhs[0]);
    sm = mxGetM(prhs[0]);
    if (sn != 1 && sm != 1) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "First input must be a vector.");
    }
    sl = sn > 1 ? sn : sm;
    
    /*  create a pointer to the input vector s */
    s = mxGetPr(prhs[0]);
    
    /* ARGUMENT 2: scalar, sample rate */
    
    /*  get the sample rate */
    fs = getScalar(prhs[1], "MATLAB:ccc:invalidInput", "Sample rate must be a scalar.");
    
    /* PREP WORK */
    
    /* setup */
    CCCSetup ccc_setup = createCCCSetup(1024, 1005, fs, true);
    
    /* figure out contour size */
    const struct ConsensusContourSize dim = cccSize(ccc_setup, (unsigned long)sl);
    if (dim.rows == 0) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Input vector must be at longer than FFT window.");
    }
    
    /* OUTPUT 1: matrix, consensus contours */
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix((size_t)dim.rows, (size_t)dim.cols, mxREAL);
    
    /*  create a C pointer to the output matrix */
    t = mxGetPr(plhs[0]);
    
    /*  call the C subroutine */
    ccc(ccc_setup, dim, s, t);
    
    /* clean up */
    destroyCCCSetup(ccc_setup);
}
