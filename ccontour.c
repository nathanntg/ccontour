/**
 * 
 */

#include "mex.h"
#include "consensus_contour.h"

static double getScalar(const mxArray *in, const char *err_id, const char *err_str) {
    /* check scalar */
    if (!mxIsDouble(in) || mxIsComplex(in) || mxGetN(in) != 1 || mxGetM(in) != 1 || mxGetNumberOfDimensions(in) != 2) {
        mexErrMsgIdAndTxt(err_id, err_str);
    }
    
    /* get the scalar input */
    return mxGetScalar(in);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    void *s, *t;
    size_t sn, sm, sl;
    double fs;
    bool is_single, is_double;
    
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
    if ((sn != 1 && sm != 1) || mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "First input must be a vector.");
    }
    sl = sn > 1 ? sn : sm;
    
    /* type */
    is_single = mxIsSingle(prhs[0]);
    is_double = mxIsDouble(prhs[0]);
    if (!is_single && !is_double) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "First input must be a single or double.");
    }
    
    /*  create a pointer to the input vector s */
    s = mxGetPr(prhs[0]);
    
    /* ARGUMENT 2: scalar, sample rate */
    
    /*  get the sample rate */
    fs = getScalar(prhs[1], "MATLAB:ccc:invalidInput", "Sample rate must be a scalar.");
    
    /* PREP WORK */
    
    if (is_double) {
        /* setup */
        CCCConfigD ccc_config = createCCCConfigD();
        cccConfigSetSampleRateD(ccc_config, fs);
        CCCSetupD ccc_setup = createCCCSetupD(ccc_config);
        destroyCCCConfigD(ccc_config);
        
        /* figure out contour size */
        const struct ConsensusContourSize dim = cccSizeD(ccc_setup, (unsigned long)sl);
        if (dim.rows == 0) {
            mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Input vector must be at longer than FFT window.");
        }
        
        /* OUTPUT 1: matrix, consensus contours */
        
        /*  set the output pointer to the output matrix */
        plhs[0] = mxCreateDoubleMatrix((size_t)dim.rows, (size_t)dim.cols, mxREAL);
        
        /*  create a C pointer to the output matrix */
        t = mxGetPr(plhs[0]);
        
        /*  call the C subroutine */
        cccSpectrogramD(ccc_setup, dim, s, t);
        
        /* clean up */
        destroyCCCSetupD(ccc_setup);
    }
    else {
        /* setup */
        CCCConfig ccc_config = createCCCConfig();
        cccConfigSetSampleRate(ccc_config, (float)fs);
        CCCSetup ccc_setup = createCCCSetup(ccc_config);
        destroyCCCConfig(ccc_config);
        
        /* figure out contour size */
        const struct ConsensusContourSize dim = cccSize(ccc_setup, (unsigned long)sl);
        if (dim.rows == 0) {
            mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Input vector must be at longer than FFT window.");
        }
        
        /* OUTPUT 1: matrix, consensus contours */
        
        /*  set the output pointer to the output matrix */
        plhs[0] = mxCreateNumericMatrix((size_t)dim.rows, (size_t)dim.cols, mxSINGLE_CLASS, mxREAL);
        
        /*  create a C pointer to the output matrix */
        t = mxGetPr(plhs[0]);
        
        /*  call the C subroutine */
        cccSpectrogram(ccc_setup, dim, s, t);
        
        /* clean up */
        destroyCCCSetup(ccc_setup);
    }
    
}
