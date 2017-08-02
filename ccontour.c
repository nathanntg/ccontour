/**
 * 
 */

#include <string.h>
#include "mex.h"
#include "consensus_contour.h"

#define MX_TEST(TEST, ERR_ID, ERR_STR) if (!(TEST)) { mexErrMsgIdAndTxt(ERR_ID, ERR_STR); }

static bool testScalar(const mxArray *in) {
    // type
    if (!mxIsDouble(in) && !mxIsSingle(in)) {
        return false;
    }
    
    // real
    if (mxIsComplex(in)) {
        return false;
    }
    
    // dimensions
    if (mxGetNumberOfDimensions(in) != 2 || mxGetM(in) != 1 || mxGetN(in) != 1) {
        return false;
    }
    
    return true;
}

static double getScalar(const mxArray *in, const char *err_id, const char *err_str) {
    /* check scalar */
    if (!testScalar(in)) {
        mexErrMsgIdAndTxt(err_id, err_str);
    }
    
    /* get the scalar input */
    return mxGetScalar(in);
}

static bool testVector(const mxArray *in) {
    // type
    if (!mxIsDouble(in) && !mxIsSingle(in)) {
        return false;
    }
    
    // real
    if (mxIsComplex(in)) {
        return false;
    }
    
    // dimensions
    if (mxGetNumberOfDimensions(in) != 2) {
        return false;
    }
    
    // make sure there is at least one singleton dimension
    if (mxGetM(in) != 1 && mxGetN(in) != 1) {
        return false;
    }
    
    // make sure there is at least one non-zero dimension
    if (mxGetM(in) < 2 && mxGetN(in) < 2) {
        return false;
    }
    
    return true;
}

static bool isPowerOfTwo(unsigned long x) {
    return ((x != 0) && ((x & (~x + 1)) == x));
}

static mxArray *processSingle(int nrhs, const mxArray *prhs[]) {
    mxArray *ret;
    void *s, *t;
    size_t sn, sm, sl;
    double fs;
    
    /* ARGUMENT 1: vector, signal */
    if (!testVector(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Signal must be a single or double real vector.");
    }
    
    /* get dimensions */
    sn = mxGetN(prhs[0]);
    sm = mxGetM(prhs[0]);
    sl = sn > 1 ? sn : sm;
    
    /*  create a pointer to the input vector s */
    s = mxGetPr(prhs[0]);
    
    /* ARGUMENT 2: scalar, sample rate */
    fs = getScalar(prhs[1], "MATLAB:ccc:invalidInput", "Sample rate must be a scalar.");
    
    /* ARGUMENT 3+: configuration */
    unsigned long fft_length = 1024;
    unsigned long fft_overlap = 1005;
    bool pow_weight = true;
    unsigned long num_timescales = 0;
    float *timescales = NULL;
    unsigned long num_angles = 0;
    float *angles = NULL;
    
    MX_TEST((nrhs - 2) % 2 == 0, "MATLAB:ccc:invalidNumInputs", "Expected an even number of input configurations.");
    for (int i = 2; i < nrhs; i += 2) {
        // parameter name
        char *nm = mxArrayToString(prhs[i]);
        
        // handle
        if (0 == strcmp(nm, "fft_length")) {
            if (!testScalar(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT length must be a scalar.");
            }
            
            double v = mxGetScalar(prhs[i + 1]);
            unsigned long tmp = (unsigned long)v;
            if (tmp < 8 || !isPowerOfTwo(tmp)) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT length must be a power of two greater than or equal to 8.");
            }
            
            fft_length = tmp;
        }
        else if (0 == strcmp(nm, "fft_overlap")) {
            if (!testScalar(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT overlap must be a scalar.");
            }
            
            double v = mxGetScalar(prhs[i + 1]);
            fft_overlap = (unsigned long)v;
        }
        else if (0 == strcmp(nm, "pow_weight")) {
            if (!mxIsLogicalScalar(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Weight by power must be a logical value.");
            }
            
            pow_weight = mxIsLogicalScalarTrue(prhs[i + 1]);
        }
        else if (0 == strcmp(nm, "timescales")) {
            if (!testVector(prhs[i + 1]) || !mxIsSingle(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Timescales must be a vector of the same type as the signal.");
            }
            
            num_timescales = mxGetM(prhs[i + 1]) * mxGetN(prhs[i + 1]);
            timescales = (float *)mxGetPr(prhs[i + 1]);
            
            if (num_timescales > 16) {
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "A maximum of 16 timescales are supported.");
            }
        }
        else if (0 == strcmp(nm, "angles")) {
            if (!testVector(prhs[i + 1]) ||! mxIsSingle(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Angles must be a vector of the same type as the signal.");
            }
            
            num_angles = mxGetM(prhs[i + 1]) * mxGetN(prhs[i + 1]);
            angles = (float *)mxGetPr(prhs[i + 1]);
            
            if (num_angles > 16) {
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "A maximum of 16 angles are supported.");
            }
        }
        else {
            mxFree(nm);
            mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Unknown parameter.");
        }
        
        // free
        mxFree(nm);
    }
    
    // validate
    if (fft_overlap >= fft_length) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT overlap must be less than the FFT length.");
    }
    
    /* PREPARE */
    /* configuration */
    CCCConfig ccc_config = createCCCConfig();
    cccConfigSetSampleRate(ccc_config, (float)fs);
    cccConfigSetFFTLength(ccc_config, fft_length);
    cccConfigSetFFTOverlap(ccc_config, fft_overlap);
    cccConfigSetWeightByPower(ccc_config, pow_weight);
    if (0 < num_timescales) {
        cccConfigSetTimescales(ccc_config, num_timescales, timescales);
    }
    if (0 < num_angles) {
        cccConfigSetAngles(ccc_config, num_angles, angles);
    }
    
    /* setup */
    CCCSetup ccc_setup = createCCCSetup(ccc_config);
    destroyCCCConfig(ccc_config);
    
    /* figure out contour size */
    const struct ConsensusContourSize dim = cccSize(ccc_setup, (unsigned long)sl);
    if (dim.rows == 0) {
        destroyCCCSetup(ccc_setup);
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Signal vector must be at longer than FFT window.");
    }
    
    /* OUTPUT 1: matrix, consensus contours */
    
    /*  set the output pointer to the output matrix */
    ret = mxCreateNumericMatrix((size_t)dim.rows, (size_t)dim.cols, mxSINGLE_CLASS, mxREAL);
    
    /*  create a C pointer to the output matrix */
    t = mxGetPr(ret);
    
    /*  call the C subroutine */
    cccSpectrogram(ccc_setup, dim, s, t);
    
    /* clean up */
    destroyCCCSetup(ccc_setup);
    
    return ret;
}

static mxArray *processDouble(int nrhs, const mxArray *prhs[]) {
    mxArray *ret;
    void *s, *t;
    size_t sn, sm, sl;
    double fs;
    
    /* ARGUMENT 1: vector, signal */
    if (!testVector(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Signal must be a single or double real vector.");
    }
    
    /* get dimensions */
    sn = mxGetN(prhs[0]);
    sm = mxGetM(prhs[0]);
    sl = sn > 1 ? sn : sm;
    
    /*  create a pointer to the input vector s */
    s = mxGetPr(prhs[0]);
    
    /* ARGUMENT 2: scalar, sample rate */
    fs = getScalar(prhs[1], "MATLAB:ccc:invalidInput", "Sample rate must be a scalar.");
    
    /* ARGUMENT 3+: configuration */
    unsigned long fft_length = 1024;
    unsigned long fft_overlap = 1005;
    bool pow_weight = true;
    unsigned long num_timescales = 0;
    double *timescales = NULL;
    unsigned long num_angles = 0;
    double *angles = NULL;
    
    MX_TEST((nrhs - 2) % 2 == 0, "MATLAB:ccc:invalidNumInputs", "Expected an even number of input configurations.");
    for (int i = 2; i < nrhs; i += 2) {
        // parameter name
        char *nm = mxArrayToString(prhs[i]);
        
        // handle
        if (0 == strcmp(nm, "fft_length")) {
            if (!testScalar(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT length must be a scalar.");
            }
            
            double v = mxGetScalar(prhs[i + 1]);
            unsigned long tmp = (unsigned long)v;
            if (tmp < 8 || !isPowerOfTwo(tmp)) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT length must be a power of two greater than or equal to 8.");
            }
            
            fft_length = tmp;
        }
        else if (0 == strcmp(nm, "fft_overlap")) {
            if (!testScalar(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT overlap must be a scalar.");
            }
            
            double v = mxGetScalar(prhs[i + 1]);
            fft_overlap = (unsigned long)v;
        }
        else if (0 == strcmp(nm, "pow_weight")) {
            if (!mxIsLogicalScalar(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Weight by power must be a logical value.");
            }
            
            pow_weight = mxIsLogicalScalarTrue(prhs[i + 1]);
        }
        else if (0 == strcmp(nm, "timescales")) {
            if (!testVector(prhs[i + 1]) || !mxIsDouble(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Timescales must be a vector of the same type as the signal.");
            }
            
            num_timescales = mxGetM(prhs[i + 1]) * mxGetN(prhs[i + 1]);
            timescales = mxGetPr(prhs[i + 1]);
            
            if (num_timescales > 16) {
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "A maximum of 16 timescales are supported.");
            }
        }
        else if (0 == strcmp(nm, "angles")) {
            if (!testVector(prhs[i + 1]) ||! mxIsDouble(prhs[i + 1])) {
                mxFree(nm);
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Angles must be a vector of the same type as the signal.");
            }
            
            num_angles = mxGetM(prhs[i + 1]) * mxGetN(prhs[i + 1]);
            angles = mxGetPr(prhs[i + 1]);
            
            if (num_angles > 16) {
                mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "A maximum of 16 angles are supported.");
            }
        }
        else {
            mxFree(nm);
            mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Unknown parameter.");
        }
        
        // free
        mxFree(nm);
    }
    
    // validate
    if (fft_overlap >= fft_length) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "FFT overlap must be less than the FFT length.");
    }
    
    /* PREPARE */
    /* configuration */
    CCCConfigD ccc_config = createCCCConfigD();
    cccConfigSetSampleRateD(ccc_config, fs);
    cccConfigSetFFTLengthD(ccc_config, fft_length);
    cccConfigSetFFTOverlapD(ccc_config, fft_overlap);
    cccConfigSetWeightByPowerD(ccc_config, pow_weight);
    if (0 < num_timescales) {
        cccConfigSetTimescalesD(ccc_config, num_timescales, timescales);
    }
    if (0 < num_angles) {
        cccConfigSetAnglesD(ccc_config, num_angles, angles);
    }
    
    /* setup */
    CCCSetupD ccc_setup = createCCCSetupD(ccc_config);
    destroyCCCConfigD(ccc_config);
    
    /* figure out contour size */
    const struct ConsensusContourSize dim = cccSizeD(ccc_setup, (unsigned long)sl);
    if (dim.rows == 0) {
        destroyCCCSetupD(ccc_setup);
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Signal vector must be at longer than FFT window.");
    }
    
    /* OUTPUT 1: matrix, consensus contours */
    
    /*  set the output pointer to the output matrix */
    ret = mxCreateDoubleMatrix((size_t)dim.rows, (size_t)dim.cols, mxREAL);
    
    /*  create a C pointer to the output matrix */
    t = mxGetPr(ret);
    
    /*  call the C subroutine */
    cccSpectrogramD(ccc_setup, dim, s, t);
    
    /* clean up */
    destroyCCCSetupD(ccc_setup);
    
    return ret;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*  check for proper number of arguments */
    MX_TEST(nrhs >= 2, "MATLAB:ccc:invalidNumInputs", "Signal and sampling rate are required (two inputs).");
    MX_TEST(nlhs == 1, "MATLAB:ccc:invalidNumOutputs", "Spectrogram output is required (one output).");
    
    if (mxIsSingle(prhs[0])) {
        plhs[0] = processSingle(nrhs, prhs);
    }
    else if (mxIsDouble(prhs[0])) {
        plhs[0] = processDouble(nrhs, prhs);
    }
    else {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Signal must be a single or double.");
    }
}
