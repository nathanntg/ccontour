/**
 * 
 */

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <Accelerate/Accelerate.h>

#define N_ANGLES 8
#define N_TIMESCALES 9

typedef vDSP_Stride t_stride;
typedef vDSP_Length t_len;

static t_len fft_length = 1024; /* must be power of two */
static t_len fft_overlap = 1005;

static bool pow_weight = true;

static double angles[N_ANGLES]; /* in radians */
static double timescales[N_TIMESCALES]; /* in ms */

struct CCTimescaleAngle
{
    t_len fft_length_half;
    
    double timescale;
    double angle;
    
    double mu_real;
    double mu_imag;
    
    double grad_x;
    double grad_y;
    
    // TODO: could potentially be int8_t
    double *sign_last;
    double *sign_cur;
};

static void setupCCTimescaleAngle(struct CCTimescaleAngle *ccta, const t_len fft_length_half, const double timescale, const double angle) {
    ccta->fft_length_half = fft_length_half;
    ccta->timescale = timescale;
    ccta->angle = angle;
    ccta->mu_real = cos(angle);
    ccta->mu_imag = sin(angle);
    ccta->grad_x = 0 - cos(angle + M_PI_2) / 2;
    ccta->grad_y = sin(angle + M_PI_2) / 2;
    ccta->sign_last = calloc(fft_length_half * 2, sizeof(double));
    ccta->sign_cur = calloc(fft_length_half * 2, sizeof(double));
}

static void destroyCCTimescaleAngle(struct CCTimescaleAngle *ccta) {
    free(ccta->sign_last);
    free(ccta->sign_cur);
}

static void ingestCCTimescaleAngle(struct CCTimescaleAngle *ccta, const double *real, const double *imag, double *out) {
    t_len i;
    double f, gx, gy;
    double *last;
    double *cur;
    
    // swap cur and last
    last = ccta->sign_cur;
    cur = ccta->sign_last;
    ccta->sign_last = last;
    ccta->sign_cur = cur;
    
    // iterate
    for (i = 0; i < ccta->fft_length_half; ++i) {
        // calculate the imaginary portion of the angle (FOIL)
        f = real[i] * ccta->mu_imag + imag[i] * ccta->mu_real;
        if (f > 0) {
            cur[i] = 1;
        }
        else if (f < 0) {
            cur[i] = -1;
        }
        else {
            cur[i] = 0;
        }
        
        // adjust by angle
        // TODO: optimize via switch case and
        gx = (cur[i] - last[i]) * ccta->grad_x;
        gy = (i > 0 ? (cur[i] - cur[i - 1]) * ccta->grad_y : 0);
        
        out[i] = (gx + gy > 0.001 ? 1 : 0);
    }
}

static void fillFftWindow(double *fft_window, const double signal_fs) {
	int i;
	t_len j;
	t_len index = 0;
	double timescale_samples;
	double w;
	double twin, twin_offset = 0.5 + (0 - (double)fft_length) / 2;
	const t_len index_offset = N_TIMESCALES * fft_length;
	
	for (i = 0; i < N_TIMESCALES; ++i) {
		timescale_samples = signal_fs * timescales[i] / 1000;
		
		/* TODO: could be heabily optimized using vvexp and such */
		for (j = 0; j < fft_length; ++j) {
			twin = twin_offset + (double)j;
			w = exp(-pow(twin / timescale_samples, 2));
			
			// store exponential
			fft_window[index + j] = w;
			
			// store derivative
			fft_window[index_offset + index + j] = -2 * w * twin / pow(timescale_samples, 2);
		}
		
		index += fft_length;
	}
}

static void windowSamples(const double *signal, const double *window, double *signal_windowed) {
	t_len i = 0;
	t_len signals = 2 * N_TIMESCALES;

	for (i = 0; i < signals; ++i) {
		vDSP_vmulD(signal, 1, window + i * fft_length, 1, signal_windowed + i * fft_length, 1, fft_length);
	}
}

struct ConsensusContourSize
{
    t_len rows;
    t_len cols;
    t_len bytes;
};

static struct ConsensusContourSize cccSize(const t_len signal_len) {
    const t_len ccn = fft_length / 2;
    const t_len ccm = 1 + (signal_len - fft_length) / (fft_length - fft_overlap);
    struct ConsensusContourSize ret;
    
    ret.rows = ccn;
    ret.cols = ccm;
    ret.bytes = sizeof(double) * ccn * ccm;
    
    return ret;
}

static void addConsensusToContours(const t_len fft_length_half, const double *consensus, double *consensus_contours) {
    for (t_len i = 0; i < fft_length_half; ++i) {
        if (consensus[i] > 1.5) {
            consensus_contours[i] += 1;
        }
    }
}

static void ccc(const double *signal, const t_len signal_len, const double signal_fs, double *consensus_contours) {
    FFTSetupD fft_setup;
    double *fft_window, *signal_windowed, *power;
    double *consensus, *consensus_cur, *consensus_pow;
    double power_scale = 2.0, two_pi = 2 * M_PI;
    DSPDoubleSplitComplex fft_temporary, fft_output;
    DSPDoubleSplitComplex p_exp, p_der, p_ratio; // pointers into the fft_output for convience
    struct CCTimescaleAngle ccta[N_TIMESCALES * N_ANGLES];
    t_len i, j, k;
    
    /* output dimensions */
    const struct ConsensusContourSize dim = cccSize(signal_len);
    
    /* PREP WORK: fill timescales and angles */
    for (i = 0; i < N_ANGLES; ++i) {
    	angles[i] = M_PI * (double)(2 + i) / N_ANGLES;
    }
    for (i = 0; i < N_TIMESCALES; ++i) {
    	timescales[i] = 0.5 + 0.2 * (double)i;
    }
    
    /* clear output */
    vDSP_vclrD(consensus_contours, 1, dim.cols * dim.rows);
    
    /* STEP 1: spectrogram, build windows and fft */
    const t_len fft_length_half = fft_length / 2;
    const t_len fft_size = (t_len)ceil(log2((double)fft_length));
    assert((1 << fft_size) == fft_length); // should be power of 2
    
    fft_setup = vDSP_create_fftsetupD(fft_size, kFFTRadix2);
    
    /* allocate window */
    fft_window = malloc(sizeof(double) * fft_length * N_TIMESCALES * 2);
    signal_windowed = malloc(sizeof(double) * fft_length * N_TIMESCALES * 2);
    
    /* allocate temporary and output */
    fft_temporary.realp = malloc(sizeof(double) * fft_length_half);
    fft_temporary.imagp = malloc(sizeof(double) * fft_length_half);
    fft_output.realp = malloc(sizeof(double) * fft_length_half * N_TIMESCALES * 2);
    fft_output.imagp = malloc(sizeof(double) * fft_length_half * N_TIMESCALES * 2);
    
    /* make pointers */
    p_exp = fft_output;
    p_der.realp = fft_output.realp + fft_length_half * N_TIMESCALES;
    p_der.imagp = fft_output.imagp + fft_length_half * N_TIMESCALES;
    p_ratio = p_der;
    
    /* power array */
    power = malloc(sizeof(double) * fft_length_half * N_TIMESCALES);
    
    /* conensus array */
    /* might be able to use singles? */
    consensus = malloc(sizeof(double) * fft_length_half * N_TIMESCALES * N_ANGLES);
    consensus_cur = malloc(sizeof(double) * fft_length_half);
    if (pow_weight) {
        consensus_pow = malloc(sizeof(double) * fft_length_half);
    }
    
    /* fill window */
    fillFftWindow(fft_window, signal_fs);
    
    /* setup timescale angle structures */
    for (i = 0; i < N_TIMESCALES; ++i) {
        for (j = 0; j < N_ANGLES; ++j) {
            setupCCTimescaleAngle(&ccta[i * N_ANGLES + j], fft_length_half, timescales[i], angles[j]);
        }
    }
    
    clock_t begin, end;
    
    /* STEP 2: spectrogram columns */
    for (i = 0; i < dim.cols; ++i) {
        // begin = clock();
        
    	// window
    	windowSamples(signal + i * (fft_length - fft_overlap), fft_window, signal_windowed);
    	
    	// pack samples
    	vDSP_ctozD((DSPDoubleComplex *)signal_windowed, 2, &fft_output, 1, fft_length_half * N_TIMESCALES * 2);
    	
    	// calculate
    	vDSP_fftm_zriptD(fft_setup, &fft_output, 1, (t_stride)fft_length_half, &fft_temporary, fft_size, N_TIMESCALES * 2, FFT_FORWARD);
    	
    	// calculate power
    	vDSP_zvabsD(&fft_output, 1, power, 1, fft_length_half * N_TIMESCALES);
    	vDSP_vsdivD(power, 1, &power_scale, power, 1, fft_length_half * N_TIMESCALES);
    	
    	// calculate ratio
    	vDSP_zvdivD(&p_exp, 1, &p_der, 1, &p_ratio, 1, fft_length_half * N_TIMESCALES);
        
        // scale / (2 * pi)
        // TODO: can probably remove, won't affect sign?
        vDSP_vsdivD(p_ratio.realp, 1, &two_pi, p_ratio.realp, 1, fft_length_half * N_TIMESCALES);
        vDSP_vsdivD(p_ratio.imagp, 1, &two_pi, p_ratio.imagp, 1, fft_length_half * N_TIMESCALES);
        
        // calculate contours for each timescale and angle
        for (j = 0; j < N_TIMESCALES; ++j) {
            for (k = 0; k < N_ANGLES; ++k) {
                ingestCCTimescaleAngle(&ccta[j * N_ANGLES + k], p_ratio.realp + j * fft_length_half, p_ratio.imagp + j * fft_length_half, consensus + (j * N_ANGLES + k) * fft_length_half);
            }
        }
        
        // look for consensus
        for (j = 0; j < (N_TIMESCALES - 1); ++j) {
            if (pow_weight) {
                vDSP_vclrD(consensus_pow, 1, fft_length_half);
            }
            
            for (k = 0; k < N_ANGLES; ++k) {
                // start with consensus / angle
                memcpy(consensus_cur, consensus + (j * N_ANGLES + k) * fft_length_half, sizeof(double) * fft_length_half);
                // add next sigma, same angle
                vDSP_vaddD(consensus_cur, 1, consensus + ((j + 1) * N_ANGLES + k) * fft_length_half, 1, consensus_cur, 1, fft_length_half);
                // add same sigma, previous angle (wrap around)
                vDSP_vaddD(consensus_cur, 1, consensus + (j * N_ANGLES + (0 == k ? N_ANGLES - 1 : k - 1)) * fft_length_half, 1, consensus_cur, 1, fft_length_half);
                
                // add to contours
                if (pow_weight) {
                    addConsensusToContours(fft_length_half, consensus_cur, consensus_pow);
                }
                else {
                    addConsensusToContours(fft_length_half, consensus_cur, consensus_contours + i * fft_length_half);
                }
            }
            
            if (pow_weight) {
                // scale by power
                vDSP_vmulD(consensus_pow, 1, power + j * fft_length_half, 1, consensus_pow, 1, fft_length_half);
                
                // add to output
                vDSP_vaddD(consensus_contours + i * fft_length_half, 1, consensus_pow, 1, consensus_contours + i * fft_length_half, 1, fft_length_half);
            }
        }
        
        // end = clock();
        // printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    }
    
    /* CLEAN UP */
    free(fft_window);
    free(signal_windowed);
    free(fft_temporary.realp);
    free(fft_temporary.imagp);
    free(fft_output.realp);
    free(fft_output.imagp);
    free(power);
    free(consensus);
    free(consensus_cur);
    if (pow_weight) {
        free(consensus_pow);
    }
    vDSP_destroy_fftsetupD(fft_setup);
    
    /* destroy timescale angle structures */
    for (i = 0; i < N_TIMESCALES; ++i) {
        for (j = 0; j < N_ANGLES; ++j) {
            destroyCCTimescaleAngle(&ccta[i * N_ANGLES + j]);
        }
    }
}

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
    if (nrhs != 2 && nrhs != 3) {
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
    if (sl < fft_length) {
        mexErrMsgIdAndTxt("MATLAB:ccc:invalidInput", "Input vector must be at longer than FFT window.");
    }
    
    /*  create a pointer to the input vector s */
    s = mxGetPr(prhs[0]);
    
    /* ARGUMENT 2: scalar, sample rate */
    
    /*  get the sample rate */
    fs = getScalar(prhs[1], "MATLAB:ccc:invalidInput", "Sample rate must be a scalar.");
    
    /* OUTPUT 1: matrix, consensus contours */
    const struct ConsensusContourSize dim = cccSize((t_len)sl);
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix((size_t)dim.rows, (size_t)dim.cols, mxREAL);
    
    /*  create a C pointer to the output matrix */
    t = mxGetPr(plhs[0]);
    
    /*  call the C subroutine */
    ccc(s, (t_len)sl, fs, t);
}
