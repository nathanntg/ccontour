/**
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "consensus_contour.h"

#define N_ANGLES 8
#define N_TIMESCALES 9

typedef vDSP_Stride t_stride;
typedef vDSP_Length t_len;

#if PRECISION == 2
#define vDSP(FN) vDSP_ ## FN ## D
#define CMATH(FN) FN
#define NUM(N) N
typedef DSPDoubleSplitComplex RealSplitComplex;
typedef DSPDoubleComplex RealComplex;
typedef FFTSetupD RealFFTSetup;
#else
#define vDSP(FN) vDSP_ ## FN
#define NUM(N) N ## f
#define CMATH(FN) FN ## f
typedef DSPSplitComplex RealSplitComplex;
typedef DSPComplex RealComplex;
typedef FFTSetup RealFFTSetup;
#endif

static t_real s_zero = 0.0;
static t_real s_one = 1.0;
static t_real s_power_scale = 2.0;

struct CCTimescaleAngle
{
    t_len fft_length_half;
    
    t_real timescale;
    t_real angle;
    
    t_real mu_real;
    t_real mu_imag;
    
    t_real grad_x;
    t_real grad_y;
    
    // TODO: could potentially be int8_t
    t_real *sign_last;
    t_real *sign_cur;
};

static void setupCCTimescaleAngle(struct CCTimescaleAngle *ccta, const t_len fft_length_half, const t_real timescale, const t_real angle) {
    ccta->fft_length_half = fft_length_half;
    ccta->timescale = timescale;
    ccta->angle = angle;
    ccta->mu_real = cos(angle);
    ccta->mu_imag = sin(angle);
    ccta->grad_x = 0 - cos(angle + M_PI_2) / NUM(2.0);
    ccta->grad_y = sin(angle + M_PI_2) / NUM(2.0);
    ccta->sign_last = calloc(fft_length_half * 2, sizeof(t_real));
    ccta->sign_cur = calloc(fft_length_half * 2, sizeof(t_real));
}

static void destroyCCTimescaleAngle(struct CCTimescaleAngle *ccta) {
    free(ccta->sign_last);
    free(ccta->sign_cur);
}

static void ingestCCTimescaleAngle(struct CCTimescaleAngle *ccta, const t_real *real, const t_real *imag, t_real *out) {
    t_len i = 0;
    t_real gx, gy;
    t_real *last;
    t_real *cur;
    
    // swap cur and last
    last = ccta->sign_cur;
    cur = ccta->sign_last;
    ccta->sign_last = last;
    ccta->sign_cur = cur;
    
    // calculate the imaginary portion of the angle (FOIL)
    vDSP(vsmsma)(real, 1, &ccta->mu_imag, imag, 1, &ccta->mu_real, cur, 1, ccta->fft_length_half);
    vDSP(vlim)(cur, 1, &s_zero, &s_one, cur, 1, ccta->fft_length_half); // limit to +1 / -1
    
    // special case: i = 0
    gx = (cur[i] - last[i]) * ccta->grad_x;
    out[i] = (gx > NUM(0.001) ? NUM(1.0) : NUM(0.0));
    
    // iterate
    for (i = 1; i < ccta->fft_length_half; ++i) {
        // adjust by angle
        // TODO: optimize via switch case and
        gx = (cur[i] - last[i]) * ccta->grad_x;
        gy = (cur[i] - cur[i - 1]) * ccta->grad_y;
        
        out[i] = (gx + gy > NUM(0.001) ? NUM(1.0) : NUM(0.0));
    }
}

static void addConsensusToContours(const t_len fft_length_half, const t_real *consensus, t_real *consensus_contours) {
    for (t_len i = 0; i < fft_length_half; ++i) {
        if (consensus[i] > NUM(1.5)) {
            consensus_contours[i] += NUM(1.0);
        }
    }
}

static t_len fftSize(const t_len length) {
    double size = ceil(log2((double)length));
    return (t_len)size;
}

struct OpaqueCCCSetup
{
    /* parameters */
    // cannot be changed
    t_len fft_length; /* must be power of two */
    t_len fft_overlap;
    
    // can be changed
    bool pow_weight;
    t_real fs;
    
    // currently not customizable
    t_real angles[N_ANGLES]; /* radians */
    t_real timescales[N_TIMESCALES]; /* ms */
    
    /* internal data */
    t_len fft_length_half;
    t_len fft_size;
    
    /* internal memory */
    RealFFTSetup fft_setup;
    t_real *fft_window, *signal_windowed, *power;
    t_real *consensus, *consensus_cur, *consensus_pow;
    RealSplitComplex fft_temporary, fft_output;
    RealSplitComplex p_exp, p_der, p_ratio; // pointers into the fft_output for convience
    struct CCTimescaleAngle ta[N_TIMESCALES * N_ANGLES];
};

static void fillFftWindow(const CCCSetup setup) {
    t_len i, j;
    t_len index = 0;
    t_real timescale_samples;
    t_real w;
    t_real twin, twin_offset = NUM(0.5) + (NUM(0.0) - (t_real)setup->fft_length) / NUM(2.);
    t_len index_offset = N_TIMESCALES * setup->fft_length;
    
    for (i = 0; i < N_TIMESCALES; ++i) {
        timescale_samples = setup->fs * setup->timescales[i] / NUM(1000.0);
        
        /* TODO: could be heabily optimized using vvexp and such */
        for (j = 0; j < setup->fft_length; ++j) {
            twin = twin_offset + (t_real)j;
            w = CMATH(exp)(-CMATH(pow)(twin / timescale_samples, 2));
            
            // store exponential
            setup->fft_window[index + j] = w;
            
            // store derivative
            setup->fft_window[index_offset + index + j] = NUM(-2.0) * w * twin / CMATH(pow)(timescale_samples, 2);
        }
        
        index += setup->fft_length;
    }
}

CCCSetup createCCCSetup(t_len fft_length, t_len fft_overlap, t_real fs, bool pow_weight) {
    t_len i, j;
    
    /* allocate memory */
    CCCSetup ret = malloc(sizeof(struct OpaqueCCCSetup));
    
    /* store parameters */
    ret->fft_length = fft_length;
    ret->fft_overlap = fft_overlap;
    ret->pow_weight = pow_weight;
    ret->fs = fs;
    
    /* fill timescales and angles */
    for (i = 0; i < N_ANGLES; ++i) {
        ret->angles[i] = M_PI * (t_real)(2 + i) / N_ANGLES;
    }
    for (i = 0; i < N_TIMESCALES; ++i) {
        ret->timescales[i] = NUM(0.5) + NUM(0.2) * (t_real)i;
    }
    
    /* STEP 1: spectrogram, build windows and fft */
    ret->fft_length_half = fft_length / 2;
    ret->fft_size = fftSize(fft_length);
    assert((1 << ret->fft_size) == ret->fft_length); // should be power of 2
    
    ret->fft_setup = vDSP(create_fftsetup)(ret->fft_size, kFFTRadix2);
    
    /* allocate window */
    ret->fft_window = malloc(sizeof(t_real) * fft_length * N_TIMESCALES * 2);
    ret->signal_windowed = malloc(sizeof(t_real) * fft_length * N_TIMESCALES * 2);
    
    /* allocate temporary and output */
    ret->fft_temporary.realp = malloc(sizeof(t_real) * ret->fft_length_half);
    ret->fft_temporary.imagp = malloc(sizeof(t_real) * ret->fft_length_half);
    ret->fft_output.realp = malloc(sizeof(t_real) * ret->fft_length_half * N_TIMESCALES * 2);
    ret->fft_output.imagp = malloc(sizeof(t_real) * ret->fft_length_half * N_TIMESCALES * 2);
    
    /* make pointers */
    ret->p_exp = ret->fft_output;
    ret->p_der.realp = ret->fft_output.realp + ret->fft_length_half * N_TIMESCALES;
    ret->p_der.imagp = ret->fft_output.imagp + ret->fft_length_half * N_TIMESCALES;
    ret->p_ratio = ret->p_der;
    
    /* power array */
    ret->power = malloc(sizeof(t_real) * ret->fft_length_half * N_TIMESCALES);
    
    /* conensus array */
    /* might be able to use singles? */
    ret->consensus = malloc(sizeof(t_real) * ret->fft_length_half * N_TIMESCALES * N_ANGLES);
    ret->consensus_cur = malloc(sizeof(t_real) * ret->fft_length_half);
    ret->consensus_pow = malloc(sizeof(t_real) * ret->fft_length_half);
    
    /* fill window */
    fillFftWindow(ret);
    
    /* setup timescale angle structures */
    for (i = 0; i < N_TIMESCALES; ++i) {
        for (j = 0; j < N_ANGLES; ++j) {
            setupCCTimescaleAngle(&ret->ta[i * N_ANGLES + j], ret->fft_length_half, ret->timescales[i], ret->angles[j]);
        }
    }
    
    return ret;
}

void destroyCCCSetup(CCCSetup setup) {
    t_len i, j;
    
    /* free memory */
    free(setup->fft_window);
    free(setup->signal_windowed);
    free(setup->fft_temporary.realp);
    free(setup->fft_temporary.imagp);
    free(setup->fft_output.realp);
    free(setup->fft_output.imagp);
    free(setup->power);
    free(setup->consensus);
    free(setup->consensus_cur);
    free(setup->consensus_pow);
    
    /* clean up FFT setup */
    vDSP(destroy_fftsetup)(setup->fft_setup);
    
    /* destroy timescale angle structures */
    for (i = 0; i < N_TIMESCALES; ++i) {
        for (j = 0; j < N_ANGLES; ++j) {
            destroyCCTimescaleAngle(&setup->ta[i * N_ANGLES + j]);
        }
    }
    
    /* free memory for whole setup */
    free(setup);
}

static void windowSamples(const CCCSetup setup, const t_real *signal) {
    t_len i = 0;
    t_len signals = 2 * N_TIMESCALES;
    
    for (i = 0; i < signals; ++i) {
        vDSP(vmul)(signal, 1, setup->fft_window + i * setup->fft_length, 1, setup->signal_windowed + i * setup->fft_length, 1, setup->fft_length);
    }
}

struct ConsensusContourSize cccSize(const CCCSetup setup, const t_len signal_len) {
    struct ConsensusContourSize ret;
    
    // store signal length
    ret.signal_len = signal_len;
    
    // insufficient signal?
    if (signal_len < setup->fft_length) {
        ret.rows = 0;
        ret.cols = 0;
        ret.bytes = 0;
        return ret;
    }
    
    const t_len ccn = setup->fft_length / 2;
    const t_len ccm = 1 + (signal_len - setup->fft_length) / (setup->fft_length - setup->fft_overlap);
    
    ret.rows = ccn;
    ret.cols = ccm;
    ret.bytes = sizeof(t_real) * ccn * ccm;
    
    return ret;
}

/* signal must be setup->fft_length long */
static void buildColumn(const CCCSetup setup, const t_real *signal, t_real *consensus_contour) {
    t_len j, k;
    
#ifdef TIMING
    clock_t begin, end;
    begin = clock();
#endif
    
    // window
    windowSamples(setup, signal);
    
    // pack samples
    vDSP(ctoz)((RealComplex *)setup->signal_windowed, 2, &setup->fft_output, 1, setup->fft_length_half * N_TIMESCALES * 2);
    
    // calculate
    vDSP(fftm_zript)(setup->fft_setup, &setup->fft_output, 1, (t_stride)setup->fft_length_half, &setup->fft_temporary, setup->fft_size, N_TIMESCALES * 2, FFT_FORWARD);
    
    // calculate power
    vDSP(zvabs)(&setup->fft_output, 1, setup->power, 1, setup->fft_length_half * N_TIMESCALES);
    vDSP(vsdiv)(setup->power, 1, &s_power_scale, setup->power, 1, setup->fft_length_half * N_TIMESCALES);
    
    // calculate ratio
    vDSP(zvdiv)(&setup->p_exp, 1, &setup->p_der, 1, &setup->p_ratio, 1, setup->fft_length_half * N_TIMESCALES);
    
    // calculate contours for each timescale and angle
    for (j = 0; j < N_TIMESCALES; ++j) {
        for (k = 0; k < N_ANGLES; ++k) {
            ingestCCTimescaleAngle(&setup->ta[j * N_ANGLES + k], setup->p_ratio.realp + j * setup->fft_length_half, setup->p_ratio.imagp + j * setup->fft_length_half, setup->consensus + (j * N_ANGLES + k) * setup->fft_length_half);
        }
    }
    
    // look for consensus
    for (j = 0; j < (N_TIMESCALES - 1); ++j) {
        if (setup->pow_weight) {
            vDSP(vclr)(setup->consensus_pow, 1, setup->fft_length_half);
        }
        
        for (k = 0; k < N_ANGLES; ++k) {
            // start with consensus / angle
            memcpy(setup->consensus_cur, setup->consensus + (j * N_ANGLES + k) * setup->fft_length_half, sizeof(t_real) * setup->fft_length_half);
            // add next sigma, same angle
            vDSP(vadd)(setup->consensus_cur, 1, setup->consensus + ((j + 1) * N_ANGLES + k) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            // add same sigma, previous angle (wrap around)
            vDSP(vadd)(setup->consensus_cur, 1, setup->consensus + (j * N_ANGLES + (0 == k ? N_ANGLES - 1 : k - 1)) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            
            // add to contours
            if (setup->pow_weight) {
                addConsensusToContours(setup->fft_length_half, setup->consensus_cur, setup->consensus_pow);
            }
            else {
                addConsensusToContours(setup->fft_length_half, setup->consensus_cur, consensus_contour);
            }
        }
        
        if (setup->pow_weight) {
            // scale by power
            vDSP(vmul)(setup->consensus_pow, 1, setup->power + j * setup->fft_length_half, 1, setup->consensus_pow, 1, setup->fft_length_half);
            
            // add to output
            vDSP(vadd)(consensus_contour, 1, setup->consensus_pow, 1, consensus_contour, 1, setup->fft_length_half);
        }
    }
    
#ifdef TIMING
    end = clock();
    printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);
#endif
}

void cccColumn(const CCCSetup setup, const t_real *signal, t_real *consensus_contour) {
    /* clear output */
    vDSP(vclr)(consensus_contour, 1, setup->fft_length_half);
    
    /* call internal function */
    buildColumn(setup, signal, consensus_contour);
}

void ccc(const CCCSetup setup, const struct ConsensusContourSize dim, const t_real *signal, t_real *consensus_contours) {
    t_len i;
    
    /* clear output */
    vDSP(vclr)(consensus_contours, 1, dim.cols * dim.rows);
    
    /* STEP 2: spectrogram columns */
    for (i = 0; i < dim.cols; ++i) {
        buildColumn(setup, signal + i * (setup->fft_length - setup->fft_overlap), consensus_contours + i * setup->fft_length_half);
    }
}
