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

static double s_power_scale = 2.0;

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

static void addConsensusToContours(const t_len fft_length_half, const double *consensus, double *consensus_contours) {
    for (t_len i = 0; i < fft_length_half; ++i) {
        if (consensus[i] > 1.5) {
            consensus_contours[i] += 1;
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
    double fs;
    
    // currently not customizable
    double angles[N_ANGLES]; /* radians */
    double timescales[N_TIMESCALES]; /* ms */
    
    /* internal data */
    t_len fft_length_half;
    t_len fft_size;
    
    /* internal memory */
    FFTSetupD fft_setup;
    double *fft_window, *signal_windowed, *power;
    double *consensus, *consensus_cur, *consensus_pow;
    DSPDoubleSplitComplex fft_temporary, fft_output;
    DSPDoubleSplitComplex p_exp, p_der, p_ratio; // pointers into the fft_output for convience
    struct CCTimescaleAngle ta[N_TIMESCALES * N_ANGLES];
};

static void fillFftWindow(const CCCSetup setup) {
    int i;
    t_len j;
    t_len index = 0;
    double timescale_samples;
    double w;
    double twin, twin_offset = 0.5 + (0 - (double)setup->fft_length) / 2;
    const t_len index_offset = N_TIMESCALES * setup->fft_length;
    
    for (i = 0; i < N_TIMESCALES; ++i) {
        timescale_samples = setup->fs * setup->timescales[i] / 1000;
        
        /* TODO: could be heabily optimized using vvexp and such */
        for (j = 0; j < setup->fft_length; ++j) {
            twin = twin_offset + (double)j;
            w = exp(-pow(twin / timescale_samples, 2));
            
            // store exponential
            setup->fft_window[index + j] = w;
            
            // store derivative
            setup->fft_window[index_offset + index + j] = -2 * w * twin / pow(timescale_samples, 2);
        }
        
        index += setup->fft_length;
    }
}

CCCSetup createCCCSetup(t_len fft_length, t_len fft_overlap, double fs, bool pow_weight) {
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
        ret->angles[i] = M_PI * (double)(2 + i) / N_ANGLES;
    }
    for (i = 0; i < N_TIMESCALES; ++i) {
        ret->timescales[i] = 0.5 + 0.2 * (double)i;
    }
    
    /* STEP 1: spectrogram, build windows and fft */
    ret->fft_length_half = fft_length / 2;
    ret->fft_size = fftSize(fft_length);
    assert((1 << ret->fft_size) == ret->fft_length); // should be power of 2
    
    ret->fft_setup = vDSP_create_fftsetupD(ret->fft_size, kFFTRadix2);
    
    /* allocate window */
    ret->fft_window = malloc(sizeof(double) * fft_length * N_TIMESCALES * 2);
    ret->signal_windowed = malloc(sizeof(double) * fft_length * N_TIMESCALES * 2);
    
    /* allocate temporary and output */
    ret->fft_temporary.realp = malloc(sizeof(double) * ret->fft_length_half);
    ret->fft_temporary.imagp = malloc(sizeof(double) * ret->fft_length_half);
    ret->fft_output.realp = malloc(sizeof(double) * ret->fft_length_half * N_TIMESCALES * 2);
    ret->fft_output.imagp = malloc(sizeof(double) * ret->fft_length_half * N_TIMESCALES * 2);
    
    /* make pointers */
    ret->p_exp = ret->fft_output;
    ret->p_der.realp = ret->fft_output.realp + ret->fft_length_half * N_TIMESCALES;
    ret->p_der.imagp = ret->fft_output.imagp + ret->fft_length_half * N_TIMESCALES;
    ret->p_ratio = ret->p_der;
    
    /* power array */
    ret->power = malloc(sizeof(double) * ret->fft_length_half * N_TIMESCALES);
    
    /* conensus array */
    /* might be able to use singles? */
    ret->consensus = malloc(sizeof(double) * ret->fft_length_half * N_TIMESCALES * N_ANGLES);
    ret->consensus_cur = malloc(sizeof(double) * ret->fft_length_half);
    ret->consensus_pow = malloc(sizeof(double) * ret->fft_length_half);
    
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
    vDSP_destroy_fftsetupD(setup->fft_setup);
    
    /* destroy timescale angle structures */
    for (i = 0; i < N_TIMESCALES; ++i) {
        for (j = 0; j < N_ANGLES; ++j) {
            destroyCCTimescaleAngle(&setup->ta[i * N_ANGLES + j]);
        }
    }
    
    /* free memory for whole setup */
    free(setup);
}

static void windowSamples(const CCCSetup setup, const double *signal) {
    t_len i = 0;
    t_len signals = 2 * N_TIMESCALES;
    
    for (i = 0; i < signals; ++i) {
        vDSP_vmulD(signal, 1, setup->fft_window + i * setup->fft_length, 1, setup->signal_windowed + i * setup->fft_length, 1, setup->fft_length);
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
    ret.bytes = sizeof(double) * ccn * ccm;
    
    return ret;
}

/* signal must be setup->fft_length long */
static void buildColumn(const CCCSetup setup, const double *signal, double *consensus_contour) {
    t_len j, k;
    
#ifdef TIMING
    clock_t begin, end;
    begin = clock();
#endif
    
    // window
    windowSamples(setup, signal);
    
    // pack samples
    vDSP_ctozD((DSPDoubleComplex *)setup->signal_windowed, 2, &setup->fft_output, 1, setup->fft_length_half * N_TIMESCALES * 2);
    
    // calculate
    vDSP_fftm_zriptD(setup->fft_setup, &setup->fft_output, 1, (t_stride)setup->fft_length_half, &setup->fft_temporary, setup->fft_size, N_TIMESCALES * 2, FFT_FORWARD);
    
    // calculate power
    vDSP_zvabsD(&setup->fft_output, 1, setup->power, 1, setup->fft_length_half * N_TIMESCALES);
    vDSP_vsdivD(setup->power, 1, &s_power_scale, setup->power, 1, setup->fft_length_half * N_TIMESCALES);
    
    // calculate ratio
    vDSP_zvdivD(&setup->p_exp, 1, &setup->p_der, 1, &setup->p_ratio, 1, setup->fft_length_half * N_TIMESCALES);
    
    // calculate contours for each timescale and angle
    for (j = 0; j < N_TIMESCALES; ++j) {
        for (k = 0; k < N_ANGLES; ++k) {
            ingestCCTimescaleAngle(&setup->ta[j * N_ANGLES + k], setup->p_ratio.realp + j * setup->fft_length_half, setup->p_ratio.imagp + j * setup->fft_length_half, setup->consensus + (j * N_ANGLES + k) * setup->fft_length_half);
        }
    }
    
    // look for consensus
    for (j = 0; j < (N_TIMESCALES - 1); ++j) {
        if (setup->pow_weight) {
            vDSP_vclrD(setup->consensus_pow, 1, setup->fft_length_half);
        }
        
        for (k = 0; k < N_ANGLES; ++k) {
            // start with consensus / angle
            memcpy(setup->consensus_cur, setup->consensus + (j * N_ANGLES + k) * setup->fft_length_half, sizeof(double) * setup->fft_length_half);
            // add next sigma, same angle
            vDSP_vaddD(setup->consensus_cur, 1, setup->consensus + ((j + 1) * N_ANGLES + k) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            // add same sigma, previous angle (wrap around)
            vDSP_vaddD(setup->consensus_cur, 1, setup->consensus + (j * N_ANGLES + (0 == k ? N_ANGLES - 1 : k - 1)) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            
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
            vDSP_vmulD(setup->consensus_pow, 1, setup->power + j * setup->fft_length_half, 1, setup->consensus_pow, 1, setup->fft_length_half);
            
            // add to output
            vDSP_vaddD(consensus_contour, 1, setup->consensus_pow, 1, consensus_contour, 1, setup->fft_length_half);
        }
    }
    
#ifdef TIMING
    end = clock();
    printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);
#endif
}

void cccColumn(const CCCSetup setup, const double *signal, double *consensus_contour) {
    /* clear output */
    vDSP_vclrD(consensus_contour, 1, setup->fft_length_half);
    
    /* call internal function */
    buildColumn(setup, signal, consensus_contour);
}

void ccc(const CCCSetup setup, const struct ConsensusContourSize dim, const double *signal, double *consensus_contours) {
    t_len i;
    
    /* clear output */
    vDSP_vclrD(consensus_contours, 1, dim.cols * dim.rows);
    
    /* STEP 2: spectrogram columns */
    for (i = 0; i < dim.cols; ++i) {
        buildColumn(setup, signal + i * (setup->fft_length - setup->fft_overlap), consensus_contours + i * setup->fft_length_half);
    }
}
