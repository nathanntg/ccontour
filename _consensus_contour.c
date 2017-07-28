static REAL TYPE(s_zero) = 0.0;
static REAL TYPE(s_one) = 1.0;
static REAL TYPE(s_two) = 2.0;

struct TYPE(CCTimescaleAngle)
{
    t_len fft_length_half;
    
    REAL timescale;
    REAL angle;
    
    REAL mu_real;
    REAL mu_imag;
    
    REAL grad_x;
    REAL grad_y;
    
    // TODO: could potentially be int8_t
    REAL *sign_last;
    REAL *sign_cur;
};

static void TYPE(setupCCTimescaleAngle)(struct TYPE(CCTimescaleAngle) *ccta, const t_len fft_length_half, const REAL timescale, const REAL angle) {
    ccta->fft_length_half = fft_length_half;
    ccta->timescale = timescale;
    ccta->angle = angle;
    ccta->mu_real = CMATH(cos)(angle);
    ccta->mu_imag = CMATH(sin)(angle);
    ccta->grad_x = 0 - CMATH(cos)(angle + (REAL)M_PI_2) / CMATH(2.0);
    ccta->grad_y = CMATH(sin)(angle + (REAL)M_PI_2) / CMATH(2.0);
    ccta->sign_last = calloc(fft_length_half * 2, sizeof(REAL));
    ccta->sign_cur = calloc(fft_length_half * 2, sizeof(REAL));
}

static void TYPE(destroyCCTimescaleAngle)(struct TYPE(CCTimescaleAngle) *ccta) {
    free(ccta->sign_last);
    free(ccta->sign_cur);
}

static void TYPE(ingestCCTimescaleAngle)(struct TYPE(CCTimescaleAngle) *ccta, const REAL *real, const REAL *imag, REAL *out) {
    t_len i = 0;
    REAL gx, gy;
    REAL *last;
    REAL *cur;
    
    // swap cur and last
    last = ccta->sign_cur;
    cur = ccta->sign_last;
    ccta->sign_last = last;
    ccta->sign_cur = cur;
    
    // calculate the imaginary portion of the angle (FOIL)
    vDSP(vsmsma)(real, 1, &ccta->mu_imag, imag, 1, &ccta->mu_real, cur, 1, ccta->fft_length_half);
    vDSP(vlim)(cur, 1, &TYPE(s_zero), &TYPE(s_one), cur, 1, ccta->fft_length_half); // limit to +1 / -1
    
    // special case: i = 0
    gx = (cur[i] - last[i]) * ccta->grad_x;
    out[i] = (gx > CMATH(0.001) ? CMATH(1.0) : CMATH(0.0));
    
    // iterate
    for (i = 1; i < ccta->fft_length_half; ++i) {
        // adjust by angle
        // TODO: optimize via switch case and
        gx = (cur[i] - last[i]) * ccta->grad_x;
        gy = (cur[i] - cur[i - 1]) * ccta->grad_y;
        
        out[i] = (gx + gy > CMATH(0.001) ? CMATH(1.0) : CMATH(0.0));
    }
}

static void TYPE(addConsensusToContours)(const t_len fft_length_half, const REAL *consensus, REAL *consensus_contours) {
    for (t_len i = 0; i < fft_length_half; ++i) {
        if (consensus[i] > CMATH(1.5)) {
            consensus_contours[i] += CMATH(1.0);
        }
    }
}


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"

struct TYPE(OpaqueCCCSetup)
{
    /* parameters */
    // cannot be changed
    t_len fft_length; /* must be power of two */
    t_len fft_overlap;
    
    // can be changed
    bool pow_weight;
    REAL fs;
    
    // currently not customizable
    REAL angles[N_ANGLES]; /* radians */
    REAL timescales[N_TIMESCALES]; /* ms */
    
    /* internal data */
    t_len fft_length_half;
    t_len fft_size;
    
    /* internal memory */
    REAL_FFT_SETUP fft_setup;
    REAL *fft_window, *signal_windowed, *power;
    REAL *consensus, *consensus_cur, *consensus_pow;
    REAL_SPLIT_COMPLEX fft_temporary, fft_output;
    REAL_SPLIT_COMPLEX p_exp, p_der, p_ratio; // pointers into the fft_output for convience
    struct TYPE(CCTimescaleAngle) ta[N_TIMESCALES * N_ANGLES];
};

#pragma clang diagnostic pop

static void TYPE(fillFftWindow)(const TYPE(CCCSetup) setup) {
    t_len i, j;
    t_len index = 0;
    REAL timescale_samples;
    REAL w;
    REAL twin, twin_offset = CMATH(0.5) + (CMATH(0.0) - (REAL)setup->fft_length) / CMATH(2.);
    t_len index_offset = N_TIMESCALES * setup->fft_length;
    
    for (i = 0; i < N_TIMESCALES; ++i) {
        timescale_samples = setup->fs * setup->timescales[i] / CMATH(1000.0);
        
        /* TODO: could be heabily optimized using vvexp and such */
        for (j = 0; j < setup->fft_length; ++j) {
            twin = twin_offset + (REAL)j;
            w = CMATH(exp)(-CMATH(pow)(twin / timescale_samples, 2));
            
            // store exponential
            setup->fft_window[index + j] = w;
            
            // store derivative
            setup->fft_window[index_offset + index + j] = CMATH(-2.0) * w * twin / CMATH(pow)(timescale_samples, 2);
        }
        
        index += setup->fft_length;
    }
}

TYPE(CCCSetup) TYPE(createCCCSetup)(t_len fft_length, t_len fft_overlap, REAL fs, bool pow_weight) {
    t_len i, j;
    
    /* allocate memory */
    TYPE(CCCSetup) ret = malloc(sizeof(struct TYPE(OpaqueCCCSetup)));
    
    /* store parameters */
    ret->fft_length = fft_length;
    ret->fft_overlap = fft_overlap;
    ret->pow_weight = pow_weight;
    ret->fs = fs;
    
    /* fill timescales and angles */
    for (i = 0; i < N_ANGLES; ++i) {
        ret->angles[i] = (REAL)M_PI * (REAL)(2 + i) / N_ANGLES;
    }
    for (i = 0; i < N_TIMESCALES; ++i) {
        ret->timescales[i] = CMATH(0.5) + CMATH(0.2) * (REAL)i;
    }
    
    /* STEP 1: spectrogram, build windows and fft */
    ret->fft_length_half = fft_length / 2;
    ret->fft_size = fftSize(fft_length);
    assert((1 << ret->fft_size) == ret->fft_length); // should be power of 2
    
    ret->fft_setup = vDSP(create_fftsetup)(ret->fft_size, kFFTRadix2);
    
    /* allocate window */
    ret->fft_window = malloc(sizeof(REAL) * fft_length * N_TIMESCALES * 2);
    ret->signal_windowed = malloc(sizeof(REAL) * fft_length * N_TIMESCALES * 2);
    
    /* allocate temporary and output */
    ret->fft_temporary.realp = malloc(sizeof(REAL) * ret->fft_length_half);
    ret->fft_temporary.imagp = malloc(sizeof(REAL) * ret->fft_length_half);
    ret->fft_output.realp = malloc(sizeof(REAL) * ret->fft_length_half * N_TIMESCALES * 2);
    ret->fft_output.imagp = malloc(sizeof(REAL) * ret->fft_length_half * N_TIMESCALES * 2);
    
    /* make pointers */
    ret->p_exp = ret->fft_output;
    ret->p_der.realp = ret->fft_output.realp + ret->fft_length_half * N_TIMESCALES;
    ret->p_der.imagp = ret->fft_output.imagp + ret->fft_length_half * N_TIMESCALES;
    ret->p_ratio = ret->p_der;
    
    /* power array */
    ret->power = malloc(sizeof(REAL) * ret->fft_length_half * N_TIMESCALES);
    
    /* conensus array */
    /* might be able to use singles? */
    ret->consensus = malloc(sizeof(REAL) * ret->fft_length_half * N_TIMESCALES * N_ANGLES);
    ret->consensus_cur = malloc(sizeof(REAL) * ret->fft_length_half);
    ret->consensus_pow = malloc(sizeof(REAL) * ret->fft_length_half);
    
    /* fill window */
    TYPE(fillFftWindow)(ret);
    
    /* setup timescale angle structures */
    for (i = 0; i < N_TIMESCALES; ++i) {
        for (j = 0; j < N_ANGLES; ++j) {
            TYPE(setupCCTimescaleAngle)(&ret->ta[i * N_ANGLES + j], ret->fft_length_half, ret->timescales[i], ret->angles[j]);
        }
    }
    
    return ret;
}

void TYPE(destroyCCCSetup)(TYPE(CCCSetup) setup) {
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
            TYPE(destroyCCTimescaleAngle)(&setup->ta[i * N_ANGLES + j]);
        }
    }
    
    /* free memory for whole setup */
    free(setup);
}

static void TYPE(windowSamples)(const TYPE(CCCSetup) setup, const REAL *signal) {
    t_len i = 0;
    t_len signals = 2 * N_TIMESCALES;
    
    for (i = 0; i < signals; ++i) {
        vDSP(vmul)(signal, 1, setup->fft_window + i * setup->fft_length, 1, setup->signal_windowed + i * setup->fft_length, 1, setup->fft_length);
    }
}

struct ConsensusContourSize TYPE(cccSize)(const TYPE(CCCSetup) setup, const t_len signal_len) {
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
    ret.bytes = sizeof(REAL) * ccn * ccm;
    
    return ret;
}

/* signal must be setup->fft_length long */
static void TYPE(buildColumn)(const TYPE(CCCSetup) setup, const REAL *signal, REAL *consensus_contour) {
    t_len j, k;
    
#ifdef TIMING
    clock_t begin, end;
    begin = clock();
#endif
    
    // window
    TYPE(windowSamples)(setup, signal);
    
    // pack samples
    vDSP(ctoz)((REAL_COMPLEX *)setup->signal_windowed, 2, &setup->fft_output, 1, setup->fft_length_half * N_TIMESCALES * 2);
    
    // calculate
    vDSP(fftm_zript)(setup->fft_setup, &setup->fft_output, 1, (t_stride)setup->fft_length_half, &setup->fft_temporary, setup->fft_size, N_TIMESCALES * 2, FFT_FORWARD);
    
    // calculate power
    vDSP(zvabs)(&setup->fft_output, 1, setup->power, 1, setup->fft_length_half * N_TIMESCALES);
    vDSP(vsdiv)(setup->power, 1, &TYPE(s_two), setup->power, 1, setup->fft_length_half * N_TIMESCALES);
    
    // calculate ratio
    vDSP(zvdiv)(&setup->p_exp, 1, &setup->p_der, 1, &setup->p_ratio, 1, setup->fft_length_half * N_TIMESCALES);
    
    // calculate contours for each timescale and angle
    for (j = 0; j < N_TIMESCALES; ++j) {
        for (k = 0; k < N_ANGLES; ++k) {
            TYPE(ingestCCTimescaleAngle)(&setup->ta[j * N_ANGLES + k], setup->p_ratio.realp + j * setup->fft_length_half, setup->p_ratio.imagp + j * setup->fft_length_half, setup->consensus + (j * N_ANGLES + k) * setup->fft_length_half);
        }
    }
    
    // look for consensus
    for (j = 0; j < (N_TIMESCALES - 1); ++j) {
        if (setup->pow_weight) {
            vDSP(vclr)(setup->consensus_pow, 1, setup->fft_length_half);
        }
        
        for (k = 0; k < N_ANGLES; ++k) {
            // start with consensus / angle
            memcpy(setup->consensus_cur, setup->consensus + (j * N_ANGLES + k) * setup->fft_length_half, sizeof(REAL) * setup->fft_length_half);
            // add next sigma, same angle
            vDSP(vadd)(setup->consensus_cur, 1, setup->consensus + ((j + 1) * N_ANGLES + k) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            // add same sigma, previous angle (wrap around)
            vDSP(vadd)(setup->consensus_cur, 1, setup->consensus + (j * N_ANGLES + (0 == k ? N_ANGLES - 1 : k - 1)) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            
            // add to contours
            if (setup->pow_weight) {
                TYPE(addConsensusToContours)(setup->fft_length_half, setup->consensus_cur, setup->consensus_pow);
            }
            else {
                TYPE(addConsensusToContours)(setup->fft_length_half, setup->consensus_cur, consensus_contour);
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

void TYPE(cccColumn)(const TYPE(CCCSetup) setup, const REAL *signal, REAL *consensus_contour) {
    /* clear output */
    vDSP(vclr)(consensus_contour, 1, setup->fft_length_half);
    
    /* call internal function */
    TYPE(buildColumn)(setup, signal, consensus_contour);
}

void TYPE(cccSpectrogram)(const TYPE(CCCSetup) setup, const struct ConsensusContourSize dim, const REAL *signal, REAL *consensus_contours) {
    t_len i;
    
    /* clear output */
    vDSP(vclr)(consensus_contours, 1, dim.cols * dim.rows);
    
    /* STEP 2: spectrogram columns */
    for (i = 0; i < dim.cols; ++i) {
        TYPE(buildColumn)(setup, signal + i * (setup->fft_length - setup->fft_overlap), consensus_contours + i * setup->fft_length_half);
    }
}

