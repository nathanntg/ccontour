static REAL TYPE(s_zero) = 0.0;
static REAL TYPE(s_one) = 1.0;
static REAL TYPE(s_two) = 2.0;

/* configuration */

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"

struct TYPE(OpaqueCCCConfig) {
    t_len fft_length;
    t_len fft_shift;
    
    bool pow_weight;
    
    REAL fs;
    
    t_len num_timescales;
    REAL *timescales;
    
    t_len num_angles;
    REAL *angles;
};

#pragma clang diagnostic pop

TYPE(CCCConfig) TYPE(createCCCConfig)() {
    t_len i;
    
    // allocate memory
    TYPE(CCCConfig) ret = (TYPE(struct OpaqueCCCConfig) *)malloc(sizeof(TYPE(struct OpaqueCCCConfig)));
    
    // default fft parameters
    ret->fft_length = 1024;
    ret->fft_shift = 19;
    ret->pow_weight = true;
    ret->fs = 44100;
    
    // default timescales
    ret->num_timescales = 9;
    ret->timescales = (REAL *)malloc(sizeof(REAL) * ret->num_timescales);
    for (i = 0; i < ret->num_timescales; ++i) {
        ret->timescales[i] = CMATH(0.5) + CMATH(0.2) * (REAL)i;
    }
    
    // default angles
    ret->num_angles = 8;
    ret->angles = (REAL *)malloc(sizeof(REAL) * ret->num_angles);
    for (i = 0; i < ret->num_angles; ++i) {
        ret->angles[i] = (REAL)M_PI * (REAL)(2 + i) / (REAL)ret->num_angles;
    }
    
    return ret;
}

void TYPE(cccConfigSetFFTLength)(TYPE(CCCConfig) config, const t_len fft_length) {
    config->fft_length = fft_length;
}

void TYPE(cccConfigSetFFTShift)(TYPE(CCCConfig) config, const t_len fft_shift) {
    config->fft_shift = fft_shift;
}

void TYPE(cccConfigSetWeightByPower)(TYPE(CCCConfig) config, const bool pow_weight) {
    config->pow_weight = pow_weight;
}

void TYPE(cccConfigSetSampleRate)(TYPE(CCCConfig) config, const REAL fs) {
    config->fs = fs;
}

void TYPE(cccConfigSetTimescales)(TYPE(CCCConfig) config, const t_len num_timescales, const REAL timescales[]) {
    // free old timescales
    free(config->timescales);
    
    // allocate new timescales
    config->num_timescales = num_timescales;
    config->timescales = (REAL *)malloc(sizeof(REAL) * config->num_timescales);
    memcpy(config->timescales, &timescales[0], sizeof(REAL) * config->num_timescales);
}

void TYPE(cccConfigSetAngles)(TYPE(CCCConfig) config, const t_len num_angles, const REAL angles[]) {
    // free old timescales
    free(config->angles);
    
    // allocate new timescales
    config->num_angles = num_angles;
    config->angles = (REAL *)malloc(sizeof(REAL) * config->num_angles);
    memcpy(config->angles, &angles[0], sizeof(REAL) * config->num_angles);
}

void TYPE(destroyCCCConfig)(TYPE(CCCConfig) config) {
    /* free memory */
    free(config->timescales);
    free(config->angles);
    
    /* free memory for whole config */
    free(config);
}

/* timescale angle */

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
    ccta->sign_last = (REAL *)calloc(fft_length_half, sizeof(REAL));
    ccta->sign_cur = (REAL *)calloc(fft_length_half, sizeof(REAL));
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
    /* config */
    TYPE(struct OpaqueCCCConfig) config;
    
    /* internal data */
    t_len fft_length_half;
    t_len fft_size;
    
    /* internal memory */
    REAL_FFT_SETUP fft_setup;
    REAL *fft_window, *signal_windowed, *power;
    REAL *consensus, *consensus_cur, *consensus_pow;
    REAL_SPLIT_COMPLEX fft_temporary, fft_output;
    REAL_SPLIT_COMPLEX p_exp, p_der, p_ratio; // pointers into the fft_output for convience
    struct TYPE(CCTimescaleAngle) *ta;
};

#pragma clang diagnostic pop

static void TYPE(fillFftWindow)(const TYPE(CCCSetup) setup) {
    t_len i, j;
    t_len index = 0;
    REAL timescale_samples;
    REAL w;
    REAL twin, twin_offset = CMATH(0.5) + (CMATH(0.0) - (REAL)setup->config.fft_length) / CMATH(2.);
    t_len index_offset = setup->config.num_timescales * setup->config.fft_length;
    
    for (i = 0; i < setup->config.num_timescales; ++i) {
        timescale_samples = setup->config.fs * setup->config.timescales[i] / CMATH(1000.0);
        
        /* TODO: could be heabily optimized using vvexp and such */
        for (j = 0; j < setup->config.fft_length; ++j) {
            twin = twin_offset + (REAL)j;
            w = CMATH(exp)(-CMATH(pow)(twin / timescale_samples, 2));
            
            // store exponential
            setup->fft_window[index + j] = w;
            
            // store derivative
            setup->fft_window[index_offset + index + j] = CMATH(-2.0) * w * twin / CMATH(pow)(timescale_samples, 2);
        }
        
        index += setup->config.fft_length;
    }
}

TYPE(CCCSetup) TYPE(createCCCSetup)(const TYPE(CCCConfig) config) {
    t_len i, j;
    
    /* allocate memory */
    TYPE(CCCSetup) ret = (TYPE(struct OpaqueCCCSetup) *)malloc(sizeof(TYPE(struct OpaqueCCCSetup)));
    
    /* store configuration, and copy pointers */
    memcpy(&ret->config, config, sizeof(TYPE(struct OpaqueCCCConfig)));
    ret->config.timescales = (REAL *)allocAndCopy(ret->config.timescales, sizeof(REAL) * ret->config.num_timescales);
    ret->config.angles = (REAL *)allocAndCopy(ret->config.angles, sizeof(REAL) * ret->config.num_angles);
    
    /* STEP 1: spectrogram, build windows and fft */
    ret->fft_length_half = ret->config.fft_length / 2;
    ret->fft_size = fftSize(ret->config.fft_length);
    assert((1 << ret->fft_size) == ret->config.fft_length); // should be power of 2
    
    ret->fft_setup = vDSP(create_fftsetup)(ret->fft_size, kFFTRadix2);
    
    /* allocate window */
    ret->fft_window = (REAL *)malloc(sizeof(REAL) * ret->config.fft_length * ret->config.num_timescales * 2);
    ret->signal_windowed = (REAL *)malloc(sizeof(REAL) * ret->config.fft_length * ret->config.num_timescales * 2);
    
    /* allocate temporary and output */
    ret->fft_temporary.realp = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half);
    ret->fft_temporary.imagp = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half);
    ret->fft_output.realp = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half * ret->config.num_timescales * 2);
    ret->fft_output.imagp = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half * ret->config.num_timescales * 2);
    
    /* make pointers */
    ret->p_exp = ret->fft_output;
    ret->p_der.realp = ret->fft_output.realp + ret->fft_length_half * ret->config.num_timescales;
    ret->p_der.imagp = ret->fft_output.imagp + ret->fft_length_half * ret->config.num_timescales;
    ret->p_ratio = ret->p_der;
    
    /* power array */
    ret->power = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half * ret->config.num_timescales);
    
    /* conensus array */
    /* might be able to use singles? */
    ret->consensus = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half * ret->config.num_timescales * ret->config.num_angles);
    ret->consensus_cur = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half);
    ret->consensus_pow = (REAL *)malloc(sizeof(REAL) * ret->fft_length_half);
    
    /* fill window */
    TYPE(fillFftWindow)(ret);
    
    /* setup timescale angle structures */
    ret->ta = (TYPE(struct CCTimescaleAngle) *)malloc(sizeof(TYPE(struct CCTimescaleAngle)) * ret->config.num_timescales * ret->config.num_angles);
    for (i = 0; i < ret->config.num_timescales; ++i) {
        for (j = 0; j < ret->config.num_angles; ++j) {
            TYPE(setupCCTimescaleAngle)(ret->ta + i * ret->config.num_angles + j, ret->fft_length_half, ret->config.timescales[i], ret->config.angles[j]);
        }
    }
    
    return ret;
}

void TYPE(destroyCCCSetup)(TYPE(CCCSetup) setup) {
    t_len i, j;
    
    /* free config memory */
    free(setup->config.timescales);
    free(setup->config.angles);
    
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
    for (i = 0; i < setup->config.num_timescales; ++i) {
        for (j = 0; j < setup->config.num_angles; ++j) {
            TYPE(destroyCCTimescaleAngle)(setup->ta + i * setup->config.num_angles + j);
        }
    }
    free(setup->ta);
    
    /* free memory for whole setup */
    free(setup);
}

static void TYPE(windowSamples)(const TYPE(CCCSetup) setup, const REAL *signal) {
    t_len i = 0;
    t_len signals = 2 * setup->config.num_timescales;
    
    for (i = 0; i < signals; ++i) {
        vDSP(vmul)(signal, 1, setup->fft_window + i * setup->config.fft_length, 1, setup->signal_windowed + i * setup->config.fft_length, 1, setup->config.fft_length);
    }
}

struct ConsensusContourSize TYPE(cccSizeConfig)(const TYPE(CCCConfig) config, const t_len signal_len) {
    struct ConsensusContourSize ret;
    
    // store signal length
    ret.signal_len = signal_len;
    
    // insufficient signal?
    if (signal_len < config->fft_length) {
        ret.rows = 0;
        ret.cols = 0;
        ret.bytes = 0;
        return ret;
    }
    
    const t_len ccn = config->fft_length / 2;
    const t_len ccm = 1 + (signal_len - config->fft_length) / config->fft_shift;
    
    ret.rows = ccn;
    ret.cols = ccm;
    ret.bytes = sizeof(REAL) * ccn * ccm;
    
    return ret;
}

struct ConsensusContourSize TYPE(cccSizeSetup)(const TYPE(CCCSetup) setup, const t_len signal_len) {
    struct ConsensusContourSize ret;
    
    // store signal length
    ret.signal_len = signal_len;
    
    // insufficient signal?
    if (signal_len < setup->config.fft_length) {
        ret.rows = 0;
        ret.cols = 0;
        ret.bytes = 0;
        return ret;
    }
    
    const t_len ccn = setup->config.fft_length / 2;
    const t_len ccm = 1 + (signal_len - setup->config.fft_length) / setup->config.fft_shift;
    
    ret.rows = ccn;
    ret.cols = ccm;
    ret.bytes = sizeof(REAL) * ccn * ccm;
    
    return ret;
}

/* freqs must be at least dim.rows long, times must be at least dims.cols long */
/* either can be NULL if not needed */
void TYPE(cccBins)(const TYPE(CCCSetup) setup, const struct ConsensusContourSize dim, REAL *freqs, REAL *times) {
    if (freqs) {
        REAL f = setup->config.fs / (REAL)setup->config.fft_length;
        for (t_len i = 0; i < dim.rows; ++i) {
            freqs[i] = (REAL)(i + 1) * f;
        }
    }
    
    if (times) {
        for (t_len j = 0; j < dim.cols; ++j) {
            times[j] = (REAL)(setup->fft_length_half + j * setup->config.fft_shift) / setup->config.fs;
        }
    }
}

/* signal must be setup->fft_length long */
static void TYPE(buildColumn)(const TYPE(CCCSetup) setup, const REAL *signal, REAL *consensus_contour) {
    t_len j, k;
    t_len len;
    
#ifdef TIMING
    clock_t begin, end;
    begin = clock();
#endif
    
    // window
    TYPE(windowSamples)(setup, signal);
    
    // pack samples
    vDSP(ctoz)((REAL_COMPLEX *)setup->signal_windowed, 2, &setup->fft_output, 1, setup->fft_length_half * setup->config.num_timescales * 2);
    
    // calculate
    vDSP(fftm_zript)(setup->fft_setup, &setup->fft_output, 1, (t_stride)setup->fft_length_half, &setup->fft_temporary, setup->fft_size, setup->config.num_timescales * 2, FFT_FORWARD);
    
    // length
    len = setup->fft_length_half * setup->config.num_timescales;
    
    // calculate power
    vDSP(zvabs)(&setup->fft_output, 1, setup->power, 1, len);
    vDSP(vsdiv)(setup->power, 1, &TYPE(s_two), setup->power, 1, len);
    
    // calculate ratio
    vDSP(zvdiv)(&setup->p_exp, 1, &setup->p_der, 1, &setup->p_ratio, 1, len);
    
    // calculate contours for each timescale and angle
    for (j = 0; j < setup->config.num_timescales; ++j) {
        for (k = 0; k < setup->config.num_angles; ++k) {
            TYPE(ingestCCTimescaleAngle)(&setup->ta[j * setup->config.num_angles + k], setup->p_ratio.realp + j * setup->fft_length_half, setup->p_ratio.imagp + j * setup->fft_length_half, setup->consensus + (j * setup->config.num_angles + k) * setup->fft_length_half);
        }
    }
    
    // look for consensus
    for (j = 0; j < (setup->config.num_timescales - 1); ++j) {
        if (setup->config.pow_weight) {
            vDSP(vclr)(setup->consensus_pow, 1, setup->fft_length_half);
        }
        
        for (k = 0; k < setup->config.num_angles; ++k) {
            // start with consensus / angle
            memcpy(setup->consensus_cur, setup->consensus + (j * setup->config.num_angles + k) * setup->fft_length_half, sizeof(REAL) * setup->fft_length_half);
            // add next sigma, same angle
            vDSP(vadd)(setup->consensus_cur, 1, setup->consensus + ((j + 1) * setup->config.num_angles + k) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            // add same sigma, previous angle (wrap around)
            vDSP(vadd)(setup->consensus_cur, 1, setup->consensus + (j * setup->config.num_angles + (0 == k ? setup->config.num_angles - 1 : k - 1)) * setup->fft_length_half, 1, setup->consensus_cur, 1, setup->fft_length_half);
            
            // add to contours
            if (setup->config.pow_weight) {
                TYPE(addConsensusToContours)(setup->fft_length_half, setup->consensus_cur, setup->consensus_pow);
            }
            else {
                TYPE(addConsensusToContours)(setup->fft_length_half, setup->consensus_cur, consensus_contour);
            }
        }
        
        if (setup->config.pow_weight) {
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
        TYPE(buildColumn)(setup, signal + i * setup->config.fft_shift, consensus_contours + i * setup->fft_length_half);
    }
}

