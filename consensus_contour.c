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

static t_len fftSize(const t_len length) {
    double size = ceil(log2((double)length));
    return (t_len)size;
}

/* float */

#define REAL float
#define REAL_SPLIT_COMPLEX DSPSplitComplex
#define REAL_COMPLEX DSPComplex
#define REAL_FFT_SETUP FFTSetup

#define vDSP(FN) vDSP_ ## FN
#define CMATH(FN) FN ## f
#define TYPE(FN) FN

#include "_consensus_contour.c"

#undef REAL
#undef REAL_SPLIT_COMPLEX
#undef REAL_COMPLEX
#undef REAL_FFT_SETUP
#undef vDSP
#undef CMATH
#undef TYPE

/* double */

#define REAL double
#define REAL_SPLIT_COMPLEX DSPDoubleSplitComplex
#define REAL_COMPLEX DSPDoubleComplex
#define REAL_FFT_SETUP FFTSetupD

#define vDSP(FN) vDSP_ ## FN ## D
#define CMATH(FN) FN
#define TYPE(FN) FN ## D

#include "_consensus_contour.c"

#undef REAL
#undef REAL_SPLIT_COMPLEX
#undef REAL_COMPLEX
#undef REAL_FFT_SETUP
#undef vDSP
#undef CMATH
#undef TYPE
