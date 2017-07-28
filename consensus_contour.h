/**
 * 
 */

#include <Accelerate/Accelerate.h>

#ifndef PRECISION
#define PRECISION 1
#endif

#if PRECISION == 2
typedef double t_real;
#else
typedef float t_real;
#endif

/* setup */
typedef struct OpaqueCCCSetup *CCCSetup;
CCCSetup createCCCSetup(unsigned long fft_length, unsigned long fft_overlap, t_real fs, bool pow_weight);
void destroyCCCSetup(CCCSetup setup);

/* sizing */
struct ConsensusContourSize
{
    unsigned long signal_len;
    unsigned long rows;
    unsigned long cols;
    unsigned long bytes;
};
struct ConsensusContourSize cccSize(const CCCSetup setup, const unsigned long signal_len);

/* actual algorithm */
void cccColumn(const CCCSetup setup, const t_real *signal, t_real *contour);
void ccc(const CCCSetup setup, const struct ConsensusContourSize dim, const t_real *signal, t_real *consensus_contours);

