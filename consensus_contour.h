/**
 * 
 */

#include <Accelerate/Accelerate.h>

/* structure to hold dimensions */
struct ConsensusContourSize
{
    unsigned long signal_len;
    unsigned long rows;
    unsigned long cols;
    unsigned long bytes;
};

/* FLOAT */

/* setup */
typedef struct OpaqueCCCSetup *CCCSetup;
CCCSetup createCCCSetup(unsigned long fft_length, unsigned long fft_overlap, float fs, bool pow_weight);
void destroyCCCSetup(CCCSetup setup);

/* sizing */
struct ConsensusContourSize cccSize(const CCCSetup setup, const unsigned long signal_len);

/* actual algorithm */
void cccColumn(const CCCSetup setup, const float *signal, float *contour);
void cccSpectogram(const CCCSetup setup, const struct ConsensusContourSize dim, const float *signal, float *consensus_contours);

/* DOUBLE */

/* setup */
typedef struct OpaqueCCCSetupD *CCCSetupD;
CCCSetupD createCCCSetupD(unsigned long fft_length, unsigned long fft_overlap, double fs, bool pow_weight);
void destroyCCCSetupD(CCCSetupD setup);

/* sizing */
struct ConsensusContourSize cccSizeD(const CCCSetupD setup, const unsigned long signal_len);

/* actual algorithm */
void cccColumnD(const CCCSetupD setup, const double *signal, double *contour);
void cccSpectogramD(const CCCSetupD setup, const struct ConsensusContourSize dim, const double *signal, double *consensus_contours);
