/**
 * 
 */

#include <Accelerate/Accelerate.h>

/* setup */
typedef struct OpaqueCCCSetup *CCCSetup;
CCCSetup createCCCSetup(unsigned long fft_length, unsigned long fft_overlap, double fs, bool pow_weight);
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
void cccColumn(const CCCSetup setup, const double *signal, double *contour);
void ccc(const CCCSetup setup, const struct ConsensusContourSize dim, const double *signal, double *consensus_contours);

