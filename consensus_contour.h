/**
 * Consensus contour
 */

#include <stdbool.h>

/* COMMON */

/* structure to hold dimensions */
struct ConsensusContourSize
{
    unsigned long signal_len;
    unsigned long rows;
    unsigned long cols;
    unsigned long bytes;
};

/* FLOAT */

/* configuration */
typedef struct OpaqueCCCConfig *CCCConfig;
CCCConfig createCCCConfig(void); // creates a default configuration
void cccConfigSetFFTLength(CCCConfig config, unsigned long fft_length);
void cccConfigSetFFTOverlap(CCCConfig config, unsigned long fft_overlap);
void cccConfigSetWeightByPower(CCCConfig config, bool pow_weight);
void cccConfigSetSampleRate(CCCConfig config, float fs);
void cccConfigSetTimescales(CCCConfig config, unsigned long num_timescales, const float timescales[]);
void cccConfigSetAngles(CCCConfig config, unsigned long num_angles, const float angles[]);
void destroyCCCConfig(CCCConfig config);

/* setup */
typedef struct OpaqueCCCSetup *CCCSetup;
CCCSetup createCCCSetup(const CCCConfig config);
void destroyCCCSetup(CCCSetup setup);

/* sizing */
struct ConsensusContourSize cccSize(const CCCSetup setup, const unsigned long signal_len);

/* actual algorithm */
void cccColumn(const CCCSetup setup, const float *signal, float *contour);
void cccSpectrogram(const CCCSetup setup, const struct ConsensusContourSize dim, const float *signal, float *consensus_contours);

/* DOUBLE */

/* configuration */
typedef struct OpaqueCCCConfigD *CCCConfigD;
CCCConfigD createCCCConfigD(void); // creates a default configuration
void cccConfigSetFFTLengthD(CCCConfigD config, unsigned long fft_length);
void cccConfigSetFFTOverlapD(CCCConfigD config, unsigned long fft_overlap);
void cccConfigSetWeightByPowerD(CCCConfigD config, bool pow_weight);
void cccConfigSetSampleRateD(CCCConfigD config, double fs);
void cccConfigSetTimescalesD(CCCConfigD config, unsigned long num_timescales, const double timescales[]);
void cccConfigSetAnglesD(CCCConfigD config, unsigned long num_angles, const double angles[]);
void destroyCCCConfigD(CCCConfigD config);

/* setup */
typedef struct OpaqueCCCSetupD *CCCSetupD;
CCCSetupD createCCCSetupD(const CCCConfigD config);
void destroyCCCSetupD(CCCSetupD setup);

/* sizing */
struct ConsensusContourSize cccSizeD(const CCCSetupD setup, const unsigned long signal_len);

/* actual algorithm */
void cccColumnD(const CCCSetupD setup, const double *signal, double *contour);
void cccSpectrogramD(const CCCSetupD setup, const struct ConsensusContourSize dim, const double *signal, double *consensus_contours);
