/**
 * 
 */

#include <Accelerate/Accelerate.h>

typedef vDSP_Stride t_stride;
typedef vDSP_Length t_len;

struct ConsensusContourSize
{
    t_len rows;
    t_len cols;
    t_len bytes;
};

struct ConsensusContourSize cccSize(const t_len signal_len);
void ccc(const double *signal, const t_len signal_len, const double signal_fs, double *consensus_contours);
