/**
 * 
 */
 
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "consensus_contour.h"

#define SIGNAL_LEN 40000
#define SIGNAL_FS 44100

static void profileSingle(const int iter) {
    int i, j, r;
    float signal[SIGNAL_LEN];
    float *out;
    double tm = 0;
    
    clock_t begin, end;
    
    CCCSetup ccc_setup = createCCCSetup(1024, 1005, SIGNAL_FS, true);
    
    const struct ConsensusContourSize dim = cccSize(ccc_setup, SIGNAL_LEN);
    out = calloc(dim.rows * dim.cols, sizeof(float));
    
    for (i = 0; i < iter; ++i) {
        for (j = 0; j < SIGNAL_LEN; ++j) {
            r = rand();
            signal[j] = (float)r / RAND_MAX;
        }
        
        begin = clock();
        cccSpectogram(ccc_setup, dim, &signal[0], out);
        end = clock();
        
        tm += (double)(end - begin) / CLOCKS_PER_SEC;
    }
    
    printf("Avg Time: %.3fms\n", tm * 1000 / (double)iter);
    
    free(out);
    destroyCCCSetup(ccc_setup);
}

static void profileDouble(const int iter) {
    int i, j, r;
    double signal[SIGNAL_LEN];
    double *out;
    double tm = 0;
    
    clock_t begin, end;
    
    CCCSetupD ccc_setup = createCCCSetupD(1024, 1005, SIGNAL_FS, true);
    
    const struct ConsensusContourSize dim = cccSizeD(ccc_setup, SIGNAL_LEN);
    out = calloc(dim.rows * dim.cols, sizeof(double));
    
    for (i = 0; i < iter; ++i) {
        for (j = 0; j < SIGNAL_LEN; ++j) {
            r = rand();
            signal[j] = (double)r / RAND_MAX;
        }
        
        begin = clock();
        cccSpectogramD(ccc_setup, dim, &signal[0], out);
        end = clock();
        
        tm += (double)(end - begin) / CLOCKS_PER_SEC;
    }
    
    printf("Avg Time: %.3fms\n", tm * 1000 / (double)iter);
    
    free(out);
    destroyCCCSetupD(ccc_setup);
}

// int argc, const char* argv[]
int main() {
	int iter = 10;
    printf("** FLOAT **\n");
    profileSingle(iter);
    printf("** DOUBLE **\n");
    profileDouble(iter);
	return 0;
}

