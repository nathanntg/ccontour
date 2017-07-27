/**
 * 
 */
 
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "consensus_contour.h"

#define SIGNAL_LEN 40000
#define SIGNAL_FS 44100

int main(int argc, const char* argv[]) {
	int iter = 10;
	int i, j;
	double signal[SIGNAL_LEN];
	double *out;
	double tm = 0;
	
	clock_t begin, end;
	
	CCCSetup ccc_setup = createCCCSetup(1024, 1005, SIGNAL_FS, true);
	
	const struct ConsensusContourSize dim = cccSize(ccc_setup, SIGNAL_LEN);
	out = calloc(dim.rows * dim.cols, sizeof(double));
	
	for (i = 0; i < iter; ++i) {
		for (j = 0; j < SIGNAL_LEN; ++j) {
			signal[j] = (double)rand() / RAND_MAX;
		}
		
		begin = clock();
		ccc(ccc_setup, dim, &signal[0], out);
		end = clock();
		
		tm += (double)(end - begin) / CLOCKS_PER_SEC;
	}
	
	printf("Avg Time: %.3fms\n", tm * 1000 / (double)iter);
	
	free(out);
	destroyCCCSetup(ccc_setup);
	
	return 0;
}

