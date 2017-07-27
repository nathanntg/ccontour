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
	
	const struct ConsensusContourSize dim = cccSize(SIGNAL_LEN);
	out = malloc(dim.bytes);
	
	for (i = 0; i < iter; ++i) {
		for (j = 0; j < SIGNAL_LEN; ++j) {
			signal[j] = (double)rand() / RAND_MAX;
		}
		
		begin = clock();
		ccc(&signal[0], SIGNAL_LEN, SIGNAL_FS, out);
		end = clock();
		
		tm += (double)(end - begin) / CLOCKS_PER_SEC;
	}
	
	printf("Avg Time: %.3fms\n", tm * 1000 / (double)iter);
	
	free(out);
	
	return 0;
}

