#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/count.h>

int main() {
    /**
     * \brief armando2D v2.0
     *
     * An SPH code for non stationary fluid dynamics.
     * This is the reviewed and improved C version of Armando v1.0
     * developed at CERN in 2008
     *
     * \date May 2, 2012
     * \author Luca Massidda
     */
    
    int i, pout;
    int hhash[10];
    int *dhash;
    
    cudaMalloc((void**) &dhash, (10 * sizeof(int)));
	thrust::device_ptr<int> thash(dhash);
    thrust::sequence(thash, thash +10, 0);
    
	pout = thrust::count(thash, thash + 5, 3);
	printf("%d\n\n", pout);
	
    thrust::copy(thash, thash +10, hhash);
    
    for (i = 0; i < 10; i++) {
		printf("%d %d \n", i, hhash[i]);
	}
	
    return 0;
}
