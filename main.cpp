#include "fem.h"
#include <iostream>
#include <omp.h>

using std::cout;
//--- driver program 
int main()
{
	double start_time = omp_get_wtime();
	Fem *femPtr = new Fem;
	femPtr->input->file();
	femPtr->Loop();
	fprintf( femPtr->logFile, "Loop time of %g s for %d steps with %d nodes\n", omp_get_wtime() - start_time, femPtr->run->nstep, femPtr->node->npoin );
	delete femPtr;
	return 0;
}

