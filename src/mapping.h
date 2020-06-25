#ifndef MAPPING_H
#define MAPPING_H
#include "memory.h"
#include "matrix.h"

class Wrap
{
	public:
		Wrap( double * );
		~Wrap();
	
		void GetDimensionlessCords(  unsigned int,  double *, double * );
		void GetCords(  unsigned int,  double *, double * );
		void CoordWrapper(  unsigned int, double ** );
	private:
		int ndime;
		double **h;
		double **h_inv;
		double *lohiPtr;
		void set_h();
		Memory *memPtr;
		Matrix *matPtr;
};
#endif 
