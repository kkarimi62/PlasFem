

#ifndef COMPUTE_X_H
#define COMPUTE_X_H

#include "compute_sxy.h"
#include "pure_shear.h"
class Fem;
class ComputeX
{
	public:
		ComputeX( Fem * );
		~ComputeX();
		void GetX();
		FILE *outFile;		
		bool cx;

	private:

		unsigned int kount, init, duration;
		double a, b, c, ds;
		PureShear *ps;
		ComputeSxy *cxy;	
		
};

#endif

