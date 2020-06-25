

#ifndef COMPUTE_DS_H
#define COMPUTE_DS_H
//#include "thermo_style.h"
#include "compute_sxy.h"
#include "dump.h"
#include "quasi_static.h"
#include "tracer.h"
//#include "compute_xl.h"
class Fem;
class ComputeDS
{
	public:
		ComputeDS( Fem * );
		~ComputeDS();
		void GetAvalanch( unsigned int );
		FILE *outFile;		
		bool cds;
		unsigned int n, navch;

	private:

		unsigned int kount, init, duration, intDur;
		double a, b, c, ds;
		unsigned int ta, tb, tc;
		ComputeSxy *cxy;	
		PureShear *ps;
//		ComputeXL *cxl;	
		Fem *femPtr;	
		Dump *dp;
		QuasiStatic *qs;
		Tracer *tr;
};

#endif

