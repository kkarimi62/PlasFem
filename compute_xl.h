//--- local distance to sigm_y
#ifndef COMPUTE_XL_H
#define COMPUTE_XL_H
#include "elem.h"
class Fem;
class ComputeXL
{
	public:
		ComputeXL( Fem * );
		~ComputeXL();
		void GetXL( bool );
		FILE *outFile;		
		bool cxl;
		unsigned int nfreq;
	private:
		Elem *elem;	
		double Invart( double * );	
		inline double Smean(  double * );

};

#endif

