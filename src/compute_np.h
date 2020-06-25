#ifndef COMPUTE_NP_H
#define COMPUTE_NP_H
#include "stdio.h"

class Fem;
class ComputeNP
{
	public:
		ComputeNP( Fem * );
		~ComputeNP();

		void GetActMap( unsigned int );
		bool cnp;
		FILE *outFile;
	private:
		 Fem *femPtr;
};

#endif
