#ifndef COMPUTE_DT_H
#define COMPUTE_DT_H
#include "stdio.h"
#include "node.h"

class Fem;
class ComputeDT
{
	public:
		ComputeDT( Fem * );
		~ComputeDT();

		void GetDuration( unsigned int );
		bool cdt;
		FILE *outFile;
		unsigned int tmin, tmax, nfrq, count;
		double sigma0, *dsigma, *epdot, *fx, *fy, *kin, svelo, sigmp, *velos;
		unsigned int *time;
	private:
		Node *node;
};

#endif
