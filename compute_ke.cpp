#include "compute_ke.h"

ComputeKE::ComputeKE()// double *mass, double *v, unsigned int nsvab ):
//xmass( mass ), veloc( v ), n( nsvab ), ndime( 2 )
{}

void ComputeKE::Init(  double *mass,  double *v,  unsigned int nsvab,  int ndim )
{
	xmass = mass;
	veloc = v;
	n = nsvab;
	ndime = ndim;
}

double ComputeKE::GetKE()
{
	ke = 0.0;
//	#pragma omp parallel
	{
//	#pragma omp for
	for( int isvab = 0; isvab < n; isvab++ )
	{
		ipoin = isvab / ndime;
//		#pragma omp atomic 
		ke += xmass[ ipoin ] * veloc[ isvab ] * veloc[ isvab ];
	}
	}
	return 0.5 * ke;
}
