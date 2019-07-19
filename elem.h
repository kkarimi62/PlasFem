#ifndef ELEM_H
#define ELEM_H
#include "memory.h"
#include "moduli.h" // for material parameters
#include "rand_numbers.h"
#include "domain.h"
#include "matrix.h"


class Fem;
class GaussCoord;
class Elem
{
	public:
		Elem( Fem * );
		~Elem();

		int nnode, ngaus, nstre;
		unsigned int nelem, ntype, nevab;
		double meanGmodu, meanKmodu, meanRecoveryStrain;
		unsigned int *ltype;
//		double *mass;
		unsigned int **lnods;
		double *kmodu, *gmodu, *uniax,  *uniax0, *strny, *frict0, *frict, *hards, *infrc, *inunx;
		double *rsidu, *rsidf, *reduc, *redup, *reduf, *rsidp; //--- residual uniax and damage parameter
		double *pcrit, *pcrit0, *mnunx;
		double *mufst, *musnd, *tprin; //--- anisotroic moduli in damage model
		double **stran, **strni, **strst, **strsi, *epstn, **stranRate, **plasticStrain, **plasticStrain0, *vol;
		double **strnt; //--- store total strain
		double **strrs;
		double *fyild;
		double *fyild0;
		double **cordg;
		double *presd, *sresd, *sxyrd;
		bool *resid;
		int *actField, *actField0, *actFieldAccumulated, *actFieldAccumulated0;

		void Init();
		void set(  unsigned int, const double, const double ); // sets shear modulus etc.		
		void GAUSSQ( double **, double * ); // quadrature points
		void Jacob2( double **, double &,  double,  double,  double **,  unsigned int );
		void Shape2( double *, double **,  double,  double );
		void Blarge( double **,  double ** );
		void ksiTOxy( double *,  double *,  double ** );
	private:
		Fem *femPtr;
		Domain *domain;
		Memory *memPtr;
		Matrix *matPtr;
		inline double det(  double ** );
		inline void inv( double **,  double ** );
		RandNumbers *randOBJ, *randOBJ2nd; // for random numbers	
//		double *shape, **deriv;
//		double **XJACM, **XJACI;
		double minE;
		GaussCoord *gcObj;
		
};

#endif
