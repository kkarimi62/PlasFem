#ifndef FORCE_H
#define FORCE_H
#include "memory.h"
#include "matrix.h"
#include "node.h"
#include "elem.h"
#include "domain.h"
#include "fix_viscous.h"
#include "fix_deform.h"
#include "time_step.h"
#include "mapping.h"
#include "yield_crit.h"
#include "yield_rule.h"
#include "weak_rule.h"
#include "recovery.h"
#include "pure_shear.h"
#include <Eigen/Sparse> //--- store hessian
#include "nthreads.h" 


using Eigen::SparseMatrix;

class Fem;
class Force
{
	public:
		Force( Fem * );
		~Force();
		void Init();
		void ComputeForce( int, unsigned int );
		double Invart(  double * ); //, int );
		double Smean( double * );
	private:
		void GetDisp( double *,  unsigned int, int );//, int );
		void GetVeloc( double *,  unsigned int ); //, int );
		void GetNodalCords( double *,  unsigned int, int pbc = 1, int initial = 0 );//, int = 0 ); // nodal coordinates
		void Distance(  double *,  int, double &, double &, int & ); //, int ); // find true distance
		inline void Modulus( double *, double, double, int );
		void AnisoMod( double *, double, double, double, double );
		void StrPrinDir( double *, double *, double * );
		inline double DotPr( double *, double *, unsigned int );
		void ReturnStress( double *, double *, double *, int *, unsigned int, unsigned int, double *, double *, double * );
		inline void Devia( double *,  double * ); //, int );
		inline void DeviaStrain( double *,  double * );
		inline double EffectiveStrain(  double * );//, int );
		Fem *femPtr;
		Memory *memObj;
		Matrix *matPtr;
		Node *node;
		Elem *elem;
		Domain *domain;
		FixViscous *fv;  // viscosity 
		FixDeform *fd;
		TimeStep *ts;   // time step
		YieldCrit *yc;
		YieldRule *yr;
		WeakRule *wr;
		Recovery *rv;
		PureShear *ps;
		

		SparseMatrix<double> *stiff;	
		Wrap *mpObj[ 2 ];	
		double **cordPtr[ 2 ],
		       **lohi;	
		double ***bmatx[ nThreads ], ***bt[ nThreads ]; //--- b-matrix, transpose of b-matrix
		double visc_coeff, second_visc_coeff;
		double *val;
		FILE *outpt;			  
//		int xxx;
};

#endif
