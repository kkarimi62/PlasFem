#ifndef TRACER_H
#define TRACER_H
#include "memory.h"
#include "node.h"
#include "elem.h"
#include "domain.h"
#include "mapping.h"
#include "time_step.h"
#include "fix_deform.h"
#include "pure_shear.h"

class Fem;
class Tracer
{
	public:
		Tracer( Fem * );
		~Tracer();	
		void TracerUpdate();
		void Init();
		void UpdateCords();
		void ComputeMSD( unsigned int );

		double **coord, *veloi, *velot, *dispt, *tdisp, *tdsnf, *xdisp;
		bool msd;
		unsigned int nfreq;
		int mrows, ncols; //--- don't declare as unsigned!
		unsigned int *cellTriangCount;
		unsigned int **lelms;
		FILE *outFile, *ofile, *ofile2, *ofile3;		
	private:
		bool PointInsideTriag( double *, double *, double &, double & );
		void GetDisp( double *,  unsigned int, bool );
		void GetNodalCords( double *,  unsigned int, int pbc = 1, int initial = 0 ); // nodal coordinates
		void GetVeloc( double *, double *,  unsigned int );
		void Distance(  double *,  int, double &, double &, int & ); // find true distance
		void FindElem( unsigned int, unsigned int &, double &, double & );
		void GetExispEtasp( double *, double *, double &, double & );
		void Nod2Gas( double *, double **, double &, double & );
		void Nod2Gas( double *, double *, double &, double & );
		void ApplyHomo();
		void inline MeanSq( double*, double* );
		Fem *femPtr;
		Memory *memObj;
		Node *node;
		Elem *elem;
		Domain *domain;
		TimeStep *ts;   // time step
		FixDeform *fd;  // shear 
		PureShear *ps;
//---
		Wrap *mpObj[ 2 ];	
//---
		double **cordPtr[ 2 ];
		double **lohi;
//---
//		unsigned int kount, init, duration, telas;
//		const double rsqTl = 1.0e-08;
		RandNumbers *randOBJ;
};

#endif
