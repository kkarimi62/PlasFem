#ifndef FEM_H
#define FEM_H
#include "node.h"
#include "stdio.h"
#include "elem.h"
#include "domain.h"
#include "time_step.h"
#include "fix_deform.h"
#include "fix_viscous.h"
#include "run.h"
#include "thermo_style.h"
#include "dump.h"
#include "dump_modify.h"
#include "input.h"
#include "force.h"
#include "mapping.h"
#include "restart.h"
#include "fix_veloc.h"
#include "fix_force.h"
#include "quasi_static.h"
#include "yield_crit.h"
#include "yield_rule.h"
#include "weak_rule.h"
#include "recovery.h"
#include "pure_shear.h"
#include "tracer.h"
#include "compute_x.h"
#include "compute_xl.h"
#include "compute_ds.h"
#include "compute_dt.h"
#include "compute_np.h"
#include <Eigen/Sparse> //--- store hessian
#include "nthreads.h"

using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;
using Eigen::VectorXd;
class Fem
{
	public:
		Fem();
		~Fem();

		FILE *infile; // standard input
		FILE *logFile; // standard input
		Node *node; // node-based quantities
		Elem *elem; // element-based variables
		Domain *domain; // simulation box
		TimeStep *ts;   // time step
		FixDeform *fd;  // shear 
		FixViscous *fv;  // viscosity 
		Run *run;
		ThermoStyle *ths;	
		Dump *dp;
		DumpModify *dpmd;
		Input *input; // loading input from hard disk
		Force *force;  // inter-particle forces
		Restart *rs;
		FixVeloc *fvl;  // shear 
		FixForce *ff;
		QuasiStatic *qs;
		YieldCrit *yc;
		YieldRule *yr;
		WeakRule *wr;
		Recovery *rv;
		PureShear *ps;
		ComputeDS *cds;
		ComputeDT *cdt;
		ComputeNP *cnp;
		ComputeX *cx;
		ComputeXL *cxl;
		Tracer *tr;
		ConjugateGradient< SparseMatrix <double> > *cg;
		void Loop();
	
	private:
		void Verlet( double *, double *, double, double );
		void ApplyHOMO();
		void Switch();
		void Update();
		void Init();
		void SetVeloc();
		void SetForce();
		bool Iterate( unsigned int & );
		double fmax, *effForce, cdrag0, *velol;
		unsigned int kount, nblck, nblck0, duration;
		VectorXd *b, *x; //--- dynamic	
};

#endif
