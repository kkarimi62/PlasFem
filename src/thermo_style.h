#ifndef THERMO_STYLE_H
#define THERMO_STYLE_H
#include "memory.h"
#include "node.h"
#include "domain.h"
#include "elem.h"
#include "time_step.h"
#include "fix_deform.h"
#include "compute_ke.h"
#include "compute_sxy.h"


class Fem; //forward declaration of parent class
class ThermoStyle
{
	public:
//		ThermoStyle( char **, Fem *, unsigned int, int );
		ThermoStyle( Fem * );
		~ThermoStyle();
		void Init();
		void Process( unsigned int );

		ComputeSxy *cSxyObj;
		ComputeKE *cKEobj;
		unsigned int n;
		int narg;
//		Fem *femPtr;
		char **args;		

	private:
//		Memory *memPtr;
		Node *node;
		Elem *elem;		
		Domain *domain;
		TimeStep *ts;
		FixDeform *fd;  // shear 
		unsigned int kount;
		FILE *output;
};

#endif
