#ifndef DUMP_MODIFY_H
#define DUMP_MODIFY_H
#include "stdio.h"
//#include "memory.h"
//#include "node.h"
//#include "domain.h"
//#include "elem.h"
//#include "time_step.h"
//#include "fix_deform.h"
//#include "tracer.h"

class Fem;
class DumpModify
{
	public:
		DumpModify( Fem * );
		~DumpModify();

		int nDump;
//		int *n;
		int *narg;
		int **args;
		bool *exist;
//		char **output;
//		Fem *femPtr;
//		void Process(  unsigned int,  int );
 //		void Init();
	private:
//		FILE **outptFile;
//		Memory *memPtr;
//		Node *node;
//		Elem *elem;		
//		Domain *domain;
//		TimeStep *ts;
//		FixDeform *fd;  // shear 
//		double **data;
//		unsigned int N;
//		Tracer *trace;
		Fem *femPtr;
//		double **cordg;
//		GaussCoord *gcObj;
};

#endif

