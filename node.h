
#ifndef NODE_H
#define NODE_H
#include "memory.h"

class Node
{
	public:
		Node();
		~Node();

		unsigned int npoin, nsvab;
		double **coord, **cordi, *dispt, *dispt0, *dldis, **coordUnWrapped, **cordiUnWrapped;
		double *velot, *veloi, *presVeloc, *extForce;
		double *fintl, *absfl;
		double *mass;
		int *fixID, *rigidID;
		double effMass, fx, fy;	
	private:
		Memory* memPtr;
};
#endif
