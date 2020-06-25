#ifndef GAUSS_COORD_H
#define GAUSS_COORD_H
#include "memory.h"
#include "node.h"
#include "domain.h"
#include "elem.h"
#include "mapping.h"

class Fem;
class GaussCoord
{
	public:
		GaussCoord( Fem * );
		~GaussCoord();
		void GetCoord( double ** );
		void Init();
	private:
//		Fem *femPtr;
		void GetNodalCords( double**,  unsigned int,  int );
		void Distance(  double **,  double *, double *, double *, int * ); // find true distance
		Memory *memObj;
		Node *node;
		Elem *elem;		
		Domain *domain;
		unsigned int kgaus;
		double **posgp; // position of gauss points
		double *weigp; // quadrature points: weight
		double **lcord;  // xyz in the initial frame
		double lohi[ 6 ];
		Wrap *mpObj;
		unsigned int nodei;
		double **Cords; // array to store xyz of each pair of particles
		double *r, *theta;//distance between two nodes and angle of the connecting line
		int *cross_pbc; // cross PBC?
		double xi, yi, xj, yj, xij, yij, rsq; // xyz
		double s0i, s1i, s0j, s1j, s0_ij, s1_ij;
		double Xj, Yj;
		double **dimensioless_cords; // declare an array to hold s0, s1 
		unsigned int initCount;
};

#endif
