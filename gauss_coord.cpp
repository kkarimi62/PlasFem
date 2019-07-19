//#include "Python.h"
#include "gauss_coord.h"
#include "fem.h"

GaussCoord::GaussCoord( Fem *femObj ):
memObj( new Memory ), elem( femObj->elem ), domain( femObj->domain ),
node( femObj->node ), initCount( 0 )
{}


GaussCoord::~GaussCoord()
{
	if( initCount > 0 )
	{
		delete [] weigp;
		memObj->Delete2ndMat( posgp, elem->ngaus );
		memObj->Delete2ndMat( lcord, elem->nnode );
		memObj->Delete2ndMat( Cords, 2 );
//		memObj->Delete2ndMat( dimensioless_cords, 2 );
		delete r;
		delete theta;
		delete cross_pbc;
		delete mpObj;
	}
	delete memObj;
}

void GaussCoord::Init()
{
	initCount++;
	memObj->doubleVec( weigp, elem->ngaus );
	memObj->doubleMat( posgp, elem->ngaus, domain->ndime );
	memObj->doubleMat( lcord, elem->nnode, domain->ndime ); // allocation
	memObj->doubleMat( Cords, 2, domain->ndime );
//	memObj->doubleMat( dimensioless_cords, 2, domain->ndime ); // allocate memory
	r = new double; //distance between two nodes and angle of the connecting line
	theta = new double;
	cross_pbc = new int; // cross PBC?

	lohi[ 0 ] = domain->xlo0;
	lohi[ 1 ] = domain->xhi0;
	lohi[ 2 ] = domain->ylo0;
	lohi[ 3 ] = domain->yhi0;
	lohi[ 4 ] = domain->xy0;
	mpObj = new Wrap( lohi ); // create Wrap obj
}

void GaussCoord::GetCoord( double **output )
{
	elem->GAUSSQ( posgp, weigp ); // gauss points: pos and weight

	// --- element loop starts
	kgaus = 0;
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{ 
		GetNodalCords( lcord, ielem, 1 ); // pbc( 1 )
		// loop over quad points starts
		for( int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			elem->ksiTOxy( output[ kgaus ], posgp[ igaus ], lcord );
			kgaus++;
		}	
	}
	//--- wrap cordg
	mpObj = new Wrap( lohi );
	mpObj->CoordWrapper( elem->ngaus * elem->nelem, output );
}
//-----------------------------------------------------------------------------------------
void GaussCoord::GetNodalCords( double **elcod,  unsigned int ielem, int pbc )
{
// --- returns elcod: nodal coordinates of the element
	for( int inode = 0; inode < elem->nnode; inode++ ) // loop over nodes
	{
		nodei = elem->lnods[ ielem ][ inode ]; // node global id
		for( int idime = 0; idime < domain->ndime; idime++ ) // assign coordantes
			elcod[ inode ][ idime ] = node->coord[ nodei ][ idime ];
		if( pbc && inode != 0 ) // --- returns xyz relative to node0 
		{
			Xj = elcod[ 0 ][ 0 ]; 
			Yj = elcod[ 0 ][ 1 ];
			for( int idime = 0; idime < domain->ndime; idime++ )
			{
				Cords[ 0 ][ idime ] = elcod[ inode ][ idime ];
				Cords[ 1 ][ idime ] = elcod[ 0 ][ idime ];
			}
			Distance( Cords, lohi, r, theta, cross_pbc ); // find the true distance
			Xj += ( *r ) * cos( *theta );
			Yj += ( *r ) * sin( *theta );
			elcod[ inode ][ 0 ] = Xj;
			elcod[ inode ][ 1 ] = Yj;
		}
	}
}
//-----------------------------------------------------------------------------------------------
void GaussCoord::Distance(  double **coord,  double *box_dim, double *r, double *theta, int *cross_pbc )
{
	// --- get s0, s1

//	printf( "xhi = %g\n", box_dim[ 1 ] );
	double tmpCord[ 2 * domain->ndime ];
	tmpCord[ 0 ] = coord[ 0 ][ 0 ];
	tmpCord[ 1 ] = coord[ 0 ][ 1 ];
	tmpCord[ 2 ] = coord[ 1 ][ 0 ];
	tmpCord[ 3 ] = coord[ 1 ][ 1 ];

	double dimensiolessTmp[ 2 * domain->ndime ];
	mpObj->GetDimensionlessCords( 2, tmpCord, dimensiolessTmp ); // output dimension...

	s0i = dimensiolessTmp[ 0 ];
	s1i = dimensiolessTmp[ 1 ];
	s0j = dimensiolessTmp[ 1 * domain->ndime ];
	s1j = dimensiolessTmp[ 1 * domain->ndime + 1 ];
	s0_ij = s0j - s0i; // dimensionless distance
	s1_ij = s1j - s1i;

	// --- translate point j to find the proper distance rij
	*cross_pbc = 0;
        if( s0_ij > 0.5 ){ 
		s0j -= 1.0;
		*cross_pbc = 1;
	}
        else if ( s0_ij < - 0.5 ){ 
		s0j += 1.0;
		*cross_pbc = 1;
        }
	if( s1_ij > 0.5 ){
		s1j -= 1.0;
		*cross_pbc = 1;
        }
	else if( s1_ij < - 0.5 ){
		s1j += 1.0;
		*cross_pbc = 1;
	}
	// --- find new coordinates
	dimensiolessTmp[ 0 ] = s0i;
	dimensiolessTmp[ 1 ] = s1i;
	dimensiolessTmp[ 1 * domain->ndime ] = s0j;
	dimensiolessTmp[ 1 * domain->ndime + 1 ] = s1j;

	mpObj->GetCords( 2, dimensiolessTmp, tmpCord ); // --- update coord
	coord[ 0 ][ 0 ] = tmpCord[ 0 ];
	coord[ 0 ][ 1 ] = tmpCord[ 1 ];
        coord[ 1 ][ 0 ] = tmpCord[ 2 ];
        coord[ 1 ][ 1 ] = tmpCord[ 3 ];

	xi = coord[ 0 ][ 0 ];
	yi = coord[ 0 ][ 1 ];
        xj = coord[ 1 ][ 0 ];
        yj = coord[ 1 ][ 1 ];
        xij = xi - xj;
        yij = yi - yj;
        rsq = xij * xij + yij * yij;
        *r = sqrt( rsq );
	*theta = atan2( yij, xij );
}
	
