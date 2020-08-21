#include "tracer.h"
#include "fem.h"
#include "stdlib.h"
#include <cassert>
#include "matrix.h" 
//#include <random>

Tracer::Tracer( Fem *femObj ):
femPtr( femObj ), memObj( new Memory ),
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ),
ts( femObj->ts ), fd( femObj->fd ), ps( femObj->ps ),
msd( 0 ) //, randOBJ( new RandNumbers( 1.0, 1000 ) )
//kount( 0 ), init( 0 ), ds( 0.0 ), duration( 0 )
{}
//------------------------------------------------------
Tracer::~Tracer() //--- the destructors
{

//	memObj->Delete2ndMat( lohi, 2 );
	
	delete mpObj[ 0 ];
	delete mpObj[ 1 ];
	delete [] veloi;
	delete [] velot;
	delete [] dispt;
	delete [] xdisp;
	delete [] tdsnf;
	delete [] tdisp;		
	if( msd )
	{
		fclose( outFile );
		fclose( ofile );
		fclose( ofile2 );
		fclose( ofile3 );
	}
}
//------------------------------------------------------
void Tracer::Init() //--- interpolate velocities onto tracers
{
	// --- some initialization
	// --- box bounds
	double dGamma = ts->dt * fd->shearRATE; //--- shear strain
	memObj->doubleMat( lohi, 2, 5 );
	lohi[ 0 ][ 0 ] = domain->xlo0;
	lohi[ 0 ][ 1 ] = domain->xhi0;
	lohi[ 0 ][ 2 ] = domain->ylo0;
	lohi[ 0 ][ 3 ] = domain->yhi0;
	lohi[ 0 ][ 4 ] = domain->xy0;
	lohi[ 1 ][ 0 ] = domain->xlo0; //--- need it to get proper disp. 
	lohi[ 1 ][ 1 ] = domain->xhi0;
	lohi[ 1 ][ 2 ] = domain->ylo0;
	lohi[ 1 ][ 3 ] = domain->yhi0;
	lohi[ 1 ][ 4 ] = domain->xy0 + dGamma * ( domain->yhi0 - domain->ylo0 );
	if( ps->shear ) //--- pure shear
	{
		//--- update xy, xlo, xhi, ylo, yhi
		lohi[ 1 ][ 4 ] = domain->xy0;
		lohi[ 1 ][ 0 ] = domain->xlo0 + 0.5 * dGamma * ( domain->xlo0 - domain->xhi0 );
		lohi[ 1 ][ 1 ] = domain->xhi0 + 0.5 * dGamma * ( domain->xhi0 - domain->xlo0 );
		lohi[ 1 ][ 2 ] = domain->ylo0 - 0.5 * dGamma * ( domain->ylo0 - domain->yhi0 );
		lohi[ 1 ][ 3 ] = domain->yhi0 - 0.5 * dGamma * ( domain->yhi0 - domain->ylo0 );
	}
	cordPtr[ 0 ] = node->cordi;
	cordPtr[ 1 ] = node->coord;
	mpObj[ 0 ] = new Wrap( lohi[ 0 ] );
	mpObj[ 1 ] = new Wrap( lohi[ 1 ] );

	elem->Init();
	double **posgp, *weigp;
	memObj->doubleMat( posgp, elem->ngaus, domain->ndime );
	memObj->doubleVec( weigp, elem->ngaus );
	elem->GAUSSQ( posgp, weigp ); // gauss points: pos and weight
	#pragma omp parallel num_threads( nThreads )		
	{
	double elcod[ elem->nevab ];
	double exisp, etasp;	
	unsigned int kgaus = 0;
	#pragma omp for
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
//		if( ielem == 0 )
//		{
//		for( unsigned int ievab = 0; ievab < elem->nevab; ievab++ )
//			printf("x=%g\n", elcod[ ievab ] );
//		}
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			exisp = posgp[ igaus ][ 0 ]; 
			etasp = posgp[ igaus ][ 1 ];
			GetNodalCords( elcod, ielem, 1, 0 ); // xyz with pbc( 1 ) and initial frame( 0 )
			//--- update coord
			kgaus = ielem * elem->ngaus + igaus;
			Nod2Gas( coord[ kgaus ], elcod, exisp, etasp ); //--- update coord		
//			if( ielem == 0 )
//				printf("Ntracer=%g\t%g\n", coord[ kgaus ][ 0 ], coord[ kgaus ][ 1 ] );
		} //--- end of gauss loop
	} //--- end of elem. loop
	} //--- end of par. sect.
	 //--- wrap coord
	Wrap mpObjLoc = lohi[ 0 ];
 	mpObjLoc.CoordWrapper( elem->nelem * elem->ngaus, coord );

	memObj->Delete2ndMat( posgp, elem->ngaus );
	delete [] weigp;
}
//------------------------------------------------------
void Tracer::TracerUpdate() //--- interpolate velocities onto tracers
{
	unsigned int kgaus, isvab, elemi;
	double exisp, etasp, dldis[ elem->nevab ], dsnaf[ elem->nevab ]; //--- dimensionless coordinates
	double vnaff[ domain->ndime ]; //--- displacement at each iteration
//---
//	ApplyHomo(); //--- update tr->cordg
/*	double width = 0.01*0.01*40.0*40.0/4.0;
	double dx, dy;
	double mnxsq = 0.0;
	for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ )
	{
		randOBJ->Normal( &dx, 0.0,  width);
		mnxsq += dx * dx;
//		printf( "dx=%e", dx );
//		tdisp[ isvab ] += dx;
		tdisp[ isvab ] = dx;
	}		
	printf( "mnxsq=%e\n", mnxsq / node->nsvab );
	mnxsq = 0.0;
*/
	#pragma omp for
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			kgaus = ielem * elem->ngaus + igaus;
//			printf( "coord[ %d ] = %f\t%f\n", kgaus, coord[ kgaus ][ 0 ], coord[ kgaus ][ 1 ] );
			FindElem( kgaus, elemi, exisp, etasp ); //--- update elemi: triangle id, exisp, etasp 
			//--- update veloi and dispt
//			GetVeloc( tdisp, dldis, elemi ); //--- update dldis( element-based ) which contains non-affine disp
			GetDisp( dldis, elemi, 0 ); // --- output node-based total( 0 ) dldis ( from tr->tdisp )
//			GetDisp( dsnaf, elemi, 1 ); // --- output node-based nonaffine( 1 ) dldis ( from tr->tdsnf )
//			matPtr.MulMatVec( femPtr->force->bmatx0, elem->nstre, elem->nevab, 
//					   dldis, elem->nevab, strag ); // strain = b * u
			Nod2Gas( vnaff, dldis, exisp, etasp ); //--- update gauss-based tracer velocity vnaff	
//			Nod2Gas( vnaff, dsnaf, exisp, etasp ); //--- update tracer velocity (non-affine)	
			for( unsigned int idime = 0; idime < domain->ndime; idime++ )
			{
				isvab = kgaus * domain->ndime + idime;
				xdisp[ isvab ] = vnaff[ idime ]; //--- instantaneous displacments
//				randOBJ->Normal( &dx, 0.0,  width);
//				dispt[ isvab ] += dx; 
//				mnxsq += vnaff[ idime ] * vnaff[ idime ];
				dispt[ isvab ] += vnaff[ idime ]; //--- update tr->dispt: non-affine only
				coord[ kgaus ][ idime ] += vnaff[ idime ]; //--- update tr->coord: aff+non-aff
			}
		}
	}
//	printf( "mnxsq=%e\n", mnxsq / elem->nelem / domain->ndime );
	//---	map tr->coord
//	#pragma omp master
/*	{
	Wrap mpObjLoc = lohi[ 0 ];
 	mpObjLoc.CoordWrapper( elem->nelem * elem->ngaus, coord );
	}
*/
//	#pragma omp barrier
}
//-----------------------------
void Tracer::FindElem( unsigned int kgaus, unsigned int &elemPtr, double &exisp, double &etasp )
// --- compute triangular element's id that contains the tracer
{
	int indx[ domain->ndime ]; //--- must be int( unsigned int? No!!! )
//	coord[ kgaus ][ 0 ]=-7.310492278376195e-01;
//	coord[ kgaus ][ 1 ]=7.642114385951331e+00;
	indx[ 0 ] = ( int )floor( ncols * ( coord[ kgaus ][ 0 ] - domain->xlo0 ) / ( domain->xhi0 - domain->xlo0 ) ) % ncols;  //--- square cell index in a periodic box
	indx[ 0 ] += ncols * ( indx[ 0 ] < 0 ); //--- modulo of a negative number
	indx[ 1 ] = ( int )floor( mrows * ( coord[ kgaus ][ 1 ] - domain->ylo0 ) / ( domain->yhi0 - domain->ylo0 ) ) % mrows; //--- in y
	indx[ 1 ] += mrows * ( indx[ 1 ] < 0 );
//	printf( "square cell index = (%d,%d) \n", indx[ 1 ], indx[ 0 ] );
	assert(  0 <= indx[ 0 ] && indx[ 0 ] < ncols ); //--- inside box?
	assert(  0 <= indx[ 1 ] && indx[ 1 ] < mrows );
	unsigned int cells[ 9 ]; //--- store cell id's that need to be scanned! index!
	cells[ 0 ] = indx[ 1 ] * ncols + indx[ 0 ]; //--- centralsquare cell index!
	int crow, ccol; //--- must be int!
	unsigned int counter = 1; //--- count
	for( int irows = -1; irows < 2; irows++ ) //--- scan neighboring cells
	{
		crow = ( indx[ 1 ] + irows ) % mrows;
//		printf( "crow0( %d ) = %d\n", irows, crow );
		crow += mrows * ( crow < 0 );
//		printf( "crow( %d ) = %d\n", irows, crow );
		assert( 0 <= crow && crow < mrows );
		for( int jcols = -1; jcols < 2; jcols++ )
		{
			if( irows == 0 && jcols == 0 ) //--- skip the central cell
				continue;
			ccol = ( indx[ 0 ] + jcols ) % ncols; //--- periodic
			ccol += ncols * ( ccol < 0 ); //--- modulo of a negative number
			assert( 0 <= ccol && ccol < ncols );
			cells[ counter ] = crow * ncols + ccol; //--- square cell index!
			counter++;
		}
	}
//	printf( "square cell index = (%d,%d) %d\n", indx[ 1 ], indx[ 0 ], cells[ 0 ] );
//	exit( EXIT_FAILURE );
	//--- loop through tiangular elements inside a given square cell and its neighbors
	bool triagFOUND = 0;
	unsigned int triagId;
//	printf( "# of triangles = %d\n", cellTriangCount[  cells[ 0 ] ] );
//        for( int jelem = 0; jelem < cellTriangCount[ cells[ 0 ] ]; jelem++ )
//		printf("%d\t", lelms[ cells[ 0 ] ][ jelem ] );
//        printf( "\n" );
//	printf( "loop in cell with triangle id ...\t");
	unsigned int celli; //--- cell id
	double elcod[ elem->nevab ]; //--- store nodal coordinates
	for( unsigned int icell = 0; icell < 9 && ( !triagFOUND ); icell++ ) //--- scan through cells
	{ 
		celli = cells[ icell ]; //--- cell id
		for( unsigned int ielem = 0; ielem < cellTriangCount[  celli ]; ielem++ ) //--- loop over triangles
		{
			triagId = lelms[ celli ][ ielem ]; //--- triangle id
//			printf( "triang id = %d\n", triagId );
			GetNodalCords( elcod, triagId, 1, 0 ); //--- xyz with pbc( 1 ) and initial frame( 0 ) update elcod
/*	if( triagId == 450 )
	{
	for( int inode = 0; inode < elem->nnode; inode++ ) //--- print nodal xy
	{
//		printf( "inode = %d\t", inode );
		for( int idime = 0; idime < domain->ndime; idime++ )
			printf( "%f\t", elcod[ inode * domain->ndime + idime ] );
		printf( "\n" );
	}		
	}
*/
			if( PointInsideTriag( coord[ kgaus ], elcod, exisp, etasp ) ) //--- update coord[ kgaus ] too!
			{
//			printf( "inside triangle %d\t", triagId );
				triagFOUND = 1;
				break;	//--- working???
			}
		}
	}
	//--- if not found->loop through every element
	if( !triagFOUND )
	{
		printf( "element containing the tracer with id = %d not found!\n", kgaus );
		for( int idime = 0; idime < domain->ndime; idime++ )
			printf( "x[%d][%d]=%16.15e\t", kgaus, idime, coord[ kgaus ][ idime ] );
		printf( "\n" );
/*		for( unsigned int ielem = 0; ielem < elem->nelem && ( !triagFOUND ); ielem++ ) //--- scan through cells
		{
			triagId = ielem;
			GetNodalCords( elcod, triagId, 1, 0 ); //--- xyz with pbc( 1 ) and initial frame( 0 ) update elcod
			if( PointInsideTriag( coord[ kgaus ], elcod, exisp, etasp ) )
			{
				triagFOUND = 1;
				break;	//--- working???
			}
		}
*/	}
	assert( triagFOUND ); //--- assert there is a tiangle
//	printf( "\ntriangle found %d !\n", triagId );
//	exit( EXIT_FAILURE );
	elemPtr = triagId; //--- assignment
}
//-----------------------------
bool Tracer::PointInsideTriag( double *coord, double *elcod, double &exisp, double &etasp )
//bool Tracer::PointInsideTriag( double &coord, double *elcod, double &exisp, double &etasp )
{
         GetExispEtasp( coord, elcod, exisp, etasp ); //--- update exisp, etasp
//	printf( "s1, s2 = %15.13e\t%15.13e\t%15.13e\n", exisp, etasp, exisp + etasp - 1.0 );
         return ( exisp >= -1.0e-12 && etasp >= -1.0e-12 && ( exisp + etasp ) - 1.0 <= 1.0e-12 ); //--- inside triangle?
 }

//-----------------------------
void Tracer::GetExispEtasp( double *cords, double *elcod, double &exisp, double &etasp )
//void Tracer::GetExispEtasp( double &cords, double *elcod, double &exisp, double &etasp )
{
	//--- calculate dimensionless coordinates of a point with respect to a tilted frame: output exisp & etasp	

//	printf( "from get exip\n" );
//	for( int inode = 0; inode < elem->nnode; inode++ ) //--- print nodal xy
//	{
//		printf( "inode = %d\t", inode );
//		for( int idime = 0; idime < domain->ndime; idime++ )
//			printf( "%f\t", elcod[ inode * domain->ndime + idime ] );
//		printf( "\n" );
//	}		
	//--- construct h matrix
	double x0 =  elcod[ 0 ];
	double y0 =  elcod[ 1 ];
	double x1 =  elcod[ 2 ];
	double y1 =  elcod[ 3 ];
	double x2 =  elcod[ 4 ];
	double y2 =  elcod[ 5 ];
	double h[ 2 ][ 2 ] = { { x1 - x0, x2 - x0  }, { y1 - y0, y2 - y0 } }; //--- h = [ a, b ]
	double det = h[ 0 ][ 0 ] * h[ 1 ][ 1 ] - h[ 0 ][ 1 ] * h[ 1 ][ 0 ]; //--- determinant
	double h_inv[ 2 ][ 2 ] = { {  h[ 1 ][ 1 ] / det, -h[ 0 ][ 1 ] / det }, 
				   { -h[ 1 ][ 0 ] / det,  h[ 0 ][ 0 ] / det } } ; //--- h inverse

//	printf("from getexip: h_inv %f\t%f\n",h_inv[ 0 ][ 0 ] , h_inv[ 0 ][ 1 ] );
//	printf("%f\t%f\n",h_inv[ 1 ][ 0 ] , h_inv[ 1 ][ 1 ] );
	//--- find minimum distance between the point and the zero-th node in the periodic box
//	double tmpCords[ domain->ndime ]; //--- tmp array to store coords
//	for( int idime = 0; idime < domain->ndime; idime++ )
//		tmpCords[ idime ] = cords[ idime ]; //--- assignment
//	double rij[ 2 ] = { tmpCords[ 0 ] - x0, tmpCords[ 1 ] - y0 }; //--- relative to the zero-th node
	double rij[ 2 ] = { cords[ 0 ] - x0, cords[ 1 ] - y0 }; //--- relative to the zero-th node
	double lx = domain->xhi0 - domain->xlo0, ly = domain->yhi0 - domain->ylo0;
//	printf("not shifted %f\t%f\n",tmpCords[ 0 ] , tmpCords[ 1 ] );
	if( rij[ 0 ] >= 0.5 * lx )
//        	tmpCords[ 0 ] -= ( lx );
        	cords[ 0 ] -= ( lx );
	else if( rij[ 0 ] < - 0.5 * lx )
//        	tmpCords[ 0 ] += ( lx );
        	cords[ 0 ] += ( lx );
	if( rij[ 1 ] >= 0.5 * ly )
//        	tmpCords[ 1 ] -= ( ly );
        	cords[ 1 ] -= ( ly );
	else if( rij[ 1 ] < - 0.5 * ly )
//        	tmpCords[ 1 ] += ( ly );
        	cords[ 1 ] += ( ly );

//	printf("from getexip:shifted %f\t%f\n",tmpCords[ 0 ] , tmpCords[ 1 ] );
	//--- dimLessCords = h_inv * ( x - x0 )                
	double dimLessCords[ 2 ];
//	rij[ 0 ] = tmpCords[ 0 ] - x0; //--- relative to the zero-th node
//	rij[ 1 ] = tmpCords[ 1 ] - y0;
	rij[ 0 ] = cords[ 0 ] - x0; //--- relative to the zero-th node
	rij[ 1 ] = cords[ 1 ] - y0;
//	printf("from getexip:relative position %f\t%f\n",rij[ 0 ] , rij[ 1 ] );
	for( int idime = 0; idime < domain->ndime; idime++ )
	{
		dimLessCords[ idime ] = 0.0;
		for( int jdime = 0; jdime < domain->ndime; jdime++ )
			dimLessCords[ idime ] += h_inv[ idime ][ jdime ] * rij[ jdime ];
	}

	exisp = dimLessCords[ 0 ]; //--- assignment
	etasp = dimLessCords[ 1 ];
//         printf( "s1, s2 = %f\t%f\n", exisp, etasp );
}
//-----------------------------
void Tracer::Nod2Gas( double *gausBased, double *nodeBased, double &S, double &T ) //--- update gausBased
{
	double shape[ elem->nnode ];
	double **deriv;
	memObj->doubleMat( deriv, domain->ndime, elem->nnode );

	elem->Shape2( shape, deriv, S, T ); //--- output shape, deriv
//	printf( "exisp, etasp = %g\t%g\n", deriv[ 0 ][ 0 ], deriv[ 1 ][ 0 ] );
	for( int idime = 0; idime < domain->ndime; idime++ )
	{
		gausBased[ idime ] = 0.0;
		for( int inode = 0; inode < elem->nnode; inode++ )
			gausBased[ idime ] += shape[ inode ] * nodeBased[ inode * domain->ndime + idime ];
	}
	memObj->Delete2ndMat( deriv, domain->ndime );
}
//-----------------------------------------------------------------------
//-----------------------------
void Tracer::Nod2Gas( double *gausBased, double **nodeBased, double &S, double &T ) //--- update gausBased
{
	double shape[ elem->nnode ];
	double **deriv;
	memObj->doubleMat( deriv, domain->ndime, elem->nnode );
	elem->Shape2( shape, deriv, S, T ); //--- output shape, deriv
//	printf( "exisp, etasp = %g\t%g\n", deriv[ 0 ][ 0 ], deriv[ 1 ][ 0 ] );
	for( int idime = 0; idime < domain->ndime; idime++ )
	{
		gausBased[ idime ] = 0.0;
		for( int inode = 0; inode < elem->nnode; inode++ )
			gausBased[ idime ] += shape[ inode ] * nodeBased[ inode ][ idime ];
	}
	memObj->Delete2ndMat( deriv, domain->ndime );
}
//-----------------------------------------------------------------------
void Tracer::GetNodalCords( double *elcod, unsigned int ielem, int pbc, int initial )
{
	unsigned int nodei;
	double Cords[ 2 * domain->ndime ];
	double r, theta;
	int cross_pbc;
	unsigned int ievab = 0;

	for( int inode = 0; inode < elem->nnode; inode++ ) // loop over nodes
	{
		nodei = elem->lnods[ ielem ][ inode ]; // node global id
		for( int idime = 0; idime < domain->ndime; idime++ ) // assign coordantes
		{
			elcod[ ievab ] = cordPtr[ initial ][ nodei ][ idime ];
			if( inode == 0 )
				Cords[ domain->ndime + idime ] = elcod[ idime ]; //--- cords[ 1 ] = elcod[ 0 ]
			Cords[ idime ] = elcod[ ievab ]; //-- cords[ 0 ] = elcod[ inode ]
//			printf("xNod=%g\n", elcod[ ievab ] );
			ievab++;
		}
//		if( pbc && inode != 0 ) // --- returns xyz relative to node0 
		if( pbc ) // --- returns xyz relative to node0 
		{
//			Xj = elcod[ 0 ]; 
//			Yj = elcod[ 1 ];
//			for( int idime = 0; idime < domain->ndime; idime++ )
//			{
//				Cords[ idime ] = elcod[ inode * domain->ndime + idime ];
//				Cords[ domain->ndime + idime ] = elcod[ idime ];
//			}
//			printf("hello!\t%d\n", inode );
			Distance( Cords, initial, r, theta, cross_pbc ); // find the true distance
//			Xj += ( r ) * cos( theta );
//			Yj += ( r ) * sin( theta );
			elcod[ ievab - domain->ndime ] = elcod[ 0 ] + ( r ) * cos( theta );
			elcod[ ievab - 1 ]     = elcod[ 1 ] + ( r ) * sin( theta );
//			printf("xNodafter=%g\n", elcod[ ievab - 1 ] );
//			printf("xNodafter=%g\n", elcod[ ievab ] );
		}
	}
//	for( unsigned int ievab = 0; ievab < elem->nevab; ievab++ )
//		printf("xAfter=%g\n", elcod[ ievab ] );
//	memObj->Delete2ndMat( Cords, 2 );
//	delete r;
//	delete theta;
//	delete cross_pbc;
}
//-----------------------------------------------------------------------
void Tracer::GetVeloc( double *output, double *elcod, unsigned int ielem )
{
	unsigned int nodei, ievab = 0, isvab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ ) // loop over nodes
	{
		nodei = elem->lnods[ ielem ][ inode ]; // node global id
		isvab = nodei * domain->ndime;
		for( int idime = 0; idime < domain->ndime; idime++ ) // assign coordantes
		{
			elcod[ ievab ] = output[ isvab ];
			ievab++;
			isvab++;
		}
	}
}
//-----------------------------------------------------------------------
void Tracer::Distance(  double *cords, int initial, double &r, double &theta, int &cross_pbc )
{
	double dimensioless_cords[ 2 * domain->ndime ];
	
	mpObj[ initial ]->GetDimensionlessCords( 2, cords, dimensioless_cords ); // output dimension...

	double s0i = dimensioless_cords[ 0 ];
	double s1i = dimensioless_cords[ 1 ];
	double s0j = dimensioless_cords[ 2 ];
	double s1j = dimensioless_cords[ 3 ];
	double s0_ij = s0j - s0i; // dimensionless distance
	double s1_ij = s1j - s1i;

	// --- translate point j to find the proper distance rij
	cross_pbc = 0;
        if( s0_ij > 0.5 )
	{ 
		s0j -= 1.0;
		cross_pbc = 1;
	}
        else if ( s0_ij < - 0.5 )
	{ 
		s0j += 1.0;
		cross_pbc = 1;
        }
	if( s1_ij > 0.5 )
	{
		s1j -= 1.0;
		cross_pbc = 1;
        }
	else if( s1_ij < - 0.5 )
	{
		s1j += 1.0;
		cross_pbc = 1;
	}
	// --- find new coordinates
	dimensioless_cords[ 0 ] = s0i;
	dimensioless_cords[ 1 ] = s1i;
	dimensioless_cords[ 2 ] = s0j;
	dimensioless_cords[ 3 ] = s1j;
	double scord[ domain->ndime * 2 ];
	mpObj[ initial ]->GetCords( 2, dimensioless_cords, scord ); // --- update coord

	double xi = scord[ 0 ];
	double yi = scord[ 1 ];
        double xj = scord[ 2 ];
        double yj = scord[ 3 ];
        double xij = xi - xj;
        double yij = yi - yj;
        double rsq = xij * xij + yij * yij;
        r = sqrt( rsq );
	theta = atan2( yij, xij );
}

//------------------------------------------------------
void Tracer::ApplyHomo()
{
//--- subtract off homogeneous deformation: update tr->tdsnf
//---
	unsigned int kgaus; //--- gauss counter
	double xydim[ domain->ndime ]; //--- dimensionless coordinates s1, s2
//---
	#pragma omp for
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			kgaus = ielem * elem->ngaus + igaus;
			mpObj[ 0 ]->GetDimensionlessCords( 1, coord[ kgaus ], xydim ); //--- update xydim lohi must be a function of t ?
			mpObj[ 1 ]->GetCords( 1, xydim, coord[ kgaus ] ); //--- update cordg
		}
	}
//---	
}
//------------------------------------------------------
void Tracer::GetDisp( double *dldis,  unsigned int ielem, bool NONAFF )
{
// --- compute proper dldis
	double elcod[ elem->nevab ], elcod0[ elem->nevab ];
//---
	GetNodalCords( elcod0, ielem, 0, 0 ); // xyz in initial frame( 0 ): no pbc( 0 )
	unsigned int ievab = 0, isvab = 0, nodei;	
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		nodei = elem->lnods[ ielem ][ inode ];
		isvab = domain->ndime * nodei;
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			if( NONAFF )
				dldis[ ievab ] = tdsnf[ isvab ]; //--- nonaffine displacements
			else 
				dldis[ ievab ] = tdisp[ isvab ]; //--- total displacements
			elcod0[ ievab ] += dldis[ ievab ];  // update original coords
			isvab++;
			ievab++;
		}
	}
//--- proper coordinates in the current frame
	ievab = 0;
	double Cords[ 2 * domain->ndime ];
	double r, theta;
	int cross_pbc;
	for( int idime = 0; idime < domain->ndime; idime++ )
		Cords[ domain->ndime + idime ] = elcod0[ idime ]; //--- cords[ 1 ] = elcod0[ 0 ]
	for( int jnode = 0; jnode <  elem->nnode; jnode++ )
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			Cords[ idime ] = elcod0[ ievab ]; //--- cords[ 0 ] = elcod[ inode ]
			ievab++;
		}
		if( NONAFF )
			Distance( Cords, 0, r, theta, cross_pbc ); //--- find coords in the initial frame
		else
			Distance( Cords, 1, r, theta, cross_pbc ); //--- find coords in the current frame
		elcod0[ ievab - domain->ndime ] = elcod0[ 0 ] + ( r ) * cos( theta );
		elcod0[ ievab - 1 ]     = elcod0[ 1 ] + ( r ) * sin( theta );
	}	
	// --- proper xyz respecting pbc in the initial frame
	GetNodalCords( elcod, ielem, 1, 0 );
	ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			dldis[ ievab ] = elcod0[ ievab ] - elcod[ ievab ];
			ievab++;
		}
	}
	//--- dldis might include +- L that needs to be taken out!
	Cords[ 0 ] = elcod0[ 0 ];
	Cords[ 1 ] = elcod0[ 1 ];
	Cords[ 2 ] = elcod[ 0 ];
	Cords[ 3 ] = elcod[ 1 ];
	double dimensioless_cords[ 2 * domain->ndime ];
	mpObj[ 1 ]->GetDimensionlessCords( 2, Cords, dimensioless_cords );
	double s0i = dimensioless_cords[ 0 ];
	double s1i = dimensioless_cords[ 1 ];
	double s0j = dimensioless_cords[ 2 ];
	double s1j = dimensioless_cords[ 3 ];
	double s0_ij = s0j - s0i; // dimensionless distance
	double s1_ij = s1j - s1i;
	assert( fabs( s0_ij ) < 0.5 );
	assert( fabs( s1_ij ) < 0.5 );
//	memObj->Delete2ndMat( elcod, elem->nnode );
//	memObj->Delete2ndMat( Cords, 2 );
//	delete r;
//	delete theta;
//	delete cross_pbc;
}

//-----------------------------------------------------------
void inline Tracer::MeanSq( double *rsqrd, double *mndis )
{
        unsigned int isvab, jsvab;
 	unsigned int ngaus = elem->nelem * elem->ngaus;
	//--- initialize ...
	for( int idofn = 0; idofn < domain->ndime; idofn++ )
	{
		mndis[ idofn ] = 0.0;
		for( int jdofn = 0; jdofn < domain->ndime; jdofn++ )
			rsqrd[ idofn * domain->ndime + jdofn ] = 0.0;
	}
//---
        for( unsigned int ipoin = 0; ipoin < ngaus; ipoin++ ) //--- msd tensor
        {
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			isvab = ipoin * domain->ndime + idofn;
			mndis[ idofn ] += dispt[ isvab ]; //--- average displacements
			for( int jdofn = 0; jdofn < domain->ndime; jdofn++ )
			{
				jsvab = ipoin * domain->ndime + jdofn;
				rsqrd[ idofn * domain->ndime + jdofn ] += dispt[ isvab ] * dispt[ jsvab ]; //--- u_\alpha u_\beta 
			}
		}
        }
	for( int idofn = 0; idofn < domain->ndime; idofn++ )
		mndis[ idofn ] /= ngaus;
	//--- average
	for( int idofn = 0; idofn < domain->ndime; idofn++ )
	{
		for( int jdofn = 0; jdofn < domain->ndime; jdofn++ )
			rsqrd[ idofn * domain->ndime + jdofn ] = rsqrd[ idofn * domain->ndime + jdofn ] / ngaus; // - mndis[ idofn ] * mndis[ jdofn ];
	}
}
//-----------------------------------------------------------
void Tracer::ComputeMSD( unsigned int itime ) // unsigned int itime )
{
	double msdis[ domain->ndime * domain->ndime ]; //--- msd
	double mndis[ domain->ndime ]; //--- mean displacements
	unsigned int isvab = 0;
 	unsigned int ngaus = elem->nelem * elem->ngaus;
//---
	MeanSq( msdis, mndis ); //--- update msd
	fprintf( outFile, "%d\t", itime );
	for( int idofn = 0; idofn < domain->ndime; idofn++ ) //--- print
	{
		for( int jdofn = 0; jdofn < domain->ndime; jdofn++ )
			fprintf( outFile, "%8.7e\t", msdis[ idofn * domain->ndime + jdofn ] );
	}
	fprintf( outFile, "\n" );
	fflush( outFile );
//---
	//--- van-hove function
	fprintf( ofile, "ITIME\n%d\n", itime );
	fprintf( ofile, "ID\tDX\tDY\n" );
        for( unsigned int ipoin = 0; ipoin < ngaus; ipoin++ ) //--- msd tensor
        {
		fprintf( ofile, "%d\t", ipoin );
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			isvab = ipoin * domain->ndime + idofn; 
			fprintf( ofile, "%4.3e\t", dispt[ isvab ] ); // - mndis[ idofn ] );
		}
		fprintf( ofile, "\n" );
        }
	fflush( ofile );
//---
	//--- positions
	fprintf( ofile3, "ITIME\n%d\n", itime );
	fprintf( ofile3, "ID\tX\tY\n" );
        for( unsigned int ipoin = 0; ipoin < ngaus; ipoin++ ) //--- msd tensor
	{
		fprintf( ofile3, "%d\t", ipoin );
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
			fprintf( ofile3, "%4.3e\t", coord[ ipoin ][ idofn ] );
		fprintf( ofile3, "\n" );
	}
	fflush( ofile3 );
	//--- instantaneous displacements
	fprintf( ofile2, "ITIME\n%d\n", itime );
	fprintf( ofile2, "ID\tDX\tDY\n" );
        for( unsigned int ipoin = 0; ipoin < ngaus; ipoin++ ) //--- msd tensor
//        for( unsigned int ipoin = 0; ipoin < node->npoin; ipoin++ ) //--- msd tensor
        {
		fprintf( ofile2, "%d\t", ipoin );
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			isvab = ipoin * domain->ndime + idofn;
//			fprintf( ofile2, "%8.7e\t", tdisp[ isvab ] );
			fprintf( ofile2, "%4.3e\t", xdisp[ isvab ] );
		}
		fprintf( ofile2, "\n" );
        }
	fflush( ofile2 );
/*
	double sdr;
	if( kount == 0 )
	{
		MeanSq( ax, ay );
		fprintf( outFile, "#DURATION\tDY^2\n" );
	}
	else
	{
		MeanSq( bx, by );
		sdr = by - ay; //--- derivative
		if( fabs( sdr ) < rsqTl ) //--- elastic loading
		{
			init = kount + 1; //--- initialize init
			telas++; 
		}
		else //--- avalanche startes!		
		{
			ds += sdr;
			duration++;
		}
		if( init - kount == 1 && duration != 0 ) //--- avalanche ends!
		{
			fprintf( outFile, "%d\t%8.7e\n", telas, ds );
			fflush( outFile );
			ds = 0.0;
			duration = 0;
			telas = 0;
		}
		ay = by;
	}
	kount++;
*/
}
