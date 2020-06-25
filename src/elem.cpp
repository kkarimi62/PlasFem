#include "elem.h"
#include "fem.h"
#include "gauss_coord.h"
#include "stdlib.h"
#include "matrix.h"
#include <cassert>
#include <iostream>

using std::cout;

Elem::Elem( Fem *femObj ):
femPtr( femObj ), 
memPtr( new Memory ),
domain( femObj->domain )
{}

Elem::~Elem()
{
	delete [] ltype;
	memPtr->Delete2ndMat( lnods, nelem );
	delete [] kmodu;
	delete [] gmodu;
	delete [] rsidu;
	delete [] rsidf;
	delete [] rsidp;
	delete [] reduc;
	delete [] redup;
	delete [] reduf;
	delete [] pcrit;
	delete [] pcrit0;
	delete [] mnunx;
	delete [] hards;
	delete [] uniax;
	delete [] uniax0;
	delete [] frict0;
	delete [] infrc;
	delete [] inunx;
	delete [] frict;
	delete [] strny;
	delete [] mufst;
	delete [] musnd;
	delete [] tprin;
	memPtr->Delete2ndMat( stran, nelem * ngaus );
	memPtr->Delete2ndMat( strni, nelem * ngaus );
	memPtr->Delete2ndMat( strst, nelem * ngaus );
	memPtr->Delete2ndMat( strsi, nelem * ngaus );
	memPtr->Delete2ndMat( strnt, nelem * ngaus );
	memPtr->Delete2ndMat( strrs, nelem * ngaus );
	delete [] epstn;
	memPtr->Delete2ndMat( stranRate, nelem * ngaus );
	memPtr->Delete2ndMat( plasticStrain, nelem * ngaus );
	memPtr->Delete2ndMat( plasticStrain0, nelem * ngaus );
	delete [] vol;
	delete [] actField;
	delete [] actField0;
	delete [] actFieldAccumulated;
	delete [] actFieldAccumulated0;
	delete [] resid;
	delete [] presd;
	delete [] sresd;
	delete [] sxyrd;
	delete [] fyild;	
	delete [] fyild0;	
//	memPtr->Delete2ndMat( deriv, domain->ndime );
//	memPtr->Delete2ndMat( XJACM, domain->ndime );
//	memPtr->Delete2ndMat( XJACI, domain->ndime );
//	delete [] shape;
	delete memPtr;
	delete randOBJ;
	delete randOBJ2nd;

	memPtr->Delete2ndMat( cordg, nelem * ngaus );
}

void Elem::Init()
{
	//--- instantiate random number generator object
	unsigned int seed;
	FILE *urandom = fopen( "/dev/urandom", "r" );
	fread( &seed, sizeof ( seed ), 1, urandom );
	randOBJ = new RandNumbers( seed );
//	randOBJ2nd = new RandNumbers( lamda2nd, seed );
	fclose( urandom );
//	printf( "ndime=%d\n", femPtr->domain->ndime );
//	printf( "nnode=%d\n", nnode );
//	memPtr->doubleVec( shape, nnode );
//	memPtr->doubleMat( deriv, domain->ndime, nnode );
//	memPtr->doubleMat( XJACM, domain->ndime, domain->ndime );
//	memPtr->doubleMat( XJACI, domain->ndime, domain->ndime );
	memPtr->doubleMat( cordg, nelem * ngaus, domain->ndime ); // allocate cordg	
	gcObj = new GaussCoord( femPtr );
	gcObj->Init();
	gcObj->GetCoord( cordg );
}
//------------------------------------------------------------------
void Elem::ksiTOxy( double *output,  double *ksi,  double **nodalXY )
{
	double *shape, **deriv;
	memPtr->doubleVec( shape, nnode );
	memPtr->doubleMat( deriv, domain->ndime, nnode );
	Shape2( shape, deriv, ksi[ 0 ], ksi[ 1 ] ); // output shape & deriv functions 
	for( int idime = 0; idime < domain->ndime; idime++ )
		for( int inode = 0; inode < nnode; inode++ )
			output[ idime ] += shape[ inode ] * nodalXY[ inode ][ idime ];
	delete [] shape;
	for( int idime = 0; idime < domain->ndime; idime++ ){
		delete [] deriv[ idime ];
	}
	delete [] deriv;

}
//------------------------------------------------------------------
void Elem::set(  unsigned int ielem, const double mnunx, const double minsy )
{
	// set gmodu and kmodu
	//--- set yield strain	
	randOBJ->Expon( &( uniax[ ielem ] ), 1.0 / ( mnunx - minsy ) ); //--- lambda = (mean - min )^-1
	uniax[ ielem ] += minsy;
	assert( uniax[ ielem ] - minsy >= 0.0 );
//	randOBJ2nd->Expon( &( frict[ ielem ] ) ); //--- weakening in friction
//	uniax[ ielem ] += minE;
//	hards[ ielem ] = HARDS;
//	uniax[ ielem ] = sqrt( 2.0 * gmodu[ ielem ] * uniax[ ielem ] ); // energy to stress
//	uniax[ ielem ] = 2.0 * sqrt( 1.0 * uniax[ ielem ] ); // energy to stress
//	strny[ ielem ] = meanRecoveryStrain; //femPtr->elem->meanRecoveryStrain;
}
//------------------------------------------------------------------
void Elem::Shape2( double *shape,  /*output*/
		   double **deriv, /*output*/
		    double S, 
		    double T )
{
// --- ( S, T ) material coordinates ranging between -1 and 1
// --- SHAPE FUNCTIONS AND DERIVATIVES FOR LINEAR TRIANGLE
	double p;
	double S2, T2, SS, TT, ST, SST, STT, ST2;

	if( nnode == 3 )
	{
		p = 1.0 - S - T;
		shape[ 0 ] = p; 
		shape[ 1 ] = S;
		shape[ 2 ] = T;

		deriv[ 0 ][ 0 ] = -1.0;
		deriv[ 0 ][ 1 ] = 1.0;
		deriv[ 0 ][ 2 ] = 0.0;
		deriv[ 1 ][ 0 ] = -1.0;
		deriv[ 1 ][ 1 ] = 0.0;
		deriv[ 1 ][ 2 ] = 1.0;
	}
	else if( nnode == 4 ){
		S2 = S * 2.0;
		T2 = T * 2.0;
		SS = S * S;
		TT = T * T;
		ST = S * T;
		SST = S * S * T;
		STT = S * T * T;
		ST2 = S * T * 2.0;
//--- SHAPE FUNCTIONS FOR 4 NODED ELEMENT
		shape[ 0 ] = ( 1.0 - S - T + ST ) * 0.25;
		shape[ 1 ] = ( 1.0 + S - T - ST ) * 0.25;
		shape[ 2 ] = ( 1.0 + S + T + ST ) * 0.25;
		shape[ 3 ] = ( 1.0 - S + T -ST ) * 0.25;
//--- AND DERIVATIVES
		deriv[ 0 ][ 0 ] = ( -1.0 + T ) * 0.25;
		deriv[ 0 ][ 1 ] = -deriv[ 0 ][ 0 ];
		deriv[ 0 ][ 2 ] = ( 1.0 + T ) * 0.25;
		deriv[ 0 ][ 3 ] = -deriv[ 0 ][ 2 ];
		deriv[ 1 ][ 0 ] = ( -1.0 + S ) * 0.25;
		deriv[ 1 ][ 1 ] = ( -1.0 - S ) * 0.25;
		deriv[ 1 ][ 2 ] = -deriv[ 1 ][ 1 ];
		deriv[ 1 ][ 3 ] = -deriv[ 1 ][ 0 ];
	}
}
//------------------------------------------------------------------------
void Elem::Blarge( double **bmatx, /*output*/
		    double **cartd )
{
	int ngash, mgash;
	ngash = -1;
	for( int inode = 0; inode < nnode; inode++ )
	{
		mgash = ngash + 1;
		ngash = mgash + 1;
		bmatx[ 0 ][ mgash ] = cartd[ 0 ][ inode ];
		bmatx[ 0 ][ ngash ] = 0.0;
		bmatx[ 1 ][ mgash ] = 0.0;
		bmatx[ 1 ][ ngash ] = cartd[ 1 ][ inode ];
		bmatx[ 2 ][ mgash ] = cartd[ 1 ][ inode ]; //* 0.5
		bmatx[ 2 ][ ngash ] = cartd[ 0 ][ inode ]; //* 0.5
	}	
}
//------------------------------------------------------------------------
void Elem::Jacob2( double **cartd, /*output*/
		   double &djacb,  /*output*/
		    double S, 
		    double T, 
		    double **elcod, 
		    unsigned int ielem )
{
	double *shape, **deriv;
	double **XJACM, **XJACI;

	memPtr->doubleVec( shape, nnode );
	memPtr->doubleMat( deriv, domain->ndime, nnode );
	memPtr->doubleMat( XJACM, domain->ndime, domain->ndime );
	memPtr->doubleMat( XJACI, domain->ndime, domain->ndime );
		
	Shape2( shape, deriv, S, T ); // output shape, deriv
	matPtr->MulMatMat( deriv, domain->ndime, nnode, elcod, nnode, domain->ndime, XJACM ); //XJACM = deriv * elcod
//	}
	djacb = det( XJACM ); 
//	cout << *djacb << '\n';
//	exit( EXIT_FAILURE );
	if( djacb <= 0.0 )
	{
//		for( int i = 0; i < 3; i++ ){
//			for( int j = 0; j < 2; j++ )
//				cout << elcod[ i ][ j ] << '\t';
//			cout << '\n';
//		}
		printf( "STOP IN JACOB2: element %d\t%f\n", ielem, djacb );
		exit( EXIT_FAILURE );
	}
	inv( XJACI, XJACM ); //inverse
	matPtr->MulMatMat( XJACI, domain->ndime, domain->ndime, deriv, domain->ndime, nnode, cartd ); //output cartd = XJACI * deriv

	delete [] shape;
	for( int idime = 0; idime < domain->ndime; idime++ )
	{
		delete [] deriv[ idime ];
		delete [] XJACM[ idime ];
		delete [] XJACI[ idime ];
	}
	delete [] deriv;
	delete [] XJACM;
	delete [] XJACI;
}

//-----------------------------------------------------------------------------
inline void Elem::inv( double **XJACI,  double **XJACM )
{
	double determinant;
	determinant = XJACM[ 0 ][ 0 ] * XJACM[ 1 ][ 1 ] - XJACM[ 0 ][ 1 ] * XJACM[ 1 ][ 0 ];
	assert( determinant != 0.0 );
	XJACI[ 0 ][ 0 ] = XJACM[ 1 ][ 1 ] / determinant;
	XJACI[ 1 ][ 1 ] = XJACM[ 0 ][ 0 ] / determinant;
	XJACI[ 0 ][ 1 ] = -XJACM[ 0 ][ 1 ] / determinant;
	XJACI[ 1 ][ 0 ] = -XJACM[ 1 ][ 0 ] / determinant;
}
//-----------------------------------------------------------------------------
inline double Elem::det(  double ** XJACM )
{
	double determinant;
	determinant = XJACM[ 0 ][ 0 ] * XJACM[ 1 ][ 1 ] - XJACM[ 0 ][ 1 ] * XJACM[ 1 ][ 0 ];
	return determinant;
}
//-----------------------------------------------------------------------------
void Elem::GAUSSQ( double **posgp, double *weigp )
{
//	nnode = femPtr->elem->nnode;
//	ngaus = femPtr->elem->ngaus;
//	int ndime = femPtr->domain->ndime;
	if( nnode == 3 )
	{
		posgp[ 0 ][ 0 ] = 1.0 / 3.0;
		posgp[ 0 ][ 1 ] = 1.0 / 3.0;
		weigp[ 0 ] = 0.5;
	}
	else if( nnode == 4 && ngaus == 4 )
	{
		double G = 0.577350269189626;
		posgp[ 0 ][ 0 ] = -1.0;
		posgp[ 0 ][ 1 ] = -1.0;
		posgp[ 1 ][ 0 ] = -1.0;
		posgp[ 1 ][ 1 ] = 1.0;
		posgp[ 2 ][ 0 ] = 1.0;
		posgp[ 2 ][ 1 ] = -1.0;
		posgp[ 3 ][ 0 ] = 1.0;
		posgp[ 3 ][ 1 ] = 1.0;
		for( int igaus = 0; igaus < ngaus; igaus++ )
		{
			weigp[ igaus ] = 1.0;
			for( int idime = 0; idime < domain->ndime; idime++ )
				posgp[ igaus ][ idime ] *= G;
		}
	}
	else if( nnode == 4 && ngaus == 1 )
	{
	        posgp[ 0 ][ 0 ] = 0.0;
	        posgp[ 0 ][ 1 ] = 0.0;
	        weigp[ 0 ] = 4.0;
	}
	
}
