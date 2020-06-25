//#include "Python.h"
#include "dump.h"
#include "fem.h"
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <string>

#define MAXARGLEN 100

using std::cout;
using std::string;

Dump::Dump( Fem *femObj ):
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ), 
ts( femObj->ts ), fd( femObj->fd ), femPtr( femObj ), dpmd( femObj->dpmd ),
//gcObj( new GaussCoord( femObj ) ),
nDump( 0 ), memPtr( new Memory ) 
{}

Dump::~Dump()
{
	size_t found;
	for( int iDump = 0; iDump < nDump; iDump++ )
	{
		string stdStr( output[ iDump ] ); // --- standard str
		found = stdStr.find( '*' ); 
		if( found == std::string::npos ) // --- only close single file
			fclose( outptFile[ iDump ] );
	}
	if( nDump != 0 )
	{
		delete [] outptFile;
//		memPtr->Delete2ndMat( cordg, elem->nelem * elem->ngaus );
	}
	delete memPtr;
//	delete gcObj;
}

void insert( char *str1, size_t pos, char *str2 )
{
	string stdStr1( str1 );
	stdStr1.erase( pos, 1 );
	string stdStr2( str2 );
	stdStr1.insert( pos, stdStr2 );
	strcpy( str1, stdStr1.c_str() );
	
}
 void Dump::Init()
{
//	memPtr->doubleMat( cordg, elem->nelem * elem->ngaus, domain->ndime ); // allocate cordg	
//	gcObj->Init();
	outptFile = new FILE*[ nDump ];
}
void Dump::Process(  unsigned int itime,  int iDump )
{
	if( itime == 0 && iDump == 0 ) //--- define file array 
		Init();
	char tmp0[ 64 ], tmp1[ 64 ];
	sprintf( tmp0, "%d", itime ); // --- int -> string
	strcpy( tmp1, output[ iDump ] ); // --- assign tmp1
	string stdStr( tmp1 ); // --- standard str
	size_t found = stdStr.find( '*' ); // --- dump.*.xyz

	if( found != std::string::npos ) //--- separate files need to be opened at every call!
	{
		insert( tmp1, found, tmp0 ); //--- format tmp1
		outptFile[ iDump ] = fopen( tmp1, "w" );
	}
	else // --- single file
	{
		if( itime == 0 ) // open at first call
			outptFile[ iDump ] = fopen( output[ iDump ], "w" );
	}


	data = new double*[ narg[ iDump ] ];
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
	{
		if ( !strcmp( args[ iDump ][ iarg ], "act" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
//			unsigned int kounter = 0;
//			for( unsigned int kgaus = 0; kgaus < elem->nelem * elem->ngaus; kgaus++ )
//				kounter += elem->actFieldAccumulated0[ kgaus ];
//			printf( "kount=%d\n", kounter );
			if( femPtr->qs->min )
			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = elem->actFieldAccumulated[ ipoin ]; 
			}
			else
			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = elem->actField[ ipoin ]; 
			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vi_x" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->veloi[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vi_y" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->veloi[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "dx" ) )
		{ 
			N = node->npoin;
//			if( femPtr->tr->msd )
//				N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
/*			if( femPtr->tr->msd )
			{
				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = femPtr->tr->dispt[ ipoin * domain->ndime ];	
			}
			else
			{
*/				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->dispt[ ipoin * domain->ndime ];
//			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "dy0" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->dispt0[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "dx0" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->dispt0[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "dy" ) )
		{ 
			N = node->npoin;
//			if( femPtr->tr->msd )
//				N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
/*			if( femPtr->tr->msd )
			{
				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = femPtr->tr->dispt[ ipoin * domain->ndime + 1 ];	
			}
			else
			{
*/				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->dispt[ ipoin * domain->ndime + 1 ];
//			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vx" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->velot[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vy" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->velot[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "fx" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->fintl[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "fy" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->fintl[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "sxy0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strsi[ igaus ][ 2 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "sxy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strst[ igaus ][ 2 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "sxx0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strsi[ igaus ][ 0 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "sxx" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strst[ igaus ][ 0 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "syy0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strsi[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "syy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strst[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "exx0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strni[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "exx" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strnt[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eyy0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strni[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eyy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strnt[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "exy0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->strni[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "exy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->strnt[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eDoTxy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->stranRate[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eDoTxx" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->stranRate[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eDoTyy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->stranRate[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "vol" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->vol[ igaus ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "frict" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->frict[ igaus ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "frict0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->frict0[ igaus ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "hards" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->hards[ igaus ];
		}
/*		if( !strcmp( args[ iDump ][ iarg ], "ep0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->epstn0[ igaus ];
		}
*/		if( !strcmp( args[ iDump ][ iarg ], "ep" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->epstn[ igaus ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "epX-Y" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->plasticStrain[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "epY-X" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->plasticStrain[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "epXY" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->plasticStrain[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "x" ) )
		{ 
			N = node->npoin;
//			if( femPtr->tr->msd )
//				N = elem->nelem * elem->ngaus ;
			data[ iarg ] = new double[ N ];
//			if( femPtr->tr->msd )
//			{
//				for( int ipoin = 0; ipoin < N; ipoin++ )
//					data[ iarg ][ ipoin ] = femPtr->tr->coord[ ipoin ][ 0 ];
//			}
//			else
//			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->coord[ ipoin ][ 0 ];
//			}
		}
		if( !strcmp( args[ iDump ][ iarg ], "xu" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->coordUnWrapped[ ipoin ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "yu" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->coordUnWrapped[ ipoin ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "y" ) )
		{ 
			N = node->npoin;
//			if( femPtr->tr->msd )
//				N = elem->nelem * elem->ngaus ;
			data[ iarg ] = new double[ N ];
//			if( femPtr->tr->msd )
//			{
//				for( int ipoin = 0; ipoin < N; ipoin++ )
//					data[ iarg ][ ipoin ] = femPtr->tr->coord[ ipoin ][ 1 ];
//			}
//			else
//			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->coord[ ipoin ][ 1 ];
//			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "n0" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( unsigned int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->lnods[ ielem ][ 0 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "n1" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( unsigned int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->lnods[ ielem ][ 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "n2" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( unsigned int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->lnods[ ielem ][ 2 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "yf0" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->fyild0[ igaus ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "yf." ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = ( elem->fyild[ igaus ] - elem->fyild0[ igaus ] ) / ts->dt;
		}
		if ( !strcmp( args[ iDump ][ iarg ], "yf" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->fyild[ igaus ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "xg" ) )
		{ 
			// computes cordg
//			gcObj->GetCoord( cordg );
			// assign data
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->cordg[ igaus ][ 0 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "yg" ) )
		{ 
//			gcObj->GetCoord( cordg );
			// assign data
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->cordg[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "K" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->kmodu[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "G" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->gmodu[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "sy0" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->uniax0[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "sy" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->uniax[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "Ey" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->strny[ ielem ];
		}
	}
	//--- write output
	fprintf( outptFile[ iDump ], "Timestep\n%d\n", itime );
	fprintf( outptFile[ iDump ], "Number of Nodes\n%d\n", N );
//	fprintf( outptFile[ iDump ], "%d\n", N );
	fprintf( outptFile[ iDump ], "Box Bounds\n%g\t%g\n%g\t%g\n%g\n",
		domain->xlo, domain->xhi, domain->ylo, 
		domain->yhi, domain->xy );
//---   print items specified in args dictionary
	fprintf( outptFile[ iDump ], "Atoms\t" );
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
		fprintf( outptFile[ iDump ], "%s\t", args[ iDump ][ iarg ] );
	fprintf( outptFile[ iDump ], "\n" );
/*
//---   set format
	char *formt[ narg[ iDump ] ];
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
	{
		formt[ iarg ] = new char[ MAXARGLEN ];
		if( dpmd->exist[ iDump ] )
			strcpy( formt[ iarg ], dpmd->args[ iDump ][ iarg ] );
		else
		{	
			printf("idump=%d\n",iDump);
			formt[ iarg ] = "%e";
		}
	}
*/
//---
	
//	printf("%d\n", iDump );
//	printf("%d\n",dpmd->args[ iDump ][ 0 ] );//, dpmd->args[ iDump ][ 1 ]);
//	printf("\n" );
	for( int ipoin = 0; ipoin < N; ipoin++ )
	{
		fprintf( outptFile[ iDump ], "%d\t", ipoin );
		for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
		{
			if( dpmd->exist[ iDump ] )
				fprintf( outptFile[ iDump ], "%*.*e\t", dpmd->args[ iDump ][ 0 ], dpmd->args[ iDump ][ 1 ], data[ iarg ][ ipoin ] );
			else
				fprintf( outptFile[ iDump ], "%g\t", data[ iarg ][ ipoin ] );
		}
		fprintf( outptFile[ iDump ], "\n" );
	}
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
	{
		delete [] data[ iarg ];
//		delete [] formt[ iarg ];
	}
	delete [] data;
//	delete [] formt;
	fflush( outptFile[ iDump ] );
	if( found != string::npos ) //--- separate files need to be closed at every call!
	{
//		printf( "bug!\n" );
		fclose( outptFile[ iDump ] );
	}
}

