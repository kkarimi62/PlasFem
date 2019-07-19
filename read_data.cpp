//#include "Python.h"
#include <iostream>
#include <stdio.h>
#include <cstring>
#include "stdlib.h"
#include "read_data.h"

//---------------------------------------------------------
ReadData::ReadData( char* dataFile, Fem *femObj ):
MAXSTRLEN( 1024 ), MAXARG( 1024 )
{
	infile = fopen( dataFile, "r" );
	line = new char[ MAXSTRLEN ]; // command string
	args = new char*[ MAXARG ]; // string for a parsed command
	for( int i = 0; i < MAXARG; i++ ) args[ i ] = new char[ MAXSTRLEN ];
	femPtr = femObj;	
	memory = new Memory();
}
//---------------------------------------------------------
ReadData::~ReadData()
{
	delete [] line;
	delete [] args;
	delete memory;
//	delete [] femPtr->node->coord;
//	delete [] femPtr->node->veloc;
//	delete [] femPtr->elem->ltype;
//	delete [] femPtr->elem->mass;
//	delete [] femPtr->elem->lnods;
}
//---------------------------------------------------------
void ReadData::process()
{
// junk
	for( int i = 0; i < 2; i++ ){
		fgets( line, MAXSTRLEN, infile );
	}
// initial time
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->ts->itime0 = atoi( args[ 0 ] );
// dimensions
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->domain->ndime = atoi( args[ 0 ] );
// npoin
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->node->npoin = atoi( args[ 0 ] );	
// nelem
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->elem->nelem = atoi( args[ 0 ] );	
// nnode
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->elem->nnode = atoi( args[ 0 ] );	
	femPtr->elem->nevab = femPtr->elem->nnode * femPtr->domain->ndime;
// ngaus
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->elem->ngaus = atoi( args[ 0 ] );	
// nstre
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->elem->nstre = atoi( args[ 0 ] );	
// ntype
	for( int i = 0; i < 2; i++ ) 
		fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->elem->ntype = atoi( args[ 0 ] );	
// lohi
	for( int i = 0; i < 2; i++ ) 
		fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->domain->xlo0 = atof( args[ 0 ] );	
	femPtr->domain->xhi0 = atof( args[ 1 ] );	
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->domain->ylo0 = atof( args[ 0 ] );	
	femPtr->domain->yhi0 = atof( args[ 1 ] );	
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->domain->xy0 = atof( args[ 0 ] );
// coord
	for( int i = 0; i < 3; i++ )
		fgets( line, MAXSTRLEN, infile );
	unsigned int npoin = femPtr->node->npoin;
	unsigned int ndime = femPtr->domain->ndime;
	memory->doubleMat( femPtr->node->cordi, npoin, ndime );	
	memory->doubleMat( femPtr->node->cordiUnWrapped, npoin, ndime );	
	femPtr->node->mass = new double[ femPtr->node->npoin ];
	for( unsigned int ipoin = 0; ipoin < npoin; ipoin++ ){ // parsing
		fgets( line, MAXSTRLEN, infile );
		parse();
		femPtr->node->mass[ ipoin ] = atof( args[ 1 ] );
		for( int idime = 0; idime < ndime; idime++ ){
			femPtr->node->cordi[ ipoin ][ idime ] = atof( args[ idime + 2 ] );
			femPtr->node->cordiUnWrapped[ ipoin ][ idime ] = atof( args[ idime + 2 ] );
		}
	}
// cordi = coord	
	memory->doubleMat( femPtr->node->coord, npoin, ndime );	
	memory->doubleMat( femPtr->tr->coord, femPtr->elem->nelem * femPtr->elem->ngaus, ndime );
	memory->doubleMat( femPtr->node->coordUnWrapped, npoin, ndime );	
//	for( unsigned int ipoin = 0; ipoin < npoin; ipoin++ ){
//		for( int idime = 0; idime < ndime; idime++ ){
//			femPtr->node->cordi[ ipoin ][ idime ] = femPtr->node->coord[ ipoin ][ idime ];
//			femPtr->node->cordiUnWrapped[ ipoin ][ idime ] = femPtr->node->coord[ ipoin ][ idime ];
//			femPtr->node->coordUnWrapped[ ipoin ][ idime ] = femPtr->node->coord[ ipoin ][ idime ];
//		}
//	}
	femPtr->node->nsvab = npoin * ndime;
// velocity
	for( int i = 0; i < 3; i++ ) 
		fgets( line, MAXSTRLEN, infile );
	memory->doubleVec( femPtr->node->velot, femPtr->node->nsvab );
	memory->doubleVec( femPtr->tr->velot, femPtr->elem->nelem * femPtr->elem->ngaus * femPtr->domain->ndime ); //--- msd

	femPtr->node->fixID = new int[ femPtr->node->nsvab ];
	femPtr->node->rigidID = new int[ femPtr->node->nsvab ];
	memory->doubleVec( femPtr->node->presVeloc, femPtr->node->nsvab );

	//parsing data
	int isvab = 0;
	for( unsigned int ipoin = 0; ipoin < npoin; ipoin++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		for( int idime = 0; idime < ndime; idime++ ){
			femPtr->node->velot[ isvab ] = atof( args[ idime + 1 ] );
			isvab++; 
		}
	}
	//midpoint velocity
	memory->doubleVec( femPtr->node->veloi, femPtr->node->nsvab );
	memory->doubleVec( femPtr->tr->veloi, femPtr->elem->nelem * femPtr->elem->ngaus * femPtr->domain->ndime ); //--- msd

//	for( int isvab = 0; isvab < femPtr->node->nsvab; isvab++ )
//		femPtr->node->veloi[ isvab ] = femPtr->node->velot[ isvab ];
// connectivity
	for( int i = 0; i < 3; i++ )
		fgets( line, MAXSTRLEN, infile );
	unsigned int nelem = femPtr->elem->nelem; 
	femPtr->elem->ltype = new unsigned int[ nelem ];
	femPtr->elem->lnods = new unsigned int*[ nelem ];
	int nnode = femPtr->elem->nnode;
	for( unsigned int ielem = 0; ielem < nelem; ielem++ )
		femPtr->elem->lnods[ ielem ] = new unsigned int[ nnode ];
	for( unsigned int ielem = 0; ielem < nelem; ielem++ ){ 		
		fgets( line, MAXSTRLEN, infile );
		parse();
		femPtr->elem->ltype[ ielem ] = atoi( args[ 1 ] );
		for( int inode = 0; inode < nnode; inode++ ) 
			femPtr->elem->lnods[ ielem ][ inode ] = atoi( args[ 2 + inode ] );
	}

//--- cells and thier corresponding triangular elements (msd)
	for( int i = 0; i < 2; i++ )
		fgets( line, MAXSTRLEN, infile );
	//--- size
	fgets( line, MAXSTRLEN, infile );
	parse();
	femPtr->tr->mrows = atoi( args[ 0 ] );
	femPtr->tr->ncols = atoi( args[ 1 ] );
//        printf( "%d\t%d\n", femPtr->tr->mrows, femPtr->tr->ncols );
	unsigned int ncels = femPtr->tr->mrows * femPtr->tr->ncols; 
	//--- header	
	fgets( line, MAXSTRLEN, infile );
	femPtr->tr->cellTriangCount = new unsigned int[ ncels ]; //--- memory allocation	
	femPtr->tr->lelms = new unsigned int*[ ncels ]; //--- store triangles' id
	unsigned int index, *tmp;
	for( unsigned int icels = 0; icels < ncels; icels++ )
	{
//        	printf( "get...\n");
		fgets( line, MAXSTRLEN, infile ); //--- read
 //       	printf( "parse...\n");
		parse();
//        	printf( "index...\n");
		index = atoi( args[ 0 ] ); //--- cell index
//        	printf( "index, count %d\t%d\n", index, atoi( args[ 1 ] ) );
		femPtr->tr->cellTriangCount[ index ] = atoi( args[ 1 ] ); //--- # of triangles
		tmp = new unsigned int[ atoi( args[ 1 ] ) ]; //--- assign memory: tmp
		for( int itrig = 0; itrig < atoi( args[ 1 ] ); itrig++ ) 
			tmp[ itrig ] = atoi( args[ 2 + itrig ] ); //--- save list of triangles' id
//        	printf( "pointer assignment ... \n");
		femPtr->tr->lelms[ index ] = tmp; //--- assignment: tmp pointer to lelms
	}
//        printf( "done\n");

//	for( int ielem =0; ielem < ncels; ielem++ )
//	{
//        printf( "%d\t", ielem );
//        for( int jelem = 0; jelem < femPtr->tr->cellTriangCount[ ielem ]; jelem++ )
//                printf("%d\t", femPtr->tr->lelms[ ielem ][ jelem ] );
//        printf( "\n" );
//	}


// initial stress
	memory->doubleVec( femPtr->elem->fyild, femPtr->elem->nelem * femPtr->elem->ngaus );
	memory->doubleVec( femPtr->elem->fyild0, femPtr->elem->nelem * femPtr->elem->ngaus );
	for( int i = 0; i < 3; i++ )
		fgets( line, MAXSTRLEN, infile );
	memory->doubleMat( femPtr->elem->strsi, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	for( unsigned int igaus = 0; igaus < femPtr->elem->nelem * femPtr->elem->ngaus; igaus++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		for( int istre = 0; istre < femPtr->elem->nstre; istre++ ) 
			femPtr->elem->strsi[ igaus ][ istre ] = atof( args[ istre + 1 ] );
	} 		
// plastic strain
	for( int i = 0; i < 3; i++ )
		fgets( line, MAXSTRLEN, infile );
	memory->doubleMat( femPtr->elem->plasticStrain0, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	for( unsigned int igaus = 0; igaus < femPtr->elem->nelem * femPtr->elem->ngaus; igaus++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		for( int istre = 0; istre < femPtr->elem->nstre; istre++ ) 
			femPtr->elem->plasticStrain0[ igaus ][ istre ] = atof( args[ istre + 1 ] );
	} 		
// displacements
//	for( int i = 0; i < 3; i++ ) 
//		fgets( line, MAXSTRLEN, infile );
	memory->doubleVec( femPtr->node->dispt, femPtr->node->nsvab );
	memory->doubleVec( femPtr->node->dispt0, femPtr->node->nsvab );
	memory->doubleVec( femPtr->tr->dispt, femPtr->elem->nelem * femPtr->elem->ngaus * femPtr->domain->ndime ); //--msd
	memory->doubleVec( femPtr->tr->xdisp, femPtr->elem->nelem * femPtr->elem->ngaus * femPtr->domain->ndime ); //--msd
	memory->doubleVec( femPtr->tr->tdsnf, femPtr->node->nsvab ); //--msd
	memory->doubleVec( femPtr->tr->tdisp, femPtr->node->nsvab ); //--msd
	memory->doubleVec( femPtr->node->dldis, femPtr->node->nsvab );
	//parsing data
//	isvab = 0;
//	for( unsigned int ipoin = 0; ipoin < npoin; ipoin++ ){
//		fgets( line, MAXSTRLEN, infile );
//		parse();
//		for( int idime = 0; idime < ndime; idime++ ){
//			femPtr->node->dispt[ isvab ] = atof( args[ idime + 1 ] );
//			isvab++; 
//		}
//	}
// forces 
	for( int i = 0; i < 3; i++ ) 
		fgets( line, MAXSTRLEN, infile );
	memory->doubleVec( femPtr->node->fintl, femPtr->node->nsvab );
	memory->doubleVec( femPtr->node->absfl, femPtr->node->nsvab );
	memory->doubleVec( femPtr->node->extForce, femPtr->node->nsvab );
	//parsing data
	isvab = 0;
	for( unsigned int ipoin = 0; ipoin < npoin; ipoin++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		for( int idime = 0; idime < ndime; idime++ ){
			femPtr->node->fintl[ isvab ] = atof( args[ idime + 1 ] );
			isvab++; 
		}
	}
// initial strain 
	for( int i = 0; i < 3; i++ )
		fgets( line, MAXSTRLEN, infile );
	memory->doubleMat( femPtr->elem->strni, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	for( unsigned int igaus = 0; igaus < femPtr->elem->nelem * femPtr->elem->ngaus; igaus++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		for( int istre = 0; istre < femPtr->elem->nstre; istre++ ) 
			femPtr->elem->strni[ igaus ][ istre ] = atof( args[ istre + 1 ] );
	} 		
// initial activity 
	for( int i = 0; i < 3; i++ )
		fgets( line, MAXSTRLEN, infile );
	femPtr->elem->actField0 = new int[ femPtr->elem->nelem * femPtr->elem->ngaus ];
	for( unsigned int igaus = 0; igaus < femPtr->elem->nelem * femPtr->elem->ngaus; igaus++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		femPtr->elem->actField0[ igaus ] = atoi( args[ 1 ] );
	} 		
// allocate memory for elem->stran
	memory->doubleMat( femPtr->elem->stran, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	memory->doubleMat( femPtr->elem->strnt, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	memory->doubleMat( femPtr->elem->strst, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	memory->doubleMat( femPtr->elem->strrs, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	memory->doubleMat( femPtr->elem->plasticStrain, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	memory->doubleMat( femPtr->elem->stranRate, femPtr->elem->nelem * femPtr->elem->ngaus, femPtr->elem->nstre );	
	memory->doubleVec( femPtr->elem->epstn, femPtr->elem->nelem * femPtr->elem->ngaus );
	memory->doubleVec( femPtr->elem->presd, femPtr->elem->nelem * femPtr->elem->ngaus );
	memory->doubleVec( femPtr->elem->sresd, femPtr->elem->nelem * femPtr->elem->ngaus );
	memory->doubleVec( femPtr->elem->sxyrd, femPtr->elem->nelem * femPtr->elem->ngaus );
	memory->doubleVec( femPtr->elem->vol, femPtr->elem->nelem * femPtr->elem->ngaus );
	femPtr->elem->actField = new int[ femPtr->elem->nelem * femPtr->elem->ngaus ];
	femPtr->elem->actFieldAccumulated= new int[ femPtr->elem->nelem * femPtr->elem->ngaus ];
	femPtr->elem->actFieldAccumulated0= new int[ femPtr->elem->nelem * femPtr->elem->ngaus ];
	femPtr->elem->resid= new bool[ femPtr->elem->nelem * femPtr->elem->ngaus ];

	fclose( infile );

}

//---------------------------------------------------------
//--- parse the string
void ReadData::parse()
{

	char *ptr;
	ptr = strtok( line, " \t" );
	int narg = 0;
	while( ptr != NULL ){
		strcpy( args[ narg ], ptr );
		ptr = strtok( NULL, " \t" );
		narg++;
	}
}
