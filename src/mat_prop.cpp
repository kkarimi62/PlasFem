//#include "Python.h"
#include <iostream>
#include <stdio.h>
#include <cstring>
#include "stdlib.h"
#include "mat_prop.h"

//---------------------------------------------------------
MatProp::MatProp( char* dataFile, Fem *femObj ):
MAXSTRLEN( 1000 ), MAXARG( 100 )
{
	infile = fopen( dataFile, "r" );
    if (infile == nullptr) {
        // Throw an error if file does not exist
        throw std::runtime_error("Error: File '" + std::string(dataFile) + "' does not exist or cannot be opened.");
    }
	line = new char[ MAXSTRLEN ]; // command string
	args = new char*[ MAXARG ]; // string for a parsed command
	for( int i = 0; i < MAXARG; i++ ) 
		args[ i ] = new char[ MAXSTRLEN ];
	femPtr = femObj;	
}
//---------------------------------------------------------
MatProp::~MatProp()
{
	delete [] line;
	delete [] args;
//	delete [] femPtr->elem->kmodu;
//	delete [] femPtr->elem->gmodu;
//	delete [] femPtr->elem->uniax;
//	delete [] femPtr->elem->strny;
}
//---------------------------------------------------------
void MatProp::process()
{
	//junk
	fgets( line, MAXSTRLEN, infile );
	// loop starts: read material properties
	unsigned int nelem = femPtr->elem->nelem;
	femPtr->elem->kmodu = new double[ nelem ];
	femPtr->elem->gmodu = new double[ nelem ];
	femPtr->elem->mufst = new double[ nelem ];
	femPtr->elem->musnd = new double[ nelem ];
	femPtr->elem->tprin = new double[ nelem ];
	femPtr->elem->uniax = new double[ nelem ];
	femPtr->elem->hards = new double[ nelem ];
	femPtr->elem->uniax0 = new double[ nelem ];
	femPtr->elem->strny = new double[ nelem ];
	femPtr->elem->frict = new double[ nelem ];
	femPtr->elem->infrc = new double[ nelem ];
	femPtr->elem->inunx = new double[ nelem ];
	femPtr->elem->frict0 = new double[ nelem ];
	femPtr->elem->rsidu = new double[ nelem ];
	femPtr->elem->rsidf = new double[ nelem ];
	femPtr->elem->rsidp = new double[ nelem ];
	femPtr->elem->reduc = new double[ nelem ];
	femPtr->elem->redup = new double[ nelem ];
	femPtr->elem->reduf = new double[ nelem ];
	femPtr->elem->pcrit = new double[ nelem ];
	femPtr->elem->pcrit0 = new double[ nelem ];
	femPtr->elem->mnunx = new double[ nelem ];
	for( int ielem = 0; ielem < nelem; ielem++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
//		if( ielem == 0 )
//			printf( "%e\n", atof( args[ 3 ] ) ); 
		femPtr->elem->kmodu[ ielem ] = atof( args[ 1 ] );
		femPtr->elem->gmodu[ ielem ] = atof( args[ 2 ] );
		femPtr->elem->mufst[ ielem ] = atof( args[ 2 ] ); //--- assuming this is isotropic
		femPtr->elem->musnd[ ielem ] = atof( args[ 2 ] );
		femPtr->elem->tprin[ ielem ] = 0.0;
		femPtr->elem->uniax0[ ielem ] = atof( args[ 3 ] );
		femPtr->elem->strny[ ielem ] = atof( args[ 4 ] );
		femPtr->elem->frict0[ ielem ] = atof( args[ 5 ] );
		femPtr->elem->hards[ ielem ] = atof( args[ 6 ] );
		femPtr->elem->rsidu[ ielem ] = atof( args[ 7 ] );
		femPtr->elem->rsidf[ ielem ] = atof( args[ 8 ] );
		femPtr->elem->reduc[ ielem ] = atof( args[ 9 ] );
		femPtr->elem->redup[ ielem ] = atof( args[ 10 ] );
		femPtr->elem->reduf[ ielem ] = atof( args[ 11 ] );
		femPtr->elem->mnunx[ ielem ] = atof( args[ 12 ] );
		femPtr->elem->pcrit[ ielem ] = atof( args[ 13 ] );
		femPtr->elem->rsidp[ ielem ] = atof( args[ 14 ] );
		}
	fclose( infile ); 
//	printf( "uniax = %e\t", femPtr->elem->uniax0[ 0 ] );
//	femPtr->elem->tprin[ 64 ] = ( 0.25 * PI - 0.5 * femPtr->elem->frict0[ 64 ] ); //--- princ direction of D
//	femPtr->elem->musnd[ 64 ] = 0.0; //--- set modulus to zero
	
	femPtr->elem->meanGmodu = femPtr->elem->gmodu[ 0 ];	
	femPtr->elem->meanKmodu = femPtr->elem->kmodu[ 0 ];	
	femPtr->elem->meanRecoveryStrain = femPtr->elem->strny[ 0 ];
}

//---------------------------------------------------------
//--- parse the string
void MatProp::parse()
{
	char *ptr;
	ptr = strtok( line, " " );
	int narg = 0;
	while( ptr != NULL ){
		strcpy( args[ narg ], ptr );
		ptr = strtok( NULL, " " );
		narg++;
	}
}
