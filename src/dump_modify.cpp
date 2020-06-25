//#include "Python.h"
#include "dump_modify.h"
#include "fem.h"
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <string>

#define MAXDUMP 100

using std::cout;
using std::string;

DumpModify::DumpModify( Fem *femObj ):
//node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ), 
//ts( femObj->ts ), fd( femObj->fd ), femPtr( femObj ),
//gcObj( new GaussCoord( femObj ) ),
nDump( 0 )//, memPtr( new Memory ) 
{
	exist = new bool[ MAXDUMP ];
	for( int i =0; i < MAXDUMP;i++)
		exist[ i ] = false; 

}

DumpModify::~DumpModify()
{}
