//#include "Python.h"
#include "thermo_style.h"
#include "fem.h"
#include <cstring>
#include <cassert>
#include "compute_ke.h"
#include "compute_sxy.h"
#include "iostream"
#include "gauss_coord.h"

using std::cout;

ThermoStyle::ThermoStyle( Fem *femObj ):
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ), 
ts( femObj->ts ), fd( femObj->fd ), output( femObj->logFile ),
cSxyObj( new ComputeSxy ), cKEobj( new ComputeKE )
{}

ThermoStyle::~ThermoStyle()
{
	delete cSxyObj;
	delete cKEobj;
}

void ThermoStyle::Init()
{
	cSxyObj->Init( elem->stran, elem->strst, elem->plasticStrain, elem->vol, 
		       elem->nelem, elem->ngaus );
	cKEobj->Init( node->mass, node->velot, node->nsvab,
		      domain->ndime );
}


void ThermoStyle::Process( unsigned int itime )
{

	if( itime == ts->itime0 )
	{
		for( int iarg = 0; iarg < narg; iarg++ )
		{
			if ( !strcmp( args[ iarg ], "pe" ) ) 
				fprintf( output, "PotEng\t" );
			if ( !strcmp( args[ iarg ], "step" ) )
				fprintf( output, "Step\t" );
			if ( !strcmp( args[ iarg ], "ke" ) ) 
				fprintf( output, "KinEng\t" );
			if ( !strcmp( args[ iarg ], "lohi" ) )
				fprintf( output, "Box\t" );
			if ( !strcmp( args[ iarg ], "sxy" ) )
				fprintf( output, "Sxy\t" );
			if ( !strcmp( args[ iarg ], "sxx" ) )
				fprintf( output, "Sxx\t" );
			if ( !strcmp( args[ iarg ], "sigma" ) )
				fprintf( output, "sigma\t" );
			if ( !strcmp( args[ iarg ], "syy" ) )
				fprintf( output, "Syy\t" );
			if ( !strcmp( args[ iarg ], "p" ) )
				fprintf( output, "Press\t" );
			if ( !strcmp( args[ iarg ], "ev" ) )
				fprintf( output, "Eps\t" );
			if ( !strcmp( args[ iarg ], "n" ) )
				fprintf( output, "n\t" );
		}
	fprintf( output, "\n" );
	}

	// --- write to std out
	Init();
	for( int iarg = 0; iarg < narg; iarg++ )
	{
		if ( !strcmp( args[ iarg ], "pe" ) )
		{ 
//				cout << femObj->get_potnt() << '\t';
		}
		if ( !strcmp( args[ iarg ], "step" ) ) 
				fprintf( output, "%d\t", itime );
		if ( !strcmp( args[ iarg ], "sxy" ) )
		{ 
			fprintf( output, "%e\t", cSxyObj->GetSxy() );
		}
		if ( !strcmp( args[ iarg ], "sxx" ) )
		{ 
			fprintf( output, "%e\t", cSxyObj->GetSxx() );
		}
		if ( !strcmp( args[ iarg ], "syy" ) )
		{ 
			fprintf( output, "%e\t", cSxyObj->GetSyy() );
		}
		if ( !strcmp( args[ iarg ], "sigma" ) )
		{ 
			fprintf( output, "%e\t", cSxyObj->GetSigma() );
		}
		if ( !strcmp( args[ iarg ], "p" ) )
		{ 
			fprintf( output, "%e\t", cSxyObj->GetP() );
		}
		if ( !strcmp( args[ iarg ], "ev" ) )
		{ 
			fprintf( output, "%e\t", cSxyObj->GetEv() );
		}
		if ( !strcmp( args[ iarg ], "n" ) )
		{ 
			kount = 0;
//			#pragma omp parallel
			{
//			#pragma omp for
			for( unsigned int kgaus = 0; kgaus < elem->nelem * elem->ngaus; kgaus++ )
//				#pragma omp atomic
				kount += elem->actField[ kgaus ];
			}
			fprintf( output, "%g\t", ( 1.0 * kount ) / ( elem->nelem * elem->ngaus ) );
		}
		if ( !strcmp( args[ iarg ], "ke" ) )
		{ 
			fprintf( output, "%8.7e\t", cKEobj->GetKE() / node->nsvab );
		}
		if ( !strcmp( args[ iarg ], "lohi" ) )
			fprintf( output, "%e\t", itime * ts->dt * fd->shearRATE );
	}
	fprintf( output, "\n" );
}
