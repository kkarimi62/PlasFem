#include "compute_dt.h"
#include "fem.h"

//--- constructor
ComputeDT::ComputeDT( Fem *femPtr ):
count( 0 ), node( femPtr->node )
{}
//--- deconstructor
ComputeDT::~ComputeDT()
{
	if( cdt )
	{
		fclose( outFile );
		delete [] epdot;
		delete [] dsigma;
		delete [] time;
		delete [] fx;
		delete [] fy;
		delete [] kin;
		delete [] velos;
	}
}
//--- member funcs
void ComputeDT::GetDuration( unsigned int duration )
{
	if( duration >= tmin && duration < tmax )
	{
		fprintf( outFile, "Duration\n%d\n", duration );
		fprintf( outFile, "t\tS(t)\tV(t)\tep(t)\tk(t)\n" ); //--- shape and dissipation rate
		for( unsigned int itime = 0; itime < count; itime++ )
			fprintf( outFile, "%d\t%e\t%e\t%e\t%e\n", time[ itime ], 
			dsigma[ itime ], epdot[ itime ], velos[ itime ], kin[ itime ] );
		fflush( outFile );
	}
	//--- initialize
	for( unsigned int itime = 0; itime < count; itime++ )
	{
		epdot[ itime ] = 0.0;
		dsigma[ itime ] = 0.0;
		velos[ itime ] = 0.0;
//		fx[ itime ] = 0.0;
//		fy[ itime ] = 0.0;
		kin[ itime ] = 0.0;
		time[ itime ] = 0;
	}
	count = 0;
}

