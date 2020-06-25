#include "fem.h"
#include "compute_x.h"


ComputeX::ComputeX( Fem *femObj ):
kount( 0 ), init( 0 ), ds( 0.0 ), duration( 0 ),
cx( 0 ), cxy( femObj->ths->cSxyObj ),
ps( femObj->ps )
{}
//-----------------------------------------------------------
ComputeX::~ComputeX()
{
	if( cx )
		fclose( outFile );
}
//-----------------------------------------------------------
void ComputeX::GetX()
{
	double sdr;
	if( kount == 0 )
	{
		fprintf( outFile, "#DURATION\tSIZE\n" );
		a = ( ps->shear ? 0.5 * ( cxy->GetSxx() - cxy->GetSyy() ) : cxy->GetSxy() );
	}
	else if( kount == 1 )
		b = ( ps->shear ? 0.5 * ( cxy->GetSxx() - cxy->GetSyy() ) : cxy->GetSxy() );
	else
	{
		c = ( ps->shear ? 0.5 * ( cxy->GetSxx() - cxy->GetSyy() ) : cxy->GetSxy() );
		sdr = 0.5 * ( c - a ); //--- derivative
		if( sdr <= 0.0 ) //--- plasticity
			 init = kount + 1; //--- initialize init
		else //--- elastic loading!		
		{
			ds += sdr;
			duration++;
		}
		if( init - kount == 1 && duration != 0 ) //--- avalanche ends!
		{
			//--- output t, s
			fprintf( outFile, "%d\t%8.7e\n", duration, ds );
			fflush( outFile );
//			printf( "avalanch: ds = %g\n", -ds );
			ds = 0.0;
			duration = 0;
		}
		a = b;
		b = c;
	}
	kount++;
}
