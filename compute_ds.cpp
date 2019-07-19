#include "fem.h"
#include "compute_ds.h"


ComputeDS::ComputeDS( Fem *femObj ):
kount( 0 ), init( 0 ), ds( 0.0 ), duration( 0 ), intDur( 0 ),
cxy( femObj->ths->cSxyObj ), cds( 0 ), femPtr( femObj ), navch( 0 ), n( 1 ), ps( femObj->ps ),
dp( femObj->dp ), qs( femObj->qs ), tr( femObj->tr )
{}
//-----------------------------------------------------------
ComputeDS::~ComputeDS()
{
	if( cds )
		fclose( outFile );
}
//-----------------------------------------------------------
void ComputeDS::GetAvalanch( unsigned int internalDuration )
{
	double sdr;
	if( kount == 0 )
	{
		fprintf( outFile, "#DURATION\tSIZE\n" );
		a = ( ps->shear ? 0.5 * ( cxy->GetSxx() - cxy->GetSyy() ) : cxy->GetSxy() );
		ta = internalDuration;	
//---
		if ( tr->msd ) //--- print header
			fprintf( tr->outFile, "#ITIME\tDX^2\tDXDY\tDYDX\tDY^2\n" );
	}
	else if( kount == 1 )
	{
		b = ( ps->shear ? 0.5 * ( cxy->GetSxx() - cxy->GetSyy() ) : cxy->GetSxy() );
		tb = internalDuration;	
	}
	else
	{
		c = ( ps->shear ? 0.5 * ( cxy->GetSxx() - cxy->GetSyy() ) : cxy->GetSxy() );
		tc = internalDuration;	
		sdr = 0.5 * ( c - a ); //--- derivative
		if( sdr > 0.0 ) //--- elastic loading
			 init = kount + 1; //--- initialize init
		else //--- avalanche starts!		
		{
			if( duration == 0 )
			{
				for( int iDump = 0; iDump < dp->nDump and qs->min; iDump++ ) //--- if quasi-static
//				for( int iDump = 0; iDump < dp->nDump and cds; iDump++ )
				{
					if( ( kount + 1 ) % dp->n[ iDump ] == 0 )
						dp->Process( kount + 1, iDump );
				}
//---
				if( femPtr->cxl->cxl )	//--- compute xl (start of an avalanch)
					femPtr->cxl->GetXL( 1 );
//---				
				if ( tr->msd )
					tr->ComputeMSD( kount ); //--- msd
			}
			ds += sdr;
			intDur += tb;
			duration++;
		}
		if( init - kount == 1 && duration != 0 ) //--- avalanche ends!
		{
			if( femPtr->cxl->cxl )	//--- compute xl (end of an avalanch)
				femPtr->cxl->GetXL( 0 );
			if( femPtr->cnp->cnp ) //--- accumulated activity
				femPtr->cnp->GetActMap( navch );
			//--- output t, s
			if( cds )
			{
//				fprintf( outFile, "%d\t%8.7e\n", intDur, - ds ); //--- internal duration!
				fprintf( outFile, "%d\t%8.7e\n", duration, - ds );
				fflush( outFile );
				navch++;
			}
			ds = 0.0;
			duration = 0;
			intDur = 0;
			for( unsigned int kgaus = 0; kgaus < femPtr->elem->nelem * femPtr->elem->ngaus; kgaus++ )
				femPtr->elem->actFieldAccumulated[ kgaus ] = 0; //--- initialize
		}
		a = b;
		b = c;
		ta = tb;
		tb = tc;
	}
	kount++;
}
