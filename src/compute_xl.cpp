#include "fem.h"
#include "compute_x.h"
#include "math.h"

ComputeXL::ComputeXL( Fem *femObj ):
elem( femObj->elem ), cxl( 0 )
{}
//-----------------------------------------------------------
ComputeXL::~ComputeXL()
{
	if( cxl )
		fclose( outFile );
}
//-----------------------------------------------------------
void ComputeXL::GetXL( bool state )
{

	unsigned int kgaus = 0, itype;
	double uniax2, effectiveStress, yfunc, theta;//, s[ elem->nstre ];
	if( state )
		fprintf( outFile, "TIME\nstart\n" );
	else
		fprintf( outFile, "TIME\nend\n" );
	fprintf( outFile, "NELEM\n%d\n", elem->nelem * elem->ngaus );
	fprintf( outFile, "ID\tX\n" );
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		itype = elem->ltype[ ielem ];
//		uniax2 = elem->uniax0[ itype ];// * elem->uniax0[ itype ];
		uniax2 = elem->uniax0[ itype ] * elem->uniax0[ itype ];
		for( int igaus = 0; igaus < elem->ngaus; igaus ++ )
		{
			effectiveStress = Invart( elem->strsi[ kgaus ] );	
			theta = atan2( elem->strsi[ kgaus ][ 2 ], 0.5 * ( elem->strsi[ kgaus ][0] - elem->strsi[ kgaus ][ 1 ] ) );
			assert(  uniax2 - effectiveStress * effectiveStress * pow( cos( theta ), 2 ) >= 0.0 );
			yfunc = sqrt( uniax2 - effectiveStress * effectiveStress * pow( cos( theta ), 2 ) ) - effectiveStress * sin( theta );
/*			s[ 0 ] = elem->strsi[ kgaus ][ 0 ];
			s[ 1 ] = elem->strsi[ kgaus ][ 1 ];
			s[ 2 ] = elem->strsi[ kgaus ][ 2 ] + yfunc;
			assert( fabs( pow( Invart( s ), 2 )  - uniax2 ) < 1.0e-03 );
*/
//			yfunc = effectiveStress - uniax2;
//			fprintf( outFile, "%8.7e\n", - yfunc ); // yfunc );// / uniax );
			fprintf( outFile, "%d\t%8.7e\n", kgaus, yfunc ); // yfunc );// / uniax );
			fflush( outFile );
			kgaus++;
		}
	}
}
//-----------------------------------------------------------
double ComputeXL::Invart(  double *u )
{
// --- invariants of tensor T
	double smean, dev[ 2 ][ 2 ], varj2;
	smean = Smean( u );
	dev[ 0 ][ 0 ] = u[ 0 ] - smean;
	dev[ 0 ][ 1 ] = u[ 2 ];
	dev[ 1 ][ 0 ] = u[ 2 ];
	dev[ 1 ][ 1 ] = u[ 1 ] - smean;
	varj2 = dev[ 0 ][ 1 ] * dev[ 1 ][ 0 ] + \
	0.5 * ( dev[ 0 ][ 0 ] * dev[ 0 ][ 0 ] + dev[ 1 ][ 1 ] * dev[ 1 ][ 1 ] );
	return sqrt( varj2 );
}
//-----------------------------------------------------------------------
inline double ComputeXL::Smean(  double *stress )
{
        return 0.5 * ( stress[ 0 ] + stress[ 1 ] );
}

