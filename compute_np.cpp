#include "compute_np.h"
#include "fem.h"

//--- constructor
ComputeNP::ComputeNP( Fem *femObj ):
femPtr( femObj )

{}
//--- deconstructor
ComputeNP::~ComputeNP()
{
	if( cnp )
	{
		fclose( outFile );
	}
}
//--- member funcs
void ComputeNP::GetActMap( unsigned int avlID )
{
	fprintf( outFile, "AVALANCHE ID\n%d\n", avlID );
	fprintf( outFile, "ELEMENT ID\tACTIVITY\n" ); //--- shape and dissipation rate
	for( unsigned int kgaus = 0; kgaus < femPtr->elem->nelem * femPtr->elem->ngaus; kgaus++ )
//		if( femPtr->elem->actFieldAccumulated[ kgaus ] > 0 )
//			fprintf( outFile, "%d\t%g\n", kgaus, ceil( femPtr->elem->actFieldAccumulated[ kgaus ] * femPtr->ts->dt / femPtr->fv->cdrag1 ) );
		fprintf( outFile, "%d\t%d\n", kgaus, femPtr->elem->actFieldAccumulated[ kgaus ] );
	fflush( outFile );
}

