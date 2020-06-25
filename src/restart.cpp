#include "restart.h"
#include "fem.h"
//#include "std"

Restart::Restart( Fem *femObj ):
n( 0 ), femPtr( femObj )
{}

Restart::~Restart(){}

void Restart::output( unsigned int itime )
{
	someFile = fopen( "restart_data.txt", "w" );

	// --- wirte data.txt
	fprintf( someFile, "#fem description\n\n" );
	fprintf( someFile, "%d\tinitial time step\n", itime );
	fprintf( someFile, "%d\tdimensions\n", femPtr->domain->ndime );
	fprintf( someFile, "%d\tpoints\n", femPtr->node->npoin );
	fprintf( someFile, "%d\telements\n", femPtr->elem->nelem );
	fprintf( someFile, "%d\t3-noded elements\n", femPtr->elem->nnode );
	fprintf( someFile, "%d\tgauss points\n", femPtr->elem->ngaus );
	fprintf( someFile, "%d\tsize of stress/strain vector\n\n", femPtr->elem->nstre );
	fprintf( someFile, "%d\telement types (elastic constants\n\n", femPtr->elem->nelem );
	fprintf( someFile, "%15.10g\t%15.10g\txlo\txhi\n", femPtr->domain->xlo, femPtr->domain->xhi ); 
	fprintf( someFile, "%15.10g\t%15.10g\tylo\tyhi\n", femPtr->domain->ylo, femPtr->domain->yhi ); 
	fprintf( someFile, "%15.10g\txy\n\n", femPtr->domain->xy0 ); // initial xy!
	// --- write xyz
	fprintf( someFile, "Nodes: id mass x y\n\n" );
	for( unsigned int ipoin = 0; ipoin < femPtr->node->npoin; ipoin++ ) 
	{
        	fprintf( someFile, "%d\t%15.10g\t", ipoin, femPtr->node->mass[ ipoin ] );
        	for( int idime = 0; idime < femPtr->domain->ndime; idime++ )
               		fprintf( someFile, "%15.10g\t", femPtr->node->coord[ ipoin ][ idime ] );
		fprintf( someFile, "\n" );
	}               	
	fprintf( someFile, "\n" );

	// --- write velocities
//	fprintf( someFile, "\n\n" );
	fprintf( someFile, "Velocities\n\n" );
	isvab = 0;
	for( unsigned int ipoin = 0; ipoin < femPtr->node->npoin; ipoin++ ) 
	{
        	fprintf( someFile, "%d\t", ipoin );
        	for( int idime = 0; idime < femPtr->domain->ndime; idime++ )
		{
               		fprintf( someFile, "%15.10g\t", femPtr->node->velot[ isvab ] );
			isvab++;
		}
		fprintf( someFile, "\n" );
	}               	
	fprintf( someFile, "\n" );

	// --- write lnods
//	fprintf( someFile, "\n\n" );
	fprintf( someFile, "Connectivity: id type nodes\n\n" );
	for( unsigned int ielem = 0; ielem < femPtr->elem->nelem; ielem++ ) 
	{
        	fprintf( someFile, "%d\t%d\t", ielem, ielem );
        	for( int inode = 0; inode < femPtr->elem->nnode; inode++ )
               		fprintf( someFile, "%d\t", femPtr->elem->lnods[ ielem ][ inode ] );
		fprintf( someFile, "\n" );
	}
	fprintf( someFile, "\n" );
	// --- write lelms
//	fprintf( someFile, "\n\n" );
	fprintf( someFile, "# Cells: rows x cols\n" );
	fprintf( someFile, "%d\t%d\n", femPtr->tr->mrows, femPtr->tr->ncols );
	fprintf( someFile, "# INDEX NUMBER ID\n" );
	unsigned int ncels = femPtr->tr->mrows * femPtr->tr->ncols; 
	for( unsigned int icels = 0; icels < ncels; icels++ )
	{
        	fprintf( someFile, "%d\t%d\t", icels, femPtr->tr->cellTriangCount[ icels ] );
		for( int itrig = 0; itrig < femPtr->tr->cellTriangCount[ icels ]; itrig++ ) 
               		fprintf( someFile, "%d\t", femPtr->tr->lelms[ icels ][ itrig ] );
		fprintf( someFile, "\n" );
	}
	fprintf( someFile, "\n" );

	// --- write stress
	fprintf( someFile, "initial stress\n\n" );
	kgaus = 0;
	for( unsigned int ielem = 0; ielem < femPtr->elem->nelem; ielem++ ) 
	{
		for( int igaus = 0; igaus < femPtr->elem->ngaus; igaus++ ) 
		{
        		fprintf( someFile, "%d\t", kgaus );
        		for( int istre = 0; istre < femPtr->elem->nstre; istre++ )
               			fprintf( someFile, "%15.14e\t", femPtr->elem->strst[ kgaus ][ istre ] );
			fprintf( someFile, "\n" );
			kgaus++;
		}
	}
	fprintf( someFile, "\n" );

	// --- write plastic strain
	fprintf( someFile, "plastic strain\n\n" );
	kgaus = 0;
	for( unsigned int ielem = 0; ielem < femPtr->elem->nelem; ielem++ ) 
	{
		for( int igaus = 0; igaus < femPtr->elem->ngaus; igaus++ )
		{ 
        		fprintf( someFile, "%d\t", kgaus );
        		for( int istre = 0; istre < femPtr->elem->nstre; istre++ )
               			fprintf( someFile, "%15.10g\t", femPtr->elem->plasticStrain[ kgaus ][ istre ] );
			fprintf( someFile, "\n" );
			kgaus++;
		}
	}
	fprintf( someFile, "\n" );

	// --- write forces
	fprintf( someFile, "internal force\n\n" );
	isvab = 0;
	for( unsigned int ipoin = 0; ipoin < femPtr->node->npoin; ipoin++ ) 
	{
        	fprintf( someFile, "%d\t", ipoin );
        	for( int idime = 0; idime < femPtr->domain->ndime; idime++ )
		{
               		fprintf( someFile, "%15.10g\t", femPtr->node->fintl[ isvab ] );
			isvab++;
		}
		fprintf( someFile, "\n" );
	}               	
	fprintf( someFile, "\n" );

	// --- write strain
	fprintf( someFile, "strain\n\n" );
	kgaus = 0;
	for( unsigned int ielem = 0; ielem < femPtr->elem->nelem; ielem++ ) 
	{
		for( int igaus = 0; igaus < femPtr->elem->ngaus; igaus++ ) 
		{
        		fprintf( someFile, "%d\t", kgaus );
        		for( int istre = 0; istre < femPtr->elem->nstre; istre++ )
               			fprintf( someFile, "%15.10g\t", femPtr->elem->stran[ kgaus ][ istre ] );
			fprintf( someFile, "\n" );
			kgaus++;
		}
	}

	// --- write act
	fprintf( someFile, "activity\n\n" );
	kgaus = 0;
	for( unsigned int ielem = 0; ielem < femPtr->elem->nelem; ielem++ ) 
	{
		for( int igaus = 0; igaus < femPtr->elem->ngaus; igaus++ ) 
		{
        		fprintf( someFile, "%d\t", kgaus );
               		fprintf( someFile, "%d\t", femPtr->elem->actField[ kgaus ] );
			fprintf( someFile, "\n" );
			kgaus++;
		}
	}

	fclose( someFile );

	// --- write moduli
	someFile = fopen( "restart_moduli.txt", "w" );
	fprintf( someFile, "#type\tbulk mod\tshear mod\tyield stress\trecovery strain\n" );
	for( unsigned int ielem = 0; ielem < femPtr->elem->nelem; ielem++ ) 
		fprintf( someFile, "%d\t%15.10g\t%15.10g\t%15.10g\t%15.10g\n", 
			 ielem, femPtr->elem->kmodu[ ielem ], femPtr->elem->gmodu[ ielem ], 
				femPtr->elem->uniax[ ielem ], femPtr->elem->strny[ ielem ] );
	fclose( someFile );

}
