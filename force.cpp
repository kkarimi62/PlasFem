#include <iostream>
#include "force.h"
#include "fem.h"
#include "domain.h"
#include "math.h"
#include "mapping.h"
#include "stdlib.h" // exit, rand,srand
#include "matrix.h"
#include "time.h" // time
#include "cassert" // assert
#include <omp.h>
#include <Eigen/Sparse> //--- store hessian

using std::cout;
using Eigen::SparseMatrix;
using Eigen::VectorXi;
using Eigen::VectorXd;

Force::Force( Fem *femObj ):
femPtr( femObj ), memObj( new Memory ), 
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ), 
fv( femObj->fv ), ts( femObj->ts ), fd( femObj->fd ), yc( femObj->yc ),
yr( femObj->yr ), ps( femObj->ps ), rv( femObj->rv ), wr( femObj->wr )
//,outpt( fopen( "epsilonPlastic.txt", "w" ) ) //--- to measure plastic strain during active phase
{
//	fprintf( outpt, "ID\texx\teyy\texy\n" );
}
//--------------------------------------------------------------------
Force::~Force()
{
         memObj->Delete2ndMat( lohi, 2 );
         for( int iThread = 0; iThread < nThreads; iThread++ )
         {
                 memObj->Delete3rdMat( bmatx[ iThread ], elem->nelem * elem->ngaus, elem->nstre  );
                 memObj->Delete3rdMat( bt[ iThread ], elem->nelem * elem->ngaus, elem->nevab );
         }
         delete mpObj[ 0 ];
         delete mpObj[ 1 ];
         delete memObj;
	 delete stiff;
//	 fclose( outpt );
}
//--------------------------------------------------------------------
void Force::Init()
{
//	val = new double[ elem->ngaus * elem->nelem * elem->nstre ]; //---store plastic strain
	//--- update bmatx, bt, elem->volu
	double **posgp, weigp[ elem->ngaus ]; //--- position of gauss points, quadrature points: weight
//	static SparseMatrix<double> stiff( node->nsvab, node->nsvab );
	stiff = new SparseMatrix<double>( node->nsvab, node->nsvab );
	( *stiff ).reserve( VectorXi::Constant( node->nsvab, 10 * domain->ndime ) ); //--- max. number of neighbors

	memObj->doubleMat( posgp, elem->ngaus, domain->ndime );
	memObj->doubleMat( lohi, 2, 5 );
	for( int iThread = 0; iThread < nThreads; iThread++ ) // --- l-values for this class
	{
		bmatx[ iThread ] = new double**[ elem->nelem * elem->ngaus ];
		bt[ iThread ]    = new double**[ elem->nelem * elem->ngaus ];
		for( int igaus = 0; igaus < elem->nelem * elem->ngaus; igaus++ )
		{
			memObj->doubleMat( bmatx[ iThread ][ igaus ], elem->nstre, elem->nevab ); 
			memObj->doubleMat( bt[ iThread ][ igaus ], elem->nevab, elem->nstre ); 
		}	
	}
	// --- box bounds
	double dGamma = ts->dt * fd->shearRATE; //--- shear strain
	double halfDGamma = 0.5 * dGamma; //--- shear strain
	lohi[ 0 ][ 0 ] = domain->xlo0;
	lohi[ 0 ][ 1 ] = domain->xhi0;
	lohi[ 0 ][ 2 ] = domain->ylo0;
	lohi[ 0 ][ 3 ] = domain->yhi0;
	lohi[ 0 ][ 4 ] = domain->xy0;
//---
	lohi[ 1 ][ 0 ] = domain->xlo0;
	lohi[ 1 ][ 1 ] = domain->xhi0;
	lohi[ 1 ][ 2 ] = domain->ylo0;
	lohi[ 1 ][ 3 ] = domain->yhi0;
	lohi[ 1 ][ 4 ] = domain->xy0 + dGamma * ( domain->yhi0 - domain->ylo0 );
	if( ps->shear ) //--- pure shear
	{
		//--- update xy, xlo, xhi, ylo, yhi
		lohi[ 1 ][ 4 ] = domain->xy0;
		lohi[ 1 ][ 0 ] = domain->xlo0 + 0.5 * halfDGamma * ( domain->xlo0 - domain->xhi0 );
		lohi[ 1 ][ 1 ] = domain->xhi0 + 0.5 * halfDGamma * ( domain->xhi0 - domain->xlo0 );
		lohi[ 1 ][ 2 ] = domain->ylo0 - 0.5 * halfDGamma * ( domain->ylo0 - domain->yhi0 );
		lohi[ 1 ][ 3 ] = domain->yhi0 - 0.5 * halfDGamma * ( domain->yhi0 - domain->ylo0 );
	}
	// --- some initialization
	cordPtr[ 0 ] = node->cordi;
	cordPtr[ 1 ] = node->coord;
	mpObj[ 0 ] = new Wrap( lohi[ 0 ] );
	mpObj[ 1 ] = new Wrap( lohi[ 1 ] );
		
	// --- set up element obj
	elem->Init();
	elem->GAUSSQ( posgp, weigp ); // gauss points: pos and weight
	#pragma omp parallel num_threads( nThreads )
	{
	// --- these are local to each thread!
	int tid = omp_get_thread_num();
	double **elcod, //--- nodal xyz
	       **cartd, //--- shape funcs derivatives
	       elcod0[ elem->nevab ];	
	memObj->doubleMat( cartd, domain->ndime, elem->nnode ); 
	memObj->doubleMat( elcod, elem->nnode, domain->ndime ); 
	double exisp, etasp; //--- gauss points positions
	unsigned int kgaus = 0, ievab; //--- counter for gauss points
	double djacb, dvolu; //--- determinant of jacobian
	unsigned itype;
	double kmodu, gmodu;
	double Dmatx[ elem->nstre * elem->nstre ], **estif;
	double **Db, **btDb;
	//--- allocate memory
	memObj->doubleMat( estif, elem->nevab, elem->nevab );
	memObj->doubleMat( Db, elem->nstre, elem->nevab );
	memObj->doubleMat( btDb, elem->nevab, elem->nevab );
	// --- element loop starts: update vol[ kgaus ], bmat[ kgaus ] , bt[ kgaus ]
	// --- don't put omp for here!
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{ 
		// --- computes elasticity matrix
		itype = elem->ltype[ ielem ];
		kmodu = elem->kmodu[ itype ];
		gmodu = elem->gmodu[ itype ];
		Modulus( Dmatx, kmodu, gmodu, elem->nstre );
		GetNodalCords( elcod0, ielem, 1, 0 );// , tid ); // xyz with pbc( 1 ) and initial frame( 0 )
		ievab = 0;
		for( int inode = 0; inode < elem->nnode; inode++ )
		{
			for( int idime = 0; idime < domain->ndime; idime++ )
			{
				elcod[ inode ][ idime ] = elcod0[ ievab ];
				ievab++;
			}
		}
//		if( ielem == 0 )
//		{
//			for( int i = 0; i < elem->nnode; i++ )
//			{
//				unsigned int nodei = elem->lnods[ ielem ][ i ];
//				printf( "x[ %d ]\t", i );
//				for( int j = 0; j < domain->ndime; j++ )
//					printf( "%g\t", elcod0[ i * 2 + j ] );
//				for( int j = 0; j < domain->ndime; j++ )
//					printf( "coord = %g\t", node->coord[ nodei ][ j ] );
//				printf( "\n" );
//			}
//		}
		//--- initialize estif( local stiff. mat.)
		for( int ievab = 0; ievab < elem->nevab; ievab++ )
		{
			for( int jevab = 0; jevab < elem->nevab; jevab++ )
				estif[ ievab ][ jevab ] = 0.0;
		}
		//--- loop over quad points starts
		for( int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			kgaus = ielem * elem->ngaus + igaus;
			exisp = posgp[ igaus ][ 0 ]; 
			etasp = posgp[ igaus ][ 1 ];
			elem->Jacob2( cartd, djacb, exisp, etasp, elcod, ielem ); //--- output cartd, djacb
			dvolu = ( djacb ) * weigp[ igaus ];
			elem->vol[ kgaus ] = dvolu;
			elem->Blarge( bmatx[ tid ][ kgaus ], cartd ); //--- output bmatx
			matPtr->Transpose( bmatx[ tid ][ kgaus ], elem->nstre, elem->nevab, bt[ tid ][ kgaus ] ); //--- output bt
			matPtr->MulMatMat( Dmatx, elem->nstre, elem->nstre, bmatx[ tid ][ kgaus ], elem->nstre, elem->nevab, Db ); //--- D . b
			matPtr->MulMatMat( bt[ tid ][ kgaus ], elem->nevab, elem->nstre, Db, elem->nstre, elem->nevab, btDb ); //--- bt * D * b
			for( int ievab = 0; ievab < elem->nevab; ievab++ )
			{
				for( int jevab = 0; jevab < elem->nevab; jevab++ )
					estif[ ievab ][ jevab ] += btDb[ ievab ][ jevab ] * dvolu;
			}
		} //--- end of quad. loop
		//--- assembly
		unsigned int ipoin, isvab, jpoin, jsvab, jevab;
		ievab = 0;
		for( int inode = 0; inode < elem->nnode; inode++ )
		{
			ipoin = elem->lnods[ ielem ][ inode ];
			isvab = ipoin * domain->ndime;
			for( int idofn = 0; idofn < domain->ndime; idofn++ )
			{
				jevab = 0;
				for( int jnode = 0; jnode < elem->nnode; jnode++ )
				{
					jpoin = elem->lnods[ ielem ][ jnode ];
					jsvab = jpoin * domain->ndime;
					for( int jdofn = 0; jdofn < domain->ndime; jdofn++ )
					{
//						stiff[ isvab ][ jsvab ] += estif[ ievab ][ jevab ]; //--- global hessian
						if( tid == 0 )
							( *stiff ).coeffRef( isvab, jsvab ) += estif[ ievab ][ jevab ];
						jevab++;
						jsvab++;
					}
				}
				ievab++;
				isvab++;
			}
		}
	} // --- end of element loop
//	printf( "%d\thello!\n", tid );
	if( tid == 0 )
	{
		( *stiff ).makeCompressed();
		( femPtr->cg )->compute( *stiff );
//		printf( "k=%g\n", stiff.coeffRef( 0, 0 ) );
	}
	memObj->Delete2ndMat( cartd, domain->ndime );
	memObj->Delete2ndMat( elcod, elem->nnode );
	memObj->Delete2ndMat( estif, elem->nevab );
	memObj->Delete2ndMat( Db, elem->nstre );
	memObj->Delete2ndMat( btDb, elem->nevab );
	} // --- end of par. section
	memObj->Delete2ndMat( posgp, elem->ngaus );
//	printf( "vol[10] = %g\n", elem->vol[ 10 ] );
//	exit( EXIT_FAILURE );
}
//--------------------------------------------------------------------
void Force::ComputeForce( int counter, unsigned int itime )
{
//	xxx = counter; 
        int deform = 0; //--- q-st
        if( counter == 1 )
                deform = 1;
	//--- local to each thread!
	double	bt_sigma[ elem->nevab ], eload[ elem->nevab ], //--- b^t * sigma & element load vector
		dldis[ elem->nevab ], //--- element dipslacements
		veloc[ elem->nevab ];
	double kmodu, //--- healing strain, K
	       gmodu, //--- G, yield stress
	       mufst, musnd, tprin; //--- anisotropic shear moduli in damage & eigen-direction
	double Dmatx[ elem->nstre * elem->nstre ], Vmatx[ elem->nstre * elem->nstre ]; 
	double strag[ elem->nstre ],//--- strain
	       dsigm[ elem->nstre ], //--- delta sigma
	       sgtri[ elem->nstre ], //--- trial stress = s_n+D * strain
	       sigmV[ elem->nstre ]; //--- viscous stress
	unsigned int kgaus, nodei, isvab, //--- gauss point counter, node id, dof's
		     ievab, itype; //--- material id
	int tid = omp_get_thread_num();
	//--- initialize node->fintl
	unsigned int chunck = node->nsvab / nThreads;
	#pragma omp for schedule( dynamic, chunck )
	for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ ) 
	{
		node->fintl[ isvab ] = 0.0;
		node->absfl[ isvab ] = 0.0;
	}
	// --- element loop starts
	chunck = elem->nelem / nThreads;
	visc_coeff = fv->cdrag0; //--- not fixed!
	second_visc_coeff = 1.0 / fv->cdrag1;
	#pragma omp for schedule( dynamic, chunck )
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		// --- initialize eload
		for( int iEvab = 0; iEvab < elem->nevab; iEvab++ )
			eload[ iEvab ] = 0.0;
		if( deform ) //---
			GetDisp( dldis, ielem, deform ); // --- output dldis
		else
		{
		ievab = 0;
		for( int inode = 0; inode < elem->nnode; inode++ )
		{
			nodei = elem->lnods[ ielem ][ inode ];
			isvab = domain->ndime * nodei;
			for( int idofn = 0; idofn < domain->ndime; idofn++ )
			{
				dldis[ ievab ] = node->dldis[ isvab ];  // update original coords
				isvab++;
				ievab++;
			}
		}
		}
		GetVeloc( veloc, ielem ); // --- output veloc
		// --- computes elasticity matrix
		itype = elem->ltype[ ielem ];
		kmodu = elem->kmodu[ itype ];
		gmodu = elem->gmodu[ itype ];
		Modulus( Dmatx, kmodu, gmodu, elem->nstre );
/*		if( yr->nlaps == 4 and yc->ncrit == 1 ) //--- damage model: anisotropic moduli for mohr-col plasticity model
		{
			mufst = elem->mufst[ itype ]; //--- 1st shear mod
			musnd = elem->musnd[ itype ]; //--- 2nd shear mod
			tprin = elem->tprin[ itype ]; //--- princ. dir
			AnisoMod( Dmatx, kmodu, mufst, musnd, tprin );
		}
*/
		Modulus( Vmatx, visc_coeff, visc_coeff, elem->nstre ); //--- visc. matrix 
		// loop over quad points starts
		for( int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			kgaus = ielem * elem->ngaus + igaus ;
			// --- stress and plastic strain update
			matPtr->MulMatVec( bmatx[ tid ][ kgaus ], elem->nstre, elem->nevab, 
					   dldis, elem->nevab, strag ); // strain = b * u
//			printf( "strag = %g\n", strag[ 2 ] );
			matPtr->MulMatVec( bmatx[ tid ][ kgaus ], elem->nstre, elem->nevab, 
					   veloc, elem->nevab, elem->stranRate[ kgaus ] ); // strain_rate = b * v
			matPtr->MulMatVec( Dmatx, elem->nstre, elem->nstre, 
					   strag, elem->nstre, dsigm ); // dsigma = D * strain
			matPtr->MulMatVec( Vmatx, elem->nstre, elem->nstre, 
					   elem->stranRate[ kgaus ], elem->nstre, sigmV ); //--- viscous sigma
			for( int istre = 0; istre < elem->nstre; istre++ )
			{
//				elem->strst[ kgaus ][ istre ] = elem->strsi[ kgaus ][ istre ]; // initialize stress
				sgtri[ istre ] = dsigm[ istre ] + elem->strst[ kgaus ][ istre ]; // trial stress
				elem->strnt[ kgaus ][ istre ] += strag[ istre ]; //--- total strain 
//				elem->plasticStrain[ kgaus ][ istre ] = elem->plasticStrain0[ kgaus ][ istre ]; // initialize plastic strain
			}
/*			if( yr->nlaps == 4 ) //--- damage model: no update incrementally
			{
				matPtr->MulMatVec( Dmatx, elem->nstre, elem->nstre, 
						   elem->strnt[ kgaus ], elem->nstre, sgtri ); //--- sigma = D * strain
				for( int istre = 0; istre < elem->nstre; istre++ )
					sgtri[ istre ] += elem->strrs[ kgaus ][ istre ];
			}
*/
			ReturnStress( elem->strst[ kgaus ], elem->stran[ kgaus ], //--- update???
				      elem->plasticStrain[ kgaus ], &( elem->actField[ kgaus ] ), 
				      ielem, kgaus, sgtri, strag, Dmatx ); //--- update strst, stran, plasticStrain 

		//--- CALCULATE THE EQUIVALENT NODAL FORCES
			for( int istre = 0; istre < elem->nstre; istre++ ) 
				sigmV[ istre ] = sigmV[ istre ] * ( 1 - elem->actField[ kgaus ] ) + elem->strst[ kgaus ][ istre ]; //--- only elastic elemnts (kelvin type)
//			matPtr->MulMatVec( bt[ tid ][ kgaus ], elem->nevab, elem->nstre, 
//					   elem->strst[ kgaus ], elem->nstre, bt_sigma ); // ---bt * sigma
			matPtr->MulMatVec( bt[ tid ][ kgaus ], elem->nevab, elem->nstre, 
					   sigmV, elem->nstre, bt_sigma ); // ---bt * sigma
			for( int ievab = 0; ievab < elem->nevab; ievab++ )
				eload[ ievab ] += bt_sigma[ ievab ] *  elem->vol[ kgaus ];
		} // gauss points loop ends
		ievab = 0;
		for( int inode = 0; inode < elem->nnode; inode++ )
		{
			nodei = elem->lnods[ ielem ][ inode ];
			isvab = domain->ndime * nodei;
			for( int idofn = 0; idofn < domain->ndime; idofn++ )
			{
				#pragma omp atomic
				node->fintl[ isvab ] += eload[ ievab ];
				#pragma omp atomic
				node->absfl[ isvab ] += fabs( eload[ ievab ] ); //--- force mag.
//				fintl[ isvab ] += eload[ ievab ];
				isvab++;
				ievab++;
			}
		}
	} //--- element loop ends
	#pragma omp for
	for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ ) //--- net force / force 
	{
		if( node->absfl[ isvab ] != 0.0 )
			node->absfl[ isvab ] = node->fintl[ isvab ] / node->absfl[ isvab ];
	}

/*	if( tid == 0 ) //--- compute plastic strain
	{
		unsigned int kounter0 = 0;
		double strain0[ elem->nstre ];
		for( unsigned int kgaus = 0; kgaus < elem->ngaus * elem->nelem; kgaus++ )
		{
			for( int istre = 0; istre < elem->nstre; istre++ )
				strain0[ istre ] = val[ kgaus * elem->nstre + istre ];
			if( EffectiveStrain( strain0 ) > 1.0e-08 )
				kounter0++;
		}
		if( kounter0 > 0 )
		{
//		fprintf( outpt, "ITIME\n%d\n", itime );
//		fprintf( outpt, "NELEM\n%d\n", kounter0 );
//		fprintf( outpt, "ID\texx\teyy\texy\n" );
		for( unsigned int kgaus = 0; kgaus < elem->ngaus * elem->nelem; kgaus++ )
		{
			for( int istre = 0; istre < elem->nstre; istre++ )
				strain0[ istre ] = val[ kgaus * elem->nstre + istre ];
			if( EffectiveStrain( strain0 ) > 1.0e-08 )
			{
				fprintf( outpt, "%d\t", kgaus );
				for( int istre = 0; istre < elem->nstre; istre++ )
				{
					fprintf( outpt, "%6.5e\t", val[ kgaus * elem->nstre + istre ] ); //--- distribution of ep?
					val[ kgaus * elem->nstre + istre ] = 0.0;
				}
				fprintf( outpt, "\n" );
			}
		}
		}
	}
*/

//	exit( EXIT_FAILURE );
}
//-----------------------------------------------------------------------
double Force::Smean(  double *stress )
{
	return 0.5 * ( stress[ 0 ] + stress[ 1 ] );
}
//-----------------------------------------------------------------------
inline double Force::DotPr( double *avect, double *bvect, unsigned int nsize )
{
	double sum = 0.0;
	for( unsigned int i = 0; i < nsize; i++ )
		sum += avect[ i ] * bvect[ i ];
	return sum;
}
//-----------------------------------------------------------------------
inline void Force::DeviaStrain( double *dev /*output*/,  double *T )
{
	dev[ 0 ] = T[ 0 ] - T[ 1 ];
	dev[ 1 ] = - dev[ 0 ];
	dev[ 2 ] = T[ 2 ]; // [ ( exx - eyy ), ( eyy - exx ), 2.0 * exy ]
}
//-----------------------------------------------------------------------
inline double Force::EffectiveStrain(  double *strain ) //, int tid )
{
	double exx_yy, exy;
	exx_yy = strain[ 0 ];
	exy = 0.5 * strain[ 2 ];
	return sqrt( exx_yy * exx_yy + 4.0 * exy * exy ); // gamma^2 = 2 e_dev : e_dev
}
//-----------------------------------------------------------------------
inline void Force::Devia( double *dev /*output*/,  double *T ) //, int tid )
{
	double smean;
	smean = Smean( T );
	dev[ 0 ] = T[ 0 ] - smean;
	dev[ 1 ] = T[ 1 ] - smean;
	dev[ 2 ] = T[ 2 ];
}
//-----------------------------------------------------------------------
double Force::Invart(  double *u ) //, int tid )
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
void Force::ReturnStress( double *stres, double *stran, double *plasticStrain, int *actField /*output*/,  unsigned int ielem,
			  unsigned int kgaus, double *sgtri, double *strag, double *Dmatx )
{
	double uniax = elem->uniax[ ielem ];
	double uniax0 = elem->inunx[ ielem ]; //--- initial yield stress
	double FRICT = elem->frict[ ielem ]; //--- friction angle
	double FRICT0 = elem->infrc[ ielem ]; //--- initial friction angle
	double gmodu = elem->gmodu[ ielem ];
	double kmodu = elem->kmodu[ ielem ];
	double strny = elem->strny[ ielem ];
	double rsidu = elem->rsidu[ ielem ]; //--- residual yield stress for weakening
	double rsidf = elem->rsidf[ ielem ]; //--- residual friction for weakening
	double rsidp = elem->rsidp[ ielem ]; //--- residual pressure cap
	double reduc = elem->reduc[ ielem ]; //--- reduction in uniax
	double redup = elem->redup[ ielem ]; //--- reduction in cap press
	double reduf = elem->reduf[ ielem ]; //--- reduction in friction
	double mnunx = elem->mnunx[ ielem ]; //--- average yield stress
	double pcrit = elem->pcrit[ ielem ]; //--- critical pressure
	double pcrit0 = elem->pcrit0[ ielem ]; //--- initial critical pressure
	double epdev[ elem->nstre ]; //--- deviatoric strain rate
	double srdev[ elem->nstre ]; //--- deviatoric strain
	double stdev[ elem->nstre ]; //--- deviatoric stress
	double dldev[ elem->nstre ]; //--- deviatoric stress
	double steff = Invart( sgtri ); //--- max shear stress: trial stress	
	double pstre = Smean( sgtri ); //--- trial pressure
	double yfunc; //--- yield criterion von-mises
	double epbar = elem->epstn[ kgaus ]; //--- plastic strain
	double hards = elem->hards[ ielem ]; //--- hardening modulus
	double dlstr[ elem->nstre ]; //--- stress decrement
	double vectb[ elem->nstre ]; //--- D * df/ds
	double stidv[ elem->nstre ]; //--- deviatoic stress at tn
	double dlprs; //--- pressure increment
	double vecta[ elem->nstre ]; //--- df/dsigma f: yield function
	bool recvy = 0; //--- bool variable: recovery
	double tprfr, tprsc; //--- principal stress dirs
	double *tpptr = &( elem->tprin[ ielem ] ); //--- ptr to modulus princ dir
	double *muptr0 = &( elem->mufst[ ielem ] ); //--- ptr to 2nd shear mod
	double *muptr = &( elem->musnd[ ielem ] ); //--- ptr to 2nd shear mod
	double inuni; //--- initial yield stress
	double infri; //--- initial friction
	double inpcr; //--- initial pressure cap
	double presd;
	double sigrs = 0.25 * uniax0 * cos( FRICT0 ); //kam ( 0.25 * uniax * cos( FRICT ) >= 1.0e-02 ? 0.25 * uniax * cos( FRICT ) : 1.0e-02 );
//---	
//	double trelx = visc_coeff / gmodu; //--- relaxation time: plastic elements have the same viscous properties as elastic elements
//	double viscm = visc_coeff; //--- viscosity
	double trelx = 1.0 / second_visc_coeff; //--- relaxation time == tau_p
	double viscm = gmodu * trelx; //--- viscosity
//---
	double preys = uniax + epbar * hards; //--- yield stress
	if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-coulomb
	{
		presd = ( sigrs - uniax * cos( FRICT ) ) / sin( FRICT ); //--- pressure at residual shear
		preys = - pstre * sin( FRICT ) + uniax * cos( FRICT ) + epbar * hards; //--- yield stress //---pstre??? trial or at n
//		if( pstre > uniax / tan( FRICT ) ) //--- bi-linear yield function
//			preys = 1.0e-08; //--- a residual ???
		if( - pstre < presd ) //---kam bi-linear yield function
			preys = sigrs; //--- a residual ???
	}
	yfunc = steff - preys; //--- yield function: current stress - yield stress
	elem->fyild[ kgaus ] = yfunc; //--- fy rate of change	
//---
	bool yield = 0;
	yield = ( yfunc >= 0.0 ) or *actField; //--- yfunc >= 0.0 ( elastic to plastic ) or actField == 1 (plastic-plastic)
//	if( kgaus == 0 )
//		printf("kgaus=%d\tuniax=%e\n", kgaus, uniax );
	if( yield and !fv->switched ) //--- fv->switched = 0 (linear solver off)
	{
//		printf("counter = %d, kgaus = %d\n", xxx, kgaus );
//		if( kgaus == 64 )
//			printf("kgaus=%d\tyfunc=%e\tuniax=%e\n", kgaus, yfunc, uniax );
//		exit( EXIT_FAILURE );
		*actField = 1; //--- elastic element becomes plastic!
		DeviaStrain( epdev, elem->stranRate[ kgaus ] ); //--- deviatoric strain rate to be used for the stress update in plastic phase
		DeviaStrain( srdev, strag ); //--- deviatoric strain 
		Devia( stdev, sgtri ); //--- trial deviatoic stress
		Devia( stidv, stres ); //--- deviatoic stress at tn
		double sigma = 0.5 * ( stres[ 0 ] - stres[ 1 ] ); //--- shear stress @ tn
		double sigxy = stres[ 2 ]; //--- 2nd shear stress
		double taumx = sqrt( sigma * sigma + sigxy * sigxy ); //--- max shear stress
		if( yr->nlaps == 0 ) //--- plastic rule
		{
			//--- compute dlmda: plastic strain magnitude
			vecta[ 0 ] = 0.5 * sigma / taumx;
			vecta[ 1 ] = - 0.5 * sigma / taumx;
			vecta[ 2 ] = sigxy / taumx; //--- df/dsigma f: yield function sxy == 0.5(sxy+syx)
			if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-coulomb
			{
				vecta[ 0 ] += 0.5 * sin( FRICT ); //--- add pressure contribution
				vecta[ 1 ] += 0.5 * sin( FRICT );
			}
			matPtr->MulMatVec( Dmatx, elem->nstre, elem->nstre, 
					   vecta, elem->nstre, vectb ); //-- vectb = D * vecta
			double dlmda = yfunc / ( hards + DotPr( vecta, vectb, elem->nstre ) ); //--- plastic strain magnitude
//			assert( dlmda >= 0.0 ); //--- inertia!!!
			epbar += dlmda; //--- total plastic strain	
			//--- compute change in stress
			for( int istre = 0; istre < elem->nstre; istre++ )
				dlstr[ istre ] = dlmda * vectb[ istre ]; //--- stress decrement
		}
		else if( yr->nlaps == 1 ) //--- visco-elastic maxwell
		{
			for( int istre = 0; istre < elem->nstre; istre++ )
				dlstr[ istre ] = stidv[ istre ] * ts->dt / trelx; //--- maxwel fluid
		}	
		//--- compute trace and deviatoric part in ds
		dlprs = Smean( dlstr ); //--- pressure decrement: dlprs
		Devia( dldev, dlstr ); //--- deviatoric decrement: dldev
//---
		for( int istre = 0; istre < elem->nstre; istre++ ) //--- update deviatoic stress
		{
			if( yr->nlaps == 0 or yr->nlaps == 1 ) //--- plastic rule or visco-elastic maxwell
				stdev[ istre ] -= dldev[ istre ]; //--- relax deviatoric part
			else if( yr->nlaps == 2 )
				stdev[ istre ] = viscm * epdev[ istre ]; //--- Newtonian fluid 
			else if( yr->nlaps == 3 ) // or yr->nlaps == 4 ) //--- void or damage
				stdev[ istre ] = 0.0; //--- void
//			plasticStrain[ istre ] += fabs( ts->dt * epdev[ istre ] ); //--- plastic strain rate is equal the total deviatoric strain rate
//                      plasticStrain[ istre ] += fabs( ts->dt * dev_strsi[ istre ] / viscm );
//---
			if( rv->nrecv == 0 )
				stran[ istre ] += 1.0; //--- time-based recovery
			else if ( rv->nrecv == 1 )
				stran[ istre ] += fabs( srdev[ istre ] ); //dt * dev_strsi[ istre ] / visc_coeff; //--- strain-based: dev. part of strag
		}
		if( yr->nlaps == 0 ) //--- plastic rule
			pstre -= dlprs; //--- relax p only for plastic rule
//		else if( yr->nlaps == 3 and yc->ncrit == 1 ) //--- void and mc no dilatancy
//		{
//			if( - pstre < 0.0 ) //--- negative pressure not allowed
//				pstre = 0.0;
//		}
		else if( yr->nlaps == 3 and yc->ncrit == 2 ) //--- void and mc with dilatancy
		{
			if( elem->resid[ kgaus ] == 0 )
			{
//				pstre -= taumx * kmodu / gmodu ; //-- pressure increases! 
				pstre -= taumx * sin( FRICT ) ; //-- pressure increases! normality rule
				if( -pstre > pcrit )
					pstre = -pcrit; //--- pressure cap
//				if( - pstre < 0.0 ) //--- negative pressure not allowed
//					pstre = 0.0;
				elem->presd[ kgaus ] = pstre;
				elem->resid[ kgaus ] = 1;
			}
			else
				pstre = elem->presd[ kgaus ]; //--- retain initial pressure
		}
//---
		stres[ 0 ] = stdev[ 0 ] + pstre; //--- total stress: dev + hydro-static
		stres[ 1 ] = stdev[ 1 ] + pstre;
		stres[ 2 ] = stdev[ 2 ];
		if( yr->nlaps == 4 ) //--- damage model
		{	
			if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-col plasticity: anisotropy
			{
				StrPrinDir( stres, &tprfr, &tprsc );
				*tpptr = tprfr + ( 0.25 * PI - 0.5 * FRICT ); //--- princ direction of D
//				*tpptr = ( 0.25 * PI - 0.5 * FRICT ); //--- princ direction of D
				*muptr = 0.0; //--- set modulus to zero
//				elem->gmodu[ ielem ] *= 0.9;
//				*muptr0 = 2.0 * elem->gmodu[ ielem ];
			}
			else if( yc->ncrit == 0 ) //--- von-mises: isotropic damage
				elem->gmodu[ ielem ] *= 0.9;	
		}
	}
	else //--- elastic loading 
	{
		for( int istre = 0; istre < elem->nstre; istre++ )
			stres[ istre ] = sgtri[ istre ]; //--- when element remains elastic
	}
	//--- recovery criterion
	if( rv->nrecv == 0 ) //--- time-based
		recvy = *actField and //--- if the element is plastic
			stran[ 0 ] * ts->dt * second_visc_coeff >= 1.0;// and //--- elapsed time is longer than the threshold
//			yfunc < 0.0; //--- stres less than threshold
	else if( rv->nrecv == 1 ) //--- strain-based
		recvy = *actField and EffectiveStrain( stran ) >= strny; //--- strain-based
	else if( rv->nrecv == 2 ) //--- stress-based recovery
	{
		steff = Invart( stres );
		if( yc->ncrit == 0 ) //--- von-mises
			preys = 0.1 * uniax; //--- yield stress: 0.1 hard-coded
		else if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-coulomb
			preys = 0.1 * ( - pstre * sin( FRICT ) + uniax * cos( FRICT ) ); //--- yield stress 
		yfunc = steff - preys; //--- yield function: current stress - yield stress
		recvy = *actField and yfunc < 0.0; //--- stress is less than some threshold
	}
	else if( rv->nrecv == 3 ) //--- damage model
	{
		assert( yr->nlaps == 4 );
		recvy = *actField;
	}
	if ( recvy ) 
	{
//		if( kgaus == 157 )
//		printf("kgaus=%d\tyfunc=%e\n", kgaus, yfunc );
		*actField = 0;
		epbar = 0.0;
		elem->resid[ kgaus ] = 0;
		for( int istre = 0; istre < elem->nstre; istre++ )
		{
//			val[ kgaus * elem->nstre + istre ] = plasticStrain[ istre ]; //--- plastic strain
			plasticStrain[ istre ] = 0.0; //--- becomes elastic again!
			stran[ istre ] = 0.0;
		}
//		do //--- while loop! be careful (might run forever!)
//		{
			if( wr->nweak == 0 ) //--- no weakening
			{
				elem->set( ielem, mnunx, rsidu ); //--- new material parameters: update uniax, frict
				uniax = elem->uniax[ ielem ];
				FRICT = elem->frict[ ielem ];
			}
			else if ( wr->nweak == 1 or wr->nweak == 2 )//--- weakening is on! 1 or 2 -> proper tilt of the initial failure envelope
			{
				inuni = uniax0; //--- initial values: non-accumulated weakening
				infri = FRICT0;
				inpcr = pcrit0;
				if( wr->nweak == 1 ) //--- accumulated weakening
				{
					inuni = uniax; //--- current values
					infri = FRICT;
					inpcr = pcrit;
				}
				if( yc->ncrit == 0 ) //--- von-mises
				{
					uniax = ( 1.0 - reduc ) * inuni; //--- no do-while
					if( uniax < rsidu ) //--- residual yield stress
						uniax = inuni;
				}
				else if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-coulomb
				{
					FRICT = ( 1.0 - reduf ) * infri; 
					if( FRICT < rsidf ) //--- residual yield stress
						FRICT = infri;
//kam					uniax = uniax0 * tan( FRICT ) / tan( FRICT0 ); //--- (p,0): fixed point of the failure envelope
//
					uniax = ( 1.0 - reduc ) * inuni; //--- (p,0): fixed point of the failure envelope
					if( uniax < rsidu ) //--- residual yield stress
						uniax = inuni;
//
					pcrit = ( 1.0 - redup ) * inpcr; 
					if( pcrit < rsidp ) //--- residual pressure cap
						pcrit = inpcr;
				}
			}		
			else if ( ( yc->ncrit == 1 or yc->ncrit == 2 ) and wr->nweak == 3 )//--- weakening is on! decrease non-accumulated cohesion only!
			{
				inuni = uniax0; //--- initial values: non-accumulated weakening
				uniax = ( 1.0 - reduc ) * inuni; //--- no do-while
				if( uniax < rsidu ) //--- residual yield stress
					uniax = inuni;
			}
//			else if( wr->nweak == 3 ) //--- non-accumulative weakening: weakened initial distribution ( only cohesion )
//			{
//				elem->set( ielem, ( 1.0 - reduc ) * mnunx, rsidu ); //--- new material parameters: update uniax, frict
//				uniax = elem->uniax[ ielem ];
//				FRICT = elem->frict[ ielem ];
//			}
			steff = Invart( stres );
			yfunc = steff - uniax;
			if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-coulomb
				yfunc = steff - ( - pstre * sin( FRICT ) + uniax * cos( FRICT ) ); //--- assert yfunc < 0
//			printf("counter = %d, kgaus = %d\t%f\t%f\t%f\n", xxx, kgaus, - pstre, uniax, yfunc );
//		} 
//		while( yfunc >= 0.0 );
//		assert( yfunc < 0.0 ); //--- below yield surface
		elem->uniax[ ielem ] = uniax; //--- assignment
		elem->frict[ ielem ] = FRICT;
		elem->pcrit[ ielem ] = pcrit;
	}
//---	 
	elem->epstn[ kgaus ] = epbar; //--- update plastic strain
//	if( kgaus == 157 )
//		printf("end of force!\n" );
}
//-----------------------------------------------------------------------
void Force::StrPrinDir( double *stres, double *tprfr, double *tprsc )
{
	double press = -0.5 * ( stres[ 0 ] + stres[ 1 ] ); //--- pressure
	double sigma = 0.5 * ( stres[ 0 ] - stres[ 1 ] ); //--- axial shear
	double sigxy = stres[ 2 ]; //--- diagonal shear
	double sprfr, sprsc; //--- princ stresses
	double shtfr, shtsc; //--- shear tractions
//---
	*tprfr = 0.5 * atan2( sigxy, sigma ); //--- stress princ. dir
	*tprsc = 0.5 * atan2( -sigxy, -sigma );
//	printf( "tprfr = %e\n", *tprfr * 180.0 / pi );
//	printf( "tprsc = %e\n", *tprsc * 180.0 / pi );
	sprfr = press + sigma * cos( 2.0 * *tprfr ) + sigxy * sin( 2.0 * *tprfr ); //--- 1st princ stress
	sprsc = press + sigma * cos( 2.0 * *tprsc ) + sigxy * sin( 2.0 * *tprsc ); //--- 2nd princ stress
	assert( sprfr >= sprsc );
//	printf( "sprfr = %e\n", sprfr );
//	printf( "sprsc = %e\n", sprsc );
//	shtfr = - sigma * sin( 2.0 * tprfr ) + sigxy * cos( 2.0 * tprfr ); //--- shear tract on princ. dir
//	shtsc = - sigma * sin( 2.0 * tprsc ) + sigxy * cos( 2.0 * tprsc );
//	printf( "shtfr = %e\n", shtfr );
//	printf( "shtsc = %e\n", shtsc );
}
//-----------------------------------------------------------------------
inline void Force::Modulus( double *Dmatx, 
		      double kmod, 
		      double gmod, 
		      int    nstre )
{	
	Dmatx[  0 ] = kmod + gmod;
	Dmatx[  1 ] = kmod - gmod;
	Dmatx[  2 ] = 0.0;

	Dmatx[ 3 ] = kmod - gmod;
	Dmatx[ 4 ] = kmod + gmod;
	Dmatx[ 5 ] = 0.0;
		
	Dmatx[ 6 ] = 0.0;
	Dmatx[ 7 ] = 0.0;
	Dmatx[ 8 ] = 1.0 * gmod;
}
//-----------------------------------------------------------------------
void Force::AnisoMod( double *Dmatx, double kmodu, double mufst, 
		      double musnd, double tprin )
//---- build anisotropic elasticity matrix 
{	
	const int ndime = domain->ndime;
	const int nevab = ndime * ndime;
	double uvect[ ndime ], uvecp[ ndime ]; //--- unit vectors at tprin
	double ndilt[ nevab ], naxsh[ nevab ], ndgsh[ nevab ]; //--- dilation + 2 shears
	double dmodu[ nevab * nevab ]; //--- moduli
	unsigned int ievab = 0, jevab = 0, isvab = 0;

//--- compute eigen-directions
	uvect[ 0 ] = cos( tprin );
	uvect[ 1 ] = sin( tprin ); //--- eigen-direction
	uvecp[ 0 ] = -sin( tprin );
	uvecp[ 1 ] = cos( tprin ); //--- eigen-direction
	for( int idime = 0; idime < ndime; idime++ )
	{
        	for( int jdime = 0; jdime < ndime; jdime++ )
        	{
                	ndilt[ ievab ] = uvect[ idime ] * uvect[ jdime ] + uvecp[ idime ] * uvecp[ jdime ]; //--- dilation
                	naxsh[ ievab ] = uvect[ idime ] * uvect[ jdime ] - uvecp[ idime ] * uvecp[ jdime ]; //--- axial shear: 1st eigen direction
                	ndgsh[ ievab ] = - ( uvect[ idime ] * uvecp[ jdime ] + uvecp[ idime ] * uvect[ jdime ] ); //--- diagl shear: 2nd eigen direction
                	ievab++;
        	}
	}	
//--- construct D
	ievab = 0;
	jevab = 0;
	isvab = 0;
	for( int idime = 0; idime < ndime; idime++ )
	{
		for( int jdime = 0; jdime < ndime; jdime++ )
		{
			ievab = idime * ndime + jdime;
			for( int kdime = 0; kdime < ndime; kdime++ )
			{
				for( int ldime = 0; ldime < ndime; ldime++ )
				{
					jevab = kdime * ndime +  ldime;
					dmodu[ isvab ] = kmodu * ndilt[ ievab ] * ndilt[ jevab ] + 
						         mufst * naxsh[ ievab ] * naxsh[ jevab ] + 
						         musnd * ndgsh[ ievab ] * ndgsh[ jevab ];
/*					if( idime == 0 and jdime == 0 and kdime == 0 and ldime == 0 )
						printf("0000 = %d\n", isvab );
					if( idime == 0 and jdime == 0 and kdime == 1 and ldime == 1 )
						printf("0011 = %d\n", isvab );
					if( idime == 0 and jdime == 0 and kdime == 0 and ldime == 1 )
						printf("0001 = %d\n", isvab );

					if( idime == 1 and jdime == 1 and kdime == 0 and ldime == 0 )
						printf("1100 = %d\n", isvab );
					if( idime == 1 and jdime == 1 and kdime == 1 and ldime == 1 )
						printf("1111 = %d\n", isvab );
					if( idime == 1 and jdime == 1 and kdime == 0 and ldime == 1 )
						printf("1101 = %d\n", isvab );

					if( idime == 0 and jdime == 1 and kdime == 0 and ldime == 0 )
						printf("0100 = %d\n", isvab );
					if( idime == 0 and jdime == 1 and kdime == 1 and ldime == 1 )
						printf("0111 = %d\n", isvab );
					if( idime == 0 and jdime == 1 and kdime == 0 and ldime == 1 )
						printf("0101 = %d\n", isvab );
*/
					isvab++;
				}
			}
		}
	}
//--- assign to Dmatx
	Dmatx[ 0 ] = dmodu[ 0 ]; //kmod + gmod;
	Dmatx[ 1 ] = dmodu[ 3 ]; //kmod - gmod;
	Dmatx[ 2 ] = dmodu[ 1 ];
//---
	Dmatx[ 3 ] = dmodu[ 12 ]; //kmod - gmod;
	Dmatx[ 4 ] = dmodu[ 15 ]; //kmod + gmod;
	Dmatx[ 5 ] = dmodu[ 13 ]; //0.0;
//---		
	Dmatx[ 6 ] = dmodu[ 4 ]; //0.0;
	Dmatx[ 7 ] = dmodu[ 7 ]; //0.0;
	Dmatx[ 8 ] = dmodu[ 5 ]; //1.0 * gmod;
}
//-----------------------------------------------------------------------
void Force::GetVeloc( double *veloc,  unsigned int ielem )//, int tid )
{
	unsigned int nodei, isvab, ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		nodei = elem->lnods[ ielem ][ inode ];
		isvab = domain->ndime * nodei;
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			veloc[ ievab ] = node->velot[ isvab ]; // displacements
			isvab++;
			ievab++;
		}
	}
}
//-----------------------------------------------------------------------
void Force::GetDisp( double *dldis,  unsigned int ielem, int deform )
{
// --- find proper displacements for stress update
	double elcod0[ elem->nevab ], elcod[ elem->nevab ];
// --- compute proper dldis
	GetNodalCords( elcod0, ielem, 0, 0 ); // xyz in initial frame( 0 ): no pbc( 0 )
	unsigned int ievab = 0, isvab = 0, nodei;	
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		nodei = elem->lnods[ ielem ][ inode ];
		isvab = domain->ndime * nodei;
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			dldis[ ievab ] = node->dldis[ isvab ]; // displacements
			elcod0[ ievab ] += dldis[ ievab ];  // update original coords
			isvab++;
			ievab++;
		}
	}
//--- current frame
	ievab = 0;
	double Cords[ 2 * domain->ndime ];
	double r, theta;
	int cross_pbc;
	for( int idime = 0; idime < domain->ndime; idime++ )
		Cords[ domain->ndime + idime ] = elcod0[ idime ]; //--- cords[ 1 ] = elcod0[ 0 ]
	for( int jnode = 0; jnode <  elem->nnode; jnode++ )
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			Cords[ idime ] = elcod0[ ievab ]; //--- cords[ 0 ] = elcod[ inode ]
			ievab++;
		}
		Distance( Cords, deform, r, theta, cross_pbc ); // find coords in the current frame
		elcod0[ ievab - domain->ndime ] = elcod0[ 0 ] + ( r ) * cos( theta );
		elcod0[ ievab - 1 ]     = elcod0[ 1 ] + ( r ) * sin( theta );
	}	
	// --- proper xyz respecting pbc
	GetNodalCords( elcod, ielem, 1, 0 );
	ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			dldis[ ievab ] = elcod0[ ievab ] - elcod[ ievab ];
			ievab++;
		}
	}
}
//-----------------------------------------------------------------------
void Force::GetNodalCords( double *elcod,  unsigned int ielem, int pbc, int initial )//, int tid )
{
// --- returns elcod: nodal coordinates of the element

	unsigned int nodei;
	double Cords[ 2 * domain->ndime ];
	double r, theta;
	int cross_pbc;
	unsigned int ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ ) // loop over nodes
	{
		nodei = elem->lnods[ ielem ][ inode ]; // node global id
		for( int idime = 0; idime < domain->ndime; idime++ ) // assign coordantes
		{
			elcod[ ievab ] = cordPtr[ initial ][ nodei ][ idime ];
			if( inode == 0 )
				Cords[ domain->ndime + idime ] = elcod[ idime ]; //--- cords[ 1 ] = elcod[ 0 ]
			Cords[ idime ] = elcod[ ievab ]; //-- cords[ 0 ] = elcod[ inode ]
//			printf("xNod=%g\n", elcod[ ievab ] );
			ievab++;
		}
		if( pbc ) // --- returns xyz relative to node0 
		{
			Distance( Cords, initial, r, theta, cross_pbc ); // find the true distance
			elcod[ ievab - domain->ndime ] = elcod[ 0 ] + ( r ) * cos( theta );
			elcod[ ievab - 1 ]     = elcod[ 1 ] + ( r ) * sin( theta );
		}
	}

}
//-----------------------------------------------------------------------
void Force::Distance(  double *coord,  int initial, double &r, double &theta, int &cross_pbc )//, int tid )
{
	// --- get s0, s1
	double dimensioless_cords[ 2 * domain->ndime ];
	double s0i, s1i, s0j, s1j, s0_ij, s1_ij;
	double xi, yi, xj, yj,
	       xij, yij,
	       rsq;

	mpObj[ initial ]->GetDimensionlessCords( 2, coord, dimensioless_cords ); // output dimension...

	s0i = dimensioless_cords[ 0 ];
	s1i = dimensioless_cords[ 1 ];
	s0j = dimensioless_cords[ 1 * domain->ndime ];
	s1j = dimensioless_cords[ 1 * domain->ndime + 1 ];
	s0_ij = s0j - s0i; // dimensionless distance
	s1_ij = s1j - s1i;

	// --- translate point j to find the proper distance rij
	cross_pbc = 0;
        if( s0_ij > 0.5 )
	{ 
		s0j -= 1.0;
		cross_pbc = 1;
	}
        else if ( s0_ij < - 0.5 )
	{ 
		s0j += 1.0;
		cross_pbc = 1;
        }
	if( s1_ij > 0.5 )
	{
		s1j -= 1.0;
		cross_pbc = 1;
        }
	else if( s1_ij < - 0.5 )
	{
		s1j += 1.0;
		cross_pbc = 1;
	}
	// --- find new coordinates
	dimensioless_cords[ 0 ] = s0i;
	dimensioless_cords[ 1 ] = s1i;
	dimensioless_cords[ 1 * domain->ndime ] = s0j;
	dimensioless_cords[ 1 * domain->ndime + 1 ] = s1j;
	double scord[ domain->ndime * 2 ];
	mpObj[ initial ]->GetCords( 2, dimensioless_cords, scord ); // --- update coord

	xi = scord[ 0 ];
	yi = scord[ 1 ];
        xj = scord[ 1 * domain->ndime + 0 ];
        yj = scord[ 1 * domain->ndime + 1 ];
        xij = xi - xj;
        yij = yi - yj;
        rsq = xij * xij + yij * yij;
        r = sqrt( rsq );
	theta = atan2( yij, xij );
}

