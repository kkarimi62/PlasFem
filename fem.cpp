#include <iostream>
#include "fem.h"
#include "stdlib.h"
#include <cmath>
#include <omp.h>
#include <cassert>
#include <algorithm> //--- min/max

using std::max_element;
using Eigen::VectorXd;
using Eigen::ConjugateGradient;
using Eigen::SparseMatrix;
static bool absCompare( double a, double b ){ return fabs( a ) < fabs( b ); }

Fem::Fem()
{
	infile = stdin;
	logFile = stdout; //fopen( "log.fem", "w" );	
	domain = new Domain();
	node   = new Node();
	srand( time( NULL ) );
	ts	   = new TimeStep();
	fd	   = new FixDeform();
	fvl	   = new FixVeloc();
	ff	   = new FixForce();
	qs	   = new QuasiStatic();
	yc	   = new YieldCrit();
	yr	   = new YieldRule();
	wr	   = new WeakRule();
	fv	   = new FixViscous();
	run    = new Run();
	ps	   = new PureShear();
	rv	   = new Recovery();
//---
	elem   = new Elem( this );
	ths	   = new ThermoStyle( this );
	dpmd	   = new DumpModify( this );
	dp	   = new Dump( this );
	input  = new Input( this );
	force  = new Force( this );
	rs     = new Restart( this );
	tr     = new Tracer( this );
	cds     = new ComputeDS( this ); //--- must be allocated after ths
	cdt     = new ComputeDT( this );
	cnp     = new ComputeNP( this );
	cx     = new ComputeX( this );
	cxl     = new ComputeXL( this ); //--- must be allocated after ths
	cg	= new ConjugateGradient< SparseMatrix <double> >;
	
}
//---
Fem::~Fem()
{
	delete cg;
	delete cxl;
	delete cx;
	delete cnp;
	delete cdt;
	delete cds;
	delete tr;
	delete rs;
	delete force;
	delete input;
	delete dp;
	delete dpmd;
	delete ths;
	delete elem;
	delete rv;
	delete ps;
	delete run;
	delete fv;
	delete wr;
	delete yr;
	delete yc;
	delete qs;
	delete ff;
	delete fvl;
	delete fd;
	delete ts;
	delete node;
	delete domain;
	fclose( logFile );

	delete [] velol;
}
//------------------------------------------------------
//--- perfom a loop for numerical integration
//------------------------------------------------------
void Fem::Loop()
{
	unsigned int ntime = run->nstep;
	Init(); // initialize state variables			
	if( fvl->fix_vel )
		SetVeloc(); // set velocities for "wall" nodes
	if( ff->fix_force )
		SetForce(); // set velocities for "wall" nodes
	// ---  td output!
	ths->Process( ts->itime0 );
	// --- output data
	for( int iDump = 0; iDump < dp->nDump; iDump++ )
		dp->Process( 0, iDump );
	#pragma omp parallel num_threads( nThreads )
	{
	unsigned int counter = 0; //--- min. steps 
	for( int itime = 0; itime < ntime && cds->navch < cds->n; itime++ ) // time integration starts
	{ 
//		if( ! fvl->fix_vel && ! ff->fix_force ) 
			ApplyHOMO(); // --- apply homogeneous deformation
		while( Iterate( counter ) ) //--- q-st
		{
			Verlet( node->veloi, node->velot, ts->dt, ts->dt ); //---output node->veloi, node->dldis: nonaffine piece
		// --- computes nodal internal forces at x( n + 1 )
			force->ComputeForce( counter, itime ); //--- if counter == 1-> update boundaries
		// --- calculates v( n + 1 )
			Verlet( node->velot, node->veloi, ts->dt, 0.0 ); //output node->velot
			Update(); //--- update state variables
		}
		if( tr->msd ) //--- tracer particle
			tr->TracerUpdate(); //--- update tr->dispt ( use tr->tdisp ) and tr->coord 
		#pragma omp master
		{
			// --- std output!
			if ( ( itime + 1 ) % ths->n == 0 )
				ths->Process( ts->itime0 + itime + 1 );
			// --- output data
			for( int iDump = 0; iDump < dp->nDump and !qs->min; iDump++ )
//			for( int iDump = 0; iDump < dp->nDump and !cds->cds; iDump++ )
			{
				if( ( itime + 1 ) % dp->n[ iDump ] == 0 )
					dp->Process( itime + 1, iDump );
			}
			if(  rs->n != 0 ) //--- write restart
			{
				if( ( itime + 1 ) % rs->n == 0 ) 
					rs->output( ts->itime0 + itime + 1 );
			}
			if( cds->cds || cxl->cxl ) //--- avalanche
				cds->GetAvalanch( duration );
			if( cx->cx )	//--- compute x
				cx->GetX();
			 if( cdt->cdt )  //--- compute dt
				cdt->GetDuration( duration );
//			if( cxl->cxl && ( itime + 1 ) % cxl->nfreq == 0 )	//--- compute xl
//				cxl->GetXL();
			if ( tr->msd && ( itime + 1 ) % tr->nfreq == 0  and !cds->cds  )
				tr->ComputeMSD( itime + 1 ); //--- msd
		}
//		printf("rate=%e\n",((elem->strst[0][0]-elem->strst[0][1])-(elem->strsi[0][0]-elem->strsi[0][1]) )/ts->dt);
		#pragma omp barrier
		Switch(); //--- switch state variables
	} // end of time loop
	} // --- end of the par. sect.
}
//------------------------------------------------------
void Fem::Init()
{
	domain->xlo = domain->xlo0;
	domain->xhi = domain->xhi0;
	domain->ylo = domain->ylo0;
	domain->yhi = domain->yhi0;
	domain->xy = domain->xy0;
	#pragma omp parallel num_threads( nThreads ) 
	{
	//--- initialize
	unsigned int kgaus = 0;
	#pragma omp for
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		elem->uniax[ ielem ] = elem->uniax0[ ielem ];
		elem->frict[ ielem ] = elem->frict0[ ielem ];
		elem->infrc[ ielem ] = elem->frict0[ ielem ]; //--- saves initial friction angle (constant)
		elem->inunx[ ielem ] = elem->uniax0[ ielem ]; //--- saves initial yield stress (constant)
		elem->pcrit0[ ielem ] = elem->pcrit[ ielem ]; //--- saves initial yield stress (constant)
		for( int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			kgaus = ielem * elem->ngaus + igaus;
			elem->actField[ kgaus ] = elem->actField0[ kgaus ];
			for( int istre = 0; istre < elem->nstre; istre++ )
			{
				elem->strst[ kgaus ][ istre ] = elem->strsi[ kgaus ][ istre ];
				elem->strrs[ kgaus ][ istre ] = elem->strsi[ kgaus ][ istre ]; //--- residual stress (damage model)
				elem->strnt[ kgaus ][ istre ] = elem->strni[ kgaus ][ istre ];
				elem->fyild[ kgaus ] = elem->fyild0[ kgaus ]; //--- residual stress (damage model)
				elem->plasticStrain[ kgaus ][ istre ] = elem->plasticStrain0[ kgaus ][ istre ];
			}
		}
	}
	unsigned int isvab = 0;
	#pragma omp for
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{ 
		isvab = ipoin * domain->ndime;
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
	                node->coord[ ipoin ][ idime ] = node->cordi[ ipoin ][ idime ];
	                node->coordUnWrapped[ ipoin ][ idime ] = node->cordiUnWrapped[ ipoin ][ idime ];
			isvab = ipoin * domain->ndime + idime;
			node->dldis[ isvab ] = 0.0;
			node->dispt[ isvab ] = node->dispt0[ isvab ];
		}
	}
	}
	force->Init(); // --- initialize "force" member
	cg->setMaxIterations( 10000 ); //--- set max. iter and tol.
	cg->setTolerance( 1.0e-08 );
	if( tr->msd )
		tr->Init(); //--- init "tracer"
	effForce = new double[ domain->ndime ];
	velol    = new double[ node->nsvab ];
	b = new VectorXd( node->nsvab );
	x = new VectorXd( node->nsvab );
}
//------------------------------------------------------
void Fem::Update()
{
//	unsigned int isvab;
	#pragma omp for
	for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ ) // cordi = coord
	{
//		for( int idime = 0; idime < domain->ndime; idime++ )
//		{
//			isvab = ipoin * domain->ndime + idime;
			node->dispt[ isvab ] += node->dldis[ isvab ];
			tr->tdisp[ isvab ] += node->dldis[ isvab ]; //--- set zero in ApplyHomo
			node->dldis[ isvab ] = 0.0;
//	                node->cordi[ ipoin ][ idime ] = node->coord[ ipoin ][ idime ];
//	                node->cordiUnWrapped[ ipoin ][ idime ] = node->coordUnWrapped[ ipoin ][ idime ];
//		}
	}

	//--- accumulated activity
	#pragma omp for
	for( unsigned int kgaus = 0; kgaus < elem->nelem * elem->ngaus; kgaus++ )
        	elem->actFieldAccumulated[ kgaus ] +=  elem->actField[ kgaus ];

}
//------------------------------------------------------
void Fem::Switch()
{
	unsigned int kgaus = 0;
	#pragma omp for
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		elem->uniax0[ ielem ] = elem->uniax[ ielem ];
		elem->frict0[ ielem ] = elem->frict[ ielem ];
		for( int igaus = 0; igaus < elem->ngaus; igaus++ ) //strsi = strst
		{
			kgaus = ielem * elem->ngaus + igaus;
			elem->actField0[ kgaus ] = elem->actField[ kgaus ];
			for( int istre = 0; istre < elem->nstre; istre++ )
			{
				elem->strsi[ kgaus ][ istre ] = elem->strst[ kgaus ][ istre ];
				elem->strni[ kgaus ][ istre ] = elem->strnt[ kgaus ][ istre ];
				elem->fyild0[ kgaus ] = elem->fyild[ kgaus ]; 
				elem->plasticStrain0[ kgaus ][ istre ] = elem->plasticStrain[ kgaus ][ istre ];
			}
		}
	}
//--- total displacements
	#pragma omp for
	for( int isvab = 0; isvab < node->nsvab; isvab++ )
		node->dispt0[ isvab ] = node->dispt[ isvab ];
}
//------------------------------------------------------
void Fem::ApplyHOMO()
{
	unsigned int iposn, kgaus; //--- dof index
	double xcent, ycent; //--- center of the box
	double dGamma = ts->dt * fd->shearRATE; //--- shear strain
	double halfDGamma = 0.5 * dGamma; //--- shear strain
//---	
	domain->xy = domain->xy0 + dGamma * ( domain->yhi - domain->ylo ); //---simple shear
	if( ps->shear ) //--- pure shear
	{
		xcent = 0.5 * ( domain->xlo0 + domain->xhi0 );
		ycent = 0.5 * ( domain->ylo0 + domain->yhi0 );
		//--- update xy, xlo, xhi, ylo, yhi
		domain->xy = domain->xy0;
		domain->xlo = domain->xlo0 + 0.5 * halfDGamma * ( domain->xlo0 - domain->xhi0 );
		domain->xhi = domain->xhi0 + 0.5 * halfDGamma * ( domain->xhi0 - domain->xlo0 );
		domain->ylo = domain->ylo0 - 0.5 * halfDGamma * ( domain->ylo0 - domain->yhi0 );
		domain->yhi = domain->yhi0 - 0.5 * halfDGamma * ( domain->yhi0 - domain->ylo0 );
	}
	// --- update displacements
	// --- initialize dispt
	#pragma omp for
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
		iposn = ipoin * domain->ndime;
		node->dldis[ iposn ] = dGamma * ( node->coord[ ipoin ][ 1 ] - domain->ylo );
		node->dldis[ iposn + 1 ] = 0.0;
		if( ps->shear )
		{
			node->dldis[ iposn ] = halfDGamma * ( node->coord[ ipoin ][ 0 ] - xcent ); //--- ux ~ gamma (x-xc)
			node->dldis[ iposn + 1 ] = - halfDGamma * ( node->coord[ ipoin ][ 1 ] - ycent ); //--- uy ~ -gamma (y-yc)
		}
	}
	#pragma omp for
	for( int isvab = 0; isvab < node->nsvab; isvab++ )
		tr->tdisp[ isvab ] = 0.0; //-node->dldis[ isvab ]; //--- affine piece
//--- homogeneous defromation for tracers
/*	if( tr->msd )
	{
	#pragma omp for
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			kgaus = ielem * elem->ngaus + igaus;
			iposn = kgaus * domain->ndime;
			tr->disph[ iposn ] = dGamma * ( tr->coord[ kgaus ][ 1 ] - domain->ylo );
			tr->disph[ iposn + 1 ] = 0.0;
			if( ps->shear )
			{
				tr->disph[ iposn ] = halfDGamma * ( tr->coord[ kgaus ][ 0 ] - xcent ); //--- ux ~ gamma (x-xc)
				tr->disph[ iposn + 1 ] = - halfDGamma * ( tr->coord[ kgaus ][ 1 ] - ycent ); //--- uy ~ -gamma (y-yc)
			}
		}
	}
	}
*/
//--- end of "ApplyHomo"
}
//------------------------------------------------------
void Fem::Verlet( double *output, double *v, double dtf, double dtv )
{
	double dtfm;
	unsigned int chunck = node->nsvab / nThreads;
/*	if( ff->fix_force )
	{
	unsigned int isvab = 0;
        for( int idofn = 0; idofn < domain->ndime; idofn++ )
		effForce[ idofn ] = 0.0;
	#pragma omp barrier //--- because of effForce = 0
	//--- effective force
	#pragma omp for
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
        	for( int idofn = 0; idofn < domain->ndime; idofn++ )
        	{
			isvab = ipoin * domain->ndime + idofn;
        	        if( node->rigidID[ isvab ] == 1 )
			{
				#pragma omp atomic
        	                effForce[ idofn ] += node->fintl[ isvab ];
			}
        	}
	}	
	}
*/
	double cdrag = fv->critDrag;
//	if( fv->switched ) //--- only in critical damping limit! 
//		cdrag = fv->critDrag;
	#pragma omp for schedule( dynamic, chunck )
	for( int isvab = 0; isvab < node->nsvab; isvab++ )
	{
		dtfm = dtf / node->mass[ isvab / domain->ndime ];
		velol[ isvab ] = cdrag * node->veloi[ isvab ]; //--- drag force
		output[ isvab ] = v[ isvab ] + 0.5 * dtfm * ( - node->fintl[ isvab ] - velol[ isvab ] );
//		if( ( fvl->fix_vel || ff->fix_force ) && node->fixID[ isvab ] == 1 ) //--- set veloc
//			output[ isvab ] = node->presVeloc[ isvab ];
//		if( ff->fix_force && node->rigidID[ isvab ] == 1 ) // --- nodes with external force
//			output[ isvab ] = v[ isvab ] + 0.5 * dtf * ( - effForce[ isvab % domain->ndime ] + node->extForce[ isvab ] ) / node->effMass;
		if( dtv != 0.0 )
			node->dldis[ isvab ] += dtv * output[ isvab ];
	}
}
//------------------------------------------------------
void Fem::SetVeloc() 
{                       
        double y = 0.0;        
        unsigned int isvab = 0, count = 0;
        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
        {
                for( int idime = 0; idime < domain->ndime; idime++ )
                        node->fixID[ isvab ] = 0;
                        isvab++;
        }

        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
        {
                isvab = ipoin * domain->ndime;
                y = node->coord[ ipoin ][ 1 ];
                if( fabs( y - domain->yhi ) < 1.0e-08 )
                {
                        node->fixID[ isvab ] = 1;
                        node->fixID[ isvab + 1 ] = 1;
                        node->presVeloc[ isvab ] = fd->shearRATE * ( domain->yhi - domain->ylo ); //--- unit length
                        node->presVeloc[ isvab + 1 ] = 0.0;
                }
                if( fabs( y - domain->ylo ) < 1.0e-08 )
                {
                        isvab = ipoin * domain->ndime;
                        node->fixID[ isvab ] = 1;
                        node->fixID[ isvab + 1 ] = 1;
                        node->presVeloc[ isvab ] = 0.0;
                        node->presVeloc[ isvab + 1 ] = 0.0;
                }
        }
}
//------------------------------------------------------
void Fem::SetForce()
{
        double y = 0.0;
        unsigned int isvab = 0, count = 0;
	for( int isvab = 0; isvab < node->nsvab; isvab++ )
		node->extForce[ isvab ] = 0.0; //--- initialize external force

	//--- initialize fixID
        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
                for( int idime = 0; idime < domain->ndime; idime++ )
                        node->fixID[ isvab ] = 0;
                        node->rigidID[ isvab ] = 0;
                        isvab++;
        }

	// --- assign fixID and extForce
        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
		isvab = ipoin * domain->ndime;
                y = node->coord[ ipoin ][ 1 ];
                if( fabs( y - domain->yhi ) < 1.0e-08 )
		{
			node->effMass += node->mass[ ipoin ]; //--- effective mass
                        node->rigidID[ isvab ] = 1; // --- nodes with fext
			node->extForce[ isvab ] = node->fx;
                        node->fixID[ isvab + 1 ] = 1; //--- constrained in y
                        node->presVeloc[ isvab + 1 ] = 0.0; // --- vy
//                        node->rigidID[ isvab + 1 ] = 1; // --- nodes with fext
//			node->extForce[ isvab + 1 ] = node->fy;
                }
                if( fabs( y - domain->ylo ) < 1.0e-08 )
		{
                        node->fixID[ isvab ] = 1;
                        node->fixID[ isvab + 1 ] = 1;
                        node->presVeloc[ isvab ] = 0.0;
                        node->presVeloc[ isvab + 1 ] = 0.0;
                }
        }
}
//------------------------------------------------------
bool Fem::Iterate( unsigned int &counter )
{

	// --- output data 
/*	for( int iDump = 0; iDump < dp->nDump; iDump++ )
	{
		if( ( counter + 1 ) % dp->n[ iDump ] == 0 )
			dp->Process( counter + 1, iDump );
	}
*/
//---
	double strsc;
	unsigned int nFreq = ( unsigned int ) ( fv->taud / ts->dt );
        //--- residual force
	fmax = fabs( *max_element( node->absfl, node->absfl + node->nsvab, absCompare ) );
//	fprintf( logFile, "%d\t%e\t%e\t%d\n", counter,  0.5*(ths->cSxyObj->GetSxx() - ths->cSxyObj->GetSyy()) , fmax, nblck ); 
	//---- number of active sites 
	kount = 0;
	#pragma omp barrier //--- because of kount = 0
	#pragma omp for //schedule( dynamic, chunck )
	for( unsigned int kgaus = 0; kgaus < elem->nelem * elem->ngaus; kgaus++ )
		#pragma omp atomic
		kount += elem->actField[ kgaus ];
	if( omp_get_thread_num() == 0 )
	{
		if( kount != 0 || counter == 0 ) //--- avalanche duration
		{
			duration = counter;
			if( cdt->cdt ) //--- s vs t
			{
				strsc = ths->cSxyObj->GetSxy(); //--- current stress
        			if( counter == 0 )
                		{
					cdt->sigma0 = strsc; //--- t = 0
					cdt->sigmp = strsc; //--- initial stress
				}		
				else //--- only at counter = 1, 2, ...
        			{
					cdt->svelo = ( strsc - cdt->sigmp ) / ts->dt; //--- stress rate (velocity)
					cdt->sigmp = strsc; //--- stress at t_n-1
        				if( duration % cdt->nfrq == 0 && duration < cdt->tmax ) //--- only every nfrq steps
					{
						cdt->velos[ cdt->count ] = cdt->svelo; //--- time derivative of stress: ds / dt
	                			cdt->dsigma[ cdt->count ] = cdt->sigma0 - strsc;
						cdt->epdot[ cdt->count ] = ths->cSxyObj->GetGammaDot(); //--- plastic strain rate
	                			cdt->time[ cdt->count ] = duration;
//						cdt->fx[ cdt->count ] = node->fintl[ 0 ]; //--- net force
//						cdt->fy[ cdt->count ] = node->fintl[ 1 ];
						cdt->kin[ cdt->count ] = ths->cKEobj->GetKE(); //--- kin energy
//	    	    				for( unsigned int kgaus = 0; kgaus < elem->nelem * elem->ngaus; kgaus++ )
//	                                		cdt->ep[ kgaus ] += ( elem->stranRate[ kgaus ][ 2 ] * ts->dt ) * ( elem->activity[ kgaus ] ); //--- plastic strain
//	                			fprintf( logFile, "count=%d\t%d\t%d\n", cdt->count, counter, duration );
                				( cdt->count )++;
					}
        			}
			}
		}
		nblck += kount; //--- accumulated activity
		if( counter % nFreq == 0 ) //--- stronger dissipation 
		{
//			fprintf( logFile, "nblock = %d\t%d\n", nblck0, nblck );
			nblck0 = nblck; //--- store in nblck0
		}
	}
        counter++; //--- counter is local
        if( counter == 1 ) //--- run at least once!
                return 1;
	if( ( counter == 2 ) && ( ! qs->min ) ) //--- no minimization
	{
                counter = 0;
		return 0;
	}
/*		if( omp_get_thread_num() == 0 )
		{
			if( counter > 2 )
				fprintf( logFile, "%d\t%e\t%e\t%d\n", counter,  0.5*(ths->cSxyObj->GetSxx() - ths->cSxyObj->GetSyy()) , fmax, nblck );
		}
*/        //--- meet the min. criterion 
	assert( qs->min );
        if( ( ( fmax < qs->ftol ) && counter >= qs->nmin && kount == 0 ) || counter > qs->nmax ) //|| //--- qs->nmin
//        if( ( ( fmax < qs->ftol ) && counter >= qs->nmin ) || counter > qs->nmax ) //--- kount coulb be non-zero
//                ( counter == 2 && fmax < 1.0e+03 * qs->ftol && kount == 0 ) ) //--- avoid further min. in the elastic regime (underdamped) 
        {
//		fprintf( logFile, "hello!\n" );
/*		#pragma omp master
		{ 
		if( fv->switched ) // --- initial dissipation 
		{
			fv->cdrag0 = cdrag0;
			fv->switched = 0;
		}
		}
		#pragma omp barrier
*/
		nblck = 0;
		nblck0 = 0;
                counter = 0;
                return 0; //--- success
        }
        else if( qs->lnsol )//--- failure then linear solver
	{
//		#pragma omp master
//		{ 
		if( ( counter - 1 ) % nFreq == nFreq - 1 && //--- time-lag 
		    nblck == nblck0 && //--- no change in accumulated activity
		    yr->nlaps != 4 ) //--- turn off linear solver for the damge model
		{
//			double xmax, xabs[ node->nsvab ];
//			VectorXd b( node->nsvab ), x( node->nsvab ); //--- dynamic
			fv->switched = 1; //--- no plasticity
			#pragma omp for //schedule( dynamic, chunck )
			for( unsigned int kgaus = 0; kgaus < elem->nelem * elem->ngaus; kgaus++ ) //--- set zero strain rates
			{
				for( unsigned int istre = 0; istre < elem->nstre; istre++ )
					elem->stranRate[ kgaus ][ istre ] = 0.0;
			}
			//--- set zero velocities
			#pragma omp for //schedule( dynamic, chunck )
			for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ )
			{
				node->veloi[ isvab ] = 0.0;
				node->velot[ isvab ] = 0.0;
			}
			while( fmax >= qs->ftol ) //--- starts linear solver
			{
				#pragma omp for //schedule( dynamic, chunck )
				for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ )
					(*b)( isvab ) = -node->fintl[ isvab ]; //--- rhs
				#pragma omp master
				{ 
				(*x) = cg->solve( (*b) ); //--- Ax=b
//				std::cout << "#iterations:     " << cg->iterations() << std::endl;
//				std::cout << "estimated error: " << cg->error()      << std::endl;
//			for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ )
//				xabs[ isvab ] = fabs( x( isvab ) );
//			xmax = fabs( *max_element( xabs, xabs + node->nsvab, absCompare ) );
				}
				#pragma omp barrier
			//--- update node->dldis
				#pragma omp for //schedule( dynamic, chunck )
				for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ )
					node->dldis[ isvab ] = (*x)( isvab );
				force->ComputeForce( counter - 1, 0 ); //--- update fintl
				Update();
				fmax = fabs( *max_element( node->absfl, node->absfl + node->nsvab, absCompare ) );
//				fprintf( logFile, "xmax=%8.7e\n", xmax );
//				fprintf( logFile, "%d\t%e\n", counter, fmax );
        			counter++; //--- counter is local
			}
			//--- assert s < sy!
			unsigned int kgaus = 0, itype;
			double yfunc, uniax, FRICT, pstre, steff, *stres, preys;
			#pragma omp for //schedule( dynamic, chunck )
			for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
			{ 
				itype = elem->ltype[ ielem ];
				uniax = elem->uniax[ itype ];
				FRICT = elem->frict[ itype ];
				for( int igaus = 0; igaus < elem->ngaus; igaus++ )
				{
					kgaus = ielem * elem->ngaus + igaus;
					assert( ! elem->actField[ kgaus ] );
					stres = elem->strst[ kgaus ];
					pstre = force->Smean( stres ); //--- pressure
					steff = force->Invart( stres ); //--- max shear stress: trial stress	
					preys = uniax;
					if( yc->ncrit == 1 or yc->ncrit == 2 ) //--- mohr-coulomb
					{
						preys = - pstre * sin( FRICT ) + uniax * cos( FRICT );
						if( pstre > uniax / tan( FRICT ) ) //--- bi-linear yield function
							preys = 1.0e-08; //--- a residual ???
					}
					yfunc = steff - preys; //--- yield criterion von-mises
//					if( yfunc >= 0.0 )
//						fprintf( logFile, "p=%8.7e\tsy=%8.7e\ts=%8.7e\tyfunc=%8.7e\tfrict=%8.7e\tpreys=%8.7e\n", - pstre, uniax, steff, yfunc, FRICT, preys );
					assert( yfunc < 0.0 ); 
				}
			}
			//--- initialize	
			fv->switched = 0;
                	counter = 0;
			nblck = 0;
			nblck0 = 0;
			return 0;
//			fv->switched = 1;
//			cdrag0 = fv->cdrag0;
//			fv->cdrag0 = fv->critDrag;
		}
                return 1;
	}
}
//------------------------------------------------------ 

