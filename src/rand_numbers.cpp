//#include "Python.h"
#include "rand_numbers.h"
//#include <random>
//#include <stdio.h>

RandNumbers::RandNumbers( /*const double mean,*/ unsigned int seed ):
rng( seed )/*, rndm( rng, boost::exponential_distribution<double>( mean ) )*/, generator( new std::default_random_engine )
{
//	printf( "constructor!\n" );
//	std::random_device rd;
}

RandNumbers::~RandNumbers()
{
//	delete generator;
}

void RandNumbers::Normal( double *output, const double mean, const double width )
{
	std::normal_distribution< double > distribution( mean, width );
	*output = distribution( *generator );
}

void RandNumbers::Expon( double *output, const double mean )
{
//	std::exponential_distribution< double > distribution( mean );
	boost::variate_generator< boost::mt19937&, boost::exponential_distribution<> > rndm( rng, boost::exponential_distribution<double>( mean ) );
	*output = rndm();
}

/*
int main()
{
	double *x = new double;
	srand( time( NULL ) );
	RandNumbers *randObj = new RandNumbers ( 1.0, rand() );
	FILE *someFile = fopen( "out.txt", "w" );
	for( int i = 0; i < 10000 ; i++ ){
		randObj->Expon( x );
		fprintf( someFile, "%15.10g\n", *x );
	}
	fclose( someFile );
}
*/
