#ifndef RAND_NUMBERS_H
#define RAND_NUMBERS_H
#include <random>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random.hpp>

class RandNumbers
{
	public:
	RandNumbers( /*const double,*/ unsigned int );
	~RandNumbers();

	void Normal( double * /*output*/, const double /*mean value*/, const double /*width*/ );
	void Expon( double * /*output*/, const double /*mean value*/ );

	private:
		boost::mt19937 rng; 
//  		boost::variate_generator< boost::mt19937&, boost::exponential_distribution<> > rndm;
		std::default_random_engine *generator;
};

#endif
