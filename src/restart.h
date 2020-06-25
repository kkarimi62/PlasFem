#ifndef RESTART_H 
#define RESTART_H

#include "stdio.h"

class Fem;
class Restart
{
	public:
		Restart( Fem * );
		~Restart();
		unsigned int n;
		void output( unsigned int );
	private:
		Fem *femPtr;
		FILE *someFile;
		unsigned int isvab, kgaus;
		
};

#endif 
