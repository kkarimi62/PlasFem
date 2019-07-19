#ifndef MEMORY_H
#define MEMORY_H


class Memory
{
	public:
		Memory();
		~Memory();

		void doubleMat( double **&, const unsigned int, const unsigned int );
		void doubleVec( double *&, const unsigned int );
		void Delete2ndMat( double **, const unsigned int );
		void Delete2ndMat( int **, const unsigned int );
		void Delete2ndMat( unsigned int **, const unsigned int );
		void Delete3rdMat( double ***, const unsigned int, const unsigned int );
		void Delete4thMat( double ****, const unsigned int, const unsigned int, const unsigned int );
};

#endif
