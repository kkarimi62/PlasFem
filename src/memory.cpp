#include "memory.h"

Memory::Memory(){}
Memory::~Memory(){}

void Memory::doubleMat( double **& matrix, const unsigned int m, const unsigned int n )
{
	matrix = new double*[ m ];
	for( int i = 0; i < m; i++ )
		matrix[ i ] = new double[ n ];
	// initialize
	for( unsigned int i = 0; i < m; i++ ) {
		for( unsigned int j = 0; j < n; j++ )
			matrix[ i ][ j ] = 0.0;
	}
}

void Memory::doubleVec( double *& vector, const unsigned int m )
{
	vector = new double[ m ];
	for( int i = 0; i < m; i++ )
		vector[ i ] = 0.0;
	
}

void Memory::Delete2ndMat( double ** matrix, const unsigned int m )
{
	for( int i = 0; i < m; i++ )
		delete [] matrix[ i ];
	delete matrix;
}
void Memory::Delete2ndMat( unsigned int ** matrix, const unsigned int m )
{
	for( int i = 0; i < m; i++ )
		delete [] matrix[ i ];
	delete matrix;
}
void Memory::Delete2ndMat( int ** matrix, const unsigned int m )
{
	for( int i = 0; i < m; i++ )
		delete [] matrix[ i ];
	delete matrix;
}

void Memory::Delete3rdMat( double *** matrix, const unsigned int m, const unsigned int n )
{
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			delete [] matrix[ i ][ j ];
	}

	for( int i = 0; i < m; i++ )
		delete [] matrix[ i ];

	delete matrix;
}

void Memory::Delete4thMat( double **** matrix, const unsigned int m, const unsigned int n, const unsigned int o )
{
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
		{
			for( int k = 0; k < n; k++ )
				delete [] matrix[ i ][ j ][ k ];
		}
	}

	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			delete [] matrix[ i ][ j ];
	}

	for( int i = 0; i < m; i++ )
		delete [] matrix[ i ];

	delete matrix;
}
