#include "matrix.h"
#include <cassert>
//#include <iostream>

//using std::cout;

Matrix::Matrix(){}
Matrix::~Matrix(){}

void Matrix::Transpose(  double **A,  unsigned int m,  unsigned int n, double **B )
{
	for( int i = 0; i < n; i++ ){ 
		for( int j = 0; j < m; j++ ){ 
		B[ i ][ j ] = A[ j ][ i ];
		}
	}
}

void Matrix::MulMatVec(  double **A,  unsigned int m,  unsigned int n,  double *u,  unsigned int o, double *v )
{
	assert( n == o );
	double sum = 0.0;
	for( int i = 0; i < m; i++ ){
		sum = 0.0;
		for( int j = 0; j < n; j++ ){
			sum += A[ i ][ j ] * u[ j ];
		}
		v[ i ] = sum;
	}
}

void Matrix::MulMatVec(  double *A,  unsigned int m,  unsigned int n,  double *u,  unsigned int o, double *v )
{
	assert( n == o );
	double sum = 0.0;
	for( int i = 0; i < m; i++ ){
		sum = 0.0;
		for( int j = 0; j < n; j++ ){
			sum += A[ i * n + j ] * u[ j ];
		}
		v[ i ] = sum;
	}
}
//--------------------------------------------------------------
void Matrix::MulMatMat(  double *A,  unsigned int m,  unsigned int n, 
			 double **B,  unsigned int o,  unsigned int p, double **C )
{
	assert( n == o );
	double sum = 0.0;
	for( int i = 0; i < m; i++ ){
		for( int k = 0; k < p; k++ ){
			sum = 0.0;
			for( int j = 0; j < n; j++ )
				sum += A[ i * n + j ] * B[ j ][ k ];
			C[ i ][ k ] = sum;
		}
	}
}
//--------------------------------------------------------------
void Matrix::MulMatMat(  double **A,  unsigned int m,  unsigned int n, 
			 double **B,  unsigned int o,  unsigned int p, double **C )
{
	assert( n == o );
	double sum = 0.0;
	for( int i = 0; i < m; i++ ){
		for( int k = 0; k < p; k++ ){
			sum = 0.0;
			for( int j = 0; j < n; j++ )
				sum += A[ i ][ j ] * B[ j ][ k ];
			C[ i ][ k ] = sum;
		}
	}
}

/*
int main()
{
	Matrix *obj;

	double **A = new double*[ 2 ];
	A[ 0 ] = new double[ 2 ];
	A[ 1 ] = new double[ 2 ];
	A[ 0 ][ 0 ] = 0.0;
	A[ 0 ][ 1 ] = 1.0;
	A[ 1 ][ 0 ] = 2.0;
	A[ 1 ][ 1 ] = 3.0;

	double *u = new double[ 2 ];
	u[ 0 ] = 1.0;
	u[ 1 ] = 2.0;

	double *v = new double[ 2 ];

	obj->MulMatVec( A, 2, 2, u, 2, v ); 
	cout << v[ 0 ] << '\n' << v[ 1 ] << '\n';

	delete [] A;
	delete [] u;
	delete [] v;
}
*/
