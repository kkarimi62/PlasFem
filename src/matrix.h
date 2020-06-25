#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
	public:	
		Matrix();
		~Matrix();

		void MulMatVec(  double **,  unsigned int,  unsigned int, 
				 double *,   unsigned int, double * );
		void MulMatVec(  double *,  unsigned int,  unsigned int, 
				 double *,   unsigned int, double * );
		void MulMatMat(  double **,  unsigned int ,  unsigned int, 
				 double **,  unsigned int,  unsigned int, double ** );
		void MulMatMat(  double *,  unsigned int ,  unsigned int, 
				 double **,  unsigned int,  unsigned int, double ** );
		void Transpose(  double **,  unsigned int,  unsigned int, double ** );
	private:
};

#endif
