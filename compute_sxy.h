#ifndef COMPUTE_SXY_H
#define COMPUTE_SXY_H

// --- calculates the incremental potential energy
class ComputeSxy
{
	public:
		ComputeSxy();
		~ComputeSxy();
		void Init(  double**,  double **, double **,  double *,  unsigned int,  unsigned int );
		double GetSxy();
		double GetSxx();
		double GetSyy();
		double GetP();
		double GetSigma();
		double GetEv();
		double GetGammaDot();

	private:
		 double **s, *area, **eps, **epsDOT;
		unsigned int nelem, ngaus;
		double sxyAvg, sum, totArea;
		unsigned int kgaus;
		
};

#endif
