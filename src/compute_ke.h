#ifndef COMPUTE_KE_H
#define COMPUTE_KE_H


// --- calculates the incremental potential energy
class ComputeKE
{
	public:
		ComputeKE();
		void Init(  double *,  double *,  unsigned int,  int );
		double GetKE();

	private:
		 double *xmass, *veloc;
		unsigned int n;
		int ndime;
		double ke;
		unsigned int ipoin;
};

#endif
