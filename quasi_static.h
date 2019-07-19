#ifndef QUASI_STATIC_H
#define QUASI_STATIC_H

class QuasiStatic
{
	public:
		QuasiStatic();
		~QuasiStatic();
		bool min;
		double ftol;
		unsigned int nmin, nmax;
		bool lnsol;
};

#endif
