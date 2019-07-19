#ifndef TIME_STEP_H
#define TIME_STEP_H

class TimeStep
{
	public:
		TimeStep();
		~TimeStep();
		double dt;
		unsigned int itime0; 
};

#endif
