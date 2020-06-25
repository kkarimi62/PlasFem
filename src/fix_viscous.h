#ifndef FIX_VISCOUS_H
#define FIX_VISCOUS_H

class FixViscous
{
	public:
		FixViscous();
		~FixViscous();
		double cdrag0, cdrag1, taud, critDrag;
		bool switched;
};

#endif
