#ifndef MAT_PROP_H 
#define MAT_PROP_H
#include "fem.h"

class MatProp
{
	public:
		MatProp( char *, Fem * );
		~MatProp();
		void process();

	private:
		Fem *femPtr;
		FILE *infile; // file name
		char *line;  // command string
		char** args; // string for a parsed command
		void parse();
		const unsigned int MAXSTRLEN;// = 100; // max string len
		const unsigned int MAXARG;// = 10; // max entry per line
};

#endif
