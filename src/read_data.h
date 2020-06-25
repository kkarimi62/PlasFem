

#ifndef READ_DATA_H
#define READ_DATA_H
#include "fem.h"
#include "memory.h"

class ReadData
{
	public:
		ReadData( char *, Fem * );
		~ReadData();
		void process();

	private:
		Fem *femPtr;
		FILE *infile; // file name
		char *line;  // command string
		char** args; // string for a parsed command
		void parse();
		Memory *memory;
		const unsigned int MAXSTRLEN;// = 100; // max string len
		const unsigned int MAXARG;// = 10; // max entry per line
};

#endif
