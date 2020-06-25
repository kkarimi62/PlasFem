#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>


class Fem;
class Input
{
	public:
		Input( Fem * );
		~Input();
		void file(); // read input file and process commands

	private:
		FILE *infile; // input script
		int narg;
		char *line, *copy, *command;
		char **arg;
		Fem *femObj;

		void parse(); // parse the command line 
		void allocate_memory();	
		void execute_command();
		void read_data(); // input script commands
		void mat_prop(); // material properties
		void time_step(); // time step
		void fix_deform(); 
		void fix_veloc();
		void fix_force();
		void min();
		void msd();
		void fix_viscous(); 
		void run();
		void thermo_style();
		void dump();
		void dump_modify();
		void restart();
		void ComputeDS();
		void ComputeDT();
		void ComputeNP();
		void ComputeX();
		void ComputeXL();
		void ncrit();
		void nweak();
		void nlaps();
		void nrecv();
		void pure_shear();
};


#endif
