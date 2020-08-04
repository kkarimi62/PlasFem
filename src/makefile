#--- Makefile for compiling cpp files in PlasFem
#--- openMP with its default compiler


#--- compiler/linker settings
#--- specify flags and libraries needed for your compiler

CC=c++

INC_DIR = /home/kkarimi/opt/eigen-eigen-07105f7124f9 
LIB_DIR = .
BOOST_INC = /home/kkarimi/opt/boost_1_59_0 

CFLAGS=-I$(INC_DIR) -I$(BOOST_INC) -fopenmp -O2 -std=c++11 -static-libstdc++ -w


DEPS = *.h
OBJ = *.o
EXE = fem_run

#--- Compilation rules

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

#--- Link target

$(EXE): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)


clean:
	rm -f *.o $(EXE)  


#---  
