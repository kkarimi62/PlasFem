#--- Makefile for compiling cpp files in PlasFem
#--- openMP with its default compiler


#--- compiler/linker settings
#--- specify flags and libraries needed for your compiler

CC=c++

INC_DIR = /home/kamran.karimi1/Project/opt/eigen-eigen-323c052e1731
LIB_DIR = .
BOOST_INC = /home/kamran.karimi1/Project/opt/boost_1_59_0

CFLAGS= -fopenmp -I$(INC_DIR) -I$(BOOST_INC) -O2 -std=c++11 -static-libstdc++ -w

#--- define the executable file
EXE = fem_run

#--- List of all source files
SRC = $(wildcard *.cpp)

#--- Object files defined from source files
OBJ = $(patsubst %.cpp, %.o, $(SRC))

#--- List of dependencies defined from list of object files
#DEPS := $(OBJ:.o=.d)
#-include $(DEPS)

#$(info VAR="$(OBJ)")

#--- Compilation rules

#%.o: %.cpp $(DEPS)
%.o: %.cpp 
	@echo Compiling $<...
	$(CC) $(CFLAGS) -c $< -o $@ 

#--- Link target
$(EXE): $(OBJ)
	$(CC) $^ -o $@  $(CFLAGS)

clean:
	rm -f *.o $(EXE)  

#---  
