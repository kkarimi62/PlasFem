# Makefile for compiling cpp files in PlasFem
#
CC=c++

INC_DIR = .
LIB_DIR = .
EXE = hellomake

CFLAGS=-I$(INC_DIR)

DEPS = map.h
OBJ = map.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o $(EXE)  


#---  
