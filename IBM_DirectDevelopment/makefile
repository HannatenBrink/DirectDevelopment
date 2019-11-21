C=g++
CFLAGS= -std=c++11 -Wall -O3 

all: RunIBM

RunIBM: main.o Individual.o Resource.o
	$(C) -L/usr/local/include main.o  Individual.o Resource.o  -o RunIBM

main.o: main.cpp main.h Individual.h Resource.h
	$(C) $(CFLAGS) -c main.cpp -o main.o




Resource.o: Resource.cpp Resource.h
	$(C) $(CFLAGS) -c Resource.cpp -o Resource.o




Individual.o: Individual.cpp Individual.h
	$(C) $(CFLAGS) -c Individual.cpp -o Individual.o

clean:
	rm *o RunIBM
