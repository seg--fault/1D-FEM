#Set up useful variables
LIBRARIES= -larmadillo -lboost_iostreams -lboost_system -lboost_filesystem -lgsl
CC_OPTS= -O3 -march=native -DARMA_NO_DEBUG

#Define primary targets to build
all: element.o main.o matrix_utilities.o wavefunction.o
	g++ $(CC_OPTS) element.o main.o matrix_utilities.o wavefunction.o $(LIBRARIES) -o FEM

#Define object files to build
element.o: element.cpp  element.h
	g++ -c $(CC_OPTS) element.cpp

main.o: main.cpp element.h wavefunction.h
	g++ -c $(CC_OPTS) main.cpp

matrix_utilities.o: matrix_utilities.cpp
	g++ -c $(CC_OPTS) matrix_utilities.cpp

wavefunction.o: element.h  wavefunction.cpp wavefunction.h
	g++ -c $(CC_OPTS) wavefunction.cpp

#Give an option to clean up intermediate files
clean:
	rm element.o main.o wavefunction.o
