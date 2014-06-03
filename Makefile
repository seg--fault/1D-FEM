all: element.o main.o wavefunction.o
	g++ -o FEM element.o main.o wavefunction.o -larmadillo -lboost_iostreams -lboost_system -lboost_filesystem -lgsl

element.o: element.cpp  element.h
	g++ -c element.cpp

main.o: main.cpp element.h
	g++ -c main.cpp 

wavefunction.o: element.h  wavefunction.cpp wavefunction.h
	g++ -c wavefunction.cpp

clean:
	rm element.o FEM main.o wavefunction.o
