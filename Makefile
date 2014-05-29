all: element.o main.o wavefunction.o
	g++ -o FEM element.o main.o wavefunction.o -larmadillo

element.o: element.cpp  element.h
	g++ -c element.cpp

main.o: main.cpp wavefunction.h
	g++ -c main.cpp

wavefunction.o: element.h wavefunction.cpp wavefunction.h
	g++ -c wavefunction.cpp

clean:
	rm wavefunction.o
