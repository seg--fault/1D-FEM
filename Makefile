all: main.o wavefunction.o
	g++ -o FEM main.o wavefunction.o -larmadillo

main.o: main.cpp wavefunction.h
	g++ -c main.cpp

wavefunction.o: wavefunction.cpp wavefunction.h
	g++ -c wavefunction.cpp
	

clean:
	rm wavefunction.o
