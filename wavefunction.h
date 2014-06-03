#ifndef __WAVEFUNCTION__
#define __WAVEFUNCTION__

#include<armadillo>
#include<vector>

#include"element.h"

using namespace arma;
using namespace std;

class wavefunction
{
	public:
	wavefunction();
	wavefunction(vector<double> endpoints, vector<double> nodes, unsigned int DOF, double (*potential)(double x, void* args), void* args);

	//Method to create elements with polynomials based on arguments
	bool create_elements(vector<double> endpoints, vector<double> nodes, unsigned int DOF);	
	
	//Method to solve for eigenfunctions and eigenvalues based on potential
	bool solve(double (*potential)(double x, void* args), void* args);			

	//Method to return the first element containing x
	unsigned int containing_element(double x);

	//Evaluate the nth wavefunction at x
	double operator()(double x, unsigned int quantum_number);

	private:
	//Degrees of Freedom
	unsigned int DOF;

	//Storage for the elments in the wavefunction sorted by lower endpoint
	vector<element*> elements;

	//Space for global solutions
	Col<double> eigenvalues;
	Mat<double> eigenvectors;
};

#endif
