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
	wavefunction(vector<double> endpoints, vector<double> nodes, unsigned int DOF, op hamiltonian, op energy,  void* ham_args, void* en_args);
	~wavefunction();

	//Creates elements with polynomials based on arguments
	bool create_elements(vector<double> endpoints, vector<double> nodes, unsigned int DOF);	
	
	//Solves for eigenfunctions and eigenvalues based on potential
	bool solve(op hamiltonian, op energy, void* ham_args, void* en_args);			

	//Evaluate the nth wavefunction at x
	double operator()(double x, unsigned int quantum_number);
	double eigenvalue(unsigned int quantum_number);	
	unsigned int n_eigenvalues();

	//Finds the limits of the wavefunction
	double upper_bound();
	double lower_bound();

	private:
	//Degrees of Freedom
	unsigned int DOF;

	//Storage for the elments in the wavefunction sorted by lower endpoint
	unsigned int n_elements;
	element* elements;

	//Space for global solutions
	Col<double> eigenvalues;
};

#endif
