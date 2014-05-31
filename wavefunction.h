#ifndef __WAVEFUNCTION__
#define __WAVEFUNCTION__

#include<armadillo>
#include<vector>

#include"element.h"
#include"potential.h"

using namespace arma;
using namespace std;

class wavefunction
{
	public:
	wavefunction();
	wavefunction(vector<double> endpoints, vector<double> nodes, unsigned int DOF, abstract_potential potential);

	//bool change_elements(vector<double> endpoints, vector<double> nodes, unsigned int DOF);
	bool create_elements(vector<double> endpoints, vector<double> nodes, unsigned int DOF);	//initialized elements and polynomials
	bool calculate_element_coefficients();

	//bool set_potential(abstract_potential potential);

	double at(double x);

	private:
	vector<element>		elements;
	vector<double>		endpoints;	//endpoints of elements in increasing order
	vector<double>		nodes;
	abstract_potential	potential;
};

#endif
