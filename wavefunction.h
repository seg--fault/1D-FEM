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

	bool create_elements(vector<double> endpoints, vector<double> nodes, unsigned int DOF);	//initialized elements and polynomials
	bool calculate_element_coefficients();

	double at(double x);

	private:
	vector<element>		elements;
	vector<double>		endpoints;	//endpoints of elements in increasing order
};

#endif
