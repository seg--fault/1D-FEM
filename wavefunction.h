#ifndef __WAVEFUNCTION__
#define __WAVEFUNCTION__

#include<armadillo>
#include<vector>

using namespace arma;
using namespace std;

class wavefunction
{
	public:
	wavefunction();

	bool initialize_polynomial(unsigned int DOF, vector<double> nodes);

	private:
	bool polynomial_initialized;
	
	Mat<double> polynomial_coefficients;
};

#endif
