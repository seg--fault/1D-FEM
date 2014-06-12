#include <armadillo>
#include <iostream>
#include <math.h>
#include <vector>

#include "gnuplot-iostream/gnuplot-iostream.h"
#include "wavefunction.h"


using namespace arma;
using namespace std;

// Hydrogen atom

double hamiltonian(double x, double (element::*wf)(double x, unsigned int m, unsigned int der), unsigned int m, unsigned int n, element* ptr, void* args)
{
	return x*x*((ptr->*wf)(x, m, 1)*(ptr->*wf)(x, n, 1) + (ptr->*wf)(x, m, 0)*-2/x*(ptr->*wf)(x, n, 0));
}

double energy(double x, double (element::*wf)(double x, unsigned int m, unsigned int der), unsigned int m, unsigned int n, element* ptr, void* args)
{
	return x*x*(ptr->*wf)(x, m, 0)*(ptr->*wf)(x, n, 0);
}

/*  Simpler Harmonic Oscillator 
double hamiltonian(double x, double (element::*wf)(double x, unsigned int m, unsigned int der), unsigned int m, unsigned int n, element* ptr, void* args)
{
	return ((ptr->*wf)(x, m, 1)*(ptr->*wf)(x, n, 1) + (ptr->*wf)(x, m, 0)*x*x*(ptr->*wf)(x, n, 0));
}

double energy(double x, double (element::*wf)(double x, unsigned int m, unsigned int der), unsigned int m, unsigned int n, element* ptr, void* args)
{
	return (ptr->*wf)(x, m, 0)*(ptr->*wf)(x, n, 0);
}*/

int main()
{
	//Create vectors for the endpoints and nodes
	vector<double> test_endpoints;
	vector<double> test_nodes;

	//Boundaries
	double upper_bound, lower_bound;
	unsigned int DOF, n_elements, n_nodes;

	cout << setprecision(15);
	
	//See what bounds to use
	cout << "Upper Bound: ";
	cin  >> upper_bound;
	cout << "Lower Bound: ";
	cin  >> lower_bound;
	cout << "Number of Elements: ";
	cin  >> n_elements;
	cout << "Number of Nodes: ";
	cin  >> n_nodes;
	cout << "Degrees of Freedom per Nodes: ";
	cin  >> DOF;
	
	//Add some endpoints
	for(double x=lower_bound; x<=upper_bound; x+=(upper_bound-lower_bound)/(double)n_elements)
	{
		test_endpoints.push_back(x);
	}

	//And add some nodes...  No more than 24 with armadillo?
	for(double x=lower_bound; x<=upper_bound; x+=(upper_bound-lower_bound)/(double)n_elements)
	{
		test_nodes.push_back(x);
	}


	//Finally, create the wavefunction
	wavefunction test_wavefunction(test_endpoints, test_nodes, DOF, &hamiltonian, &energy, NULL, NULL);

	//Space for sample points
	vector< pair< double, double > > samples;

	//Make the gnuplot iostream object
	Gnuplot gp;

	//Set up some options on gnuplot
	gp << "set xrange [" << lower_bound << ":" << upper_bound << "]\n";
	gp << "set yrange [-1:1.5]\n";

	//Holds queried quantum number
	unsigned int QN = 0;

	while(QN < test_wavefunction.n_eigenvalues())
	{
		//Query
		cout << "Quantum Number = ";
		cin >> QN;

		//Return the eigenvalue
		cout << "Eigenvalue[" << QN << "] = " << test_wavefunction.eigenvalue(QN) << endl;

		//Make sure there aren't any samples here
		samples.clear();

		for(double x = lower_bound; x < upper_bound; x+=0.01)//test_wavefunction.upper_bound(); x += 0.1)
		{
			//Make some samples
			samples.push_back(make_pair(x, test_wavefunction(x, QN)));
		}

		gp << "plot" << gp.file1d(samples) << "with lines title ''" << endl;
	}
			
	//Nothing went too badly
	return 0;
}
