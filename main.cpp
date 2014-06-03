#include <armadillo>
#include <iostream>
#include <math.h>
#include <vector>

#include "wavefunction.h"


using namespace arma;
using namespace std;

double test_potential(double x, void* args)
{
	return pow(x,2);
}


int main()
{
	//Create vectors for the endpoints and nodes
	vector<double> test_endpoints;
	vector<double> test_nodes;

	//Add some endpoints
	for(double x=-6;x<=6;x+=1)
	{
		test_endpoints.push_back(x);
	}

	//And add some nodes...  No more than 24 with armadillo?
	for(double x=-6;x<=6;x+=0.25)
	{
		test_nodes.push_back(x);
	}


	//Finally, create the wavefunction
	wavefunction test_wavefunction(test_endpoints, test_nodes, 2, &test_potential, NULL);

	//Nothing went too badly
	return 0;
}
