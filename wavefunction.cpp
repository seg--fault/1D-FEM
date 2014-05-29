#include<armadillo>
#include<math.h>

#include "wavefunction.h"

using namespace arma;

double factorial( int n )
{
	if( n <= 0 )  return 1 ; // safeguard 0 and -ve
	double res = n ;
	while( --n>1 ) res *= n ;
	return res ;
}

wavefunction::wavefunction()
{
	polynomial_initialized=false;
}

bool wavefunction::initialize_polynomial(unsigned int DOF, vector<double> nodes)
{
	for(int i=0;i<nodes.size();i++)
	{
		if(nodes.at(i)<-1 || nodes.at(i)>1)
		{
			std::cout<<"Error initializing polynomial.  Bad node value.";
			return false;
		}
	}

	Mat<double> system(DOF*nodes.size(),DOF*nodes.size());;

	for(int i=0;i<nodes.size();i++)
	{
		for(int j=0;j<DOF;j++)
		{
			for(int k=0;k<DOF*nodes.size();k++)
			{
				if(k<j)
					system(i*DOF+j,k) = 0;
				if(k==j)
					system(i*DOF+j,k) = 1;
				if(k>j)
					system(i*DOF+j,k) = pow(nodes.at(i),(double)(k-j))*factorial(k)/factorial(k-j);
			}
		}
	}
	
	polynomial_coefficients = system.i();
}


