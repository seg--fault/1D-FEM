//For matrix operations
#include <armadillo>

//Function prototypes
#include "matrix_utilities.h"

using namespace arma;

//used in calculating variable row derivatives
double factorial(int n)
{
	if(n <= 0)  return 1 ; 		//safeguard 0 and negative values
	double res = n ;
	while(--n>1) res *= n ;		//start multiplying
	return res ;
}

Row<double> variable_row(double x, unsigned int size)
{
	Row<double> ret(size);


#ifdef ERROR_CHECKING
	if(size<=0)
	{
		cout << "Error. Can't make variable row of size <= 0.";
		return ret;
	}
#endif

	double acc = 1;

	for(int i=0; i<size; i++)
	{
		ret(i) = acc;
		acc = acc*x;
	}

	return ret;
}

Row<double> variable_row(double x, unsigned int size, unsigned int n)
{
	Row<double> ret(size);

#ifdef ERROR_CHECKING
	if(size<=0)
	{
		cout << "Error.  Cannot create variable row of size <= 0";
		return ret;
	}
#endif

	double acc = 1;

	for(int i=0; i<size; i++)
	{
		if(i<n)
			ret(i) = 0;		//derivative of a constant
		else
		{
			ret(i) = acc*factorial(i)/factorial(i-n);	//add in the correct values
			acc *= x;
		}
	}

	return ret;
}
