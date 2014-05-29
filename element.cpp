#include<armadillo>
#include<iostream>
#include<math.h>

#include "element.h"

using namespace arma;
using namespace std;

double factorial(int n)
{
	if(n <= 0)  return 1 ; 		//safeguard 0 and negative values
	double res = n ;
	while(--n>1) res *= n ;		//start multiplying
	return res ;
}

element::element()
{
	polynomial_initialized=false;	//the polynomial matrix starts out uninitialized
}

bool element::initialize_polynomial(unsigned int DOF, vector<double> nodes)
{
	if(polynomial_initialized)
	{
		cout<<"Warning. The polynomial has already been initialized.  Reinitializing anyway";
	}
	else
	{
		for(int i=0; i<nodes.size(); i++)
		{
			if(nodes.at(i)<-1 || nodes.at(i)>1)
			{
				std::cout<<"Error initializing polynomial.  Bad node value.";
				return false;
			}
		}

		Mat<double> system(DOF*nodes.size(), DOF*nodes.size());;	//set up a matrix to fill with our values

		for(int i=0; i<nodes.size(); i++)
		{
			for(int j=0; j<DOF; j++)
			{
				for(int k=0; k<DOF*nodes.size(); k++)
				{
					if(k<j)
						system(i*DOF+j, k) = 0;		//derivative of a constant
					if(k==j)
						system(i*DOF+j, k) = 1;		//is the original constant
					if(k>j)
						system(i*DOF+j, k) = pow(nodes.at(i), (double)(k-j))*factorial(k)/factorial(k-j);	//add in the correct values
				}
			}
		}

		polynomial_coefficients = system.i();				//invert the matrix to find the polynomial coefficients
	}
	
	polynomial_initialized = true;		//the polynomial has been initialized
	return true;				//we succeeded
}

Mat<double> element::get_polynomial()
{
	if(!polynomial_initialized)
	{
		cout<<"Error, reading uninitialized polynomial coefficients.";
	}

	return polynomial_coefficients;
}

void element::set_polynomial(Mat<double> coefficients)
{
	if(polynomial_initialized)
	{
		cout<<"Warning. Overwriting polynomial coefficients!";
	}

	polynomial_coefficients = coefficients;
}

Col<double> element::get_coefficients()
{
	return element_coefficients;
}

void element::set_coefficients(Col<double> coefficients)
{
	element_coefficients = coefficients;
}

double element::glbl_to_lcl(double x)
{
	return (2*x-(top+bottom))/(top-bottom);
}

double element::lcl_to_glbl(double x)
{
	return (top-bottom)*x/2+(top+bottom)/2;
}

double element::at(double x)
{
	return at_lcl(glbl_to_lcl(x));
}

double element::at(double x, unsigned int n)
{
	return at_lcl(glbl_to_lcl(x), n);
}

double element::at_lcl(double x)
{
	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients;
	return ret(0,0);
}

double element::at_lcl(double x, unsigned int n)
{
	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients;
	return ret(0,0);
}

Row<double> element::variable_row(double x, unsigned int size)
{
	Row<double> ret(size);

	ret(0) = 1;

	if(size>0)
	{
		for(int i=1; i<size; i++)
		{
			ret(i) = pow(x, (double)(i));
		}
	}
	else
	{
		cout<<"You're an idiot!";
	}

	return ret;
}
Row<double> element::variable_row(double x, unsigned int size, unsigned int n)
{
	Row<double> ret(size);

	if(size>0)
	{
		for(int i=0; i<size; i++)
		{
			if(i<n)
				ret(i) = 0;		//derivative of a constant
			if(i==n)
				ret(i) = 1;		//is the original constant
			if(i>n)
				ret(i) = pow(x, (double)(i-n))*factorial(i)/factorial(i-n);	//add in the correct values
		}
	}
	else
	{
		cout<<"DBAD";
	}

	return ret;
}
