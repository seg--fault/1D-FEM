#include<algorithm>
#include<armadillo>
#include<gsl/gsl_integration.h>
#include<iostream>
#include<math.h>

#include "element.h"
#include "matrix_utilities.h"


using namespace arma;
using namespace std;

element::element()
{
	//Set up safe defaults for top and bottom until they are assigned
	top=1;
	bottom=-1;

	//The polynomial and element haven't been found yet
	polynomial_initialized=false;
	element_solved=false;
}

element::element(unsigned int DOF, vector<double> nodes, double arg_bottom, double arg_top)
{
	//Set the top and bottom bounds
	top = arg_top;
	bottom = arg_bottom;

	//The polynomial and element haven't been solved yet
	polynomial_initialized=false;
	element_solved=false;

	//Create polynomial coefficients
	initialize_polynomial(DOF, nodes);
}

//Creates polynomial coefficients based on given nodes in local coordinates and the degrees of freedom
bool element::initialize_polynomial(unsigned int DOF, vector<double> nodes)
{
#ifdef ERROR_CHECKING
	//Are we overwriting the polynomial?
	if(polynomial_initialized)
		cout<<"Warning. The polynomial has already been initialized.  Reinitializing anyway";

	//Is the value for Degrees of Freedom sane?
	if(DOF<1)
	{
		cout<<"Error. Bad value for DOF";
		return false;
	}

	//Do we have enough nodes
	if(nodes.size()<2)
	{
		cout<<"Error. Too few nodes to construct polynomial.";
		return false;
	}

	//are the nodes placed correctly
	for(int i=0; i<nodes.size(); i++)
	{
		if(nodes.at(i)<-1 || nodes.at(i)>1)
		{
			std::cout<<"Error initializing polynomial.  Bad node value.";
			return false;
		}
	}
#endif

	//The matrix that we will invert to find the polynomial coefficients
	Mat<double> system(DOF*nodes.size(), DOF*nodes.size());;

	//Go through the rows and fill them in
	for(int i=0; i<nodes.size(); i++)
	{
		for(int j=0; j<DOF; j++)
		{
			//Fill it in with the polynomial variables at the node
			system.row(i*DOF+j) = variable_row(glbl_to_lcl(nodes.at(i)), DOF*nodes.size(), j);
		}
	}

	//Invert the matrix to solve for the coefficients
	polynomial_coefficients = system.i();

	//The polynomial has now been initialized
	polynomial_initialized = true;	
	return true;		
}

//Hand back polynomial coefficients
Mat<double> element::get_polynomial()
{
#ifdef ERROR_CHECKING
	//Have the polynomials been found yet?
	if(!polynomial_initialized)
	{
		cout<<"Error, reading uninitialized polynomial coefficients.";
	}
#endif

	//Return the polynomial coefficients
	return polynomial_coefficients;
}

//Sets the polynomial coefficients
void element::set_polynomial(Mat<double> coefficients)
{
#ifdef ERROR_CHECKING
	//Are we going to overwrite the polynomial?
	if(polynomial_initialized)
	{
		cout<<"Warning. Overwriting polynomial coefficients!";
	}
#endif

	//Replace the polynomial coefficients and the polynomial has been initialized
	polynomial_coefficients = coefficients;
	polynomial_initialized=true;
}

//Hand back matrix of solution coefficients
Mat<double> element::get_coefficients()
{
#ifdef ERROR_CHECKING
	if(!element_solved)
	{
		cout << "Error. Returning incorrect coefficients.";
	}
#endif

	return element_coefficients;
}

void element::set_coefficients(Mat<double> coefficients)
{
#ifdef ERROR_CHECKING
	if(!polynomial_initialized)
	{
		cout << "Error.  Cannot set element coefficients, polynomial is unititialized.";
		return;
	}

	if(coefficients.n_rows!=polynomial_coefficients.n_rows)
	{
		cout << "Error.  Cannot set element coefficients, wrong number of elements.";
		return;
	}
#endif
	//Resize the matrix
	element_coefficients.set_size(coefficients.n_rows, coefficients.n_cols);

	//Transfer in the elements
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

int element::in(double x)
{
	//return 0 if in element, 1 if below and 2 if above
	if(x < bottom)
		return 1;
	if(x > top)
		return 2;
	return 0;
}

//Value of eigenvector at x in global coordinates
double element::at(double x, unsigned int quantum_number)
{
	return at_lcl(glbl_to_lcl(x), quantum_number);
}

//Value of derivative of eigenvector at x in global coordinates
double element::at(double x, unsigned int quantum_number, unsigned int n)
{
	return at_lcl(glbl_to_lcl(x), quantum_number, n);
}

double element::at_lcl(double x, unsigned int quantum_number)
{
	if(x<-1 || x>1)
	{
		return 0;
	}

	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients.col(quantum_number);
	return ret(0,0);
}

double element::at_lcl(double x, unsigned int quantum_number, unsigned int n)
{
	if(x<-1 || x>1)
	{	
		return 0;
	}

	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients.col(quantum_number);
	return ret(0,0);
}

double element::integration_wrapper(double x, void* params)
{
	//Set up the arguments pointer
	wrapper_args* param = (wrapper_args*)params;

	//Evaluate the operator at x with the correct wavefunctions
	return param->m_op(x, &element::polynomial_at, param->m, param->n, param->ptr, param->op_args);
}

Mat<double> element::lcl_mat( op m_op, void* args, gsl_integration_glfixed_table* table )
{
	//Create a return matrix
	Mat<double> ret(polynomial_coefficients.n_rows, polynomial_coefficients.n_rows);

	//set up the function and its arguments
	gsl_function integrate_me;
	struct wrapper_args int_args = {0, 0, m_op, args, this};

	//insert the corresponding wrapper function
	integrate_me.function = &element::integration_wrapper;
	integrate_me.params   = &int_args;

	//iterate through each element of the matrix
	for(int i=0; i<polynomial_coefficients.n_rows; i++)
	{
		int_args.m=i;
		for(int j=0; j<polynomial_coefficients.n_rows; j++)
		{
			int_args.n=j;

			//Integrate
			ret(i, j)=gsl_integration_glfixed(&integrate_me, bottom, top, table);
		}
	}

	return ret;
}

double element::polynomial_at(double x, unsigned int m, unsigned int n_der)
{
	//Return the nth derivative of the mth 
	Mat<double> ret = pow(2/(top - bottom), n_der)*( variable_row( glbl_to_lcl(x), polynomial_coefficients.n_rows, n_der )*polynomial_coefficients.col(m) );


	return ret(0,0);
}


double element::get_top()
{
	return top;
}

bool element::set_bottom(double x)
{
	bottom = x;
}

bool element::set_top(double x)
{
	top = x;
}

double element::get_bottom()
{
	return bottom;
}

unsigned int element::size()
{
	return polynomial_coefficients.n_rows;
}
