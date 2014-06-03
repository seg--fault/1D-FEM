#include<algorithm>
#include<armadillo>
#include<gsl/gsl_integration.h>
#include<iostream>
#include<math.h>

#include "element.h"


using namespace arma;
using namespace std;

//used in calculating variable row derivatives
double factorial(int n)
{
	if(n <= 0)  return 1 ; 		//safeguard 0 and negative values
	double res = n ;
	while(--n>1) res *= n ;		//start multiplying
	return res ;
}

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

	//Debugging...
	cout << "Creating element with top = " << top << " and bottom = " << bottom << "\n";

	//The polynomial and element haven't been solved yet
	polynomial_initialized=false;
	element_solved=false;

	initialize_polynomial(DOF, nodes);
}

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

Mat<double> element::get_polynomial()
{
#ifdef ERROR_CHECKING
	if(!polynomial_initialized)
	{
		cout<<"Error, reading uninitialized polynomial coefficients.";
	}
#endif

	return polynomial_coefficients;
}

void element::set_polynomial(Mat<double> coefficients)
{
#ifdef ERROR_CHECKING
	if(polynomial_initialized)
	{
		cout<<"Warning. Overwriting polynomial coefficients!";
	}
#endif

	polynomial_coefficients = coefficients;
	polynomial_initialized=true;
}

Col<double> element::get_coefficients()
{
#ifdef ERROR_CHECKING
	if(!element_solved)
	{
		cout << "Error. Returning incorrect coefficients.";
	}
#endif

	return element_coefficients;
}

void element::set_coefficients(Col<double> coefficients)
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

bool element::in(double x)
{
	return ((x>=bottom)&&(x<=top));
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
	if(x<-1 || x>1)
	{
		return 0;
	}

	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients;
	return ret(0,0);
}

double element::at_lcl(double x, unsigned int n)
{
	if(x<-1 || x>1)
	{	
		return 0;
	}

	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients;
	return ret(0,0);
}

Row<double> element::variable_row(double x, unsigned int size)
{
	Row<double> ret(size);

	ret(0) = 1;

#ifdef ERROR_CHECKING
	if(size<=0)
	{
		cout << "Error. Can't make variable row of size <= 0.";
		return ret;
	}
#endif

	for(int i=1; i<size; i++)
	{
		ret(i) = pow(x, (double)(i));
	}

	return ret;
}

Row<double> element::variable_row(double x, unsigned int size, unsigned int n)
{
	Row<double> ret(size);

#ifdef ERROR_CHECKING
	if(size<=0)
	{
		cout << "Error.  Cannot create variable row of size <= 0";
		return ret;
	}
#endif
	for(int i=0; i<size; i++)
	{
		if(i<n)
			ret(i) = 0;		//derivative of a constant
		if(i==n)
			ret(i) = 1;		//is the original constant
		if(i>n)
			ret(i) = pow(x, (double)(i-n))*factorial(i)/factorial(i-n);	//add in the correct values
	}

	return ret;
}

double element::energy_function_wrapper(double x, void* param)
{
	//Create a pointer used in accessing arguments
	function_mat_args* args = (function_mat_args*)param;
	
	//leading x^2 for calculating hydrogen atom
	return /*pow(args->ptr->lcl_to_glbl(x),2)*/(args->ptr->function_mat(x, args->m, args->n, args->n_der));
}

double element::hamiltonian_function_wrapper(double x, void* params)
{
	//Set up the arguments pointer and a scaling factor used for the transition from global to local coordinates
	hamiltonian_function_mat_args* param = (hamiltonian_function_mat_args*)params;
	double h = param->ptr->get_top()-param->ptr->get_bottom();				

	//Feed the arguments from this function into the arguments required for the energy function
	function_mat_args other_param={param->m, param->n, param->n_der+1, param->ptr};

	//Calculate the term in the action corresponding to kinetic energy
	double T = 4/h/h*element::energy_function_wrapper(x, &other_param);
	
	//Change to a first derivative and calculate the the term in the action from the potential energy
	other_param.n_der=param->n_der;
	double V = element::energy_function_wrapper(x, &other_param)*param->potential(param->ptr->lcl_to_glbl(x), param->potential_args);

	//Add the potential action to the kinetic action and return them
	return T + V;
}

Mat<double> element::lcl_hamiltonian_mat(double (*potential)(double x, void* args), void* args, gsl_integration_glfixed_table* table)
{
	//Initialize a return matrix
	Mat<double> ret(polynomial_coefficients.n_rows, polynomial_coefficients.n_rows);

#ifdef ERROR_CHECKING
	if(!polynomial_initialized)
	{
		cout << "Error. Trying to create hamiltonian without initialized polynomials";
	}
#endif

	//create the function to be integrated along with its arguments
	gsl_function integrate_me;
	struct hamiltonian_function_mat_args integrate_me_params={0, 0, 0, this, potential, args};

	//insert into the gsl function data structure
	integrate_me.function	= &element::hamiltonian_function_wrapper;
	integrate_me.params	= &integrate_me_params;

	//go through each element of the matrix and perform the integration on the corresponding pair of function
	for(int i=0; i<polynomial_coefficients.n_rows; i++)
	{
		integrate_me_params.m=i;	
		for(int j=0; j<polynomial_coefficients.n_rows; j++)
		{
			//Integrate the ith times the jth polynomials
			integrate_me_params.n=j;
			ret(i, j)=(top-bottom)/2*gsl_integration_glfixed(&integrate_me, -1, 1, table);
		}
	}

	//Debugging...
	//cout << "lcl hamiltonian mat for <" << bottom << ", " << top << ">:\n";
	//ret.print();
	//cout <<"\n";

	//Hand back the results
	return ret;
}

Mat<double> element::lcl_energy_mat(gsl_integration_glfixed_table* table)
{
	//Create a return matrix
	Mat<double> ret(polynomial_coefficients.n_rows, polynomial_coefficients.n_rows);

	//set up the function and its arguments
	gsl_function integrate_me;
	struct function_mat_args args = {0, 0, 0, this};

	//insert the corresponding wrapper function
	integrate_me.function = &element::energy_function_wrapper;
	integrate_me.params   = &args;

	//integrate
	for(int i=0; i<polynomial_coefficients.n_rows; i++)
	{
		args.m=i;
		for(int j=0; j<polynomial_coefficients.n_rows; j++)
		{
			args.n=j;
			ret(i, j)=(top-bottom)/2*gsl_integration_glfixed(&integrate_me, -1, 1, table);
		}
	}

	//Debugging...
	//cout << "lcl energy matrix for <" << bottom << ", " << top << ">:\n";
	//ret.print();
	//cout <<"\n";

	return ret;
}

double element::function_mat(double x, unsigned int m, unsigned int n, unsigned int n_der)
{
	//Return the nth derivative of the mth polynomial times the nth derivative of the nth polynomial at x
	Mat<double> ret = (variable_row(x,polynomial_coefficients.n_rows, n_der)*polynomial_coefficients.col(m))%(variable_row(x,polynomial_coefficients.n_rows, n_der)*polynomial_coefficients.col(n));

	return ret(0,0);
}


double element::get_top()
{
	return top;
}

double element::get_bottom()
{
	return bottom;
}

unsigned int element::size()
{
	return polynomial_coefficients.n_rows;
}
