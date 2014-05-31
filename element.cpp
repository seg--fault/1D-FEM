#include<algorithm>
#include<armadillo>
#include<gsl/gsl_integration.h>
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
	element_solved=false;
	nDOF=0;
	nnodes=0;
	top=1;
	bottom=-1;
}

element::element(unsigned int DOF, vector<double> nodes, double arg_bottom, double arg_top)
{
	polynomial_initialized=false;
	element_solved=false;
	top = arg_top;
	bottom = arg_bottom;

	initialize_polynomial(DOF, nodes);
}

bool element::initialize_polynomial(unsigned int DOF, vector<double> nodes)
{
	if(polynomial_initialized)
		cout<<"Warning. The polynomial has already been initialized.  Reinitializing anyway";

	if(DOF<1)
	{
		cout<<"Error. Bad value for DOF";
		return false;
	}

	//add check to remove duplicate nodes, make sure they are in range and use -1, 1
	
	if(nodes.size()<2)
	{
		cout<<"Error. Too few nodes to construct polynomial.";
		return false;
	}

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
	nDOF=DOF;
	nnodes=nodes.size();
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
	polynomial_initialized=true;
}

Col<double> element::get_coefficients()
{
	if(!element_solved)
	{
		cout << "Error. Returning incorrect coefficients.";
	}

	return element_coefficients;
}

void element::set_coefficients(Col<double> coefficients)
{
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
		return 0;
	Mat<double> ret = (variable_row(x, polynomial_coefficients.n_cols)*polynomial_coefficients)*element_coefficients;
	return ret(0,0);
}

double element::at_lcl(double x, unsigned int n)
{
	if(x<-1 || x>1)
		return 0;
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

struct function_mat_args{unsigned int m; unsigned int n; unsigned int n_der; element* ptr;};

double function_mat_for_real(double x, void* param)
{
	function_mat_args* args = (function_mat_args*)param;
	
	return args->ptr->function_mat(x, args->m, args->n, args->n_der);
}

struct hamiltonian_function_mat_args{unsigned int m; unsigned int n; unsigned int n_der; element* ptr; double (*potential)(double x, void* args); void* potential_args;};

double hamiltonian_function_mat(double x, void* params)
{
	hamiltonian_function_mat_args* param = (hamiltonian_function_mat_args*)params;

	function_mat_args other_param={param->m, param->n, param->n_der+1, param->ptr};		//set up to get the second derivative of function mat
	double h = param->ptr->get_top()-param->ptr->get_bottom();				//calculate scaling factor
	double T = 4/h/h*function_mat_for_real(x, &other_param);				//calculate kinetic energey ie. 1/h^2(p*')(p')
	
	other_param.n_der=param->n_der;
	return T + function_mat_for_real(x, &other_param)*param->potential(param->ptr->lcl_to_glbl(x), param->potential_args);		//add it to the potential energy and return
}

Mat<double> element::lcl_hamiltonian_mat(double (*potential)(double x, void* args), void* args, gsl_integration_glfixed_table* table)
{
	Mat<double> ret(nnodes*nDOF, nnodes*nDOF);

	gsl_function integrate_me;
	struct hamiltonian_function_mat_args integrate_me_params={0, 0, 0, this, potential, args};

	integrate_me.function	= &hamiltonian_function_mat;
	integrate_me.params	= &integrate_me_params;

	for(int i=0; i<nnodes*nDOF; i++)
	{
		integrate_me_params.m=i;
		for(int j=0; j<nnodes*nDOF; j++)
		{
			integrate_me_params.n=j;
			ret(i, j)=(top-bottom)/2*gsl_integration_glfixed(&integrate_me, -1, 1, table);
		}
	}

	return ret;
}

Mat<double> element::lcl_energy_mat(gsl_integration_glfixed_table* table)
{
	Mat<double> ret(nnodes*nDOF, nnodes*nDOF);	//return matrix

	gsl_function integrate_me;
	struct function_mat_args args = {0, 0, 0, this};

	integrate_me.function = &function_mat_for_real;
	integrate_me.params   = &args;

	for(int i=0; i<nnodes*nDOF; i++)
	{
		args.m=i;
		for(int j=0; j<nnodes*nDOF; j++)
		{
			args.n=j;
			ret(i, j)=(top-bottom)/2*gsl_integration_glfixed(&integrate_me, -1, 1, table);
		}
	}

	return ret;
}

double element::function_mat(double x, unsigned int m, unsigned int n, unsigned int n_der)
{
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
