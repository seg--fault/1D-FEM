#ifndef __ELEMENT__
#define __ELEMENT__

#include<armadillo>
#include<gsl/gsl_integration.h>
#include<vector>


using namespace arma;
using namespace std;

class element;

//Operator function type
typedef double (*op)(double x, double (element::*wf)(double x, unsigned int m, unsigned int der), unsigned int m, unsigned int n, element* ptr, void* args);

class element
{
	public:
	element();
	element(unsigned int DOF, vector<double> nodes, double arg_bottom, double arg_top);

	//Methods for dealing with interpolation polynomials
	bool		initialize_polynomial(unsigned int DOF, vector<double> nodes);
	Mat<double>	get_polynomial();
	void		set_polynomial(Mat<double> coefficients);

	//Methods for element coefficients
	Mat<double> 	get_coefficients();
	void		set_coefficients(Mat<double> coefficients);

	//local variable and global variable conversion methods
	double		glbl_to_lcl(double x);
	double		lcl_to_glbl(double x);

	//Whether x is in the element or not
	int		in(double x);

	//Methods for finding the value of the wavefunction within this element
	double		at(double x, unsigned int quantum_number);
	double		at(double x, unsigned int quantum_number, unsigned int n);
	double		at_lcl(double x, unsigned int quantum_number);
	double		at_lcl(double x, unsigned int quantum_number, unsigned int n);

	//Methods for constructing the local matrices from operators
	Mat<double>	lcl_mat(op m_op, void* args, gsl_integration_glfixed_table* table ); 

	//Function and wrappers for calculating matrices
	double		polynomial_at(double x, unsigned int m, unsigned int n_der);	
	static double	integration_wrapper(double x, void* params);			

	//Read only access for private variables
	double 		get_top();
	double		get_bottom();
	bool 		set_top(double x);
	bool		set_bottom(double x);
	unsigned int 	size();

	//Private:  to be fixed
	Mat<double> 	polynomial_coefficients;	//holds the coefficients of each polynomial, one per collumn
	bool 		polynomial_initialized;		//whether the polynomial coefficients exist or not

	private:
	Mat<double>	element_coefficients;		//holds coefficients for approximation
	bool		element_solved;
	
	double		top, bottom;			//upper and lower limits of the element
};

//structure to pass arguments to the wrappers
struct wrapper_args{unsigned int m; unsigned int n; op m_op; void* op_args; element* ptr;};

#endif
