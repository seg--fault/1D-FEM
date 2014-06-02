#ifndef __ELEMENT__
#define __ELEMENT__

#include<armadillo>
#include<gsl/gsl_integration.h>
#include<vector>


using namespace arma;
using namespace std;

class element
{
	public:
	element();
	element(unsigned int DOF, vector<double> nodes, double arg_bottom, double arg_top);

	//Methods for dealing with interpolation polynomials
	bool		initialize_polynomial(unsigned int DOF, vector<double> nodes);		//Warning.  Current method can only handle 16 nodes per element w/ DOF=2
	Mat<double>	get_polynomial();
	void		set_polynomial(Mat<double> coefficients);

	//Methods for element coefficients
	Col<double> 	get_coefficients();
	void		set_coefficients(Col<double> coefficients);

	//local variable and global variable conversion methods
	double		glbl_to_lcl(double x);
	double		lcl_to_glbl(double x);

	//Whether x is in the element or not
	bool		in(double x);

	//Methods for finding the value of the wavefunction within this element
	double		at(double x);
	double		at(double x, unsigned int n);
	double		at_lcl(double x);
	double		at_lcl(double x, unsigned int n);

	//Methods for creating row vectors and their derivatives where row(i) = x^i
	Row<double>	variable_row(double x, unsigned int size);
	Row<double>	variable_row(double x, unsigned int size, unsigned int n);	

	//Methods for constructing the matrices needed for the variational method
	Mat<double>	lcl_hamiltonian_mat(double (*potential)(double x, void* args), void* args, gsl_integration_glfixed_table* table); //calculates the hamiltonian matrix for a potential
	Mat<double>	lcl_energy_mat(gsl_integration_glfixed_table* table);	//calculate the energy matrix for this element

	//Function and wrappers for calculating matrices
	double		function_mat(double x, unsigned int m, unsigned int n, unsigned int n_der);	//returns the m'th times the n'th interpolation polynomial
	static double	energy_function_wrapper(double x, void* args);					//energy matrix wrapper
	static double	hamiltonian_function_wrapper(double x, void* args);				//hamiltonian wrapper (uses potential)

	//Read only access for private variables
	double 		get_top();
	double		get_bottom();

	private:
	Mat<double> 	polynomial_coefficients;	//holds the coefficients of each polynomial, one per collumn
	bool 		polynomial_initialized;		//whether the polynomial coefficients exist or not

	Col<double>	element_coefficients;		//holds coefficients for approximation
	bool		element_solved;
	
	double		top, bottom;			//upper and lower limits of the element
};

//structure to pass arguments to the energy wrapper
struct function_mat_args{unsigned int m; unsigned int n; unsigned int n_der; element* ptr;};

//structure to pass arguments to the hamiltonian wrapper
struct hamiltonian_function_mat_args{unsigned int m; unsigned int n; unsigned int n_der; element* ptr; double (*potential)(double x, void* args); void* potential_args;};

#endif
