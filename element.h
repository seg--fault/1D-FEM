#ifndef __ELEMENT__
#define __ELEMENT__

#include<armadillo>
#include<gsl/gsl_integration.h>
#include<vector>

#include"potential.h"

using namespace arma;
using namespace std;

class element
{
	public:
	element();
	element(unsigned int DOF, vector<double> nodes, double arg_bottom, double arg_top);

	bool		initialize_polynomial(unsigned int DOF, vector<double> nodes);
	Mat<double>	get_polynomial();
	void		set_polynomial(Mat<double> coefficients);

	Col<double> 	get_coefficients();
	void		set_coefficients(Col<double> coefficients);

	double		glbl_to_lcl(double x);
	double		lcl_to_glbl(double x);

	bool		in(double x);

	double		at(double x);
	double		at(double x, unsigned int n);
	double		at_lcl(double x);
	double		at_lcl(double x, unsigned int n);

	Row<double>	variable_row(double x, unsigned int size);				//returns the variable row for x	//probably should be moved elsewhere for neatness
	Row<double>	variable_row(double x, unsigned int size, unsigned int n);		//nth derivative row

	Mat<double>	lcl_hamiltonian_mat(double (*potential)(double x, void* args), void* args, gsl_integration_glfixed_table* table); //calculates the hamiltonian matrix for a potential
	Mat<double>	lcl_energy_mat(gsl_integration_glfixed_table* table);	//calculate the energy matrix for this element

	double		function_mat(double x, unsigned int m, unsigned int n, unsigned int n_der);

	double 		get_top();
	double		get_bottom();

	private:
	unsigned int	nDOF, nnodes;			//degrees of freedom per node and the number of nodes for the approximation polynomials
	Mat<double> 	polynomial_coefficients;	//holds the coefficients of each polynomial, one per collumn
	bool 		polynomial_initialized;		//whether the polynomial coefficients exist or not

	Col<double>	element_coefficients;		//holds coefficients for approximation
	bool		element_solved;
	
	double		top, bottom;			//upper and lower limits of the element
};

#endif
