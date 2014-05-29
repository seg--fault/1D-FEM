#ifndef __ELEMENT__
#define __ELEMENT__

#include<armadillo>
#include<vector>

using namespace arma;
using namespace std;

class element
{
	public:
	element();

	bool		initialize_polynomial(unsigned int DOF, vector<double> nodes);
	Mat<double>	get_polynomial();
	void		set_polynomial(Mat<double> coefficients);

	Col<double> 	get_coefficients();
	void		set_coefficients(Col<double> coefficients);

	double		glbl_to_lcl(double x);
	double		lcl_to_glbl(double x);

	double		at(double x);
	double		at(double x, unsigned int n);
	double		at_lcl(double x);
	double		at_lcl(double x, unsigned int n);

	Row<double>	variable_row(double x, unsigned int size);				//returns the variable row for x
	Row<double>	variable_row(double x, unsigned int size, unsigned int n);		//nth derivative row

	Mat<double>	lcl_matrix();

	private:
	bool 		polynomial_initialized;
	Mat<double> 	polynomial_coefficients;
	bool		element_solved;
	Col<double>	element_coefficients;		//holds coefficients for approximation
	double		top, bottom;			//upper and lower limits of the element
};

#endif
