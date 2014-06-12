#ifndef __MATRIX_UTILITIEZ__
#define __MATRIX_UTILITIEZ__

#include <armadillo>

using namespace arma;

//Return factorial of n
double factorial(int n);

//Returns a row with <size> elements of the form  (1, x, x^2, x^3 ... x^<size)
Row<double> variable_row(double x, unsigned int size);

//Returns the <n>th derivative of the above row
Row<double> variable_row(double x, unsigned int size, unsigned int n);

#endif
