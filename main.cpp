#include<armadillo>
#include<algorithm>
#include<boost/tuple/tuple.hpp>
#include<iostream>
#include<math.h>

#include"element.h"
#include"gnuplot-iostream/gnuplot-iostream.h"

using namespace arma;
using namespace std;

//potential for testing element solution
double test_potential(double x, void* args)
{
	return pow(x,2);	//change lower bound to -x_max
	return -2/x;		//change lower bound to zero
}

//comparison function to sort eigen(function, value) combination
bool eigen_order_test(pair<double, Col<double> > i, pair<double, Col<double> > j)
{
	return i.first < j.first;
}

int main()
{
	unsigned int test_DOF=0;
	vector<double> test_nodes;

	double x_max = 0;

	while(x_max<=0)
	{
	cout << "X-Max: ";
	cin >> x_max;
	}

	int num_nodes=0;
	while(num_nodes<2)
	{
		cout << "Number of nodes: ";
		cin >>num_nodes;
	}

	while(test_DOF<1)
	{
		cout <<"DOF: ";
		cin >> test_DOF;
	}
	

	for(double x=-1;x<=1;x+=2/((double)num_nodes-1))
	{
		test_nodes.push_back(x);
	}

	element test_element(test_DOF, test_nodes, -1*x_max, x_max);

	//set up an integration table
	gsl_integration_glfixed_table* integration_table =  gsl_integration_glfixed_table_alloc (32);

	cout << "Calculating Matrices\n";

	//get the energy matrix
	Mat<double> NRG_mat = test_element.lcl_energy_mat(integration_table);

	//get the hamiltonian matrix
	Mat<double> HAM_mat = test_element.lcl_hamiltonian_mat(&test_potential, &test_element, integration_table);

	//get rid of the old table
	gsl_integration_glfixed_table_free (integration_table);

	//set and solve the eigenvalue problem
	cx_vec	eigenvalues;
	cx_mat	eigenvectors;
	cout << "Solving System.\n";
	eig_pair(eigenvalues, eigenvectors, HAM_mat, NRG_mat);

	//place to store the eigen(functions, values)
	vector<pair<double, Col<double> > > eigen_order;

	//insert the eigen(functions, values) to be sorted
	for(int i=0;i<eigenvalues.n_rows;i++)
	{
		eigen_order.push_back(make_pair(real(eigenvalues(i)), real(eigenvectors.col(i))));
	}

	//sort the eigen(functions, values) by eigenvalue
	sort(eigen_order.begin(), eigen_order.end(), eigen_order_test);

	
	//plot the eigenfuncions
	Gnuplot gp;

	gp << "set xrange [" << test_element.get_bottom() << ":"<< test_element.get_top() << "]\nset yrange [-0.25:0.25]\n";

	int i=0;
	vector<pair<double, double> > samples;

	while(true)
	{
		gp << "\n";

		cout<<"Show eigenfunction #";
		cin>>i;

		if(i<0 || i>=test_DOF*test_nodes.size())
		{
			cout<<"Quitting";
			break;
		}

		cout << "Eigenfunction = " << eigen_order.at(i).first << "\n";

		test_element.set_coefficients(eigen_order.at(i).second);

		samples.clear();

		for(double x=test_element.get_bottom();x<test_element.get_top();x+=0.01)
		{
			samples.push_back(make_pair(x,test_element.at(x)));
		}

		gp << "plot " << gp.file1d(samples) << " with lines title ''" << endl;
	}


	gp << "\n";

	return 0;
}
