#include<armadillo>
#include<boost/tuple/tuple.hpp>
#include<iostream>
#include<math.h>

#include"element.h"
#include"gnuplot-iostream/gnuplot-iostream.h"

using namespace arma;
using namespace std;

double test_potential(double x, void* args)
{
	return pow(x,2);
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

	//get the energy matrix
	Mat<double> NRG_mat = test_element.lcl_energy_mat(integration_table);
	NRG_mat.print();

	cout<<"\n";

	//get the hamiltonian matrix
	Mat<double> HAM_mat = test_element.lcl_hamiltonian_mat(&test_potential, &test_element, integration_table);
	HAM_mat.print();

	//get rid of the old table
	gsl_integration_glfixed_table_free (integration_table);

	//set and solve the eigenvalue problem
	cx_vec	eigenvalues;
	cx_mat	eigenvectors;
	eig_pair(eigenvalues, eigenvectors, HAM_mat, NRG_mat);

	cout<<"\n";
	eigenvalues.print();
	cout<<"\n";
	eigenvectors.print();

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

		cout << "Eigenfunction = " << eigenvalues(i) << "\n";

		test_element.set_coefficients(real(eigenvectors.col(i)));

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
