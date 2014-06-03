#include<algorithm>
#include<armadillo>
#include<gsl/gsl_integration.h>
#include<iostream>
#include<math.h>

#include "wavefunction.h"

using namespace arma;
using namespace std;

wavefunction::wavefunction()
{
}

wavefunction::wavefunction(vector<double> endpoints, vector<double> nodes, unsigned int DOF_arg, double (*potential)(double x, void* args), void* args)
{
	//Set the corresponding degree of freedom value
	DOF = DOF_arg;

	//Create and initialize the elements based on arguments
	create_elements(endpoints, nodes, DOF);

	//Solve the system with the given potential
	solve(potential, args);
}


//Comparison function to reverse sort nodes and endpoints
bool reverse_sort(double i, double j)
{
	return i>j;
}

//Comparison function to compare doubles accurately
bool double_compare(double i, double j)
{
	//So dirty...
	return (float)i == (float)j;
}

bool wavefunction::create_elements(vector<double> endpoints, vector<double> nodes, unsigned int DOF)
{
	if(DOF == 0)
	{
		cout << "Error. Cannot create elements with DoF = 0";
		return false;
	}

	//Make sure that the endpoints are also nodes
	nodes.insert(nodes.end(), endpoints.begin(), endpoints.end());
	
	//Make sure that all inputs are in order (reverse order to use vector.pop_back() )
	sort(endpoints.begin(),endpoints.end(), reverse_sort);
	sort(nodes.begin(), nodes.end(), reverse_sort);

	//Remove duplicate entries in endpoints using the double_compare function
	vector<double>::iterator it;
	it = unique (endpoints.begin(), endpoints.end(), double_compare);  
	endpoints.resize( distance(endpoints.begin(),it) );
	
	//and then in nodes
	it = unique (nodes.begin(), nodes.end(), double_compare);  
	nodes.resize(distance(nodes.begin(),it));

	if(endpoints.size()<2)
	{
		cout << "Error.  Too few endpoints.\n";
		return false;
	}

	//Vector to hold which nodes are in the current element
	vector<double> nodes_in_crnt_element;

	//Go through the list backwards
	for(int i=endpoints.size()-1; i>0; i--)
	{
		//Debugging Test...
		//cout << "Element Creation loop.  I = " << i << "\n";
		//cout << "Endpoints are <" << endpoints.at(i) << ", " << endpoints.at(i-1) << ">\n";

		//Calculate which nodes are in the element
		//Initialize the variables
		nodes_in_crnt_element.clear();

		//Add the last (smallest) node to the current node list
		nodes_in_crnt_element.push_back(nodes.back());

		//Debugging
		//cout << "Adding node at x = " << nodes.back() << "\n";

		double node_back, endpoint_back;
		//go through each node and add it to the list of current nodes if its is in here
		do
		{
			//Get rid of the last last (smallest) nodes
			nodes.pop_back();

			//Add the last node to the list as it must be in there
			nodes_in_crnt_element.push_back(nodes.back());

			//Debugging...
			node_back = nodes.back();
			endpoint_back = endpoints.at(i-1);

			//cout << "Adding node at x = " << nodes.back() << "\n";

		//Compare after casting to floats so bad things don't happen
		} while((float)nodes.back() < (float)endpoints.at(i-1));

		//Create a new element with the given parameters
		elements.push_back(new element(DOF, nodes_in_crnt_element, endpoints.at(i), endpoints.at(i-1)));
	}
	
	return true;	//Awesome
}


//Comparison function for sorting eigenvalues
bool eigen_sort_compare(pair<double, Col<double> > i, pair<double, Col<double> > j)
{
	return i.first < j.first;
}

bool wavefunction::solve(double (*potential)(double x, void* args), void* args)
{
	//the total number of cefficients in the wavefunction
	unsigned int total_coefficients = 0;

	//calculate the total by summing
	for(int i=0; i<elements.size(); i++)
	{
		total_coefficients+=elements.at(i)->size();
	}

	//subtract out the degeneracies
	total_coefficients-=DOF*(elements.size()-1);

	//Create the global matrices
	Mat<double> glbl_energy_mat		(total_coefficients, total_coefficients, fill::zeros); 
	Mat<double> glbl_hamiltonian_mat	(total_coefficients, total_coefficients, fill::zeros);

	//Set up the current offsets for the submatrices
	unsigned int mat_offset_1 = 0;
	unsigned int mat_offset_2 = 0;

	//initialize the shared integration table
	gsl_integration_glfixed_table* integration_table =  gsl_integration_glfixed_table_alloc (32);

	//Interaction...
	cout << "Constructing Global Matrices...\n";

	//Add each sub matrix to the global matrix
	for(int i=0; i<elements.size(); i++)
	{
		//Set the correct matrix offsets
		mat_offset_1 =  mat_offset_2;
		mat_offset_2 += elements.at(i)->size()-DOF;
		
		//Debugging
		cout << "Adding element <" << elements.at(i)->get_bottom() << ", " << elements.at(i)->get_top() << ">\n";

		//Add the lcl matrices to the global matrix
		glbl_energy_mat.submat(mat_offset_1, mat_offset_1, mat_offset_2+DOF-1, mat_offset_2+DOF-1) 	+= elements.at(i)->lcl_energy_mat(integration_table);
		glbl_hamiltonian_mat.submat(mat_offset_1, mat_offset_1, mat_offset_2+DOF-1, mat_offset_2+DOF-1) += elements.at(i)->lcl_hamiltonian_mat(potential, args,  integration_table);

	}

	//Debugging...
	cout << "Energy Matrix:\n";
	glbl_energy_mat.print();
	cout << "Hamiltonian Matrix:\n";
	glbl_hamiltonian_mat.print();

	//Interaction...
	cout << "Solving Matrix system...\n";

	//Create containers to hold the results of solving the global matrix
	cx_vec cx_eigenvalues;
	cx_mat cx_eigenvectors;

	//Solve the huge system.  Oh   God..
	eig_pair(cx_eigenvalues, cx_eigenvectors, glbl_hamiltonian_mat, glbl_energy_mat);

	//Interact
	cout << "Sorting Solutions...\n";

	//Container for sorting solutions
	vector<pair<double, Col<double> > > eigen_sort;

	//Create pairs of solutions to be sorted
	for(int i=0; i<cx_eigenvalues.n_rows; i++)
	{
		eigen_sort.push_back(make_pair(real(cx_eigenvalues(i)),real(cx_eigenvectors.col(i))));
	}

	//Sort the solutions
	sort(eigen_sort.begin(), eigen_sort.end(), eigen_sort_compare);
	
	//Resize our eigenvector and eigenvalue storage
	eigenvalues.set_size(cx_eigenvalues.n_rows);
	eigenvectors.set_size(cx_eigenvectors.n_rows, cx_eigenvectors.n_cols);

	//Insert the sorted solutions into our storage
	for(int i=0; i<cx_eigenvalues.n_rows; i++)
	{
		eigenvalues(i) = eigen_sort.at(i).first;
		eigenvectors.col(i) = eigen_sort.at(i).second;

		//Test,  print out the ordered eigenvalues
		cout << "Eigenvalue[" << i << "] = " << eigenvalues(i) << "\n";
	}

	return true;
}
