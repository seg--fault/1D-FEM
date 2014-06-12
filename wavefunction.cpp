#include<algorithm>
#include<armadillo>
#include<gsl/gsl_integration.h>
#include<iostream>
#include<math.h>

#include "wavefunction.h"

using namespace arma;
using namespace std;

//Immediately create elements and then solve the problem
wavefunction::wavefunction(vector<double> endpoints, vector<double> nodes, unsigned int DOF_arg, op hamiltonian, op energy, void* ham_args, void* en_args)
{
	//Start off with no elements
	n_elements = 0;

	//Set the corresponding degree of freedom value
	DOF = DOF_arg;

	//Create and initialize the elements based on arguments
	create_elements(endpoints, nodes, DOF);

	//Solve the system with the given potential
	solve(hamiltonian, energy, ham_args, en_args);
}

//Frees up all assigned memory
wavefunction::~wavefunction()
{
	//Unallocate memory for the elements
	delete[] elements;
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

//Creates all elements needed for solving equation
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

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Allocate memory for the new elements
	elements = new element[endpoints.size() - 1];

	//Make sure we know how many endpoints there are
	n_elements = endpoints.size() - 1;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Vector to hold which nodes are in the current element
	vector<double> nodes_in_crnt_element;

	//Go through the list backwards
	for(int i = n_elements; i>0; i--)
	{
		//Add the first node (last on the stack) to the current node list
		nodes_in_crnt_element.push_back(nodes.back());

		//Go through the list and continue to add the last node until it isn't in the element
		do
		{
			//Get rid of the last last (smallest) nodes
			nodes.pop_back();

			//Add the last node to the list as it must be in there
			nodes_in_crnt_element.push_back(nodes.back());

		//Compare after casting to floats so bad things don't happen
		} while((float)nodes.back() < (float)endpoints.at(i-1));

		//Initialize the current element
		elements[n_elements - i].set_bottom(endpoints.at(i));
		elements[n_elements - i].set_top(endpoints.at(i-1));
		elements[n_elements - i].initialize_polynomial(DOF, nodes_in_crnt_element);

		//Zero out the nodes in this element
		nodes_in_crnt_element.clear();
	}
	
	return true;	
}


//Comparison function for sorting eigenvalues
bool eigen_sort_compare(pair<double, Col<double> > i, pair<double, Col<double> > j)
{
	return i.first < j.first;
}

//Find the element coefficients, sort them and insert them into the elements
bool wavefunction::solve(op hamiltonian, op energy, void* ham_args, void* en_args)
{
	//the total number of cefficients in the wavefunction
	unsigned int total_coefficients = 0;

	//calculate the total by summing
	for(int i=0; i<n_elements; i++)
	{
		total_coefficients+=elements[i].size();
	}

	//subtract out the degeneracies
	total_coefficients-=DOF*(n_elements-1);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Interaction...
	cout << "Constructing Global Matrices...\n";

	//Create the global matrices
	Mat<double> glbl_energy_mat		(total_coefficients, total_coefficients, fill::zeros); 
	Mat<double> glbl_hamiltonian_mat	(total_coefficients, total_coefficients, fill::zeros);

	//initialize the shared integration table
	gsl_integration_glfixed_table* integration_table =  gsl_integration_glfixed_table_alloc (32);

	//Set up the current offsets for the submatrices
	unsigned int mat_offset_1 = 0;
	unsigned int mat_offset_2 = 0;

	//Add each sub matrix to the global matrix
	for(int i=0; i<n_elements; i++)
	{
		//Set the correct matrix offsets
		mat_offset_1 =  mat_offset_2;
		mat_offset_2 += elements[i].size()-DOF;
		
		//Add the lcl matrices to the global matrix
		glbl_energy_mat.submat(mat_offset_1, mat_offset_1, mat_offset_2+DOF-1, mat_offset_2+DOF-1) 	+= elements[i].lcl_mat(energy, en_args, integration_table);
		glbl_hamiltonian_mat.submat(mat_offset_1, mat_offset_1, mat_offset_2+DOF-1, mat_offset_2+DOF-1) += elements[i].lcl_mat(hamiltonian, ham_args,  integration_table);

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Interaction...
	cout << "Solving Matrix system...\n";

	//Create containers to hold the results of solving the global matrix
	cx_vec cx_eigenvalues;
	cx_mat cx_eigenvectors;

	//Solve the huge system.  Oh   God..
	eig_pair(cx_eigenvalues, cx_eigenvectors, glbl_hamiltonian_mat, glbl_energy_mat);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	//Interact
	cout << "Normalizing Solutions...\n";

	//Row eigenvectors because armadillo is dumb
	cx_mat cx_eigenvectors_transpose = cx_eigenvectors.t();

	//Go through each solution
	for(int i=0; i<glbl_energy_mat.n_rows; i++)
	{
		//The norm of the solution
		cx_mat norm;

		//Perform the "integration"
		norm = cx_eigenvectors_transpose.row(i)*glbl_energy_mat*cx_eigenvectors.col(i);

		//Actually normalize that solution
		cx_eigenvectors.col(i) = (real(cx_eigenvectors(0, i)) < 0 ? -1 : 1 )*cx_eigenvectors.col(i) / real(sqrt(norm(0,0)));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
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
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Create a temporary storage for our sorted eigenvectors
	Mat<double> eigenvectors;

	//Resize our eigenvector and eigenvalue storage
	eigenvalues.set_size(cx_eigenvalues.n_rows);
	eigenvectors.set_size(cx_eigenvectors.n_rows, cx_eigenvectors.n_cols);

	//Insert the sorted solutions into our storage
	for(int i=0; i<cx_eigenvalues.n_rows; i++)
	{
		eigenvalues(i) = eigen_sort.at(i).first;
		eigenvectors.col(i) = eigen_sort.at(i).second;
	}

	//Move the eigenvectors back to the elements
	unsigned int coefficient_offset_a = 0;
	unsigned int coefficient_offset_b = 0;

	for(int i=0; i<n_elements; i++)
	{
		//Set the ending coefficient
		coefficient_offset_b = coefficient_offset_a + elements[i].size() - 1;

		//Transfer the coefficients
		mat temp = eigenvectors.rows(coefficient_offset_a, coefficient_offset_b);
		elements[i].set_coefficients(temp);

		//Set the new beginning coefficient
		coefficient_offset_a = coefficient_offset_b + 1 - DOF;
	}

	return true;
}

double wavefunction::operator()(double x, unsigned int quantum_number)
{
	//If trying to access an eigenvector which doesn't exist
	if(quantum_number >= eigenvalues.n_rows)
	{
		cout << "Error.  Eigenvector Doesn't exist.";
		return 0;
	}

	if(x < lower_bound() || x > upper_bound())
	{
		cout << "x is out of bounds\n";
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Holds information for binary search
	int binary_search_front = 0;
	int binary_search_back  = n_elements - 1;
	int binary_search_middle;
	int compare = -1;

	//Now go through the elements until we have found the correct one
	do
	{
		//Find the midpoint
		binary_search_middle = (binary_search_front + binary_search_back)/2;

		//Check if we are in the middle element
		compare = elements[binary_search_middle].in(x);

		//If x is below the element
		if(compare == 1)
			binary_search_back = binary_search_middle - 1;
		//If x is above the element
		if(compare == 2)
			binary_search_front = binary_search_middle + 1;
	}
	while(compare != 0);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Return the value of the polynomial there
	return elements[binary_search_middle].at(x, quantum_number);
}


double wavefunction::eigenvalue(unsigned int quantum_number)
{
	//Check if the quantum number is in range
	if(quantum_number >= eigenvalues.n_rows)
	{
		cout << "Error. No such solution";
		return 0;
	}
	
	//Hand back the proper eigenvalue
	return eigenvalues(quantum_number);
}

unsigned int wavefunction::n_eigenvalues()
{
	return eigenvalues.n_rows;
}

double wavefunction::lower_bound()
{
	//Return the lower limit of the first element
	return elements[0].get_bottom();
}

double wavefunction::upper_bound()
{
	//Return the upper limit of the last element
	return elements[n_elements - 1].get_top();
}
