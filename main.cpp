#include<iostream>

#include"wavefunction.h"

using namespace std;

int main()
{
	wavefunction test;

	unsigned int myDOF;
	int nnodes;
	vector<double> mynodes;

	cout<<"DOF per node: ";
	cin>>myDOF;
	cout<<"Number of nodes: ";
	cin>>nnodes;

	for(int i=0;i<nnodes;i++)
	{
		cout<<"Node #"<<i<< " position: ";
		double pos;
		cin>>pos;
		cout<<"\n";
		if(pos>1 || pos<-1)
		{
			cout<<"Invalid Postion (-1<x<1)\n";
			i--;
			continue;
		}
		mynodes.push_back(pos);
	}
	
	test.initialize_polynomial(myDOF,mynodes);
	
	return 0;
}
