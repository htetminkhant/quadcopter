#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
using namespace std;
int main()
{
	double start_time=0.0,end_time=10.0, dt=0.005;
	int N=0;
	std::vector<double>myvec;
	myvec.push_back(start_time);
	while(start_time<=end_time)
	{  
		start_time+=dt;
		myvec.push_back(start_time);
		
	}
	for(int i=0;i<myvec.size();i++)
	{   
		++N;
		std::cout << myvec[i]<<" ";

	}
	cout<<endl;
	cout<<"N"<<"="<<N;
	cout<<endl;
	return 0;
}
