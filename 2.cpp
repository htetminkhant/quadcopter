// 2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <time.h>
using namespace std;
int main()
{	srand( time(0));
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
	int x[3][1]={{0},{0},{10}};
	int xdot[3][1]={{0},{0},{0}};
	int theta[3][1]={{0},{0},{0}};
	cout<<x[0][0]<<"\t"<<x[1][0]<<"\t"<<x[2][0];
	
	cout<<endl;
	
	cout<<endl;
	double randnumber;
	int deviation =100,temp;
	double pi=3.142;
	double degree [3][1] ;
	for (int i=0;i<3;i++)
	{
		randnumber=(double)  rand()/ RAND_MAX;
		temp=(2*deviation*randnumber-deviation);
		degree[i][0]=(temp*pi)/180;
		cout<<"Radian number is "<< randnumber <<" temp answer is "<< temp <<" is equal to " << degree[i][0]<<endl;
	}
	cout<<endl;
	_getch();
	return 0;
}
