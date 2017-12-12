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
	for(int i=0;i<3;i++)
	{
		cout<<endl<<x[i][0];
	}
	cout<<endl;
	int xdot[3][1]={{0},{0},{0}};
	for(int i=0;i<3;i++)
	{
		cout<<endl<<xdot[i][0];
	}
	cout<<endl;
	int theta[3][1]={{0},{0},{0}};
	for(int i=0;i<3;i++)
	{
		cout<<endl<<theta[i][0];
	}
	cout<<endl;
	int deviation =100;
	double pi=3.142;
	double degree [3][1] ;
	for (int i=0;i<3;i++)
	{
	degree[i][0]=((rand()%200+(-100))*pi)/180;
	cout<<endl<<degree[i][0];
	}
	cout<<endl;
	return 0;
}
