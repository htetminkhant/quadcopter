#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
using namespace std;
void myvec1(vector<double>&);
double thetadot2omega(double thetadot[3][1],double theta[3][1]);
double t_;
int main()
{	srand( time(0));
std::vector<double>myvec2;
myvec1(myvec2);
int x[3][1]={{0},{0},{10}};
int xdot[3][1]={{0},{0},{0}};
double theta[3][1]={{0},{0},{0}};
double randnumber;
int deviation =100,temp;
double pi=3.142;
double thetadot [3][1] ;
for (int i=0;i<3;i++)
{
	randnumber=(double)  rand()/ RAND_MAX;
	temp=(2*deviation*randnumber-deviation);
	thetadot[i][0]=(temp*pi)/180;		
};
thetadot2omega(thetadot,theta);
cout<<endl;
_getch();
return 0;
}
void myvec1(vector<double>& myvec)
{
	double start_time=0.0,end_time=10.0, dt=0.005;
	int N=0;
	myvec.push_back(start_time);
	while(start_time<=end_time)
	{  
		start_time+=dt;
		myvec.push_back(start_time);

	}
	for(int i=0;i<myvec.size();i++)
	{   
		++N;
		t_=myvec[i];
		std::cout << t_<<" ";

	}
}
double thetadot2omega(double thetadot[3][1],double theta[3][1])
{
	double phi = theta[0][0] ,theta_ = theta[1][0], psi = theta[2][0];
	double w [3][3] = {
		{1,0,-sin(theta_)},
		{0,cos(phi),cos(theta_)*sin(phi)},
		{0,-sin(phi),cos(theta_)*cos(phi)}
	};
	double omega[3][1];
	for (int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			omega[i][j]=0;
			for(int k=0;k<3;k++)
			{
				omega[i][j]=omega[i][j]+w[i][k]*thetadot[k][j];
			}
			cout << omega[i][j]<<"\n ";
		}
	}
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			return omega[i][j];
			
		}
	}
}
