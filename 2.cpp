#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
using namespace std;
void createvector(vector<double>&myvec);
void thetadot2omega(vector<double>thetadot,vector<double>theta);
int main()
{	
	srand( time(0));
	vector<double>myvec2;
	createvector(myvec2);
	int N=0;
	for(int i=0;i<myvec2.size();i++)
	{   
		++N;
	}
	vector<double>x(2,0);
	x.push_back(10);
	vector<double>xdot(3,0);
	vector<double>theta(3,0);
	double randnumber;
	int deviation =100,temp;
	double pi=3.142;
	vector<double>thetadot;
	for(int i=0;i<3;++i)
	{
		randnumber=(double)  rand()/ RAND_MAX;
		temp=(2*deviation*randnumber-deviation);
		thetadot.push_back((temp*pi)/180);
		cout<<thetadot[i];
	};
	
	thetadot2omega(thetadot,theta);
	cout<<endl;
	_getch();
	return 0;
}
void createvector(vector<double>& myvec)
{
	double start_time=0.0,end_time=10.0, dt=0.005;
	myvec.push_back(start_time);
	while(start_time<=end_time)
	{  
		start_time+=dt;
		myvec.push_back(start_time);

	}

}
void thetadot2omega(vector<double>thetadot,vector<double>theta)
{
	double phi = theta[0],theta_ = theta[1], psi = theta[2];
	vector<vector<double>>w(3,vector<double>(3));
	
	w[0][0]=1;
	w[0][1]=0;
	w[0][2]=-sin(theta_);
	w[1][0]=0;
	w[1][1]=cos(phi);
	w[1][2]=cos(theta_)*sin(phi);
	w[2][0]=0;
	w[2][1]=-sin(phi);
	w[2][2]=cos(theta_)*cos(phi);
	vector<vector<double>>omega(3,vector<double>(3));
	for (int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			omega[i][j]=0;
			for(int k=0;k<3;k++)
			{
				omega[i][j]=omega[i][j]+w[i][k]*thetadot[k];
			}
			cout << "omega "<<omega[i][j]<<"\n ";
		}
	}
}
