#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
using namespace std;
void acceleration(double inputs,vector<double>theta,vector<double>xdot_,double m,double g,double k,double kd,vector<vector<double>>a);
void createvector(vector<double>&myvec);
vector<vector<double>> thetadot2omega(vector<double>thetadot,vector<double>angle);
int main()
{	
	srand( time(0));
	vector<double>timevector;
	createvector(timevector);
	int N=0;
	for(int i=0;i<timevector.size();i++)
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
	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25;
	vector<vector<double>>I(3,vector<double>(3));
	
	vector<vector<double>>a(3,vector<double>(1));
	I[0][0]=5e-3;
	I[0][1]=0;
	I[0][2]=0;
	I[1][0]=0;
	I[1][1]=5e-3;
	I[1][2]=0;
	I[2][0]=0;
	I[2][1]=0;
	I[2][2]=10e-3;
	
	for (int j=0;j<N;j++)
	{
		i=timevector[j];
		
		acceleration(i,theta,xdot,m,g,k,kd,a);
		
	}
	cout<<endl;
	vector<vector<double>>omega=thetadot2omega(thetadot,theta);
	for (int i=0;i<3;i++)
		for (int j=0;j<1;j++)
			cout<<omega[i][j]<<"  ";
	
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

vector<vector<double>>thetadot2omega(vector<double>thetadot,vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	vector<vector<double>>w(3,vector<double>(3));
	vector<vector<double>>omega(3,vector<double>(1));
	w[0][0]=1;
	w[0][1]=0;
	w[0][2]=-sin(theta_);
	w[1][0]=0;
	w[1][1]=cos(phi);
	w[1][2]=cos(theta_)*sin(phi);
	w[2][0]=0;
	w[2][1]=-sin(phi);
	w[2][2]=cos(theta_)*cos(phi);
	
	for (int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			omega[i][j]=0;
			for(int k=0;k<3;k++)
			{
				omega[i][j]=omega[i][j]+w[i][k]*thetadot[k];
			}
			//cout << "omega "<<omega[i][j]<<"\n ";
		}
	}
	return omega;
}

void acceleration(double inputs,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd,vector<vector<double>>a)
{
	double phi = angle[2],theta_ = angle[1], psi = angle[0];
	vector<vector<double>>R(3,vector<double>(3));
	R[0][0]=cos(phi)*cos(theta_);
	R[1][0]=cos(theta_)*sin(phi);
	R[2][0]=-sin(theta_);
	R[0][1]=cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi);
	R[1][1]=cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi);
	R[2][1]=cos(theta_)*sin(phi);
	R[0][2]=sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_);
	R[1][2]=cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi);
	R[2][2]=cos(theta_)*cos(phi);
	vector<double>timevector;
	createvector(timevector);
	double sumtime=0;
	for(int i=0;i<2001;i++)
	{   
		sumtime+=timevector[i];
	}
	vector<vector<double>>T(3,vector<double>(1));
	T[0][0]=0;
	T[1][0]=0;
	T[2][0]=k*sumtime;
	vector<vector<double>>gravity(3,vector<double>(1));
	gravity[0][0]=0;
	gravity[1][0]=0;
	gravity[2][0]=-g;
	vector<vector<double>>T_(3,vector<double>(1));
	for (int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			T_[i][j]=0;
			for(int k=0;k<3;k++)
			{
				T_[i][j]=T_[i][j]+R[i][k]*T[k][j];
			}
		}
	}
	vector<vector<double>>FD(3,vector<double>(1));
	vector<vector<double>>T2(3,vector<double>(1));
	
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			FD[i][j]=-kd*xdot_[i];
			T2[i][j]=(1/m)*T_[i][j];
			a[i][j]=gravity[i][j]+T2[i][j]+FD[i][j];
		}
	};
}
