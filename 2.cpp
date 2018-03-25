#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <fstream>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
void createvector( int &N,vector<double>&myvec,double start_time,double end_time,double dt);
vector<double> acceleration(double inputs,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
vector<double> angular_acceleration(int j,vector<double>omega,matrix<double> I,double L,double b,double k);
vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
static const double start_time = 
int main()
{	
	srand( time(0));
	double start_time=0.0,end_time=4.0, dt=0.005,randnumber,pi=3.142;
	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25;
	int N=0,deviation =100,temp,j;
	vector<double>timevector,x(3),xdot(3,0),theta(3,0),thetadot(3);
	createvector(N,timevector,start_time,end_time,dt);
	x <<= 0, 0, 10;
	for(int i=0;i<3;++i)
	{
		randnumber=(double)  rand()/ RAND_MAX;
		temp=(2*deviation*randnumber-deviation);
		thetadot[i]=((temp*pi)/180);
	};
	matrix<double> I(3,3);
	I <<= 5e-3, 0, 0,
         0, 5e-3, 0,
         0, 0, 10e-3;
	std::ofstream outfile;
	outfile.open("C:/Users/Htet Min Khant/Documents/superadvisor/mysimulation/mysimulation.csv");
	
	for ( j=0;j<timevector.size();j++)
	{
		i=timevector[j];
		vector<double>omega=thetadot2omega(thetadot,theta);
		vector<double>a=acceleration(i,theta,xdot,m,g,k,kd);
		vector<double>omegadot=angular_acceleration(j,omega,I,L,b,k);
		omega=omega+dt*omegadot;
		vector<double>thetadot=omega2thetadot(omega,theta);
		theta=theta+dt*thetadot;
		xdot=xdot+dt*a;
		x=x+dt*xdot;
		for (int h=0;h<3;h++)
		outfile<<x[h]<<"\n";
		outfile<<"\n";
	}
	outfile.close();
	std::cout<<std::endl;
	_getch();
	return 0;
}


void createvector( int &N,vector<double>&myvec,double start_time,double end_time,double dt)
{
	N=0;
	myvec[N]=(start_time);
	while(start_time<=end_time)
	{  
		++N;
		start_time+=dt;
		myvec[N]=(start_time);
		
	}
}

vector<double>thetadot2omega(vector<double>thetadot,vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
        0, cos(phi), cos(theta_)*sin(phi),
        0, -sin(phi), cos(theta_)*cos(phi);
	vector<double>omega=prod(w,thetadot );
	return omega;
}

vector<double> acceleration(double inputs,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
{
	double phi = angle[2],theta_ = angle[1], psi = angle[0];
	matrix<double>R(3,3);
	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
        cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
        -sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
	vector<double>timevector;
	double start_time=0.0,end_time=10.0, dt=0.005;
	int N=0;
	createvector(N,timevector,start_time,end_time,dt);
	double sumtime=0;
	for(int i=0;i<N;i++)
	{   
		sumtime+=timevector[i];
	}
	vector<double>T(3);
	T <<= 0, 0, k*12e+6;
	vector<double>gravity(3);
	gravity <<= 0, 0,-g;
	vector<double>T_=prod(R,T);
	vector<double>FD(3);
	vector<double>T2(3);
	vector<double>a(3);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<1;j++)
		{
			FD[i]=-kd*xdot_[i];
			T2[i]=(1/m)*T_[i];
			a[i]=gravity[i]+T2[i]+FD[i];
		}
	};
	return a;
}
vector<double> angular_acceleration(int j,vector<double>omega,matrix<double> I,double L,double b,double k)
{	
	
	double start_time=0.0,end_time=10.0, dt=0.005;
	int j_=j,N=0;
	vector<double>timevector;
	createvector(N,timevector,start_time,end_time,dt);
	vector<double>tau(3),omega_(3),omega2(3);
	tau <<= L*k*((3e+6)-(3e+6)), L*k*((3e+6)-(3e+6)),b*((3e+6)-(3e+6)+(3e+6)-(3e+6));
	omega_=prod(I,omega);
	omega2=(omega,omega_);
	matrix<double>invI=-I;
	vector<double>omegadot2=tau-omega2;
	vector<double>omegadot=prod(invI,omegadot2);
	return omegadot;
}

vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
        0, cos(phi), cos(theta_)*sin(phi),
        0, -sin(phi), cos(theta_)*cos(phi);
	matrix<double>invw=-w;
	vector<double>thetadot=prod(invw,omega);
	return thetadot;
}
