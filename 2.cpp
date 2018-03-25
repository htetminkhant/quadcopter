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
vector<double> input();
matrix<double>rotation(vector<double>angle);
vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
double sum(vector<double>in);
vector<double> thrust(vector<double>in,double k);
vector<double> torques(vector<double>in,double L,double b,double k);
vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
static const double start_time=0.0,end_time=4.0, dt=0.005;
double randnumber,pi=3.142;
int main()
{	
	srand( time(0));
	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25;
	int N=0,deviation =100,temp=0,j;
	vector<double>timevector,x(3),xdot(3,0),theta(3,0),thetadot(3),omega,a,omegadot,in(4,0);
	matrix<double> I(3,3);
	I <<= 5e-3, 0, 0,
         0, 5e-3, 0,
         0, 0, 10e-3;
	createvector(N,timevector,start_time,end_time,dt);
	x <<= 0, 0, 10;
	thetadot=deg2rad(deviation,thetadot,temp);
	std::ofstream outfile;
	outfile.open("C:/Users/Htet Min Khant/Documents/superadvisor/mysimulation/mysimulation.csv");
	for ( j=0;j<timevector.size();j++)
	{
		i=timevector[j];
		in=input();
		omega=thetadot2omega(thetadot,theta);
		a=acceleration(in,theta,xdot,m,g,k,kd);
		omegadot=angular_acceleration(in,I,omega,L,b,k);
		omega=omega+dt*omegadot;
		thetadot=omega2thetadot(omega,theta);
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
vector<double> deg2rad(int deviation,vector<double>thetadot,int temp)
{
	for(int i=0;i<3;++i)
	thetadot[i]=((((2*deviation*((double)  rand()/ RAND_MAX)-deviation))*pi)/180);
	return thetadot;
}
vector<double> input()
{
	vector<double> in(4,0);
	in<<=700+150,700,700+150,700;
	return in;
}
double sum(vector<double>in)
{
	double suminput=0;
	for(int i=0;i<4;i++)
		suminput+=in[i];
	return suminput;
}
vector<double> thrust(vector<double>in,double k)
{
	vector<double>T(3);
	T <<= 0, 0, k*sum(in);
	return T;
}
vector<double> torques(vector<double>in,double L,double b,double k)
{
	vector<double>tau(3);
	tau <<= L*k*(in[1]-in[3]), L*k*(in[2]-in[4]),b*(in[1]-in[2]+in[3]-in[4]);
	return tau;
}

matrix<double>rotation(vector<double>angle)
{
	double phi = angle[2],theta_ = angle[1], psi = angle[0];
	matrix<double>R(3,3);
	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
        cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
        -sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
	return R;
}

vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
{
	
	vector<double>gravity(3),T_,FD(3),T2(3),a(3);
	gravity <<= 0, 0,-g;
	matrix<double>R(3,3);
	R=rotation(angle);
	T_=prod(R,thrust(in,k));
	for(int i=0;i<3;i++)
	{
			FD[i]=-kd*xdot_[i];
			a[i]=gravity[i]+(1/m)*T_[i]+FD[i];
	};
	return a;
}
vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k)
{	
	vector<double>tau(3);
	tau=torques(in,L,b,k);
	vector<double>omegadot=prod(-I,tau-(omega,prod(I,omega)));
	return omegadot;
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
vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
        0, cos(phi), cos(theta_)*sin(phi),
        0, -sin(phi), cos(theta_)*cos(phi);
	vector<double>thetadot=prod(-w,omega);
	return thetadot;
}
