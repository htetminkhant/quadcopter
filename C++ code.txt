



#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iterator>
#include <matrix.h>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/iterator_range.hpp>
#include "engine.h"
#define  BUFSIZE 256
using namespace boost::numeric::ublas;

vector<double>time(vector<double>myvec);
vector<double>input();
matrix<double>rotation(vector<double>angle);
vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
double sum(vector<double>in);
vector<double> thrust(vector<double>in,double k);
vector<double> torques(vector<double>in,double L,double b,double k);
vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
vector<double>omega2thetadot(vector<double>omega,vector<double>angle);

int main()
{	

	Engine *ep;
	mxArray *X = NULL, *Xdot = NULL, *Theta = NULL, *Thetadot = NULL, *R_ = NULL; 
	char buffer[BUFSIZE+1];
    if (!(ep = engOpen(""))) 
	{
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return EXIT_FAILURE;
    }

	X = mxCreateDoubleMatrix(1, 3, mxREAL);
	Xdot = mxCreateDoubleMatrix(1, 3, mxREAL);
	Theta = mxCreateDoubleMatrix(1, 3, mxREAL);
	Thetadot = mxCreateDoubleMatrix(1, 3, mxREAL);
	R_ = mxCreateDoubleMatrix(3, 3, mxREAL);

	double i, g = 9.81 ,m = 0.5, L = 0.25, k = 3e-6, b = 1e-7, kd = 0.25, dt = 0.005, j;
	int temp = 0, deviation = 100;
	vector<double>timevec(2001), x(3), xdot(3,0), theta(3,0), thetadot(3), omega, a, omegadot, in(4,0);
	matrix<double>I(3,3), R(3,3);
	I <<= 5e-3, 0, 0,
		0, 5e-3, 0,
		0, 0, 10e-3;
	x <<= 0, 0, 10;
	thetadot = deg2rad(deviation,thetadot,temp);
	
	initialState()
	for ( j = 0; j < 1; j++)
	{
		i = timevec[j];
		in = input();

		omega = thetadot2omega(thetadot,theta);
		a = acceleration(in,theta,xdot,m,g,k,kd);
		omegadot = angular_acceleration(in,I,omega,L,b,k);
		omega = omega+dt*omegadot;
		thetadot = omega2thetadot(omega,theta);
		theta = theta+(dt*thetadot);
		xdot = xdot+dt*a;
		x = x+dt*xdot;
		//R = rotation(theta);
		visMatlab()
		
		std::copy(x.begin(), x.end(), mxGetPr(X));
		std::copy(xdot.begin(), xdot.end(), mxGetPr(Xdot));
		std::copy(theta.begin(), theta.end(), mxGetPr(Theta));
		std::copy(thetadot.begin(), thetadot.end(), mxGetPr(Thetadot));
		//std::copy(&R(0,0), &R(0,0)+3*3,mxGetPr(R_));
		
		engPutVariable(ep, "X", X);
		engPutVariable(ep, "Xdot", Xdot);
		engPutVariable(ep, "Theta", Theta);
		engPutVariable(ep, "Thetadot", Thetadot);
		//engPutVariable(ep, "r", R_);
		//engEvalString(ep, " R = r.';");

		//engEvalString(ep, "makehgtform('translate',x)");
		
		
		//std::cout<<R<<std::endl;
		//std::cout<<std::endl;
		

	}

	mxDestroyArray(X);
	mxDestroyArray(Xdot);
	mxDestroyArray(Theta);
	mxDestroyArray(Thetadot);
	mxDestroyArray(R_);
	engEvalString(ep, "close;");
	_getch();
	return 0;
}

void initialState(Engine *ep)
{
	double t = 0;
	engPutVariable(ep, "t", t);
	double dt = 0.05;
	engPutVariable(ep, "dt", dt);
	engEvalString(ep, "createPlot()");
//[model, thrusts] = quadcopter();
//setVisLimits()
}

void visMatlab(Enginae *ep, double x, double theta, double thetadot)
{
	engPutVariable(ep, "theta", theta);
	engPutVariable(ep, "thetadot", thetadot);
	engEvalString(ep, "t = t + dt");
//subplot(plots(1));
//move = makehgtform('translate',x);
//rotate = getRotation(theta);    
//set(model,'Matrix', move * rotate);
//visThrust(thrusts, move, rotate)
//[xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x);
//axis([xmin xmax ymin ymax zmin zmax]);
//drawnow;
//visVelocities(t, theta, thetadot)
}

vector<double>time(vector<double>myvec)
{
	double start_time = 0.0, end_time = 10.0, dt = 0.005;
	unsigned int N = 0;
	myvec[N] = start_time;
	while(start_time <= end_time)
	{  
		++N;
		start_time += dt;
		myvec[N] = start_time;
		return myvec;
	}
}
vector<double> deg2rad(int deviation, vector<double>thetadot, int temp)
{
	double pi = 3.142;
	for(int i = 0; i < 3; ++i)
		thetadot[i] = ((((2*deviation*((double)  rand()/ RAND_MAX)-deviation)) *pi) /180);
	return thetadot;
}
vector<double>input()
{
	vector<double>in(4);
	in <<= (700+10)*(740), 700*700, 740*(700+10), 700*700;
	return in;
}
double sum(vector<double>in)
{
	double suminput = 0;
	for (int i = 0; i < 4; i++)
		suminput += in[i];
	return suminput;
}
vector<double>thrust(vector<double>in, double k)
{
	vector<double>T(3);
	T <<= 0, 0, k*sum(in);
	return T;
}
vector<double>torques(vector<double>in, double L, double b, double k)
{
	vector<double>tau(3);
	tau <<= L * k * (in[0]-in[2]), L*k*(in[1]-in[3]), b*(in[0]-in[1]+in[2]-in[3]);
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

vector<double>acceleration(vector<double>in, vector<double>angle, vector<double>xdot_, double m, double g, double k, double kd)
{

	vector<double>gravity(3), T_, FD(3), a(3);
	matrix<double>R(3,3);
	gravity <<= 0, 0,-g;
	R = rotation(angle);
	T_ = prod(R,thrust(in,k));
	FD = -kd*xdot_;
	a = gravity+(1/m)*T_+FD;
	return a;
}
vector<double>angular_acceleration(vector<double>in, matrix<double> I, vector<double>omega, double L, double b, double k)
{	
	vector<double>tau(3), omegadot(3);
	tau = torques(in, L, b, k);
	omegadot = prod(-I, tau- ( omega, prod(I, omega)));
	return omegadot;
}
vector<double>thetadot2omega(vector<double>thetadot, vector<double>angle)
{
	double phi = angle[0], theta_ = angle[1], psi = angle[2];
	vector<double>omega;
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
		0, cos(phi), cos(theta_)*sin(phi),
		0, -sin(phi), cos(theta_)*cos(phi);
	omega = prod(w,thetadot );
	return omega;
}
vector<double>omega2thetadot(vector<double>omega, vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	vector<double>thetadot;
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
		0, cos(phi), cos(theta_)*sin(phi),
		0, -sin(phi), cos(theta_)*cos(phi);
	thetadot = prod(-w,omega);
	return thetadot;
}