//
//#include "stdio.h"
//#include <stdlib.h>
//#include <conio.h>
//#include <iostream>
//#include <vector>
//#include <time.h>
//#include <math.h>
//#include <fstream>
//#include "engine.h"
//#include <boost/numeric/ublas/assignment.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#define  BUFSIZE 256
//using namespace boost::numeric::ublas;
//void createvector(vector<double>&myvec);
//vector<double> input();
//matrix<double>rotation(vector<double>angle);
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
//double sum(vector<double>in);
//vector<double> thrust(vector<double>in,double k);
//vector<double> torques(vector<double>in,double L,double b,double k);
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
//vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
//int main()
//{	
//	Engine *ep;
//	srand( (unsigned int)time(0));
//	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25 ,dt=0.005,j;
//	mxArray *T = NULL;
//	if (!(ep = engOpen(""))) {
//		std::cout<< "\nCan't start MATLAB engine\n";
//		return EXIT_FAILURE;
//	}
//	//vector<double>x(3);
//	
//	int temp=0,deviation =100;
//	vector<double>timevector(2001),x(3),xdot(3,0),theta(3,0),thetadot(3),omega,a,omegadot,in(4,0);
//	
//	matrix<double> I(3,3);
//	I <<= 5e-3, 0, 0,
//		0, 5e-3, 0,
//		0, 0, 10e-3;
//	x <<= 0, 0, 10;
//	createvector(timevector);
//	thetadot=deg2rad(deviation,thetadot,temp);
//	//std::ofstream outfile;
//	//outfile.open("C:/Users/Htet Min Khant/Documents/superadvisor/mysimulation/mysimulation2.csv");
//	
//	for ( j=0;j<timevector.size();j++)
//	{
//		i=timevector[j];
//		in=input();
//		omega=thetadot2omega(thetadot,theta);
//		a=acceleration(in,theta,xdot,m,g,k,kd);
//		omegadot=angular_acceleration(in,I,omega,L,b,k);
//		omega=omega+dt*omegadot;
//		thetadot=omega2thetadot(omega,theta);
//		theta=theta+(dt*thetadot);
//		xdot=xdot+dt*a;
//		x=x+dt*xdot;
//		//for (int i=0;i<3;i++)
//			//outfile<<x[h]<<"\n";
//		//outfile<<"\n";
//		//std::cout<<x[i]<<std::endl;
//		//std::cout<<std::endl;
//
//		T = mxCreateDoubleMatrix(3,1,mxREAL);
//		memcpy((void *)mxGetPr(T), (void *)time, sizeof(time));
//		
//		engPutVariable(ep, "T", T);
//		engEvalString(ep, "D = T");
//		engEvalString(ep, "plot(T,D);");
//		engEvalString(ep, "title('Position vs. Time for a falling object');");
//		engEvalString(ep, "xlabel('Time (seconds)');");
//		engEvalString(ep, "ylabel('Position (meters)');");
//		std::cout<<"Hit return to continue\n\n";
//		fgetc(stdin);
//		std::cout<<"Done for Part I.\n";
//		mxDestroyArray(T);
//		engClose(ep);
//	
//	
//	}
//	//outfile.close();
//	//std::cout<<timevector.size()<<std::endl;
//	
//	return EXIT_SUCCESS;
//}
//
//
//
//void createvector(vector<double>& myvec)
//{
//	double start_time=0.0,end_time=10.0, dt=0.005;
//	unsigned int N=0;
//	myvec[N]=start_time;
//	while(start_time<=end_time)
//	{  
//		++N;
//		start_time+=dt;
//		myvec[N]=start_time;
//	}
//}
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp)
//{
//	double pi=3.142;
//	for(int i=0;i<3;++i)
//		thetadot[i]=((((2*deviation*((double)  rand()/ RAND_MAX)-deviation))*pi)/180);
//	return thetadot;
//}
//vector<double> input()
//{
//	vector<double> in(4);
//	in<<=(700+10)*(740),700*700,740*(700+10),700*700;
//	return in;
//}
//double sum(vector<double>in)
//{
//	double suminput=0;
//	for(int i=0;i<4;i++)
//		suminput+=in[i];
//	return suminput;
//}
//vector<double> thrust(vector<double>in,double k)
//{
//	vector<double>T(3);
//	T <<= 0, 0, k*sum(in);
//	return T;
//}
//vector<double> torques(vector<double>in,double L,double b,double k)
//{
//	vector<double>tau(3);
//	tau <<=L*k*(in[0]-in[2]), L*k*(in[1]-in[3]),b*(in[0]-in[1]+in[2]-in[3]);
//	return tau;
//}
//
//matrix<double>rotation(vector<double>angle)
//{
//	double phi = angle[2],theta_ = angle[1], psi = angle[0];
//	matrix<double>R(3,3);
//	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
//		cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
//		-sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
//	return R;
//}
//
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
//{
//
//	vector<double>gravity(3),T_,FD(3),a(3);
//	matrix<double>R(3,3);
//	gravity <<= 0, 0,-g;
//	R=rotation(angle);
//	T_=prod(R,thrust(in,k));
//	FD=-kd*xdot_;
//	a=gravity+(1/m)*T_+FD;
//	return a;
//}
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k)
//{	
//	vector<double>tau(3),omegadot(3);
//	tau=torques(in,L,b,k);
//	omegadot=prod(-I,tau-(omega,prod(I,omega)));
//	return omegadot;
//}
//vector<double>thetadot2omega(vector<double>thetadot,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>omega;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	omega=prod(w,thetadot );
//	return omega;
//}
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>thetadot;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	thetadot=prod(-w,omega);
//	return thetadot;
//}
//
//
//
//







//
//
//
//#include "stdio.h"
//#include <stdlib.h>
//#include <conio.h>
//#include <algorithm>
//#include <iostream>
//#include <vector>
//#include <time.h>
//#include <math.h>
//#include <fstream>
//#include <iterator>
//#include <boost/numeric/ublas/assignment.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/range/iterator_range.hpp>
//#include "engine.h"
//
//using namespace boost::numeric::ublas;
//
//void createvector(vector<double>&myvec);
//vector<double> input();
//matrix<double>rotation(vector<double>angle);
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
//double sum(vector<double>in);
//vector<double> thrust(vector<double>in,double k);
//vector<double> torques(vector<double>in,double L,double b,double k);
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
//vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
////int vistulation();
//#define  BUFSIZE 256
//int main()
//{	
//	Engine *ep;
//	mxArray *R_ = NULL;
//	mxArray *X = NULL;
//	double *xp,*rp;
//	char buffer[BUFSIZE+1];
//    if (!(ep = engOpen(""))) {
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//    }
//	//engEvalString(ep, "plots = [subplot(3, 2, 1:4), subplot(3, 2, 5), subplot(3, 2, 6)];");
//	//engEvalString(ep, "subplot(plots(1));");
//	//engEvalString(ep, "axis([-10 30 -20 20 5 30]);");
//	//engEvalString(ep, "zlabel('Height');");
//	//engEvalString(ep, "title('Quadcopter Flight Simulation'); ");
//	srand( (unsigned int)time(0));
//	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25 ,dt=0.005;
//	unsigned int j;
//	int temp=0,deviation =100;
//	vector<double>timevector(2001),x(3),xdot(3,0),theta(3,0),thetadot(3),omega,a,omegadot,in(4,0);
//	matrix<double> I(3,3),R(3,3);
//	I <<= 5e-3, 0, 0,
//		0, 5e-3, 0,
//		0, 0, 10e-3;
//	x <<= 0, 0, 10;
//	createvector(timevector);
//	thetadot=deg2rad(deviation,thetadot,temp);
//	std::ofstream outfile;
//	//vistulation();
//	//outfile.open("C:/Users/Htet Min Khant/Documents/superadvisor/mysimulation/mysimulation3.csv");
//	//for ( j=0;j<timevector.size();j++)
//	for ( j=0;j<1;j++)
//	{
//		i=timevector[j];
//		in=input();
//		omega=thetadot2omega(thetadot,theta);
//		a=acceleration(in,theta,xdot,m,g,k,kd);
//		omegadot=angular_acceleration(in,I,omega,L,b,k);
//		omega=omega+dt*omegadot;
//		thetadot=omega2thetadot(omega,theta);
//		theta=theta+(dt*thetadot);
//		xdot=xdot+dt*a;
//		x=x+dt*xdot;
//		R=rotation(theta);
//		X = mxCreateDoubleMatrix(1, 3, mxREAL);
//		R_ = mxCreateDoubleMatrix(3, 3, mxREAL);
//		std::copy(x.begin(), x.end(), mxGetPr(X));
//		std::copy(R.begin1(), 3*3*sizeof(R(0,0)), mxGetPr(R_));
//		engPutVariable(ep, "r", R_);
//		//for ( int i=0; i<3; i++ )
//		//{
//		//	memcpy((void *)mxGetPr(X)+i*sizeof(mxREAL),&x[i], sizeof(x));
//		//}
//		/*printf("%f, %f, %f", x[0],x[1],x[2]);*/
//		printf("%f, %f, %f", R(0,0),R(0,1),R(0,2));
//		printf("%f, %f, %f", R(1,0),R(1,1),R(1,2));
//		printf("%f, %f, %f", R(2,0),R(2,1),R(2,2));
//		//engEvalString(ep, " move = makehgtform('translate',x);");
//		/*for (int h=0;h<3;h++)
//			outfile<<x[h]<<"\n";
//		outfile<<"\n";*/
//
//	}
//	outfile.close();
//	std::cout<<timevector.size()<<std::endl;
//	std::cout<<std::endl;
//	//mxDestroyArray(X);
//	//mxDestroyArray(R_);
//	//engClose(ep);
//	_getch();
//	return 0;
//}
//
//void createvector(vector<double>& myvec)
//{
//	double start_time=0.0,end_time=10.0, dt=0.005;
//	unsigned int N=0;
//	myvec[N]=start_time;
//	while(start_time<=end_time)
//	{  
//		++N;
//		start_time+=dt;
//		myvec[N]=start_time;
//	}
//}
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp)
//{
//	double pi=3.142;
//	for(int i=0;i<3;++i)
//		thetadot[i]=((((2*deviation*((double)  rand()/ RAND_MAX)-deviation))*pi)/180);
//	return thetadot;
//}
//vector<double> input()
//{
//	vector<double> in(4);
//	in<<=(700+10)*(740),700*700,740*(700+10),700*700;
//	return in;
//}
//double sum(vector<double>in)
//{
//	double suminput=0;
//	for(int i=0;i<4;i++)
//		suminput+=in[i];
//	return suminput;
//}
//vector<double> thrust(vector<double>in,double k)
//{
//	vector<double>T(3);
//	T <<= 0, 0, k*sum(in);
//	return T;
//}
//vector<double> torques(vector<double>in,double L,double b,double k)
//{
//	vector<double>tau(3);
//	tau <<=L*k*(in[0]-in[2]), L*k*(in[1]-in[3]),b*(in[0]-in[1]+in[2]-in[3]);
//	return tau;
//}
//
//matrix<double>rotation(vector<double>angle)
//{
//	double phi = angle[2],theta_ = angle[1], psi = angle[0];
//	matrix<double>R(3,3);
//	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
//		cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
//		-sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
//	return R;
//}
//
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
//{
//
//	vector<double>gravity(3),T_,FD(3),a(3);
//	matrix<double>R(3,3);
//	gravity <<= 0, 0,-g;
//	R=rotation(angle);
//	T_=prod(R,thrust(in,k));
//	FD=-kd*xdot_;
//	a=gravity+(1/m)*T_+FD;
//	return a;
//}
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k)
//{	
//	vector<double>tau(3),omegadot(3);
//	tau=torques(in,L,b,k);
//	//std::cout << inner_prod(omega,prod(I,omega)) << std::endl;
//	//std::cout << std::endl;
//	omegadot=prod(-I,(tau-element_prod(omega,prod(I,omega))));
//	return omegadot;
//}
//vector<double>thetadot2omega(vector<double>thetadot,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>omega;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	omega=prod(w,thetadot );
//	return omega;
//}
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>thetadot;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	thetadot=prod(-w,omega);
//	return thetadot;
//}
/*int vistulation()
{
	Engine *ep;
	char buffer[BUFSIZE+1];
    if (!(ep = engOpen(""))) {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return EXIT_FAILURE;
    }
	engEvalString(ep, "plots = [subplot(3, 2, 1:4), subplot(3, 2, 5), subplot(3, 2, 6)];");
	engEvalString(ep, "subplot(plots(1));");
	engEvalString(ep, "axis([-10 30 -20 20 5 30]);");
	engEvalString(ep, "zlabel('Height');");
	engEvalString(ep, "title('Quadcopter Flight Simulation'); ");
}*/




//
//
//#include "stdio.h"
//#include <stdlib.h>
//#include <conio.h>
//#include <algorithm>
//#include <iostream>
//#include <vector>
//#include <time.h>
//#include <math.h>
//#include <fstream>
//#include <iterator>
//#include <boost/numeric/ublas/assignment.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/range/iterator_range.hpp>
//#include "engine.h"
//
//using namespace boost::numeric::ublas;
//
//void createvector(vector<double>&myvec);
//vector<double> input();
//matrix<double>rotation(vector<double>angle);
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
//double sum(vector<double>in);
//vector<double> thrust(vector<double>in,double k);
//vector<double> torques(vector<double>in,double L,double b,double k);
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
//vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
//#define  BUFSIZE 256
//int main()
//{	
//	Engine *ep;
//	mxArray *R_ = NULL;
//	mxArray *X = NULL,*Theta=NULL,*Xdot=NULL,*Thetadot=NULL;
//	double *xp,*rp;
//	char buffer[BUFSIZE+1];
//    if (!(ep = engOpen(""))) {
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//    }
//	srand( (unsigned int)time(0));
//	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25 ,dt=0.005;
//	unsigned int j;
//	int temp=0,deviation =100;
//	vector<double>timevector(2001),x(3),xdot(3,0),theta(3,0),thetadot(3),omega,a,omegadot,in(4,0);
//	matrix<double> I(3,3),R(3,3);
//	I <<= 5e-3, 0, 0,
//		0, 5e-3, 0,
//		0, 0, 10e-3;
//	x <<= 0, 0, 10;
//	createvector(timevector);
//	thetadot=deg2rad(deviation,thetadot,temp);
//	X = mxCreateDoubleMatrix(1, 3, mxREAL);
//	Theta = mxCreateDoubleMatrix(1, 3, mxREAL);
//	Xdot = mxCreateDoubleMatrix(1, 3, mxREAL);
//	Thetadot = mxCreateDoubleMatrix(1, 3, mxREAL);
//	R_ = mxCreateDoubleMatrix(3, 3, mxREAL);
//	std::ofstream outfile;
//	//vistulation();
//	//outfile.open("C:/Users/Htet Min Khant/Documents/superadvisor/mysimulation/mysimulation3.csv");
//	//for ( j=0;j<timevector.size();j++)
//	for ( j=0;j<1;j++)
//	{
//		i=timevector[j];
//		in=input();
//		omega=thetadot2omega(thetadot,theta);
//		a=acceleration(in,theta,xdot,m,g,k,kd);
//		omegadot=angular_acceleration(in,I,omega,L,b,k);
//		omega=omega+dt*omegadot;
//		thetadot=omega2thetadot(omega,theta);
//		theta=theta+(dt*thetadot);
//		xdot=xdot+dt*a;
//		x=x+dt*xdot;
//		R=rotation(theta);
//		std::copy(x.begin(), x.end(), mxGetPr(X));
//		engPutVariable(ep, "x", X);
//		std::copy(theta.begin(), theta.end(), mxGetPr(Theta));
//		engPutVariable(ep, "theta", Theta);
//		std::copy(xdot.begin(), xdot.end(), mxGetPr(Xdot));
//		engPutVariable(ep, "xdot", Xdot);
//		std::copy(thetadot.begin(), thetadot.end(), mxGetPr(Thetadot));
//		engPutVariable(ep, "thetadot", Thetadot);
//		std::copy(&R(0,0), &R(0,0)+3*3,mxGetPr(R_));
//		engPutVariable(ep, "r", R_);
//		engEvalString(ep, " R = r.';");
//		//std::copy(R.begin1(), 3*3*sizeof(R(0,0)), mxGetPr(R_));
//		//boost::copy(R, mxGetPr(R_));
//		//for ( int i=0; i<3; i++ )
//		//{
//		//	memcpy((void *)mxGetPr(X)+i*sizeof(mxREAL),&x[i], sizeof(x));
//		//}
//		/*printf("%f, %f, %f", x[0],x[1],x[2]);*/
//		printf("%f, %f, %f", R(0,0),R(0,1),R(0,2));
//		printf("%f, %f, %f", R(1,0),R(1,1),R(1,2));
//		printf("%f, %f, %f", R(2,0),R(2,1),R(2,2));
//		//engEvalString(ep, " move = makehgtform('translate',x);");
//		/*for (int h=0;h<3;h++)
//			outfile<<x[h]<<"\n";
//		outfile<<"\n";*/
//
//	}
//	outfile.close();
//	std::cout<<timevector.size()<<std::endl;
//	std::cout<<std::endl;
//	//mxDestroyArray(X);
//	//mxDestroyArray(R_);
//	//engClose(ep);
//	_getch();
//	return 0;
//}
//
//void createvector(vector<double>& myvec)
//{
//	double start_time=0.0,end_time=10.0, dt=0.005;
//	unsigned int N=0;
//	myvec[N]=start_time;
//	while(start_time<=end_time)
//	{  
//		++N;
//		start_time+=dt;
//		myvec[N]=start_time;
//	}
//}
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp)
//{
//	double pi=3.142;
//	for(int i=0;i<3;++i)
//		thetadot[i]=((((2*deviation*((double)  rand()/ RAND_MAX)-deviation))*pi)/180);
//	return thetadot;
//}
//vector<double> input()
//{
//	vector<double> in(4);
//	in<<=(700+10)*(740),700*700,740*(700+10),700*700;
//	return in;
//}
//double sum(vector<double>in)
//{
//	double suminput=0;
//	for(int i=0;i<4;i++)
//		suminput+=in[i];
//	return suminput;
//}
//vector<double> thrust(vector<double>in,double k)
//{
//	vector<double>T(3);
//	T <<= 0, 0, k*sum(in);
//	return T;
//}
//vector<double> torques(vector<double>in,double L,double b,double k)
//{
//	vector<double>tau(3);
//	tau <<=L*k*(in[0]-in[2]), L*k*(in[1]-in[3]),b*(in[0]-in[1]+in[2]-in[3]);
//	return tau;
//}
//
//matrix<double>rotation(vector<double>angle)
//{
//	double phi = angle[2],theta_ = angle[1], psi = angle[0];
//	matrix<double>R(3,3);
//	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
//		cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
//		-sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
//	return R;
//}
//
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
//{
//
//	vector<double>gravity(3),T_,FD(3),a(3);
//	matrix<double>R(3,3);
//	gravity <<= 0, 0,-g;
//	R=rotation(angle);
//	T_=prod(R,thrust(in,k));
//	FD=-kd*xdot_;
//	a=gravity+(1/m)*T_+FD;
//	return a;
//}
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k)
//{	
//	vector<double>tau(3),omegadot(3);
//	tau=torques(in,L,b,k);
//	//std::cout << inner_prod(omega,prod(I,omega)) << std::endl;
//	//std::cout << std::endl;
//	omegadot=prod(-I,(tau-element_prod(omega,prod(I,omega))));
//	return omegadot;
//}
//vector<double>thetadot2omega(vector<double>thetadot,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>omega;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	omega=prod(w,thetadot );
//	return omega;
//}
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>thetadot;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	thetadot=prod(-w,omega);
//	return thetadot;
//}
///*int vistulation()
//{
//	Engine *ep;
//	char buffer[BUFSIZE+1];
//    if (!(ep = engOpen(""))) {
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//    }
//	engEvalString(ep, "plots = [subplot(3, 2, 1:4), subplot(3, 2, 5), subplot(3, 2, 6)];");
//	engEvalString(ep, "subplot(plots(1));");
//	engEvalString(ep, "axis([-10 30 -20 20 5 30]);");
//	engEvalString(ep, "zlabel('Height');");
//	engEvalString(ep, "title('Quadcopter Flight Simulation'); ");
//}*/




//=======================================================
//                MY PROGRAM
////=======================================================
//
//#include "stdio.h"
//#include <stdlib.h>
//#include <conio.h>
//#include <iostream>
//#include <vector>
//#include <time.h>
//#include <math.h>
//#include <fstream>
//#include <iterator>
//#include <matrix.h>
//#include <boost/numeric/ublas/assignment.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/range/iterator_range.hpp>
//#include "engine.h"
//#define  BUFSIZE 256
//using namespace boost::numeric::ublas;
//
//vector<double>time(vector<double>myvec);
//vector<double>input();
//matrix<double>rotation(vector<double>angle);
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
//double sum(vector<double>in);
//vector<double> thrust(vector<double>in,double k);
//vector<double> torques(vector<double>in,double L,double b,double k);
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
//vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
//
//int main()
//{	
//
//	Engine *ep;
//	mxArray *X = NULL, *Xdot = NULL, *Theta = NULL, *Thetadot = NULL, *R_ = NULL; 
//	char buffer[BUFSIZE+1];
//    if (!(ep = engOpen(""))) 
//	{
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//    }
//
//	X = mxCreateDoubleMatrix(1, 3, mxREAL);
//	Xdot = mxCreateDoubleMatrix(1, 3, mxREAL);
//	Theta = mxCreateDoubleMatrix(1, 3, mxREAL);
//	Thetadot = mxCreateDoubleMatrix(1, 3, mxREAL);
//	R_ = mxCreateDoubleMatrix(3, 3, mxREAL);
//
//	double i, g = 9.81 ,m = 0.5, L = 0.25, k = 3e-6, b = 1e-7, kd = 0.25, dt = 0.005, j;
//	int temp = 0, deviation = 100;
//	vector<double>timevec(2001), x(3), xdot(3,0), theta(3,0), thetadot(3), omega, a, omegadot, in(4,0);
//	matrix<double>I(3,3), R(3,3);
//	I <<= 5e-3, 0, 0,
//		0, 5e-3, 0,
//		0, 0, 10e-3;
//	x <<= 0, 0, 10;
//	thetadot = deg2rad(deviation,thetadot,temp);
//	
//
//	for ( j = 0; j < 1; j++)
//	{
//		i = timevec[j];
//		in = input();
//
//		omega = thetadot2omega(thetadot,theta);
//		a = acceleration(in,theta,xdot,m,g,k,kd);
//		omegadot = angular_acceleration(in,I,omega,L,b,k);
//		omega = omega+dt*omegadot;
//		thetadot = omega2thetadot(omega,theta);
//		theta = theta+(dt*thetadot);
//		xdot = xdot+dt*a;
//		x = x+dt*xdot;
//		R = rotation(theta);
//		
//		
//		std::copy(x.begin(), x.end(), mxGetPr(X));
//		std::copy(xdot.begin(), xdot.end(), mxGetPr(Xdot));
//		std::copy(theta.begin(), theta.end(), mxGetPr(Theta));
//		std::copy(thetadot.begin(), thetadot.end(), mxGetPr(Thetadot));
//		std::copy(&R(0,0), &R(0,0)+3*3,mxGetPr(R_));
//		
//		engPutVariable(ep, "X", X);
//		engPutVariable(ep, "Xdot", Xdot);
//		engPutVariable(ep, "Theta", Theta);
//		engPutVariable(ep, "Thetadot", Thetadot);
//		engPutVariable(ep, "r", R_);
//		engEvalString(ep, " R = r.';");
//
//		/*engEvalString(matlab, "makehgtform('translate',x)");
//		
//		
//		std::cout<<R<<std::endl;
//		std::cout<<std::endl;
//		*/
//
//	}
//
//	mxDestroyArray(X);
//	mxDestroyArray(Xdot);
//	mxDestroyArray(Theta);
//	mxDestroyArray(Thetadot);
//	engEvalString(ep, "close;");
//	_getch();
//	return 0;
//}
//
//
//
//vector<double>time(vector<double>myvec)
//{
//	double start_time = 0.0, end_time = 10.0, dt = 0.005;
//	unsigned int N = 0;
//	myvec[N] = start_time;
//	while(start_time <= end_time)
//	{  
//		++N;
//		start_time += dt;
//		myvec[N] = start_time;
//		return myvec;
//	}
//}
//vector<double> deg2rad(int deviation, vector<double>thetadot, int temp)
//{
//	double pi = 3.142;
//	for(int i = 0; i < 3; ++i)
//		thetadot[i] = ((((2*deviation*((double)  rand()/ RAND_MAX)-deviation)) *pi) /180);
//	return thetadot;
//}
//vector<double>input()
//{
//	vector<double>in(4);
//	in <<= (700+10)*(740), 700*700, 740*(700+10), 700*700;
//	return in;
//}
//double sum(vector<double>in)
//{
//	double suminput = 0;
//	for (int i = 0; i < 4; i++)
//		suminput += in[i];
//	return suminput;
//}
//vector<double>thrust(vector<double>in, double k)
//{
//	vector<double>T(3);
//	T <<= 0, 0, k*sum(in);
//	return T;
//}
//vector<double>torques(vector<double>in, double L, double b, double k)
//{
//	vector<double>tau(3);
//	tau <<= L * k * (in[0]-in[2]), L*k*(in[1]-in[3]), b*(in[0]-in[1]+in[2]-in[3]);
//	return tau;
//}
//
//matrix<double>rotation(vector<double>angle)
//{
//	double phi = angle[2],theta_ = angle[1], psi = angle[0];
//	matrix<double>R(3,3);
//	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
//		cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
//		-sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
//	return R;
//}
//
//vector<double>acceleration(vector<double>in, vector<double>angle, vector<double>xdot_, double m, double g, double k, double kd)
//{
//
//	vector<double>gravity(3), T_, FD(3), a(3);
//	matrix<double>R(3,3);
//	gravity <<= 0, 0,-g;
//	R = rotation(angle);
//	T_ = prod(R,thrust(in,k));
//	FD = -kd*xdot_;
//	a = gravity+(1/m)*T_+FD;
//	return a;
//}
//vector<double>angular_acceleration(vector<double>in, matrix<double> I, vector<double>omega, double L, double b, double k)
//{	
//	vector<double>tau(3), omegadot(3);
//	tau = torques(in, L, b, k);
//	omegadot = prod(-I, tau- ( omega, prod(I, omega)));
//	return omegadot;
//}
//vector<double>thetadot2omega(vector<double>thetadot, vector<double>angle)
//{
//	double phi = angle[0], theta_ = angle[1], psi = angle[2];
//	vector<double>omega;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	omega = prod(w,thetadot );
//	return omega;
//}
//vector<double>omega2thetadot(vector<double>omega, vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>thetadot;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	thetadot = prod(-w,omega);
//	return thetadot;



//
//
//
//#include "stdio.h"
//#include <stdlib.h>
//#include <conio.h>
//#include <algorithm>
//#include <iostream>
//#include <vector>
//#include <time.h>
//#include <math.h>
//#include <fstream>
//#include <iterator>
//#include <cstdlib>
//#include <boost/numeric/ublas/assignment.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/range/iterator_range.hpp>
//#include "engine.h"
//
//using namespace boost::numeric::ublas;
//
//void createvector(vector<double>&myvec);
//void initialState(Engine *ep);
//void visMatlab(Engine *ep, mxArray *X, mxArray * theta,mxArray * thetadot);
//struct Data;
//Data *mxcreatedoublematrix ( mxArray *&ioX, mxArray *&ioXdot, mxArray *&ioTheta, mxArray *&ioThetadot);
//void Transmit2MATLAB(Engine *ep, vector<double>x, vector<double>xdot, vector<double>theta, vector<double>thetadot, mxArray*X,  mxArray*Xdot,  mxArray*Theta,  mxArray*Thetadot);
//vector<double> input();
//matrix<double>rotation(vector<double>angle);
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
//double sum(vector<double>in);
//vector<double> thrust(vector<double>in,double k);
//vector<double> torques(vector<double>in,double L,double b,double k);
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k);
//vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
//#define  BUFSIZE 256
//struct Data {
//mxArray*X;
//mxArray*Xdot;
//mxArray*Theta;
//mxArray*Thetadot;
//};
//int main()
//{	
//	Engine *ep;
//	mxArray *X = NULL,*Theta=NULL,*Xdot=NULL,*Thetadot=NULL;
//	char buffer[BUFSIZE+1];
//	if (!(ep = engOpen(""))) { 
//		fprintf(stderr, "\nCan't start MATLAB engine\n");
//		return EXIT_FAILURE;
//	}
//	srand( (unsigned int)time(0));
//	double i,g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25 ,dt=0.005;
//	unsigned int j;
//	int temp=0,deviation =100;
//	vector<double>timevector(2001),x(3),xdot(3,0),theta(3,0),thetadot(3),omega,a,omegadot,in(4,0);
//	matrix<double> I(3,3);
//	I <<= 5e-3, 0, 0,
//		0, 5e-3, 0,
//		0, 0, 10e-3;
//	x <<= 0, 0, 10;
//	createvector(timevector);
//	thetadot=deg2rad(deviation,thetadot,temp);
//	Data* returnData = mxcreatedoublematrix(X, Xdot, Theta, Thetadot);
//	X = returnData->X;
//	Xdot = returnData->Xdot;
//	Theta = returnData->Theta;
//	Thetadot = returnData->Thetadot;
//	/*returnData = mxcreatedoublematrix()
//	X = returnData->X
//	Xdot = returnData->XDot
//	...*/
//	initialState(ep);
//
//	for ( j = 0;j <2001; j++)
//	{
//		i=timevector[j];
//		in=input();
//		omega=thetadot2omega(thetadot,theta);
//		a=acceleration(in,theta,xdot,m,g,k,kd); 
//		omegadot=angular_acceleration(in,I,omega,L,b,k);
//		omega=omega+dt*omegadot;
//		thetadot=omega2thetadot(omega,theta);
//		theta=theta+(dt*thetadot);
//		xdot=xdot+dt*a;
//		x=x+dt*xdot;
//		
//		Transmit2MATLAB(ep, x, xdot, theta, thetadot, X, Xdot, Theta, Thetadot);
//		
//		visMatlab(ep,X, Theta,Thetadot);
//
//	}
//	
//	mxDestroyArray(X);
//	mxDestroyArray(Xdot);
//	mxDestroyArray(Theta);
//	mxDestroyArray(Thetadot);
//	engClose(ep);
//	_getch();
//	return 0;
//}
//void initialState(Engine *ep)
//{
//	mxArray *T = NULL,*Dt=NULL;
//	T = mxCreateDoubleMatrix(1, 1, mxREAL);
//	*mxGetPr(T) = 0;
//	Dt = mxCreateDoubleMatrix(1, 1, mxREAL);
//	*mxGetPr(Dt) = 0.05;
//	engPutVariable(ep, "t", T);
//	engPutVariable(ep, "dt", Dt); 
//	engEvalString(ep, "cd 'E:\\Visualizepart'");
//	engEvalString(ep, "createfigure();");
//	engEvalString(ep, "plots = createplot();");
//	engEvalString(ep, "subplot(plots(1));");
//	engEvalString(ep, "[model, thrusts] = quadcopter();");//[model, thrusts] = quadcopter();
//	engEvalString(ep, "createlabel();");
//	//engEvalString(ep, "x, xdot,theta,thetadot");//setVisLimits()
//}
//void Transmit2MATLAB(Engine *ep, vector<double>x, vector<double>xdot, vector<double>theta, vector<double>thetadot, mxArray*X,  mxArray*Xdot,  mxArray*Theta,  mxArray*Thetadot)
//{
//
//	std::copy(x.begin(), x.end(), mxGetPr(X));
//	std::copy(theta.begin(), theta.end(), mxGetPr(Theta));
//	std::copy(xdot.begin(), xdot.end(), mxGetPr(Xdot));
//	std::copy(thetadot.begin(), thetadot.end(), mxGetPr(Thetadot));
//	engPutVariable(ep, "x", X);
//	engPutVariable(ep, "xdot", Xdot);
//	engPutVariable(ep, "theta", Theta);
//	engPutVariable(ep, "thetadot", Thetadot);
//}
//
//
//void visMatlab(Engine *ep, mxArray *X, mxArray *theta, mxArray *thetadot)
//{
//	engEvalString(ep, "t = t + dt;");
//	engEvalString(ep, "plots = createplot();");
//	engEvalString(ep, "subplot(plots(1));");
//	engEvalString(ep, "move = makehgtform('translate',x);");//move = makehgtform('translate',x);
//	engEvalString(ep, "rotate = getrotate(theta);");//rotate = getRotation(theta); 
//	engEvalString(ep, "set(model,'Matrix', move * rotate);");//set(model,'Matrix', move * rotate);
//	engEvalString(ep, "getthrustscales(thrusts, move, rotate);");//visThrust(thrusts, move, rotate)
//	engEvalString(ep, "[xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x);");//[xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x);
//	engEvalString(ep, "axis([xmin xmax ymin ymax zmin zmax]);");//axis([xmin xmax ymin ymax zmin zmax]);
//	engEvalString(ep, "hold on;");//drawnow;
//	engEvalString(ep, "getveldis(t, thetadot, theta);");//visVelocities(t, theta, thetadot)
//}
//
//
//Data* mxcreatedoublematrix ( mxArray *&ioX, mxArray *&ioXdot, mxArray *&ioTheta, mxArray *&ioThetadot)
//{
//	Data* returnData = new Data;
//	returnData->X = mxCreateDoubleMatrix(1, 3, mxREAL);
//	returnData->Xdot = mxCreateDoubleMatrix(1, 3, mxREAL);
//	returnData->Theta = mxCreateDoubleMatrix(1, 3, mxREAL);
//	returnData->Thetadot = mxCreateDoubleMatrix(1, 3, mxREAL);
//	return returnData;	
//}
//
//
//void createvector(vector<double>& myvec)
//{
//	double start_time=0.0,end_time=10.0, dt=0.005;
//	unsigned int N=0;
//	myvec[N]=start_time;
//	while(start_time<=end_time)
//	{  
//		++N;
//		start_time+=dt;
//		myvec[N]=start_time;
//	}
//}
//vector<double> deg2rad(int deviation,vector<double>thetadot,int temp)
//{
//	double pi=3.142;
//	for(int i=0;i<3;++i)
//		thetadot[i]=((((2*deviation*((double)  rand()/ RAND_MAX)-deviation))*pi)/180);
//	return thetadot;
//}
//vector<double> input()
//{
//	vector<double> in(4);
//	in<<=(700+10)*(740),700*700,740*(700+10),700*700;
//	return in;
//}
//double sum(vector<double>in)
//{
//	double suminput=0;
//	for(int i=0;i<4;i++)
//		suminput+=in[i];
//	return suminput;
//}
//vector<double> thrust(vector<double>in,double k)
//{
//	vector<double>T(3);
//	T <<= 0, 0, k*sum(in);
//	return T;
//}
//vector<double> torques(vector<double>in,double L,double b,double k)
//{
//	vector<double>tau(3);
//	tau <<=L*k*(in[0]-in[2]), L*k*(in[1]-in[3]),b*(in[0]-in[1]+in[2]-in[3]);
//	return tau;
//}
//
//matrix<double>rotation(vector<double>angle)
//{
//	double phi = angle[2],theta_ = angle[1], psi = angle[0];
//	matrix<double>R(3,3);
//	R <<= cos(phi)*cos(theta_),cos(phi) * sin(theta_) * sin(psi) - cos(psi) * sin(phi),sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta_),
//		cos(theta_)*sin(phi), cos(phi) * cos(psi) + sin(phi) * sin(theta_) * sin(psi),cos(psi) * sin(phi) * sin(theta_) - cos(phi) * sin(psi),
//		-sin(theta_), cos(theta_)*sin(phi),cos(theta_)*cos(phi);
//	return R;
//}
//
//vector<double> acceleration(vector<double>in,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
//{
//
//	vector<double>gravity(3),T_,FD(3),a(3);
//	matrix<double>R(3,3);
//	gravity <<= 0, 0,-g;
//	R=rotation(angle);
//	T_=prod(R,thrust(in,k));
//	FD=-kd*xdot_;
//	a=gravity+(1/m)*T_+FD;
//	return a;
//}
//vector<double> angular_acceleration(vector<double>in,matrix<double> I,vector<double>omega,double L,double b,double k)
//{	
//	vector<double>tau(3),omegadot(3);
//	tau=torques(in,L,b,k);
//	//std::cout << inner_prod(omega,prod(I,omega)) << std::endl;
//	//std::cout << std::endl;
//	omegadot=prod(-I,(tau-element_prod(omega,prod(I,omega))));
//	return omegadot;
//}
//vector<double>thetadot2omega(vector<double>thetadot,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>omega;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	omega=prod(w,thetadot );
//	return omega;
//}
//vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
//{
//	double phi = angle[0],theta_ = angle[1], psi = angle[2];
//	vector<double>thetadot;
//	matrix<double>w(3,3);
//	w <<= 1, 0, -sin(theta_),
//		0, cos(phi), cos(theta_)*sin(phi),
//		0, -sin(phi), cos(theta_)*cos(phi);
//	thetadot=prod(-w,omega);
//	return thetadot;
//}





#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <iterator>
#include <cstdlib>
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

using namespace boost::numeric::ublas;

void createvector(vector<double>&myvec);
void initialState(Engine *ep);
struct Data;
void visMatlab(Engine *ep, mxArray *X, mxArray *Xdot, mxArray * Theta, mxArray *Thetadot);
Data* mxcreatedoublematrix(mxArray *&ioX, mxArray *&ioXdot, mxArray *&ioTheta, mxArray *&ioThetadot);
void Transmit2MATLAB(Engine *ep, vector<double>x, vector<double>xdot, vector<double>theta, vector<double>thetadot, mxArray*X,  mxArray*Xdot,  mxArray*Theta,  mxArray*Thetadot);
vector<double> input(double ii);
matrix<double>rotation(vector<double>angle);
vector<double> deg2rad(int deviation,vector<double>thetadot,int temp);
double sum(vector<double>inputs);
vector<double> thrust(vector<double>inputs,double k);
vector<double> torques(vector<double>inputs,double L,double b,double k);
vector<double> acceleration(vector<double>inputs, vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd);
vector<double> angular_acceleration(vector<double>inputs, matrix<double> I,vector<double>omega,double L,double b,double k);
vector<double> thetadot2omega(vector<double>thetadot,vector<double>angle);
vector<double>omega2thetadot(vector<double>omega,vector<double>angle);
#define  BUFSIZE 256
struct Data{
	mxArray * X;
	mxArray * Xdot;
	mxArray * Theta;
	mxArray * Thetadot;
};

int main()
{	
	Engine *ep;
	mxArray *X = NULL,*Theta=NULL,*Xdot=NULL,*Thetadot=NULL;
	/*double *xp,*rp;*/
	char buffer[BUFSIZE+1];
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}
	// TEST get rid of it!!!
	//mxArray* Tm = mxCreateDoubleMatrix(1, 1, mxREAL);
	//*mxGetPr(Tm) = 0;
	//mxArray* Ts = mxCreateDoubleScalar(0);
	//assert(*mxGetPr(Ts) == 0);
	//return
	
	srand( (unsigned int)time(0));
	double g=9.81,m=0.5,L=0.25,k=3e-6,b=1e-7,kd=0.25 ,dt=0.05;
	unsigned int j;
	int temp=0,deviation =100;
	vector<double>timevector(2001),ain(4), x(3),xdot(3,0),theta(3,0),thetadot(3),omega(3),a,omegadot,in(4,0);

	matrix<double> I(3,3);
	I <<= 5e-3, 0, 0,
		0, 5e-3, 0,
		0, 0, 10e-3;
	x <<= 0, 0, 10;

	thetadot=deg2rad(deviation,thetadot,temp);
	Data* returnData = mxcreatedoublematrix(X, Xdot, Theta, Thetadot);
	X = returnData->X;
	Xdot = returnData->Xdot;
	Theta = returnData->Theta;
	Thetadot = returnData->Thetadot;
	
	initialState(ep);

	for ( j = 0 ; j  < 50; j++)
	{
	
		double i = 0;
		i += dt; 
		ain = input(i);
		std::cout<<"ain>>"<<ain<<std::endl;
		omega=thetadot2omega(thetadot,theta);
		a=acceleration(ain,theta,xdot,m,g,k,kd); 
		omegadot=angular_acceleration(ain,I,omega,L,b,k);
		omega=omega+dt*omegadot;
		thetadot=omega2thetadot(omega,theta);
		std::cout<<"thetadot>>"<<thetadot<<std::endl;
		theta=theta+(dt*thetadot);
		std::cout<<"theta>>"<<theta<<std::endl;
		xdot=xdot+dt*a;
		std::cout<<"xdot>>"<<xdot<<std::endl;
		x=x+dt*xdot;
		std::cout<<"x>>"<<x;
		Transmit2MATLAB(ep, x, xdot, theta, thetadot, X, Xdot, Theta, Thetadot);
		visMatlab(ep, X, Xdot, Theta, Thetadot);
		std::cout<<std::endl;
	}
	
	mxDestroyArray(X);
	mxDestroyArray(Xdot);
	mxDestroyArray(Theta);
	mxDestroyArray(Thetadot);

	engClose(ep);

	_getch();
	return 0;
}
void initialState(Engine *ep)
{
	mxArray *T = NULL,*Dt=NULL;
	
	T = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(T) == 0;
	Dt = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(Dt) = 0.05;
	engPutVariable(ep, "t", T);
	engPutVariable(ep, "dt", Dt); 
	engEvalString(ep, "cd 'E:\\Matlab & C++ code\\Test private function'");
	engEvalString(ep, "createfigure();");
	engEvalString(ep, "plots=createplot();");
	engEvalString(ep, "subplot(plots(1));");
	engEvalString(ep, "[model, thrusts] = quadcopter();");//[model, thrusts] = quadcopter();
	engEvalString(ep, "createaxis();");
	engEvalString(ep, "createlabel();");
	
}


void visMatlab(Engine *ep, mxArray *X, mxArray *Xdot, mxArray *Theta, mxArray *Thetadot)
{
	
	engEvalString(ep, "t = t + dt;");
	engEvalString(ep, "subplot(plots(1));");
	engEvalString(ep, "move = makehgtform('translate',x);");//move = makehgtform('translate',x);
	engEvalString(ep, "rotate = getRotation(theta);");//rotate = getRotation(theta); 
	engEvalString(ep, "set(model,'Matrix', move * rotate);");//set(model,'Matrix', move * rotate);
	engEvalString(ep, "visThrust(thrusts, move, rotate);");//visThrust(thrusts, move, rotate)
	engEvalString(ep, "[xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x);");//[xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x);
	engEvalString(ep, "axis([xmin xmax ymin ymax zmin zmax]);");//axis([xmin xmax ymin ymax zmin zmax]);
	engEvalString(ep, "hold on;");//drawnow;
	engEvalString(ep, "drawnow;");//drawnow;
	engEvalString(ep, "visVelocities(t, theta, thetadot);");//visVelocities(t, theta, thetadot)
}


void Transmit2MATLAB(Engine *ep, vector<double>x, vector<double>xdot, vector<double>theta, vector<double>thetadot, mxArray*X,  mxArray*Xdot,  mxArray*Theta,  mxArray*Thetadot)
{
	
	std::copy(x.begin(), x.end(), mxGetPr(X));
	std::copy(theta.begin(), theta.end(), mxGetPr(Theta));
	std::copy(xdot.begin(), xdot.end(), mxGetPr(Xdot));
	std::copy(thetadot.begin(), thetadot.end(), mxGetPr(Thetadot));
	engPutVariable(ep, "x", X);
	engPutVariable(ep, "xdot", Xdot);
	engPutVariable(ep, "theta", Theta);
	engPutVariable(ep, "thetadot", Thetadot);
	
}



Data* mxcreatedoublematrix(mxArray *&ioX, mxArray *&ioXdot, mxArray *&ioTheta, mxArray *&ioThetadot)
{
	Data * returnData = new Data;
	returnData->X = mxCreateDoubleMatrix(1, 3, mxREAL);
	returnData->Xdot = mxCreateDoubleMatrix(1, 3, mxREAL);
	returnData->Theta = mxCreateDoubleMatrix(1, 3, mxREAL);
	returnData->Thetadot = mxCreateDoubleMatrix(1, 3, mxREAL);
	return returnData;

}


void createvector(vector<double>& myvec)
{
	double start_time=0.0,end_time=10.0, dt=0.005;
	unsigned int N=0;
	myvec[N]=start_time;
	while(start_time<=end_time)
	{  
		++N;
		start_time+=dt;
		myvec[N]=start_time;
	}
}
vector<double> deg2rad(int deviation,vector<double>thetadot,int temp)
{
	double pi=3.142;
	for(int i=0;i<3;++i)
		thetadot[i]=((((2*deviation*((double)  rand()/ RAND_MAX)-deviation))*pi)/180);
	return thetadot;
}
vector<double> input(double ii)
{
	vector<double> in(4);
	double value;
	in<<=(718)*2,(550)*2,(718)*2,(550)*2;
	return in;
}
double sum(vector<double>inputs)
{
	double suminput=0;
	for(int i=0;i<4;i++)
		suminput+=inputs[i];
	return suminput;
}
vector<double> thrust(vector<double>inputs, double k)
{
	
	vector<double>T(3);
	T <<= 0, 0, k*sum(inputs);
	return T;
}
vector<double> torques(vector<double>inputs,double L,double b,double k)
{
	
	vector<double>tau(3);
	tau <<=L*k*(inputs[0]-inputs[2]), L*k*(inputs[1]-inputs[3]),b*(inputs[0]-inputs[1]+inputs[2]-inputs[3]);
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

vector<double> acceleration(vector<double>inputs,vector<double>angle,vector<double>xdot_,double m,double g,double k,double kd)
{

	vector<double>gravity(3),T_,FD(3),a(3);
	matrix<double>R(3,3);
	gravity <<= 0, 0,-g;
	R=rotation(angle);
	T_=prod(R,thrust(inputs,k));
	FD=-kd*xdot_;
	a=gravity+(1/m)*T_+FD;
	return a;
}
vector<double> angular_acceleration(vector<double>inputs,matrix<double> I,vector<double>omega,double L,double b,double k)
{	
	vector<double>tau(3),omegadot(3);
	tau=torques(inputs,L,b,k);
	//std::cout << inner_prod(omega,prod(I,omega)) << std::endl;
	//std::cout << std::endl;
	omegadot=prod(-I,(tau-element_prod(omega,prod(I,omega))));
	return omegadot;
}
vector<double>thetadot2omega(vector<double>thetadot,vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	vector<double>omega;
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
		0, cos(phi), cos(theta_)*sin(phi),
		0, -sin(phi), cos(theta_)*cos(phi);
	omega=prod(w,thetadot );
	return omega;
}
vector<double>omega2thetadot(vector<double>omega,vector<double>angle)
{
	double phi = angle[0],theta_ = angle[1], psi = angle[2];
	vector<double>thetadot;
	matrix<double>w(3,3);
	w <<= 1, 0, -sin(theta_),
		0, cos(phi), cos(theta_)*sin(phi),
		0, -sin(phi), cos(theta_)*cos(phi);
	thetadot=prod(-w,omega);
	return thetadot;
}



