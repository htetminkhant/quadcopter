// 2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <vector>

int main()
{
	double start_time=0.0,end_time=10.0, dt=0.005;
	std::vector<double>myvec;
	do
	{   start_time+=dt;
		myvec.push_back(start_time);

	}
	while(start_time<=end_time);
	std::vector<double> :: iterator i;
	for(int i=0;i<myvec.size();i++)
	{
		std::cout << myvec[i]<<" ";
	}
	return 0;
}
