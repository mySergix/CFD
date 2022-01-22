#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
using namespace std;

#define PI 3.141592653589793

class ReadTXT{
	private:
		//Problems data
		int *NumericalData;
		double *PrecisionData;
		double *GeometryData;
		double *FluidData;
		double *GravityData;
		string Scheme;
		string Scheme2;

		//Boundary conditions for each problem
		double *BoundaryData;
		int *GradientData;
	
	public:
		ReadTXT();		
		void ReadFromTXT();

		int Get_NumericalData(int);
		double Get_PrecisionData(int);
		double Get_GeometryData(int);
		double Get_FluidData(int);
		double Get_GravityData(int);

		double Get_BoundaryData(int);	
		int Get_GradientData(int);
		string Get_ConvectiveScheme();
		string Get_ConvectiveScheme2();
	
};
