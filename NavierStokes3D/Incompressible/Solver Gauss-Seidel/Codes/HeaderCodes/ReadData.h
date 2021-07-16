#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

class ReadData{
	private:
		
		
	public:

		string ConvectiveScheme1;
		string ConvectiveScheme2;

		string DIRECTORIO;
		
		double *GeometryData; //Datos de la geometría del problema
		double *ProblemPhysicalData; //Datos físicos sobre las condiciones del problema
		double *ProblemData; //Datos del problema
		int *ProblemNumericalData; //Datos numéricos del problema
		
		//Constructor de la clase
		ReadData(Memory, string);

		void ReadInputs(); //Lector datos en ficheros
		void ReadArrays(string, int , double*);

		string Get_Convective1(){ return ConvectiveScheme1; } 
		string Get_Convective2(){ return ConvectiveScheme2; }
		
};
