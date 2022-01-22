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

		string TipoProblema; //Problema Tobera/Tuberia
		string TipoTobera; //Tipo de Tobera del problema
		string TypeMeshing; //Tipo de mallado del problema

		double *GeometryData; //Datos de la geometría del problema
		int *NumericalData; //Datos numéricos del problema
		double *ProblemData; //Datos del problema
		double *SpalartAllmarasData; //Datos sobre el modelo de turbulencia Spalart-Allmaras
		double *HeatData; //Datos físicos del problema sobre transferencia de calor
		double *ProblemPhysicalData; //Datos físicos sobre las condiciones del problema

		double *CombustionData; //Datos químicos sobre la combustión y sus especies presentes
		double *JannafO2; //Coeficientes función Cp(T) y H(T) del O2
		double *JannafH2; //Coeficientes función Cp(T) y H(T) del H2
		double *JannafH2O; //Coeficientes función Cp(T) y H(T) del H2O
	
		double *ViscosityData; //Datos sobre la viscosidad de las especies termoquímicas

		//double *BoundaryData; //Condiciones de contorno del problema

		//Constructor de la clase
		ReadData(Memory);

		void ReadInputs(); //Lector datos en ficheros
		void ReadArrays(string, int , double*);

		string Get_ProblemType(){ return TipoProblema; } //Enviar dato tipo de problema (Tubería o Tobera)
		string Get_NozzleType(){ return TipoTobera; } //Enviar dato tipo de tobera (Parte divergente, NASA o RAO)
		string Get_MeshingType(){ return TypeMeshing; } //Enviar dato tipo de mallado utilizado (Collocated o Staggered)

		
};
