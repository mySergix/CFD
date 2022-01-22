#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#define PI 3.141592653589793

class SpalartAllmaras{	

	public:
		//Constructor de la clase
		SpalartAllmaras();
		
		//Metodos de la clase

		//Métodos modelo de turbulencia Spalart-Allmaras
		void SA_Term1(Mesher); //Cálculo del término 1 de la ecuación de Spalart-Allmaras
		void SA_Term2(); //Cálculo del término 2 de la ecuación de Spalart-Allmaras
		void SA_Term3(Mesher); //Cálculo del término 3 de la ecuación de Spalart-Allmaras
		void SA_Term4(Mesher); //Cálculo del término 4 de la ecuación de Spalart-Allmaras
		void SA_Term5(Mesher); //Cálculo del término 5 de la ecuación de Spalart-Allmaras
		
		void FmuSA(); //Cálculo de todas las contribuciones del modelo Spalart-Allmaras
		
		void SpalartAllmarasPreparation(Mesher); //Cálculo de variables y matrices necesarias para aplicar el modelo Spalart-Allmaras
		void InitialSpalartAllmaras();
		void SpalartAllmaras(); //Cálculo de la viscosidad con el modelo Spalart-Allmaras

			
};
