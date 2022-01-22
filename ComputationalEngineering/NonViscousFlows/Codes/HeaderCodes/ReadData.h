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

		string Problema; //Problema Tobera/Tuberia

		double *GeometryData; //Datos de la geometría del problema
		int *NumericalData; //Datos numéricos del problema
		double *ProblemData; //Datos del problema
		double *ProblemPhysicalData; //Datos físicos sobre las condiciones del problema

		int *MeshStudyNX;
		int *MeshStudyNY;
		double *FactorRelax;
		
		double *AttackAngle;
		
		//Constructor de la clase
		ReadData(Memory, int);

		void ReadInputs(); //Lector datos en ficheros
		void ReadArrays(string, int , double*);
		void ReadMeshStudy(string, int);

};
