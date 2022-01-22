#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#define PI 3.141592653589793

class PostProc{	

	public:
		//Constructor de la clase
		PostProc(Memory, ReadData, Mesher, Solver);
		
		//Metodos de la clase

		void EscalarVTK3D(string, string, string, double*, double*, int, int, int);
		void EscalarVTK3DU(string, string, string, double*, double*, int, int, int);
		void VectorialVTK3D(Mesher, string, string, string, double*, double*, double*, double*, int, int, int);

		void Get_DrivenResults(Solver, Mesher, double, int);
		void Get_DifferentiallyResults(Solver, Mesher, double, int);

		void PrintTxt(Solver);

        void ExecutePostProcessing(Solver, Mesher, Memory, ParPro);

        void DeleteEverything(Solver, Mesher, ReadData);
};
