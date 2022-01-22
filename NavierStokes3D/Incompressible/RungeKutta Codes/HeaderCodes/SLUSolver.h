#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "umfpack.h"

using namespace std;

#define PI 3.141592653589793

class SLUSolver{	

	public:
		//Constructor de la clase
		SLUSolver(Solver);
		
		int *Ap;
		int *Ai;
		double *Ax;
		double *b;

		int NX;
		int NY;
		int NZ;

		//Metodos de la clase
        void Build_Ap_vector(int);

		void Build_BP_vector(double*);
		void Build_Ai_Ax_vectors(int, Solver);

		void Get_LaplacianMatrixA(Solver, double *& null, void *& Symbolic, void *& Numeric);
		void Get_Pressure(double *x, void *& Numeric, double *& null);

};
