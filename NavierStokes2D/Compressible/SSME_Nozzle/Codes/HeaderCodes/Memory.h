#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
using namespace std;

#define PI 3.141592653589793

class Memory{	

	public:
		//Constructor de la clase
		Memory();
		
		//Metodos de la clase
			//Métodos de alojamiento de memoria
			int *AllocateInt1D(int);

			double *AllocateDouble1D(int);
			double **AllocateDouble2D(int, int);
			double ***AllocateDouble3D(int, int, int, string);
			
};
