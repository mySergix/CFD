#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#define PI 3.141592653589793

class Memory{	

	public:
		//Constructor de la clase
		Memory();
		
		//Metodos de la clase

			//Métodos de alojamiento de memoria
			int *AllocateInt(int, int, int, int);
			double *AllocateDouble(int, int, int, int);
			
};
