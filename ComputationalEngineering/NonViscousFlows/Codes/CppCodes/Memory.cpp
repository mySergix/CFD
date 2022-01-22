#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Memory.h"
using namespace std;

#define DIRECTORIO "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/"

Memory::Memory(){
	
}

//Memoria dinámica matriz (double)
double *Memory::AllocateDouble(int NX, int NY, int Dim){
double *M1;

	M1 = new double [NX*NY*Dim];				
	return M1;
}

//Memoria dinámica matriz (int)
int *Memory::AllocateInt(int NX, int NY, int Dim){ 
int *M1;

	M1 = new int [NX*NY*Dim];				
	return M1;
}



