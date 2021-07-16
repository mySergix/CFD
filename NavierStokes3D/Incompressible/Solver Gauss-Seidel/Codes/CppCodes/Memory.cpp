#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

#include "../HeaderCodes/Memory.h"

using namespace std;

Memory::Memory(){
	
}

//Memoria dinámica matriz (double)
double *Memory::AllocateDouble(int NX, int NY, int NZ, int Dim){
double *M1;

	M1 = new double [NX*NY*NZ*Dim];				
	return M1;

}

//Memoria dinámica matriz (int)
int *Memory::AllocateInt(int NX, int NY, int NZ, int Dim){ 
int *M1;

	M1 = new int [NX*NY*NZ*Dim];				
	return M1;

}




