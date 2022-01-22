#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
using namespace std;

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/"

Memory::Memory(){
	
}

//Memoria dinámica matriz 3D (double)
double ***Memory::AllocateDouble3D(int NA, int NR, int Dim, string Mesh){
int i, j, NMesh;

	double ***M1;

	if(Mesh == "P"){ NMesh = 1; }
	else if(Mesh == "V"){ NMesh = 2; }
	else if(Mesh == "U"){ NMesh = 3; }
	
		if(NMesh == 1 || NMesh == 2){ M1 = new double **[NA]; }
			
		else if(NMesh == 3){ M1 = new double **[NA+1]; }	
	
		for (i = 0; i < NA+1; i++){
			if(i < NA){
				if(NMesh == 1 || NMesh == 3){
					M1[i] = new double *[NR];
				}
				else{
					M1[i] = new double *[NR+1];
				}	
			}
			else{
				if(NMesh == 3){
				M1[i] = new double *[NR];
				}
			
			}
			
		}
		for (i = 0; i < NA+1; i++){
			for (j = 0; j < NR+1; j++){
				if(i < NA && j < NR){
					M1[i][j] = new double [Dim];
				
				}
				else if(j == NR && i < NA){
					if(NMesh == 2){
						M1[i][j] = new double [Dim];
					}
				}
				else if(j < NR && i == NA){
					if(NMesh == 3){
						M1[i][j] = new double [Dim];
					}
				}		
			}
		}
	return M1;
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




