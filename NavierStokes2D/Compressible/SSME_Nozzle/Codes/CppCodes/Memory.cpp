#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "/home_nobck/sergiogus/TFG/Codes/HeaderCodes/Memory.h"
using namespace std;

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

//Memoria dinámica matriz 2D (double)
double **Memory::AllocateDouble2D(int NA, int NR){
int i, j;
double **M1;

	M1 = new double *[NA];

	for(i = 0; i < NA; i++){ M1[i] = new double [NR]; }
				
	return M1;
}

//Memoria dinámica matriz 1D (double)
double *Memory::AllocateDouble1D(int N){ 
double *M1;
		M1 = new double [N];
		
		return M1;
}

//Memoria dinámica matriz 1D (int)
int *Memory::AllocateInt1D(int N){ 
int *M1;
		M1 = new int [N];
		
		return M1;
}




