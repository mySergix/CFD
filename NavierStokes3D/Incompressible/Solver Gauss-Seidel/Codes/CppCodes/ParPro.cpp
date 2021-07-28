#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <cassert>
#include <time.h>
#include <mpi.h>

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/ParPro.h"

using namespace std;

//Global Index P Mesh
#define GP(i,j,k,Dim) (NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*((i) + HaloPressure) + ((j) + HaloPressure) + ((k) + HaloPressure)*(NY + 2*HaloPressure) + (NX + 2*HaloPressure)*(NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*(Dim)

//Global Index U Mesh
#define GU(i,j,k,Dim) (NY + 2*HaloU)*(NZ + 2*HaloU)*((i) + HaloU) + ((j) + HaloU) + ((k) + HaloU)*(NY + 2*HaloU) + (NX + 1 + 2*HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU)*(Dim)

//Global Index V Mesh
#define GV(i,j,k,Dim) (NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*((i) + HaloV) + ((j) + HaloV) + ((k) + HaloV)*(NY + 1 + 2*HaloV) + (NX + 2*HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*(Dim)

//Global Index W Mesh
#define GW(i,j,k,Dim) (NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*((i) + HaloW) + ((j) + HaloW) + ((k) + HaloW)*(NY + 2*HaloW) + (NX + 2*HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*(Dim)

//Local Index Pressure (P) Mesh
#define LP(i,j,k,dim) (NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*((i) - Ix + HaloPressure) + (j + HaloPressure) + (k + HaloPressure)*(NY + 2*HaloPressure)

//Local Index Velocity (U) Mesh
#define LU(i,j,k,dim) (NY + 2*HaloU)*(NZ + 2*HaloU)*((i) - Ix + HaloU) + ((j) + HaloU) + ((k) + HaloU)*(NY + 2*HaloU)

//Local Index Velocity (V) Mesh
#define LV(i,j,k,dim) (NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*((i) - Ix + HaloV) + ((j) + HaloV) + ((k) + HaloV)*(NY + 1 + 2*HaloV)

//Local Index Velocity (W) Mesh
#define LW(i,j,k,dim) (NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*((i) - Ix + HaloW) + ((j) + HaloW) + ((k) + HaloW)*(NY + 2*HaloW)

//Local Index Pressure Coefficients (A)
#define LA(i,j,k,dim) (NY*NZ)*((i) - Ix) + (j) + (k)*NY

ParPro::ParPro(ReadData R1){

		NX = R1.ProblemNumericalData[2];
		NY = R1.ProblemNumericalData[3];
		NZ = R1.ProblemNumericalData[4];
		Halo = 2;

		HaloPressure = 1;
		HaloU = 2;
		HaloV = 2;
		HaloW = 2;

		HP = 2;
}

void ParPro::Rango(){
	int a;
	a = MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
}

void ParPro::Processes(){
	int a;
	a = MPI_Comm_size(MPI_COMM_WORLD, &Procesos);
}

void ParPro::Initial_WorkSplit(int NX, int &Ix, int &Fx){
int Intervalo, Residuo;
int p;
	Intervalo = NX/Procesos;
	Residuo = NX%Procesos;

	if(Rank != Procesos-1){
		Ix = Rank*Intervalo;
		Fx = (Rank+1)*Intervalo;
	}
	else{
		Ix = Rank*Intervalo;
		Fx = (Rank+1)*Intervalo + Residuo;
	}
}

void ParPro::Get_Worksplit(int NX, int Procesos, int p, int &pix, int &pfx){
int Intervalo, Residuo;

	Intervalo = NX/Procesos;
	Residuo = NX%Procesos;

	if(p != Procesos-1){
		pix = p*Intervalo;
		pfx = (p+1)*Intervalo;
	}
	else{
		pix = p*Intervalo;
		pfx = (p+1)*Intervalo + Residuo;
	}
}

//Comunicar Matriz de Presión
void ParPro::CommunicateDataLP(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;	

	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LP(Fx - HaloPressure, 0, 0, 0)], (HaloPressure)*(NY)*(NZ), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LP(Ix - HaloPressure, 0, 0, 0)], (HaloPressure)*(NY)*(NZ), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LP(Ix, 0, 0, 0)], (HaloPressure)*(NY)*(NZ), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LP(Fx, 0, 0, 0)], (HaloPressure)*(NY)*(NZ), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

//Comunicar Matriz de Velocidad U
void ParPro::CommunicateDataLU(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;

	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LU(Fx - HaloU, - HaloU, - HaloU, 0)], (HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LU(Ix - HaloU, - HaloU, - HaloU,0)], (HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LU(Ix + 1, - HaloU, - HaloU,0)], (HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LU(Fx + 1,- HaloU,- HaloU,0)], (HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

//Comunicar Matriz de Velocidad V
void ParPro::CommunicateDataLV(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;

	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LV(Fx - HaloV, - HaloV, - HaloV, 0)], (HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LV(Ix - HaloV, - HaloV, - HaloV, 0)], (HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LV(Ix,- HaloV, - HaloV, 0)], (HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LV(Fx,- HaloV,- HaloV, 0)], (HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD, &ST);
	}

}

//Comunicar Matriz de Velocidad W
void ParPro::CommunicateDataLW(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;
	
	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LW(Fx - HaloW,- HaloW, - HaloW, 0)], (HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LW(Ix - HaloW, - HaloW, - HaloW,0)], (HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LW(Ix, - HaloW, - HaloW, 0)], (HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW), MPI_DOUBLE, Rank - 1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LW(Fx, - HaloW, - HaloW, 0)], (HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW), MPI_DOUBLE, Rank + 1, 0, MPI_COMM_WORLD, &ST);
	}

}

void ParPro::SendMatrixToZeroMP(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LP(Ix, - HaloPressure, - HaloPressure, 0)], (Fx - Ix)*(NY + 2*HaloPressure)*(NZ + 2*HaloPressure), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					GlobalMatrix[GP(i,j,k,0)] = LocalMatrix[LP(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GP(pix, - HaloPressure, - HaloPressure, 0)], (pfx - pix)*(NY + 2*HaloPressure)*(NZ+ 2*HaloPressure), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

	MPI_Barrier(MPI_COMM_WORLD);	
}

void ParPro::SendMatrixToZeroMU(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LU(Ix, - HaloU, - HaloU, 0)], (Fx - Ix + 1)*(NY + 2*HaloU)*(NZ + 2*HaloU), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx + 1; i++){
			for(j = - HaloU; j < NY + HaloU; j++){
				for(k = - HaloU; k < NZ + HaloU; k++){
					GlobalMatrix[GU(i,j,k,0)] = LocalMatrix[LU(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GU(pix, - HaloU, - HaloU, 0)], (pfx - pix + 1)*(NY + 2*HaloU)*(NZ + 2*HaloU), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

void ParPro::SendMatrixToZeroMV(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LV(Ix, - HaloV, - HaloV, 0)], (Fx - Ix)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = - HaloV; j < NY + 1 + HaloV; j++){
				for(k = - HaloV; k < NZ + HaloV; k++){
					GlobalMatrix[GV(i,j,k,0)] = LocalMatrix[LV(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GV(pix, - HaloV, - HaloV, 0)], (pfx - pix)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

void ParPro::SendMatrixToZeroMW(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LW(Ix, - HaloW, - HaloW, 0)], (Fx - Ix)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		
		for(i = Ix; i < Fx; i++){
			for(j = - HaloW; j < NY + HaloW; j++){
				for(k = - HaloW; k < NZ + 1 + HaloW; k++){
					GlobalMatrix[GW(i,j,k,0)] = LocalMatrix[LW(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GW(pix, - HaloW, - HaloW, 0)], (pfx - pix)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);	
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

void ParPro::SendDataToZero(double DataSent, double *DataReceived){
int i, p;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&DataSent, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);		
	}
	

	if(Rank == 0){
		DataReceived[Rank] = DataSent;
		for(p = 1; p < Procesos; p++){
			MPI_Recv(&DataReceived[p], 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}
	}
}

void ParPro::SendDataToAll(double DataSent, double &DataReceived){
int p;
MPI_Status ST;

	if(Rank == 0){
		for(p = 1; p < Procesos; p++){
			MPI_Send(&DataSent, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);		
		}
		DataReceived = DataSent;
	}
	
	if(Rank != 0){
		MPI_Recv(&DataReceived, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &ST);
	}

}
/*
void ParPro::SendMatrixToZeroBP(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LAL(Ix, 0, 0, 0)], (Fx-Ix)*(NY)*(NZ), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					GlobalMatrix[LAG(i,j,k,0)] = LocalMatrix[LAL(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[LAG(pix, 0, 0, 0)], (pfx - pix)*(NY)*(NZ), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

	MPI_Barrier(MPI_COMM_WORLD);	
}
*/
void ParPro::Execute(){

	Rango();
	Processes();
	
}
