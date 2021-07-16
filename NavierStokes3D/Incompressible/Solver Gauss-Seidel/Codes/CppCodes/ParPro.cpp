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

#define GP(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+2*HP)*(NY+2*HP)*(NZ+2*HP)*(dim)) //Global Index P Mesh
#define GU(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+1+2*HP)*(NY+2*HP)*(NZ+2*HP)*(dim)) //Global Index U Mesh
#define GV(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) + ((NX+2*HP)*(NY+1+2*HP)*(NZ+2*HP)*(dim)) //Global Index V Mesh
#define GW(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+2*HP)*(NY+2*HP)*(NZ+1+2*HP)*(dim)) //Global Index W Mesh

//Local Index P Mesh
#define LPL(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Left Core
#define LPC(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Center Cores
#define LPR(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Right Core

//Local Index U Mesh
#define LUL(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Left Core
#define LUC(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Center Cores
#define LUR(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Right Core

//Local Index V Mesh
#define LVL(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i) + HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Left Core
#define LVC(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Center Cores
#define LVR(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Right Core

//Local Index W Mesh
#define LWL(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i) + HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Left Core
#define LWC(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Center Cores
#define LWR(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Right Core

//Coeficientes A
#define LAG(i,j,k,dim) ((NY*NZ)*((i))) + ((j) + (k)*NY)
#define LAL(i,j,k,dim) ((NY*NZ)*((i) - Ix)) + ((j) + (k)*NY)

ParPro::ParPro(ReadData R1){

		NX = R1.ProblemNumericalData[2];
		NY = R1.ProblemNumericalData[3];
		NZ = R1.ProblemNumericalData[4];
		Halo = 2;

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
		MPI_Send(&LocalSend[LPC(Fx - Halo,- HP, - HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LPC(Ix - Halo,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LPC(Ix,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LPC(Fx,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

//Comunicar Matriz de Velocidad U
void ParPro::CommunicateDataLU(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;

	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LUC(Fx - Halo,- HP, - HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LUC(Ix - Halo,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LUC(Ix + 1,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LUC(Fx + 1,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

//Comunicar Matriz de Velocidad V
void ParPro::CommunicateDataLV(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;

	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LVC(Fx - Halo,- HP, - HP,0)], (Halo)*(NY + 1 + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LVC(Ix - Halo,- HP,- HP,0)], (Halo)*(NY + 1 + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LVC(Ix,- HP,- HP,0)], (Halo)*(NY + 1 + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LVC(Fx,- HP,- HP,0)], (Halo)*(NY + 1 + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD, &ST);
	}

}

//Comunicar Matriz de Velocidad W
void ParPro::CommunicateDataLW(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;
	
	if(Rank != Procesos - 1){
		MPI_Send(&LocalSend[LWC(Fx - Halo,- HP, - HP,0)], (Halo)*(NY + 2*HP)*(NZ + 1 + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD);
	}

	if(Rank != 0){
		MPI_Recv(&LocalReceive[LWC(Ix - Halo,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 1 + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rank != 0){
		MPI_Send(&LocalSend[LWC(Ix,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 1 + 2*HP), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD);
	}

	if(Rank != Procesos - 1){
		MPI_Recv(&LocalReceive[LWC(Fx,- HP,- HP,0)], (Halo)*(NY + 2*HP)*(NZ + 1 + 2*HP), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD, &ST);
	}

}

void ParPro::SendMatrixToZeroMP(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LPC(Ix,- HP,- HP,0)], (Fx-Ix)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = - HP; j < NY + HP; j++){
				for(k = - HP; k < NZ + HP; k++){
					GlobalMatrix[GP(i,j,k,0)] = LocalMatrix[LPL(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GP(pix,- HP,- HP,0)], (pfx - pix)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

	MPI_Barrier(MPI_COMM_WORLD);	
}

void ParPro::SendMatrixToZeroMU(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LUC(Ix,- HP,- HP,0)], (Fx-Ix + 1)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx + 1; i++){
			for(j = - HP; j < NY + HP; j++){
				for(k = - HP; k < NZ + HP; k++){
					GlobalMatrix[GU(i,j,k,0)] = LocalMatrix[LUL(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GU(pix,- HP,- HP,0)], (pfx - pix + 1)*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

void ParPro::SendMatrixToZeroMV(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LVC(Ix,- HP,- HP,0)], (Fx-Ix)*(NY + 1 + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = - HP; j < NY + 1 + HP; j++){
				for(k = - HP; k < NZ + HP; k++){
					GlobalMatrix[GV(i,j,k,0)] = LocalMatrix[LVL(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GV(pix,- HP,- HP,0)], (pfx - pix)*(NY + 1 + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

void ParPro::SendMatrixToZeroMW(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LWC(Ix,- HP,- HP,0)], (Fx-Ix)*(NY + 2*HP)*(NZ + 1 + 2*HP), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = - HP; j < NY + HP; j++){
				for(k = - HP; k < NZ + 1 + HP; k++){
					GlobalMatrix[GW(i,j,k,0)] = LocalMatrix[LWL(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[GW(pix,- HP,- HP,0)], (pfx - pix)*(NY + 2*HP)*(NZ + 1 + 2*HP), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
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

void ParPro::SendMatrixToAll_Pressure(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int NZ, int Procesos, int Ix, int Fx){
int i, j, k, p;
int pix, pfx;
MPI_Status ST;

	if(Rank == 0){

		for(i = Ix; i < Fx + Halo; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					LocalMatrix[LPC(i,j,k,0)] = GlobalMatrix[LAG(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos - 1; p++){
			
			Get_Worksplit(NX, Procesos, p, pix, pfx);	
			MPI_Send(&GlobalMatrix[LAG(pix - Halo, 0, 0, 0)], (pfx-pix + 2*Halo)*(NY)*(NZ), MPI_DOUBLE, p, 0, MPI_COMM_WORLD);

		}

		p = Procesos - 1;
		Get_Worksplit(NX, Procesos, p, pix, pfx);	
		MPI_Send(&GlobalMatrix[LAG(pix - Halo, 0, 0, 0)], (pfx-pix + Halo)*(NY)*(NZ), MPI_DOUBLE, p, 0, MPI_COMM_WORLD);

	}
	if(Rank != 0 && Rank != Procesos - 1){
		MPI_Recv(&LocalMatrix[LPC(Ix - Halo, 0, 0, 0)], (Fx - Ix + 2*Halo)*(NY)*(NZ), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &ST);
	}
	else if(Rank == Procesos - 1){
		MPI_Recv(&LocalMatrix[LPC(Ix - Halo, 0, 0, 0)], (Fx - Ix + Halo)*(NY)*(NZ), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &ST);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);	
	
}

void ParPro::Execute(){

	Rango();
	Processes();
	
}
