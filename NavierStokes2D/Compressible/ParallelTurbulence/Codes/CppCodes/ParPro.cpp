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

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ReadData.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ParPro.h"
using namespace std;

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/"

#define G(i,j,dim) (((j) + (i)*NR) + NA*NR*(dim)) //Global Index
#define LNH(i,j,dim) (((j) + ((i) - Ix + Halo)*NR)) //Local Index No Halo 
#define LSH(i,j,dim) (((j) + ((i) - Ix)*NY) + NY*(Fx-Ix + 2*Halo)*dim) //Local Index Si Halo

ParPro::ParPro(ReadData R1){

		NA = R1.NumericalData[4];
		NR = R1.NumericalData[5];
		Halo = 1;

}

void ParPro::Rango(){
	int a;
	a = MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
}

void ParPro::Processes(){
	int a;
	a = MPI_Comm_size(MPI_COMM_WORLD, &Procesos);
}

void ParPro::Initial_WorkSplit(int NA, int &Ix, int &Fx){
int Intervalo, Residuo;
int p;
	Intervalo = NA/Procesos;
	Residuo = NA%Procesos;

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

	Intervalo = NA/Procesos;
	Residuo = NA%Procesos;

	if(p != Procesos-1){
		pix = p*Intervalo;
		pfx = (p+1)*Intervalo;
	}
	else{
		pix = p*Intervalo;
		pfx = (p+1)*Intervalo + Residuo;
	}
}

void ParPro::SendData(double *LocalSend, int Ix, int Fx){
	
	if(Rank != 0 && Rank != Procesos-1){
		MPI_Send(&LocalSend[LNH(Ix,0,0)], Halo*NR, MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(&LocalSend[LNH(Fx - Halo,0,0)], Halo*NR, MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD);
	}
	else if(Rank == 0){ 
		MPI_Send(&LocalSend[LNH(Fx - Halo,0,0)], Halo*NR, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD); 
	}
	else{
		MPI_Send(&LocalSend[LNH(Ix,0,0)], Halo*NR, MPI_DOUBLE, Procesos-2, 0, MPI_COMM_WORLD); 
	}

}

void ParPro::ReceiveData(double *LocalReceive, int Ix, int Fx){
int p;
MPI_Status ST;

	if(Rank != 0 && Rank != Procesos-1){
		MPI_Recv(&LocalReceive[LNH(Ix - Halo,0,0)], Halo*NR, MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD, &ST);
		MPI_Recv(&LocalReceive[LNH(Fx,0,0)], Halo*NR, MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD, &ST);
	}
	else if(Rank == 0){ 
		MPI_Recv(&LocalReceive[LNH(Fx,0,0)], Halo*NR, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &ST);
	}
	else if(Rank == Procesos-1){ 
		MPI_Recv(&LocalReceive[LNH(Ix - Halo,0,0)], Halo*NR, MPI_DOUBLE, Procesos-2, 0, MPI_COMM_WORLD, &ST); 
	}
	

}

void ParPro::SendMatrixToZero(double *LocalMatrix, double *GlobalMatrix, int NX, int NY, int Procesos, int Ix, int Fx){
int i, j, p;
int pix, pfx;
MPI_Status ST;

	if(Rank != 0){
		MPI_Send(&LocalMatrix[LNH(Ix,0,0)], NR*(Fx-Ix), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rank == 0){
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NR; j++){
				GlobalMatrix[G(i,j,0)] = LocalMatrix[LNH(i,j,0)];
			}
		}
		
		for(p = 1; p < Procesos; p++){
			Get_Worksplit(NA, Procesos, p, pix, pfx);	
			MPI_Recv(&GlobalMatrix[G(pix,0,0)], NR*(pfx - pix), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

		
	}

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

void ParPro::Execute(){

	Rango();
	Processes();
	
}
