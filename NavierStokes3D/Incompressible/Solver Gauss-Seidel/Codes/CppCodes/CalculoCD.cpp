#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/ParPro.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/Solver.h"
#include "../HeaderCodes/CalculoCD.h"

#define PI 3.141592653589793

#define GP(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+2*HP)*(NY+2*HP)*(NZ+2*HP)*(dim)) //Global Index P Mesh
#define GU(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+1+2*HP)*(NY+2*HP)*(NZ+2*HP)*(dim)) //Global Index U Mesh
#define GV(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) + ((NX+2*HP)*(NY+1+2*HP)*(NZ+2*HP)*(dim)) //Global Index V Mesh
#define GW(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+2*HP)*(NY+2*HP)*(NZ+1+2*HP)*(dim)) //Global Index W Mesh

//Local Index P Mesh
#define LPL(i,j,k,dim) ((NY)*(NZ))*((i)) + (((j)) + ((k))*(NY)) //Left Core
#define LPC(i,j,k,dim) ((NY)*(NZ))*((i) - Ix + Halo) + (((j)) + ((k))*(NY)) //Center Cores
#define LPR(i,j,k,dim) ((NY)*(NZ))*((i) - Ix + Halo) + (((j)) + ((k))*(NY)) //Right Core

//Local Index U Mesh
#define LUL(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Left Core
#define LUC(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Center Cores
#define LUR(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Right Core

//Local Index V Mesh
#define LVL(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i) - Ix + HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Left Core
#define LVC(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Center Cores
#define LVR(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Right Core

//Local Index W Mesh
#define LWL(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i) + HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Left Core
#define LWC(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Center Cores
#define LWR(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) //Right Core

//Condiciones de Contorno

//Parte Izquierda
#define PLEFT(i,j,k) NY*(k) + (j)
#define ULEFT(i,j,k) NY*(k) + (j) 
#define VLEFT(i,j,k) (NY+1)*(k) + (j)
#define WLEFT(i,j,k) NY*(k) + (j)

//Parte Derecha
#define PRIGHT(i,j,k) NY*(k) + (j)
#define URIGHT(i,j,k) NY*(k) + (j)
#define VRIGHT(i,j,k) (NY+1)*(k) + (j)
#define WRIGHT(i,j,k) NY*(k) + (j)

//Parte Inferior
#define PBOT(i,j,k) NZ*(i - Ix) + k
#define UBOT(i,j,k) NZ*(i - Ix) + k
#define VBOT(i,j,k) NZ*(i - Ix + 1) + k
#define WBOT(i,j,k) (NZ+1)*(i - Ix + 1) + k

//Parte Superior
#define PTOP(i,j,k) NZ*(i - Ix) + k
#define UTOP(i,j,k) NZ*(i - Ix) + k
#define VTOP(i,j,k) NZ*(i - Ix + 1) + k
#define WTOP(i,j,k) (NZ+1)*(i - Ix + 1) + k

//Parte Here
#define PHERE(i,j,k) NY*(i - Ix) + j
#define UHERE(i,j,k) NY*(i - Ix) + j
#define VHERE(i,j,k) (NY+1)*(i - Ix + 1) + j
#define WHERE(i,j,k) NY*(i - Ix + 1) + j
 
//Parte There 
#define PTHERE(i,j,k) NY*(i - Ix) + j 
#define UTHERE(i,j,k) NY*(i - Ix) + j
#define VTHERE(i,j,k) (NY+1)*(i - Ix + 1) + j
#define WTHERE(i,j,k) NY*(i - Ix + 1) + j

#define Interpolacion(CoordObjetivo, Coord1, Valor1, Coord2, Valor2) Valor1 + ((Valor2 - Valor1)/(Coord2 - Coord1))*(CoordObjetivo - Coord1) //Interpolación Lineal

//Constructor del mallador
CalculoCD::CalculoCD(ReadData R1, ParPro MPI1, Mesher MESH, Solver S1){
		
	//Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3]; 
	NZ = R1.ProblemNumericalData[4];

	EsquemaLargo = R1.ConvectiveScheme1;
	EsquemaCorto = R1.ConvectiveScheme2;

	//Parámetros de computación paralela
	Rank = MPI1.Rank;
	Procesos = MPI1.Procesos;
	Ix = MESH.Ix;
	Fx = MESH.Fx;
	Halo = MPI1.Halo;
	
	HP = MESH.HP;
	
	//Datos Físicos del Problema
	Rho = R1.ProblemPhysicalData[0];
	Uref = R1.ProblemPhysicalData[1];
	Reynolds = R1.ProblemPhysicalData[2];

	Rayleigh = R1.ProblemPhysicalData[3];
	Cp = R1.ProblemPhysicalData[4];
	Prandtl = R1.ProblemPhysicalData[5];

	gx = R1.ProblemPhysicalData[6];
	gy = R1.ProblemPhysicalData[7];
	gz = R1.ProblemPhysicalData[8];

	Tleft = R1.ProblemPhysicalData[9]; 
	Tright = R1.ProblemPhysicalData[10]; 

	Tbot = R1.ProblemPhysicalData[11]; 
	Ttop = R1.ProblemPhysicalData[12]; 

	//Cálculos extra para cada problema
	if(Problema == 1){ //Problema Driven Cavity

		mu = (Rho*Uref*Xdominio)/Reynolds;

	}

	else if(Problema == 2){ //Problema DIfferentially Heated

		//Cálculo de To
		if(abs(0.50*(Tleft + Tright)) > 0.0){
			To = abs(0.50*(Tleft + Tright));
			Difference = abs(Tleft - Tright);
		}
		else{
			To = abs(0.50*(Tbot + Ttop));
			Difference = abs(Tbot - Ttop);
		}

		Beta = 1.0/To;
		Producto = (pow(Rho,2)*abs(gy)*pow(Xdominio,3)*Beta*Difference*Prandtl)/Rayleigh;

		mu = sqrt(Producto);

		K = (Cp*mu)/Prandtl;

	}

}

//Cálculo del término difusivo de la velocidad U
void Solver::Get_DiffusiveU(Mesher MESH, Solver S1, double& DiffusiveU, double *UFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx + 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		
		for(i = Ix; i < Fx + 1; i++){

			//Partes Inferior y Superior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Flas de las Esquinas

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}
	
	}
	else if(Rank == 0){

		//Centro
		for(i = Ix + 1; i < Fx + 1; i++){
			for(k = 1; k < NZ-1; k++){
				for(j = 1; j < NY-1; j++){
					DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix + 1; i < Fx + 1; i++){

			//Partes Inferior y Superior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Flas de las Esquinas

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}
	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix; i < Fx; i++){

			//Partes Inferior y Superior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Flas de las Esquinas

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - S1.Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - S1.Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(S1.Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(S1.Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}	
		
}
/*
//Cálculo del término difusivo de la velocidad V
void Solver::Get_DiffusiveV(Mesher MESH, double *VFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		//Partes Here y There
		for(i = Ix; i < Fx; i++){

			for(j = 1; j < NY; j++){

				//Parte Here
				k = 0;

				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*( 
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - Vhere[VHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(Vthere[VTHERE(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}
		}

	}
	else if(Rank == 0){
	
		for(i = Ix + 1; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY; j++){

				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - Vhere[VHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY; j++){
				
				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(Vthere[VTHERE(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

		}

		
		

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - Vleft[VLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}
		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){

			DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - Vleft[VLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - Vhere[VHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}
		
		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){
		
			DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - Vleft[VLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(Vthere[VTHERE(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
		}
	
	}
	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;

			for(j = 1; j < NY; j++){

				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - Vhere[VHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}

			//Parte There
			k = NZ - 1;

			for(j = 1; j < NY; j++){

				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LVC(i+1,j,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(Vthere[VTHERE(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

		}


		//Parte Derecha
		i = NX - 1;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(Vright[VRIGHT(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}
		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){

			DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(Vright[VRIGHT(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LVC(i,j,k+1,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - Vhere[VHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){
			
			DiffusiveV[LVC(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(Vright[VRIGHT(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LVC(i,j+1,k,0)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LVC(i,j,k,0)] - VFIELD[LVC(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(Vthere[VTHERE(i,j,k)] - VFIELD[LVC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LVC(i,j,k,0)] - VLPRES[LVC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}

}

//Cálculo del término difusivo de la velocidad W
void Solver::Get_DiffusiveW(Mesher MESH, double *WFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}

		}
		
	}
	else if(Rank == 0){

		//Centro
		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}

		}

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){
				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - Wleft[WLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ; k++){

			DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - Wleft[WLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){

			DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - Wleft[WLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}

	}
	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LWC(i+1,j,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

			}

		}

		//Parte Derecha
		i = NX - 1;

		//Centro
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){

				DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*( 0.0
											+ MESH.SupMW[GW(i,j,k,1)]*(Wright[WRIGHT(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;

		for(k = 1; k < NZ; k++){

			DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(Wright[WRIGHT(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LWC(i,j+1,k,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){
			
			DiffusiveW[LWC(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(Wright[WRIGHT(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LWC(i,j,k+1,0)] - WFIELD[LWC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LWC(i,j,k,0)] - WFIELD[LWC(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}
	
	}

}

//Cálculo del término difusivo de la temperatura T
void Solver::Get_DiffusiveT(Mesher MESH, double *TFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix; i < Fx; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){
				
				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}
		
	}
	else if(Rank == 0){

		//Centro
		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix + 1; i < Fx; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}


		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){
				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}
		}

		//Partes Superior e Inferior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}


		//Esquina Abajo Here
		j = 0;
		k = 0;

		DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
								+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
								- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
								+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
								- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
								+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
								- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
								);

		//Esquina Abajo There
		j = 0;
		k = NZ - 1;

		DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
								+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
								- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
								+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
								- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
								+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
								- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
								);

		//Esquina Arriba Here
		j = NY - 1;
		k = 0;

		DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
								+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
								- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
								+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
								- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
								+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
								- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
								);

		//Esquina Arriba There
		j = NY - 1;
		k = NZ - 1;

		DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
								+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
								- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
								+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
								- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
								+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
								- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
								);

	}

	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				}
			}
		}

		for(i = Ix; i < Fx - 1; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){
				
				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
		}

		//Parte Derecha
		i = NX - 1;

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				}
			}

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TFIELD[LPC(i,j+1,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TFIELD[LPC(i,j,k+1,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TFIELD[LPC(i+1,j,k,0)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TFIELD[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TFIELD[LPC(i,j,k,0)] - TFIELD[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);


	}

}

*/