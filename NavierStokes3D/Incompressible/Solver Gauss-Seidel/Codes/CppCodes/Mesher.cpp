#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <cassert>
#include <time.h>
#include <bits/stdc++.h>
#include <string> 
#include <mpi.h>

using namespace std;

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/ParPro.h"
#include "../HeaderCodes/Mesher.h"

#define PI 3.141592653589793

//Global Index P Mesh
#define GP(i,j,k,Dim) (NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*((i) + HaloPressure) + ((j) + HaloPressure) + ((k) + HaloPressure)*(NY + 2*HaloPressure) + (NX + 2*HaloPressure)*(NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*(Dim)

//Global Index U Mesh
#define GU(i,j,k,Dim) (NY + 2*HaloU)*(NZ + 2*HaloU)*((i) + HaloU) + ((j) + HaloU) + ((k) + HaloU)*(NY + 2*HaloU) + (NX + 1 + 2*HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU)*(Dim)

//Global Index V Mesh
#define GV(i,j,k,Dim) (NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*((i) + HaloV) + ((j) + HaloV) + ((k) + HaloV)*(NY + 1 + 2*HaloV) + (NX + 2*HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*(Dim)

//Global Index W Mesh
#define GW(i,j,k,Dim) (NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*((i) + HaloW) + ((j) + HaloW) + ((k) + HaloW)*(NY + 2*HaloW) + (NX + 2*HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*(Dim)

//Constructor del mallador
Mesher::Mesher(Memory M1, ReadData R1, ParPro MPI1, string InputDirectorio){
		
	//Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	ProgresionLineal = R1.ProblemNumericalData[1];

	NZ = R1.ProblemNumericalData[4];

	OptionX = R1.ProblemNumericalData[5];
	OptionY = R1.ProblemNumericalData[6];
	OptionZ = R1.ProblemNumericalData[7];

	SFX = R1.ProblemData[0];
	SFY = R1.ProblemData[1];
	SFZ = R1.ProblemData[2];

	//Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	Xcentroide = R1.GeometryData[3];
	Ycentroide = R1.GeometryData[4];

	Xcuadrado = R1.GeometryData[5];
	Ycuadrado = R1.GeometryData[6];

	DIRECTORIO = InputDirectorio;
	
	//Datos necesarios para computación paralela
	Rank = MPI1.Rank;
	Procesos = MPI1.Procesos;

	HaloPressure = 1;
	HaloU = 2;
	HaloV = 2;
	HaloW = 2;

	HP = 2;

	if(Problema == 1 || Problema == 2){
		NX = R1.ProblemNumericalData[2];
		NY = R1.ProblemNumericalData[3]; 
	}
	else if(Problema == 3){ //Problema Square Cylinder

		NX2 = 80;
		NY2 = 80;

		if(OptionX == 1){
			NX1 = (Xcentroide - 0.50*Xcuadrado)/(Xcuadrado/NX2);
			NX3 = (Xdominio - (Xcentroide + 0.50*Xcuadrado))/(Xcuadrado/NX2);
		}
		else if(OptionX == 2){
			NX1 = - SFX/( - SFX + atanh( tanh(SFX)*  ( 1.0 - Xcuadrado/(NX2*(Xcentroide - 0.50*Xcuadrado)))   ));
			NX3 = SFX/( SFX + atanh( tanh(SFX)*  (Xcuadrado/(NX2*(Xdominio - (Xcentroide + 0.50*Xcuadrado))) - 1)   ));
		}

		if(OptionY == 1){
			NY1 = (Ycentroide - 0.50*Ycuadrado)/(Ycuadrado/NX2);
			NY3 = (Ydominio - (Ycentroide + 0.50*Ycuadrado))/(Ycuadrado/NX2);
		}
		else if(OptionY == 2){
			NY1 = - SFY/( - SFY + atanh( tanh(SFY)*  ( 1.0 - Ycuadrado/(NY2*(Ycentroide - 0.50*Ycuadrado)))   ));
			NY3 = SFY/( SFY + atanh( tanh(SFY)*  (Ycuadrado/(NY2*(Ydominio - (Ycentroide + 0.50*Ycuadrado))) - 1)   ));
		}

		NX = NX1 + NX2 + NX3;
		NY = NY1 + NY2 + NY3;

	}

}

//Alojamiento de memoria para cada matriz
void Mesher::AllocateMemory(Memory M1){

	//Matrices del Mallado del problema
	MP = M1.AllocateDouble(NX + 2*HaloPressure, NY + 2*HaloPressure, NZ + 2*HaloPressure, 3); //Coordenadas matriz de presión/temperatura
	MU = M1.AllocateDouble(NX + 1 + 2*HaloU, NY + 2*HaloU, NZ + 2*HaloU, 3); //Coordenadas matriz velocidad U
	MV = M1.AllocateDouble(NX + 2*HaloV, NY + 1 + 2*HaloV, NZ + 2*HaloV, 3); //Coordenadas matriz velocidad V
	MW = M1.AllocateDouble(NX + 2*HaloW, NY + 2*HaloW, NZ + 1 + 2*HaloW, 3); //Coordenadas matriz velocidad W

	//Código:
	//0 -> Coordenada X
	//1 -> Coordenada Y
	//2 -> Coordenada Z
	
	//Matrices de distancias de volúmenes de control
	DeltasMP = M1.AllocateDouble(NX + 2*HaloPressure, NY + 2*HaloPressure, NZ + 2*HaloPressure, 3); //Deltas de la matriz de Presión/Temperatura
	DeltasMU = M1.AllocateDouble(NX + 1 + 2*HaloU, NY + 2*HaloU, NZ + 2*HaloU, 3); //Deltas de la matriz de velocidades (U)
	DeltasMV = M1.AllocateDouble(NX + 2*HaloV, NY + 1 + 2*HaloV, NZ + 2*HaloV, 3); //Deltas de la matriz de velocidades (V)
	DeltasMW = M1.AllocateDouble(NX + 2*HaloW, NY + 2*HaloW, NZ + 1 + 2*HaloW, 3); //Deltas de la matriz de velocidades (W)

	//Código:
	//0 -> Delta X
	//1 -> Delta Y
	//2 -> Delta Z

	//Matrices de superficies de volúmenes de control
	SupMP = M1.AllocateDouble(NX + 2*HaloPressure, NY + 2*HaloPressure, NZ + 2*HaloPressure, 6); //Superficies de los volúmenes de la matriz de Presión/Temperatura
	SupMU = M1.AllocateDouble(NX + 1 + 2*HaloU, NY + 2*HaloU, NZ + 2*HaloU, 6); //Superficies de los volúmenes de la matriz de velocidades (U)
	SupMV = M1.AllocateDouble(NX + 2*HaloV, NY + 1 + 2*HaloV, NZ + 2*HaloV, 6); //Superficies de los volúmenes de la matriz de velocidades (V)
	SupMW = M1.AllocateDouble(NX + 2*HaloW, NY + 2*HaloW, NZ + 1 + 2*HaloW, 6); //Superficies de los volúmenes de la matriz de velocidades (V)

	//Código:
	//0 -> West
	//1 -> East
	//2 -> South
	//3 -> North
	//4 -> Here
	//5 -> There

	//Matrices de volúmenes de los volúmenes de control
	VolMP = M1.AllocateDouble(NX + 2*HaloPressure, NY + 2*HaloPressure, NZ + 2*HaloPressure, 1); //Volúmenes de los volúmenes de la matriz de Presión/Temperatura
	VolMU = M1.AllocateDouble(NX + 1 + 2*HaloU, NY + 2*HaloU, NZ + 2*HaloU, 1); //Volúmenes de los volúmenes de la matriz de velocidades (U)
	VolMV = M1.AllocateDouble(NX + 2*HaloV, NY + 1 + 2*HaloV, NZ + 2*HaloV, 1); //Volúmenes de los volúmenes de la matriz de velocidades (V)
	VolMW = M1.AllocateDouble(NX + 2*HaloW, NY + 2*HaloW, NZ + 1 + 2*HaloW, 1); //Volúmenes de los volúmenes de la matriz de velocidades (W)

}

//Creación de todas las mallas de tipo Staggered
void Mesher::Get_Meshes(){
int i, j, k;
double I, J, K;
double nx = NX;
double ny = NY;
double nz = NZ;

	//Código:
	//0 -> Coordenada X
	//1 -> Coordenada Y
	//2 -> Coordenada Z
	
	//COORDENADAS X

	if(Problema == 1 || Problema == 2){

		if(OptionX == 1){ //Discretización Regular

			//Malla Velocidades U
			for(i = - HaloU; i < NX + 1 + HaloU; i++){
				for(j = - HaloU; j < NY + HaloU; j++){
					for(k = - HaloU; k < NZ + HaloU; k++){
						MU[GU(i,j,k,0)] = i*(Xdominio/nx);
					}
				}
			}

		}
		else if(OptionX == 2){ //Discretización Tangencial Hiperbólica

			//Malla Velocidades U
			for(i = - HaloU; i < NX + 1 + HaloU; i++){
				for(j = - HaloU; j < NY + HaloU; j++){
					for(k = - HaloU; k < NZ + HaloU; k++){
						I = i;
						MU[GU(i,j,k,0)] = (Xdominio/2.0)*(1.0 + (tanh(SFX*((2.0*I - nx)/nx)) + tanh(SFX))/tanh(SFX) - 1.0);
					}
				}
			}

		}

		//Matriz de Presiones P
		for(i = - HaloPressure; i < NX + HaloPressure; i++){
			for(j = - HaloPressure; j < NY + HaloPressure; j++){
				for(k = - HaloPressure; k < NZ + HaloPressure; k++){
					MP[GP(i,j,k,0)] = 0.50*(MU[GU(i,j,k,0)] + MU[GU(i+1,j,k,0)]);
				}
			}
		}

		//Matriz de Velocidades V
		for(i = - HaloPressure; i < NX + HaloPressure; i++){
			for(k = - HaloPressure; k < NZ + HaloPressure; k++){
				MV[GV(i,NY + HaloV,k,0)] = 0.50*(MU[GU(i,NY - 1 + HaloU,k,0)] + MU[GU(i + 1,NY-1 + HaloU,k,0)]);
				for(j = - HaloPressure; j < NY + HaloPressure; j++){
					MV[GV(i,j,k,0)] = 0.50*(MU[GU(i,j,k,0)] + MU[GU(i+1,j,k,0)]);
				}
			}
		}

		//Matriz de Velocidades W
		for(i = - HaloPressure; i < NX + HaloPressure; i++){
			for(j = - HaloPressure; j < NY + HaloPressure; j++){
				MW[GW(i,j,NZ + HaloPressure,0)] = 0.50*(MU[GU(i,j,NZ-1 + HaloPressure,0)] + MU[GU(i+1,j,NZ-1 + HaloPressure,0)]);
				for(k = - HaloPressure; k < NZ + HaloPressure; k++){
					MW[GW(i,j,k,0)] = 0.50*(MU[GU(i,j,k,0)] + MU[GU(i+1,j,k,0)]);
				}
			}
		}

	}
	else if(Problema == 3){

		//Coordenadas X

		double nx1 = NX1;
		double nx2 = NX2;
		double nx3 = NX3;

		if(OptionX == 1){ //Discretizacion regular

			//Primera parte
			for(i = - HP; i < NX1; i++){
				for(j = - HP; j < NY + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						MU[GU(i,j,k,0)] = i*((Xcentroide - 0.50*Xcuadrado)/nx1);
					}
				}
			}

			//Zona centro
			for(i = 0; i < NX2; i++){
				for(j = - HP; j < NY + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						MU[GU(NX1 + i,j,k,0)] = (Xcentroide - 0.50*Xcuadrado) + i*((Xcuadrado)/nx2);
					}
				}
			}

			//Tercera parte
			for(i = 0; i < NX3 + 1 + HP; i++){
				for(j = - HP; j < NY + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						MU[GU(NX1 + NX2 + i,j,k,0)] = (Xcentroide + 0.50*Xcuadrado) + i*((Xdominio - (Xcentroide + 0.50*Xcuadrado))/nx3);
					}
				}
			}

		}
		else if(OptionX == 2){ //Discretizacion tangencial hiperbolica

			//Primera parte
			for(i = - HP; i < NX1; i++){
				for(j = - HP; j < NY + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						I = i;
						MU[GU(i,j,k,0)] = (Xcentroide - 0.50*Xcuadrado)*((tanh(SFX*((I)/nx1)))/tanh(SFX));
					}
				}
			}

			//Zona centro
			for(i = 0; i < NX2; i++){
				for(j = - HP; j < NY + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						MU[GU(NX1 + i,j,k,0)] = (Xcentroide - 0.50*Xcuadrado) + i*((Xcuadrado)/nx2);
					}
				}
			}

			//Tercera parte
			for(i = 0; i < NX3 + 1 + HP; i++){
				for(j = - HP; j < NY + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						I = i;
						MU[GU(NX1 + NX2 + i,j,k,0)] = (Xcentroide + 0.50*Xcuadrado) + (Xdominio - (Xcentroide + 0.50*Xcuadrado))*((tanh(SFX*((I - nx3)/nx3)) + tanh(SFX))/tanh(SFX));
					}
				}
			}

		}
		
		//Matriz de Presiones P
		for(i = - HP; i < NX + HP; i++){
			for(j = - HP; j < NY + HP; j++){
				for(k = - HP; k < NZ + HP; k++){
					MP[GP(i,j,k,0)] = 0.50*(MU[GU(i,j,k,0)] + MU[GU(i+1,j,k,0)]);
				}
			}
		}

		//Matriz de Velocidades V
		for(i = - HP; i < NX + HP; i++){
			for(k = - HP; k < NZ + HP; k++){
				MV[GV(i,NY+HP,k,0)] = 0.50*(MU[GU(i,NY-1 + HP,k,0)] + MU[GU(i+1,NY-1 + HP,k,0)]);
				for(j = - HP; j < NY + HP; j++){
					MV[GV(i,j,k,0)] = 0.50*(MU[GU(i,j,k,0)] + MU[GU(i+1,j,k,0)]);
				}
			}
		}

		//Matriz de Velocidades W
		for(i = - HP; i < NX + HP; i++){
			for(j = - HP; j < NY + HP; j++){
				MW[GW(i,j,NZ + HP,0)] = 0.50*(MU[GU(i,j,NZ-1 + HP,0)] + MU[GU(i+1,j,NZ-1 + HP,0)]);
				for(k = - HP; k < NZ + HP; k++){
					MW[GW(i,j,k,0)] = 0.50*(MU[GU(i,j,k,0)] + MU[GU(i+1,j,k,0)]);
				}
			}
		}

	}

	//COORDENADAS Y

	if(Problema == 1 || Problema == 2){

		if(OptionY == 1){ //Discretización Regular

			//Malla Velocidades V
			for(i = - HaloV; i < NX + HaloV; i++){
				for(j = - HaloV; j < NY + 1 + HaloV; j++){
					for(k = - HaloV; k < NZ + HaloV; k++){
						MV[GV(i,j,k,1)] = j*(Ydominio/ny);
					}
				}
			}

		}
		else if(OptionY == 2){ //Discretización Tangencial Hiperbólica

			//Malla Velocidades V
			for(i = - HaloV; i < NX + HaloV; i++){
				for(j = - HaloV; j < NY + 1 + HaloV; j++){
					for(k = - HaloV; k < NZ + HaloV; k++){
						J = j;
						MV[GV(i,j,k,1)] = (Ydominio/2.0)*(1.0 + (tanh(SFY*((2.0*J - ny)/ny)) + tanh(SFY))/tanh(SFY) - 1.0);
					}
				}
			}

		}

		//Matriz de Presiones P
		for(i = - HaloPressure; i < NX + HaloPressure; i++){
			for(j = - HaloPressure; j < NY + HaloPressure; j++){
				for(k = - HaloPressure; k < NZ + HaloPressure; k++){
					MP[GP(i,j,k,1)] = 0.50*(MV[GV(i,j,k,1)] + MV[GV(i,j+1,k,1)]);
				}
			}
		}

		//Matriz de Velocidades U
		for(k = - HaloU; k < NZ + HaloU; k++){
			for(j = - HaloU; j < NY + HaloU; j++){
				MU[GU(NX + HaloU,j,k,1)] = 0.50*(MV[GV(NX-1 + HaloV,j,k,1)] + MV[GV(NX-1 + HaloV,j+1,k,1)]);
				for(i = - HaloU; i < NX + HaloU; i++){
					MU[GU(i,j,k,1)] = 0.50*(MV[GV(i,j,k,1)] + MV[GV(i,j+1,k,1)]);
				}
			}
		}

		//Matriz de Velocidades W
		for(i = - HaloW; i < NX + HaloW; i++){
			for(j = - HaloW; j < NY + HaloW; j++){
				MW[GW(i,j,NZ + HaloW,1)] = 0.50*(MV[GV(i,j,NZ-1 + HaloV,1)] + MV[GV(i,j+1,NZ-1 + HaloV,1)]);
				for(k = - HP; k < NZ + HaloW; k++){
					MW[GW(i,j,k,1)] = 0.50*(MV[GV(i,j,k,1)] + MV[GV(i,j+1,k,1)]);
				}
			}
		}

	}
	else if(Problema == 3){

		//Coordenadas Y

		double ny1 = NY1;
		double ny2 = NY2;
		double ny3 = NY3;

		if(OptionY == 1){ //Discretización Regular

			//Parte Inferior
			for(i = - HP; i < NX + HP; i++){
				for(j = - HP; j < NY1; j++){
					for(k = - HP; k < NZ + HP; k++){
						MV[GV(i,j,k,1)] = j*((Ycentroide - 0.50*Ycuadrado)/ny1);
					}
				}
			}

			//Parte Central
			for(i = - HP; i < NX + HP; i++){
				for(j = 0; j < NY2; j++){
					for(k = - HP; k < NZ + HP; k++){
						MV[GV(i,NY1 + j,k,1)] = (Ycentroide - 0.50*Ycuadrado) + j*((Ycuadrado)/ny2);
					}
				}
			}

			//Parte Superior
			for(i = - HP; i < NX + HP; i++){
				for(j = 0; j < NY3 + 1 + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						MV[GV(i,NY1 + NY2 + j,k,1)] = (Ycentroide + 0.50*Ycuadrado) + j*((Ydominio - (Ycentroide + 0.50*Ycuadrado))/ny3);
					}
				}
			}

		}
		else if(OptionY == 2){ //Discretización Tangencial Hiperbólica

			//Parte Inferior
			for(i = - HP; i < NX + HP; i++){
				for(j = - HP; j < NY1; j++){
					for(k = - HP; k < NZ + HP; k++){
						J = j;
						MV[GV(i,j,k,1)] = (Ycentroide - 0.50*Ycuadrado)*((tanh(SFY*((J)/ny1)))/tanh(SFY));
					}
				}
			}

			//Parte Central
			for(i = - HP; i < NX + HP; i++){
				for(j = 0; j < NY2; j++){
					for(k = - HP; k < NZ + HP; k++){
						MV[GV(i,NY1 + j,k,1)] = (Ycentroide - 0.50*Ycuadrado) + j*((Ycuadrado)/ny2);
					}
				}
			}

			//Parte Superior
			for(i = - HP; i < NX + HP; i++){
				for(j = 0; j < NY3 + 1 + HP; j++){
					for(k = - HP; k < NZ + HP; k++){
						J = j;
						MV[GV(i,NY1 + NY2 + j,k,1)] = (Ycentroide + 0.50*Ycuadrado) + (Ydominio - (Ycentroide + 0.50*Ycuadrado))*((tanh(SFY*((J - ny3)/ny3)) + tanh(SFY))/tanh(SFY));
					}
				}
			}

		}

		//Matriz de Presiones P
		for(i = - HP; i < NX + HP; i++){
			for(j = - HP; j < NY + HP; j++){
				for(k = - HP; k < NZ + HP; k++){
					MP[GP(i,j,k,1)] = 0.50*(MV[GV(i,j,k,1)] + MV[GV(i,j+1,k,1)]);
				}
			}
		}

		//Matriz de Velocidades U
		for(k = - HP; k < NZ + HP; k++){
			for(j = - HP; j < NY + HP; j++){
				MU[GU(NX + HP,j,k,1)] = 0.50*(MV[GV(NX-1 + HP,j,k,1)] + MV[GV(NX-1 + HP,j+1,k,1)]);
				for(i = - HP; i < NX + HP; i++){
					MU[GU(i,j,k,1)] = 0.50*(MV[GV(i,j,k,1)] + MV[GV(i,j+1,k,1)]);
				}
			}
		}

		//Matriz de Velocidades W
		for(i = - HP; i < NX + HP; i++){
			for(j = - HP; j < NY + HP; j++){
				MW[GW(i,j,NZ + HP,1)] = 0.50*(MV[GV(i,j,NZ-1 + HP,1)] + MV[GV(i,j+1,NZ-1 + HP,1)]);
				for(k = - HP; k < NZ + HP; k++){
					MW[GW(i,j,k,1)] = 0.50*(MV[GV(i,j,k,1)] + MV[GV(i,j+1,k,1)]);
				}
			}
		}

	}

	//COORDENADAS Z

		if(OptionZ == 1){ //Discretización Regular

			//Malla Velocidades W
			for(i = - HaloW; i < NX + HaloW; i++){
				for(j = - HaloW; j < NY + HaloW; j++){
					for(k = - HaloW; k < NZ + 1 + HaloW; k++){
						MW[GW(i,j,k,2)] = k*(Zdominio/NZ);
					}
				}
			}

		}
		else if(OptionZ == 2){ //Discretización Tangencial Hiperbólica

			//Malla Velocidades W
			for(i = - HaloW; i < NX + HaloW; i++){
				for(j = - HaloW; j < NY + HaloW; j++){
					for(k = - HaloW; k < NZ + 1 + HaloW; k++){
						K = k;
						MW[GW(i,j,k,2)] = (Zdominio/2.0)*(1.0 + (tanh(SFZ*((2.0*K - nz)/nz)) + tanh(SFZ))/tanh(SFZ) - 1.0);
					}
				}
			}

		}

		//Matriz de Presiones P
		for(i = - HaloPressure; i < NX + HaloPressure; i++){
			for(j = - HaloPressure; j < NY + HaloPressure; j++){
				for(k = - HaloPressure; k < NZ + HaloPressure; k++){
					MP[GP(i,j,k,2)] = 0.50*(MW[GW(i,j,k,2)] + MW[GW(i,j,k+1,2)]);
				}
			}
		}

		//Matriz de Velocidades U
		for(k = - HaloU; k < NZ + HaloU; k++){
			for(j = - HaloU; j < NY + HaloU; j++){
				MU[GU(NX + HaloU,j,k,2)] = 0.50*(MW[GW(NX - 1 + HaloW,j,k,2)] + MW[GW(NX - 1 + HaloW,j,k+1,2)]);
				for(i = - HaloU; i < NX + HaloU; i++){
					MU[GU(i,j,k,2)] = 0.50*(MW[GW(i,j,k,2)] + MW[GW(i,j,k+1,2)]);
				}
			}
		}

		//Matriz de Velocidades V
		for(i = - HaloV; i < NX + HaloV; i++){
			for(k = - HaloV; k < NZ + HaloV; k++){
				MV[GV(i,NY + HaloV,k,2)] = 0.50*(MW[GW(i,NY-1 + HaloV,k,2)] + MW[GW(i,NY-1 + HaloV,k+1,2)]);
				for(j = - HP; j < NY + HP; j++){
					MV[GV(i,j,k,2)] = 0.50*(MW[GW(i,j,k,2)] + MW[GW(i,j,k+1,2)]);
				}
			}
		}

}

//Pasar los resultados de las mallas a un txt
void Mesher::PrintTxt(){
int i, j, k;
string Carpeta = "GnuPlotResults/Meshes/";
ofstream file;
string FileName;
stringstream InitialNameMP;
string FinalNameMP;

	FileName = "MallaMP.txt";

	InitialNameMP<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMP = InitialNameMP.str();
    file.open(FinalNameMP.c_str());

	for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				file<<MP[GP(i,j,k,0)]<<"\t"<<MP[GP(i,j,k,1)]<<"\t"<<MP[GP(i,j,k,2)]<<endl;
			}
			file<<endl;
		}   	
    }

	file.close();

stringstream InitialNameMU;
string FinalNameMU;

	FileName = "MallaMU.txt";

	InitialNameMU<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMU = InitialNameMU.str();
    file.open(FinalNameMU.c_str());

	for(i = 0; i < NX + 1; i++){
        for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				file<<MU[GU(i,j,k,0)]<<"\t"<<MU[GU(i,j,k,1)]<<"\t"<<MU[GU(i,j,k,2)]<<endl;
			}
			file<<endl;
		}   	
    }

	file.close();

stringstream InitialNameMV;
string FinalNameMV;

	FileName = "MallaMV.txt";

	InitialNameMV<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMV = InitialNameMV.str();
    file.open(FinalNameMV.c_str());

	for(i = 0; i < NX; i++){
        for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){
				file<<MV[GV(i,j,k,0)]<<"\t"<<MV[GV(i,j,k,1)]<<"\t"<<MV[GV(i,j,k,2)]<<endl;
			}
			file<<endl;
		}   	
    }

	file.close();

stringstream InitialNameMW;
string FinalNameMW;

	FileName = "MallaMW.txt";

	InitialNameMW<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMW = InitialNameMW.str();
    file.open(FinalNameMW.c_str());

	for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){
				file<<MW[GW(i,j,k,0)]<<"\t"<<MW[GW(i,j,k,1)]<<"\t"<<MW[GW(i,j,k,2)]<<endl;
			}
			file<<endl;
		}   	
    }

	file.close();

stringstream InitialNameVolMV;
string FinalNameVolMV;

	FileName = "VolMV.txt";

	InitialNameVolMV<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameVolMV = InitialNameVolMV.str();
    file.open(FinalNameVolMV.c_str());

	for(i = 0; i < NX; i++){
        for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){
				file<<"I: "<<i<<", J: "<<j<<", K: "<<k<<", SupN: "<<SupMV[GV(i,j,k,3)]<<", SupH: "<<SupMV[GV(i,j,k,4)]<<", SupT: "<<SupMV[GV(i,j,k,5)]<<endl;
			}
			file<<endl;
		}   	
    }

	file.close();

}

//Cálculo de las distancias entre nodos en cada una de las matrices
void Mesher::Get_Deltas(){
int i, j, k;

	//Código:
	//0 -> Delta X
	//1 -> Delta Y
	//2 -> Delta Z

	//Deltas de la Matriz P
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				DeltasMP[GP(i,j,k,0)] = MU[GU(i+1,j,k,0)] - MU[GU(i,j,k,0)]; //Delta X
				DeltasMP[GP(i,j,k,1)] = MV[GV(i,j+1,k,1)] - MV[GV(i,j,k,1)]; //Delta Y
				DeltasMP[GP(i,j,k,2)] = MW[GW(i,j,k+1,2)] - MW[GW(i,j,k,2)]; //Delta Z
			}
		}
	}

	//Deltas de la Matriz U
	for(j = 0; j < NY; j++){
		for(k = 0; k < NZ; k++){

			DeltasMU[GU(0,j,k,0)] = MP[GP(0,j,k,0)] - MU[GU(0,j,k,0)]; //Delta X Parte Izquierda
			DeltasMU[GU(NX,j,k,0)] = MU[GU(NX,j,k,0)] - MP[GP(NX-1,j,k,0)]; //Delta X Parte Derecha

			DeltasMU[GU(0,j,k,1)] = MV[GV(0,j+1,k,1)] - MV[GV(0,j,k,1)]; //Delta Y Parte Izquierda
			DeltasMU[GU(NX,j,k,1)] = MV[GV(NX-1,j+1,k,1)] - MV[GV(NX-1,j,k,1)]; //Delta Y Parte Derecha

			DeltasMU[GU(0,j,k,2)] = MW[GW(0,j,k+1,2)] - MW[GW(0,j,k,2)]; //Delta Z Parte Izquierda
			DeltasMU[GU(NX,j,k,2)] = MW[GW(NX-1,j,k+1,2)] - MW[GW(NX-1,j,k,2)]; //Delta Z Parte Derecha

			for(i = 1; i < NX; i++){
				DeltasMU[GU(i,j,k,0)] = MP[GP(i,j,k,0)] - MP[GP(i-1,j,k,0)]; //Delta X
				DeltasMU[GU(i,j,k,1)] = MV[GV(i,j+1,k,1)] - MV[GV(i,j,k,1)]; //Delta Y
				DeltasMU[GU(i,j,k,2)] = MW[GW(i,j,k+1,2)] - MW[GW(i,j,k,2)]; //Delta Z
			}

		}
	}

	//Deltas de la Matriz V
	for(i = 0; i < NX; i++){
		for(k = 0; k < NZ; k++){

			DeltasMV[GV(i,0,k,0)] = MU[GU(i+1,0,k,0)] - MU[GU(i,0,k,0)]; //Delta X Parte Abajo
			DeltasMV[GV(i,NY,k,0)] = MU[GU(i+1,NY-1,k,0)] - MU[GU(i,NY-1,k,0)]; //Delta X Parte Arriba

			DeltasMV[GV(i,0,k,1)] = MP[GP(i,0,k,1)] - MV[GV(i,0,k,1)]; //Delta Y Parte Abajo
			DeltasMV[GV(i,NY,k,1)] = MV[GV(i,NY,k,1)] - MP[GP(i,NY-1,k,1)]; //Delta Y Parte Arriba

			DeltasMV[GV(i,0,k,2)] = MW[GW(i,0,k+1,2)] - MW[GW(i,0,k,2)]; //Delta Z Parte Abajo
			DeltasMV[GV(i,NY,k,2)] = MW[GW(i,NY-1,k+1,2)] - MW[GW(i,NY-1,k,2)]; //Delta Z Parte Arriba

			for(j = 1; j < NY; j++){
				DeltasMV[GV(i,j,k,0)] = MU[GU(i+1,j,k,0)] - MU[GU(i,j,k,0)]; //Delta X
				DeltasMV[GV(i,j,k,1)] = MP[GP(i,j,k,1)] - MP[GP(i,j-1,k,1)]; //Delta Y
				DeltasMV[GV(i,j,k,2)] = MW[GW(i,j,k+1,2)] - MW[GW(i,j,k,2)]; //Delta Z
			}

		}
	}

	//Deltas de la Matriz W
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){

			DeltasMW[GW(i,j,0,0)] = MU[GU(i+1,j,0,0)] - MU[GU(i,j,0,0)]; //Delta X Parte Cercana
			DeltasMW[GW(i,j,NZ,0)] = MU[GU(i+1,j,NZ-1,0)] - MU[GU(i,j,NZ-1,0)]; //Delta X Parte Lejana

			DeltasMW[GW(i,j,0,1)] = MV[GV(i,j+1,0,1)] - MV[GV(i,j,0,1)]; //Delta Y Parte Cercana
			DeltasMW[GW(i,j,NZ,1)] = MV[GV(i,j+1,NZ,1)] - MV[GV(i,j,NZ,1)]; //Delta Y Parte Lejana

			DeltasMW[GW(i,j,0,2)] = MP[GP(i,j,0,2)] - MW[GW(i,j,0,2)]; //Delta Z Parte Cercana
			DeltasMW[GW(i,j,NZ,2)] = MW[GW(i,j,NZ,2)] - MP[GP(i,j,NZ-1,2)]; //Delta Z Parte Lejana

			for(k = 1; k < NZ; k++){
				DeltasMW[GW(i,j,k,0)] = MU[GU(i+1,j,k,0)] - MU[GU(i,j,k,0)]; //Delta X
				DeltasMW[GW(i,j,k,1)] = MV[GV(i,j+1,k,1)] - MV[GV(i,j,k,1)]; //Delta Y
				DeltasMW[GW(i,j,k,2)] = MP[GP(i,j,k,2)] - MP[GP(i,j,k-1,2)]; //Delta Z
			}
		}
	}

}

//Cálculo de las superficies de cada uno de los volúmenes de control
void Mesher::Get_Surfaces(){
int i, j, k;

	//Código:
	//0 -> West
	//1 -> East
	//2 -> South
	//3 -> North
	//4 -> Here
	//5 -> There

	//Superficies Volumen Matriz P
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				SupMP[GP(i,j,k,0)] = DeltasMU[GU(i,j,k,1)]*DeltasMU[GU(i,j,k,2)]; //West
				SupMP[GP(i,j,k,1)] = DeltasMU[GU(i+1,j,k,1)]*DeltasMU[GU(i+1,j,k,2)]; //East

				SupMP[GP(i,j,k,2)] = DeltasMV[GV(i,j,k,0)]*DeltasMV[GV(i,j,k,2)]; //South
				SupMP[GP(i,j,k,3)] = DeltasMV[GV(i,j+1,k,0)]*DeltasMV[GV(i,j+1,k,2)]; //North

				SupMP[GP(i,j,k,4)] = DeltasMW[GW(i,j,k,0)]*DeltasMW[GW(i,j,k,1)]; //Here
				SupMP[GP(i,j,k,5)] = DeltasMW[GW(i,j,k+1,0)]*DeltasMW[GW(i,j,k+1,1)]; //There
			}
		}
	}

	//Superficies Volumen Matriz U
	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){

			//Parte Izquierda
			SupMU[GU(0,j,k,0)] = DeltasMU[GU(0,j,k,1)]*DeltasMU[GU(0,j,k,2)]; //West
			SupMU[GU(0,j,k,1)] = DeltasMP[GP(0,j,k,1)]*DeltasMP[GP(0,j,k,2)]; //East

			SupMU[GU(0,j,k,2)] = DeltasMU[GU(0,j,k,0)]*DeltasMU[GU(0,j,k,2)]; //South
			SupMU[GU(0,j,k,3)] = DeltasMU[GU(0,j,k,0)]*DeltasMU[GU(0,j,k,2)]; //North

			SupMU[GU(0,j,k,4)] = DeltasMU[GU(0,j,k,0)]*DeltasMU[GU(0,j,k,1)]; //Here
			SupMU[GU(0,j,k,5)] = DeltasMU[GU(0,j,k,0)]*DeltasMU[GU(0,j,k,1)]; //There

			//Parte Derecha
			SupMU[GU(NX,j,k,0)] = DeltasMP[GP(NX-1,j,k,1)]*DeltasMP[GP(NX-1,j,k,2)]; //West
			SupMU[GU(NX,j,k,1)] = DeltasMU[GU(NX,j,k,1)]*DeltasMU[GU(NX,j,k,2)]; //East

			SupMU[GU(NX,j,k,2)] = DeltasMU[GU(NX,j,k,0)]*DeltasMU[GU(NX,j,k,2)]; //South
			SupMU[GU(NX,j,k,3)] = DeltasMU[GU(NX,j,k,0)]*DeltasMU[GU(NX,j,k,2)]; //North

			SupMU[GU(NX,j,k,4)] = DeltasMU[GU(NX,j,k,0)]*DeltasMU[GU(NX,j,k,1)]; //Here
			SupMU[GU(NX,j,k,5)] = DeltasMU[GU(NX,j,k,0)]*DeltasMU[GU(NX,j,k,1)]; //There

			for(i = 1; i < NX; i++){
				SupMU[GU(i,j,k,0)] = DeltasMP[GP(i-1,j,k,1)]*DeltasMP[GP(i-1,j,k,2)]; //West
				SupMU[GU(i,j,k,1)] = DeltasMP[GP(i,j,k,1)]*DeltasMP[GP(i,j,k,2)]; //East

				SupMU[GU(i,j,k,2)] = DeltasMU[GU(i,j,k,0)]*DeltasMU[GU(i,j,k,2)]; //South
				SupMU[GU(i,j,k,3)] = DeltasMU[GU(i,j,k,0)]*DeltasMU[GU(i,j,k,2)]; //North

				SupMU[GU(i,j,k,4)] = DeltasMU[GU(i,j,k,0)]*DeltasMU[GU(i,j,k,1)]; //Here
				SupMU[GU(i,j,k,5)] = DeltasMU[GU(i,j,k,0)]*DeltasMU[GU(i,j,k,1)]; //There
			}
		}
	}

	//Superficies Volumen Matriz V
	for(i = 0; i < NX; i++){
		for(k = 0; k < NZ; k++){

			//Parte Abajo
			SupMV[GV(i,0,k,0)] = DeltasMV[GV(i,0,k,1)]*DeltasMV[GV(i,0,k,2)]; //West
			SupMV[GV(i,0,k,1)] = DeltasMV[GV(i,0,k,1)]*DeltasMV[GV(i,0,k,2)]; //East

			SupMV[GV(i,0,k,2)] = DeltasMV[GV(i,0,k,0)]*DeltasMV[GV(i,0,k,2)]; //South
			SupMV[GV(i,0,k,3)] = DeltasMP[GP(i,0,k,0)]*DeltasMP[GP(i,0,k,2)]; //North

			SupMV[GV(i,0,k,4)] = DeltasMV[GV(i,0,k,0)]*DeltasMV[GV(i,0,k,1)]; //Here
			SupMV[GV(i,0,k,5)] = DeltasMV[GV(i,0,k,0)]*DeltasMV[GV(i,0,k,1)]; //There

			//Parte Arriba
			SupMV[GV(i,NY,k,0)] = DeltasMV[GV(i,NY,k,1)]*DeltasMV[GV(i,NY,k,2)]; //West
			SupMV[GV(i,NY,k,1)] = DeltasMV[GV(i,NY,k,1)]*DeltasMV[GV(i,NY,k,2)]; //East

			SupMV[GV(i,NY,k,2)] = DeltasMP[GP(i,NY-1,k,0)]*DeltasMP[GP(i,NY-1,k,2)]; //South
			SupMV[GV(i,NY,k,3)] = DeltasMV[GV(i,NY,k,0)]*DeltasMV[GV(i,NY,k,2)]; //North

			SupMV[GV(i,NY,k,4)] = DeltasMV[GV(i,NY,k,0)]*DeltasMV[GV(i,NY,k,1)]; //Here
			SupMV[GV(i,NY,k,5)] = DeltasMV[GV(i,NY,k,0)]*DeltasMV[GV(i,NY,k,1)]; //There

			for(j = 1; j < NY; j++){
				SupMV[GV(i,j,k,0)] = DeltasMV[GV(i,j,k,1)]*DeltasMV[GV(i,j,k,2)]; //West
				SupMV[GV(i,j,k,1)] = DeltasMV[GV(i,j,k,1)]*DeltasMV[GV(i,j,k,2)]; //East

				SupMV[GV(i,j,k,2)] = DeltasMP[GP(i,j-1,k,0)]*DeltasMP[GP(i,j-1,k,2)]; //South
				SupMV[GV(i,j,k,3)] = DeltasMP[GP(i,j,k,0)]*DeltasMP[GP(i,j,k,2)]; //North

				SupMV[GV(i,j,k,4)] = DeltasMV[GV(i,j,k,0)]*DeltasMV[GV(i,j,k,1)]; //Here
				SupMV[GV(i,j,k,5)] = DeltasMV[GV(i,j,k,0)]*DeltasMV[GV(i,j,k,1)]; //There
			}
		}
	}

	//Superficie Volumen Matriz W
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){

			//Parte Here
			SupMW[GW(i,j,0,0)] = DeltasMW[GW(i,j,0,1)]*DeltasMW[GW(i,j,0,2)]; //West
			SupMW[GW(i,j,0,1)] = DeltasMW[GW(i,j,0,1)]*DeltasMW[GW(i,j,0,2)]; //East

			SupMW[GW(i,j,0,2)] = DeltasMW[GW(i,j,0,0)]*DeltasMW[GW(i,j,0,2)]; //South
			SupMW[GW(i,j,0,3)] = DeltasMW[GW(i,j,0,0)]*DeltasMW[GW(i,j,0,2)]; //North

			SupMW[GW(i,j,0,4)] = DeltasMW[GW(i,j,0,0)]*DeltasMW[GW(i,j,0,1)]; //Here
			SupMW[GW(i,j,0,5)] = DeltasMP[GP(i,j,0,0)]*DeltasMP[GP(i,j,0,1)]; //There

			//Parte There
			SupMW[GW(i,j,NZ,0)] = DeltasMW[GW(i,j,NZ,1)]*DeltasMW[GW(i,j,NZ,2)]; //West
			SupMW[GW(i,j,NZ,1)] = DeltasMW[GW(i,j,NZ,1)]*DeltasMW[GW(i,j,NZ,2)]; //East

			SupMW[GW(i,j,NZ,2)] = DeltasMW[GW(i,j,k,0)]*DeltasMW[GW(i,j,NZ,2)]; //South
			SupMW[GW(i,j,NZ,3)] = DeltasMW[GW(i,j,k,0)]*DeltasMW[GW(i,j,NZ,2)]; //North

			SupMW[GW(i,j,NZ,4)] = DeltasMP[GP(i,j,NZ-1,0)]*DeltasMP[GP(i,j,NZ-1,1)]; //Here
			SupMW[GW(i,j,NZ,5)] = DeltasMW[GW(i,j,NZ,0)]*DeltasMW[GW(i,j,NZ,1)]; //There

			for(k = 1; k < NZ; k++){
				SupMW[GW(i,j,k,0)] = DeltasMW[GW(i,j,k,1)]*DeltasMW[GW(i,j,k,2)]; //West
				SupMW[GW(i,j,k,1)] = DeltasMW[GW(i,j,k,1)]*DeltasMW[GW(i,j,k,2)]; //East

				SupMW[GW(i,j,k,2)] = DeltasMW[GW(i,j,k,0)]*DeltasMW[GW(i,j,k,2)]; //South
				SupMW[GW(i,j,k,3)] = DeltasMW[GW(i,j,k,0)]*DeltasMW[GW(i,j,k,2)]; //North

				SupMW[GW(i,j,k,4)] = DeltasMP[GP(i,j,k-1,0)]*DeltasMP[GP(i,j,k-1,1)]; //Here
				SupMW[GW(i,j,k,5)] = DeltasMP[GP(i,j,k,0)]*DeltasMP[GP(i,j,k,1)]; //There
			}
		}
	}
}

//Cálculo de los volúmenes de control de cada volúmen
void Mesher::Get_Volumes(){
int i, j, k;

	//Volúmenes Matriz P
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				VolMP[GP(i,j,k,0)] = DeltasMP[GP(i,j,k,0)]*SupMP[GP(i,j,k,0)];
			}
		}
	}

	//Volúmenes Matriz U
	for(i = 0; i < NX + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				VolMU[GU(i,j,k,0)] = DeltasMU[GU(i,j,k,0)]*SupMU[GU(i,j,k,0)];
			}
		}
	}

	//Volúmenes Matriz V
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){
				VolMV[GV(i,j,k,0)] = DeltasMV[GV(i,j,k,1)]*DeltasMV[GV(i,j,k,0)]*DeltasMV[GV(i,j,k,2)];
			}
		}
	}

	
	//Volúmenes Matriz W
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){
				VolMW[GW(i,j,k,0)] = DeltasMW[GW(i,j,k,2)]*SupMW[GW(i,j,k,4)];
			}
		}
	}

}

//Pasar los resultados a un archivo VTK en 3D
void Mesher::MallaVTK3D(string Carpeta, string Variable, string NombreFile, double *MC, int Nx, int Ny, int Nz){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<DIRECTORIO<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<Nx<<"   "<<Ny<<"   "<<Nz<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<Nx*Ny*Nz<<"   "<<"double"<<endl;
	
	for(k = 0; k < Nz; k++){
		for(j = 0; j < Ny; j++){
			for(i = 0; i < Nx; i++){
				file<<MC[GP(i,j,k,0)]<<"   "<<MC[GP(i,j,k,1)]<<"   "<<MC[GP(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
	file<<"POINT_DATA"<<"   "<<Nx*Ny*Nz<<endl;
    file<<"SCALARS "<<Variable<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
    file<<endl;
	for(i = 0; i < Nx; i++){
		for(j = 0; j < Ny; j++){
			for(k = 0; k < Nz; k++){
				file<<1.5<<" ";
			}
		}
	}

    file.close();
}

//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1, ParPro MPI1){
int i, j, k;

	MPI1.Initial_WorkSplit(NX, Ix, Fx);
	
	AllocateMemory(M1); //Alojamiento de memoria para cada matriz 
	Get_Meshes(); //Creación de todas las mallas
	Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
	Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control
	Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen

	if(Rank == 0){	
	//	PrintTxt(); //Pasar los resultados de las mallas a un txt
	//	MallaVTK3D("ParaviewResults/MeshResults/", "MallaP", "MallaMP", MP, NX, NY, NZ);
		cout<<"Mesh created."<<endl;
	}

}

