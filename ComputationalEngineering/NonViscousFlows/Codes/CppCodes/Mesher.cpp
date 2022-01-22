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
#include <time.h>
#include <mpi.h>

using namespace std;

#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Mesher.h"

#define DIRECTORIO "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/"

#define PI 3.141592653589793

#define sind(x) sin(x * (PI/180.0)) //Cálculo seno en grados
#define cosd(x) cos(x * (PI/180.0)) //Cálculo coseno en grados
#define tand(x) tan(x * (PI/180.0)) //Cálculo tangente en grados

#define Hyp(x1, x2, y1, y2) sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0)) //Cálculo hipotenusa

#define G(i,j,dim) (((j) + (i)*NY) + NX*NY*(dim)) //Global Index
#define GU(i,j,dim) (((j) + (i)*NY) + (NX+1)*NY*(dim)) //Global Index Matriz U
#define GR(i,j,dim) (((j) + (i)*(NY+1)) + NX*(NY+1)*(dim)) //Global Index Matriz R

#define MU(i,j,dim) ((j) + NY*(dim)) 
#define MR(i,j,dim) ((i) + (Fx - Ix)*(dim) + (Fx - Ix)*(4)*(j))

#define VR(i,j,dim) ((i) + (Fx - Ix)*(j))

#define LNH(i,j,dim) (((j) + ((i) - Ix + Halo)*NY)) //Local Index No Halo 
#define LSH(i,j,dim) (((j) + ((i) - Ix)*NY) + NY*(Fx-Ix + 2*Halo)*dim) //Local Index Si Halo

//Constructor del mallador
Mesher::Mesher(Memory M1, ReadData R1, ParPro MPI1, int i){
		
	//Datos sobre el problema
		Problema = R1.Problema; //Problema (Canal/Canal con Cilindro/Perfil)

	//Datos para la computación en paralelo
		Procesos = MPI1.Procesos;
		Rank = MPI1.Rank;

	//Datos sobre la geometría del problema

		//Geometría del caNAl
		ChannelLength = R1.GeometryData[0]; //Longitud del canal
		ChannelHeight = R1.GeometryData[1]; //Altura del canal
		ChannelDepth = R1.GeometryData[2]; //Profundidad del canal

		//Densidad del mallado	
		NX = R1.NumericalData[0];
		NY = R1.NumericalData[1];
		// //Número de nodos dirección X
		//NY =  //Número de nodos dirección Y
		
		//Opciones del mallado
		OpcionR = R1.NumericalData[2]; //Tipo de discretización en la dirección axial
		OpcionA = R1.NumericalData[3]; //Tipo de discretización en la dirección radial
		StretFactorX = R1.ProblemData[0]; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección axial
		StretFactorY = R1.ProblemData[1]; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección radial
		
}

//Cálculo del la posición axial de los nodos de velocidad axial
void Mesher::Get_AxialCoordinates(){
int i;
double I;
double nx1 = NX/2;
int NX2 = NX - NX/2;
double nx2 = NX2;

	if(OpcionA == 1){ //Regular
		for(i = 0; i < NX + 1; i++){
			I = i;
			AxialCoord[i] = I*(ChannelLength/NX);
		}
	}
	else if(OpcionA == 2){ //Tangencial hiperbólica
		
		for(i = 0; i < NX/2 + 1; i++){
			I = i;
			AxialCoord[i] = (0.50*ChannelLength)*(tanh(StretFactorX*(I/nx1))/tanh(StretFactorX));	
		}
		for(i = 1; i < NX2 + 1; i++){
			I = i;
			AxialCoord[i+NX2] = (0.50*ChannelLength) + (0.50*ChannelLength)*((tanh(StretFactorX*(I/nx2 - 1.0)) + tanh(StretFactorX))/tanh(StretFactorX));
		}		
	}
	

}

//Alojamiento de memoria para cada matriz
void Mesher::AllocateMemory(Memory M1){
	
	AxialCoord = M1.AllocateDouble(NX + 1, 1, 1); //Matriz coordenada axial de la distribución

	//Matrices de coordeNXdas de discretización
	MP = M1.AllocateDouble(NX, NY, 2); //CoordeNXdas matriz de presión/temperatura
	MU = M1.AllocateDouble(NX + 1, NY, 2); //CoordeNXdas matriz velocidades axiales
	MR = M1.AllocateDouble(NX, NY + 1, 2); //CoordeNXdas matriz velocidades radiales

	//Matrices de distancias de volúmenes de control
	DeltasMP = M1.AllocateDouble(NX, NY, 2); //Deltas X y R de la matriz de Presión/Temperatura
	DeltasMU = M1.AllocateDouble(NX + 1, NY, 2); //Deltas X y R de la matriz de velocidades axiales (U)
	DeltasMR = M1.AllocateDouble(NX, NY + 1, 2); //Deltas X y R de la matriz de velocidades radiales (V)

	//Matrices de superficies de los volúmenes de control
	SupMP = M1.AllocateDouble(NX, NY, 4);
						
}

//Creación de todas las mallas de tipo Staggered
void Mesher::Get_Mesh(){
int i, j;
double ny1 = NY/2;
int NY2 = NY - NY/2;
double ny2 = NY2;
double J;

	//Coordenadas Axiales

	//Coordenadas axiales matriz de velocidades axiales (U)
	for(i = 0; i < NX + 1; i++){
		for(j = 0; j < NY; j++){
			MU[GU(i,j,0)] = AxialCoord[i];	
		}
	}

	//CoordeNXdas axiales matriz de velocidades radial (V)
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY + 1; j++){
			MR[GR(i,j,0)] = 0.5*(AxialCoord[i] + AxialCoord[i+1]);
		}
	}

	//CoordeNXdas axiales matriz de Presión/Temperatura
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			MP[G(i,j,0)] = 0.5*(AxialCoord[i] + AxialCoord[i+1]);
		}
	}

	//CoordeNXdas Radiales

	//CoordeNXdas radiales de la matriz de velocidades radial (V)
	if(OpcionR == 1){ //Regular
		for(i = 0; i < NX; i++){
			for(j = 0; j < NY + 1; j++){
				J = j;
				MR[GR(i,j,1)] = J*(ChannelHeight/NY);
			}
		}
	}
	else if(OpcionR == 2){ //Tangencial hiperbólica
		for(i = 0; i < NX; i++){
			for(j = 0; j < NY/2 + 1; j++){
				J = j;
				MR[GR(i,j,1)] = (0.50*ChannelHeight)*(tanh(StretFactorY*(J/ny1))/tanh(StretFactorY));	
			}
			for(j = 1; j < NY2 + 1; j++){
				J = j;
				MR[GR(i,j+NY/2,1)] = (0.50*ChannelHeight) + (0.50*ChannelHeight)*((tanh(StretFactorY*(J/ny2 - 1.0)) + tanh(StretFactorY))/tanh(StretFactorY));
			}
		}
	}
			
	//CoordeNXdas radiales de la matriz de Presión/Temperatura
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			MP[G(i,j,1)] = 0.5*(MR[GR(i,j,1)] + MR[GR(i,j + 1,1)]);
		}
	}

	//Coordenadas radiales de la matriz de velocidades axial (U)
	for(j = 0; j < NY; j++){
		//Parte Izquierda
		MU[GU(0,j,1)] = MP[G(0,j,1)];

		//Parte Derecha
		MU[GU(NX,j,1)] = MP[G(NX - 1,j,1)];

		for(i = 1; i < NX; i++){ 	
			MU[GU(i,j,1)] = 0.5*(MP[G(i-1,j,1)] + MP[G(i,j,1)]);
		}
	}

}

//Pasar los resultados de las mallas a un txt
void Mesher::PrintTxt(){
int i, j;
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
			file<<MP[G(i,j,0)]<<"\t"<<MP[G(i,j,1)]<<"\t"<<endl;
			}
        	file<<endl;
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
			file<<MU[GU(i,j,0)]<<"\t"<<MU[GU(i,j,1)]<<"\t"<<endl;
			}
        	file<<endl;
    	}

	file.close();

stringstream InitialNameMR;
string FinalNameMR;

	FileName = "MallaMR.txt";

	InitialNameMR<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMR = InitialNameMR.str();
        file.open(FinalNameMR.c_str());

		for(i = 0; i < NX; i++){
        	for(j = 0; j < NY + 1; j++){
			file<<MR[GR(i,j,0)]<<"\t"<<MR[GR(i,j,1)]<<"\t"<<endl;
			}
        	file<<endl;
    	}

	file.close();

}

//Cálculo de las distancias entre nodos en cada uNX de las matrices
void Mesher::Get_Deltas(){
int i, j;

	//Deltas X e Y de la matriz de Presión/Temperatura
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			DeltasMP[G(i,j,0)] =  MU[GU(i + 1,j,0)] - MU[GU(i,j,0)]; //Deltas X
			DeltasMP[G(i,j,1)] =  MR[GR(i,j + 1,1)] - MR[GR(i,j,1)]; //Deltas Y
		}
	}

	//Deltas X e Y de la matriz de velocidades axiales (U)	
	for(j = 0; j < NY; j++){
		//Parte Izquierda
		DeltasMU[GU(0,j,0)] = MP[G(0,j,0)] - MU[GU(0,j,0)]; //Deltas X
		DeltasMU[GU(0,j,1)] = MR[GR(0,j + 1,1)] - MR[GR(0,j,1)]; //Deltas Y

		//Parte Derecha
		DeltasMU[GU(NX,j,0)] = MU[GU(NX,j,0)] - MP[G(NX - 1,j,0)]; //Deltas X
		DeltasMU[GU(NX,j,1)] = MR[GR(NX - 1,j + 1,1)] - MR[GR(NX - 1,j,1)]; //Deltas Y

		for(i = 1; i < NX; i++){
			DeltasMU[GU(i,j,0)] = MP[G(i,j,0)] - MP[G(i - 1,j,0)]; //Deltas X
			DeltasMU[GU(i,j,1)] = MR[GR(i,j + 1,1)] - MR[GR(i,j,1)]; //Deltas Y
		}
	}
	
	//Deltas X e Y de la matriz de velocidades Verticales (V)
	for(i = 0; i < NX; i++){
		//Parte abajo
		DeltasMR[GR(i,0,0)] = MU[GU(i + 1,0,0)] - MU[GU(i,0,0)]; //Deltas X
		DeltasMR[GR(i,0,1)] = MP[G(i,0,1)] - MR[GR(i,0,1)]; //Deltas Y

		//Parte arriba
		DeltasMR[GR(i,NY,0)] = MU[GU(i + 1,NY - 1,0)] - MU[GU(i,NY - 1,0)]; //Deltas X
		DeltasMR[GR(i,NY,1)] = MR[GR(i,NY,1)] - MP[G(i,NY - 1,1)]; //Deltas Y

		for(j = 1; j < NY; j++){
			DeltasMR[GR(i,j,0)] = MU[GU(i + 1,j,0)] - MU[GU(i,j,0)]; //Deltas X
			DeltasMR[GR(i,j,1)] = MP[G(i,j,1)] - MP[G(i,j - 1,1)]; //Deltas Y
		}
	}

}

//Cálculo de las superficies de los volúmenes de control
void Mesher::Get_Surfaces(){
int i, j;

	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			SupMP[G(i,j,0)] = DeltasMP[G(i,j,0)]*ChannelDepth;
			SupMP[G(i,j,1)] = DeltasMP[G(i,j,0)]*ChannelDepth;
			SupMP[G(i,j,2)] = DeltasMP[G(i,j,1)]*ChannelDepth;
			SupMP[G(i,j,3)] = DeltasMP[G(i,j,1)]*ChannelDepth;
		}
	}

}

//Pasar los resultados a un archivo VTK en 2D
void Mesher::MallaVTK2D(string Carpeta, string Variable, string NombreFile, double *MC, int NX, int NY){
int i, j;

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
        file<<"DIMENSIONS"<<"   "<<NX<<"   "<<NY<<"   "<<1<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<NX*NY<<"   "<<"double"<<endl;
	
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				file<<MC[G(i,j,0)]<<"   "<<MC[G(i,j,1)]<<"   "<<0.0<<endl;	
			}	
		}
        
        file<<endl;
		file<<"POINT_DATA"<<"   "<<NX*NY<<endl;
        file<<"SCALARS "<<Variable<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
        file<<endl;
        for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				file<<0.0<<" ";	
			}
		}

        file.close();
}

//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1, ParPro MPI1){
	
	MPI1.Initial_WorkSplit(NX, Ix, Fx);
	AllocateMemory(M1); //Alojamiento de memoria para cada matriz 
	Get_AxialCoordinates(); //Cálculo del la posición axial de los nodos de velocidad axial
	Get_Mesh(); //Creación de todas las mallas
	Get_Deltas(); //Cálculo de las distancias entre nodos en cada uNX de las matrices
	Get_Surfaces();

	if(Rank == 0){
		PrintTxt();
		MallaVTK2D("ParaviewResults/MeshResults/", "MallaBase", "MallaMP", MP, NX, NY);
	}

	cout<<"Mesh created."<<endl;

}

