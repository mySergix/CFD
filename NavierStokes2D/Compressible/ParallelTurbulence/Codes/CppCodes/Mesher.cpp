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

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ReadData.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ParPro.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Geometry.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Mesher.h"

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/"

#define PI 3.141592653589793

#define sind(x) sin(x * (PI/180.0)) //Cálculo seno en grados
#define cosd(x) cos(x * (PI/180.0)) //Cálculo coseno en grados
#define tand(x) tan(x * (PI/180.0)) //Cálculo tangente en grados

#define Hyp(x1, x2, y1, y2) sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0)) //Cálculo hipotenusa

#define G(i,j,dim) (((j) + (i)*NR) + NA*NR*(dim)) //Global Index
#define GU(i,j,dim) (((j) + (i)*NR) + (NA+1)*NR*(dim)) //Global Index Matriz U
#define GR(i,j,dim) (((j) + (i)*(NR+1)) + NA*(NR+1)*(dim)) //Global Index Matriz R

#define MU(i,j,dim) ((j) + NR*(dim)) 
#define MR(i,j,dim) ((i) + (Fx - Ix)*(dim) + (Fx - Ix)*(4)*(j))

#define VR(i,j,dim) ((i) + (Fx - Ix)*(j))

#define LNH(i,j,dim) (((j) + ((i) - Ix + Halo)*NR)) //Local Index No Halo 
#define LSH(i,j,dim) (((j) + ((i) - Ix)*NY) + NY*(Fx-Ix + 2*Halo)*dim) //Local Index Si Halo

//Constructor del mallador
Mesher::Mesher(Memory M1, ReadData R1, Geometry G1, ParPro MPI1){
		
	//Datos sobre el problema
		Problema = R1.Get_ProblemType(); //Problema (Tobera/Tubería/Canal)
		Tobera = R1.Get_NozzleType();	//Tipo de tobera (NASA/RAO)
	
	//Datos necesarios para computación paralela
		Rank = MPI1.Rank;
		Procesos = MPI1.Procesos;
		
	//Datos sobre la geometría del problema

		//Selección problema canal

		ChannelLength = R1.GeometryData[14]; //Longitud del canal
		ChannelHeight = R1.GeometryData[15]; //Altura del canal
		ChannelDepth = R1.GeometryData[16]; //Profundidad del canal

		//Selección problema tubería

		PipeLength = R1.GeometryData[0]; //Longitud de la tubería
		PipeDiameter = R1.GeometryData[1]; //Diametro de la tubería
	
		//Selección problema tobera

		//Datos geometría sección convergente
		L1 = R1.GeometryData[2]/1000.0; //Longitud de la cámara de combustión (Sección 1)
		D1 = R1.GeometryData[3]/1000.0; //Diámetro de la cámara de combustión (Sección 1)

		Rad1 = R1.GeometryData[4]/1000.0; //Radio Sección 2
		Theta1 = R1.GeometryData[5]; //Ángulo Sección 2

		Rad2 = R1.GeometryData[6]/1000.0; //Radio Sección 4
		Theta2 = R1.GeometryData[7]; //Ángulo Sección 4

		//Datos geometría sección divergente
		Dt = R1.GeometryData[8]/1000.0; //Diámetro de la garganta de la tobera
		Rad3 = R1.GeometryData[9]/1000.0; //Radio Seccion 5 de la tobera
		
		Theta3 = R1.GeometryData[10]; //Ángulo del punto de tangencia sección divergente
	
		xE = R1.GeometryData[11]/1000.0; //Coordenada axial de la salida de la tobera
		rE = R1.GeometryData[12]/2000.0; //Radio de salida de la tobera

		ThetaE = R1.GeometryData[13]; //Ángulo de salida de la tobera

		//Datos del mallado
		TypeMeshing = R1.Get_MeshingType(); //Tipo de mallado seleccionado (Collocated/Staggered)
		
		//Opciones del mallado
		OpcionSmooth = R1.NumericalData[0];
		OpcionSide = R1.NumericalData[1];
		OpcionA = R1.NumericalData[7]; //Tipo de discretización en la dirección axial
		OpcionR = R1.NumericalData[8]; //Tipo de discretización en la dirección radial
		StretFactorX = R1.ProblemData[0]; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección axial
		StretFactorR = R1.ProblemData[1]; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección radial

		//Densidad del mallado
		NAConv = R1.NumericalData[2]; //Dirección axial parte convergente
		NADiv =  R1.NumericalData[3]; //Dirección axial parte divergente
		NA = R1.NumericalData[4]; //Número de nodos totales
		NR = R1.NumericalData[5]; //Dirección radial
		NAng = R1.NumericalData[6]; //Dirección angular (solo para postproceso)

		//Opciones del mallado en la pared
		OpcionYplus = R1.ProblemData[2]; //Opción para poner un nodo en Y+
		Yplus = R1.ProblemData[2]; //Y+ seleccionado

		//Cálculos de datos necesarios
		m = -tand(Theta1); //Pendiente recta Sección 3
		rW = (Dt/2.0) +  Rad3*(1.0 - cosd(Theta3)); //Radio del punto de tangencia sección divergente
		Xt = L1 + Rad1*sind(Theta1) + Rad2*sind(Theta2) + ((Dt/2.0 + Rad2*(1.0 - cosd(Theta2))) - (D1/2.0 - Rad1*(1.0 - cosd(Theta1))))/m;
		//Coordenada axial de la garganta de la tobera
		xW = Rad3*sind(Theta3); //Coordenada axial del punto de tangencia sección divergente
		
		//Separaciones secciones
		P1 = L1; //Separación secciones 1 y 2
		P2 = P1 + Rad1*sind(Theta1); //Separación secciones 2 y 3
		P3 = P2 + ((Dt/2.0 + Rad2*(1.0 - cosd(Theta2))) - (D1/2.0 - Rad1*(1.0 - cosd(Theta1))))/m; //Separación secciones 3 y 4
		P4 = Xt; //Separación secciones 4 y 5
		P5 = Xt + xW; //Separación secciones 5 y 6
		P6 = Xt + xE; //Final geometría

		if(Problema == "Tuberia"){
			Rinicio = PipeDiameter/2.0; //Radio inicial de la geometría
			Rfinal = PipeDiameter/2.0; //Radio final de la geometría	

			Xinicio = 0.0; //Coordenada axial inicial de la geometría
			Xfinal = PipeLength; //Coordenada axial final de la geometría
		}
		else if(Problema == "Tobera"){
			Rinicio = D1/2.0; //Radio inicial de la geometría
			Rfinal = rE; //Radio final de la geometría

			Xinicio = 0.0; //Coordenada axial inicial de la geometría
			Xfinal = P6; //Coordenada axial final de la geometría

			ThroatArea = 0.25*PI*pow(Dt,2.0);
		}
		else if(Problema == "Canal"){
			Rinicio = ChannelHeight;
			Xinicio = 0.0;

			Rfinal = ChannelHeight;
			Xfinal = ChannelLength;
		}

		
}

//Cálculo del número de nodos totales en función de las opciones de discretización
void Mesher::GetTotalAxialNodes(){
double Delta;
double naConv = NAConv;
double naDiv = NADiv;

	if(OpcionSmooth == 0){
		NAConv = (P4/Xfinal)*NA;
		NADiv = NA - NAConv;
	}
	
		else if(OpcionSmooth == 1){
			if(OpcionA == 1){//Regular
				if(OpcionSide == 1){//Parte Convergente
					Delta = Xt/naConv;
					NADiv = (xE - Xt)/Delta;
					NA = NAConv + NADiv;
				}
				else if(OpcionSide == 2){//Parte Divergente
					Delta = (xE - Xt)/naDiv;
					NAConv = Xt/Delta;
					NA = NAConv + NADiv;
				}	
			}
			else if(OpcionA == 2){//Tangente Hiperbólica
				if(OpcionSide == 1){//Parte Convergente
					Delta = Xt*(1.0 - tanh(StretFactorX*(naConv/(naConv+1.0)))/tanh(StretFactorX));
					NADiv =  1.0/(atanh((Delta*tanh(StretFactorX))/(xE-Xt) - tanh(StretFactorX))/StretFactorX + 1.0);
					NA = NAConv + NADiv;
				}
				else if(OpcionSide == 2){//Parte Divergente
					Delta = (xE - Xt)*((tanh(StretFactorX*(1.0/(naDiv + 1.0) - 1.0)) + tanh(StretFactorX))/tanh(StretFactorX));
					NAConv = StretFactorX/(StretFactorX - atanh((1.0 - Delta/Xt)*tanh(StretFactorX)));
					NA = NAConv + NADiv;
				}
			}
			else if(OpcionA == 3){//Senoidal - Cosenoidal
				if(OpcionSide == 1){//Parte convergente
					Delta = Xt*(1.0 - sind(90.0*((naConv-1.0)/naConv)));
					NADiv = (PI/2.0)/acos(1.0 - Delta/(xE-Xt));
					NA = NAConv + NADiv;
				}
				else if(OpcionSide == 2){//Parte divergente
					Delta = (xE-Xt)*(1.0 - cosd(90.0*(1.0/naDiv)));
					NAConv = (PI/2.0)/(PI/2.0 - asin(1.0 - Delta/Xt));
					NA = NAConv + NADiv;
				}
			}
		}
			
}



//Cálculo del la posición axial de los nodos de velocidad axial
void Mesher::Get_AxialCoordinates(){
int i;
double na = NA;
double I;
double naConv = NAConv;
double naDiv = NADiv;
		
		if(Problema == "Tuberia"){ 
			if(OpcionA == 1){
				for(i = 0; i < NA+1; i++){ 
				I = i;
				AxialCoord[i] = I*(PipeLength/NA);
				}
			}
			else if(OpcionA == 2){
				for(i = 0; i < NA+1; i++){
					I = i;
					AxialCoord[i] = PipeLength*(1.0 + (tanh(StretFactorX*((I - na)/na)) + tanh(StretFactorX))/tanh(StretFactorX) - 1.0);	
				} 
			}//(1+(tanh(StretFactor*((2*I-nxtot)/nxtot)) + tanh(StretFactor))/tanh(StretFactor) - 1)
			
		}
		else if(Problema == "Tobera"){
			if(OpcionA == 1){ //Regular
				for(i = 0; i < NAConv+1; i++){
					I = i;
					AxialCoord[i] = I*(P4/naConv);	
				} 
				for(i = 1; i < NADiv+1; i++){
					I = i;
					AxialCoord[i+NAConv] = P4 + ((P6 - P4)/naDiv)*I;
				}

			}
			else if(OpcionA == 2){ //Tangencial hiperbólica
				for(i = 0; i < NAConv+1; i++){
					I = i;
					AxialCoord[i] = P4*(tanh(StretFactorX*(I/naConv))/tanh(StretFactorX));	
				} 
				for(i = 1; i < NADiv+1; i++){
					I = i;
					AxialCoord[i+NAConv] = P4 + (P6 - P4)*((tanh(StretFactorX*(I/naDiv - 1.0)) + tanh(StretFactorX))/tanh(StretFactorX));
				}

			}
			else if(OpcionA == 3){ //Senoidal - Cosenoidal
				for(i = 0; i < NAConv+1; i++){	
					I = i;
					AxialCoord[i] = P4*sind(90.0*(I/naConv));
				}
				for(i = 1; i < NADiv+1; i++){
					I = i;
					AxialCoord[i+NAConv] = P4 + (P6 - P4)*(1.0 - cosd(90.0*(I/naDiv)));
				}
			}
		}
		else if(Problema == "Canal"){
			for(i = 0; i < NA+1; i++){ 
				I = i;
				AxialCoord[i] = I*(ChannelLength/NA);
			}
		}
	

}

//Cálculo del radio de cada coordenada axial discretizada
void Mesher::Get_Radius(Memory M1, Geometry G1){
int i;

	if(Problema == "Tuberia"){
		for(i = 0; i < NA + 1; i++){ Radius[i] = PipeDiameter/2.0; }
	}
	else if(Problema == "Tobera"){
		for(i = 0; i < NA + 1; i++){
			if(AxialCoord[i] <= P1){
				Radius[i] = Rinicio;
			}
			else if(AxialCoord[i] > P1 && AxialCoord[i] <= P2){
				Radius[i] = sqrt(pow(Rad1,2.0) - pow(0.5*(AxialCoord[i] + AxialCoord[i+1]) - P1,2.0)) - (Rad1 - D1/2.0);
			}
			else if(AxialCoord[i] > P2 && AxialCoord[i] <= P3){
				Radius[i] = D1/2.0 - Rad1*(1 - cosd(Theta1)) + m*(0.5*(AxialCoord[i] + AxialCoord[i+1]) - P2);
			}
			else if(AxialCoord[i] > P3 && AxialCoord[i] <= P4){
				Radius[i] = -sqrt(pow(Rad2,2.0) - pow(0.5*(AxialCoord[i] + AxialCoord[i+1]) - P4,2.0)) + (Dt/2.0 + Rad2);
			}
			else if(AxialCoord[i] > P4 && AxialCoord[i] <= P5){
				Radius[i] = -sqrt(pow(Rad3,2.0) - pow(0.5*(AxialCoord[i] + AxialCoord[i+1]) - P4,2.0)) + (Dt/2.0 + Rad3);
			}
			else{		
				if(Tobera == "NASA"){
					Radius[i] = G1.CoeficientesExpresion[0] + sqrt(G1.CoeficientesExpresion[1] + G1.CoeficientesExpresion[2]*(0.5*(AxialCoord[i] + AxialCoord[i+1]) - P5)) + G1.CoeficientesExpresion[3]*(0.5*(AxialCoord[i] + AxialCoord[i+1]) - P5);
				}
				else if(Tobera == "RAO"){
					Radius[i] = (- G1.CoeficientesExpresion[1] + sqrt(pow(G1.CoeficientesExpresion[1],2.0) - 4.0*G1.CoeficientesExpresion[0]*(G1.CoeficientesExpresion[2] - 0.5*(AxialCoord[i] + AxialCoord[i+1]))))/(2*G1.CoeficientesExpresion[0]);
				}
			}
		}
	}
	else if(Problema == "Canal"){
		for(i = 0; i < NA + 1; i++){ Radius[i] = ChannelHeight; }
	}
	
}

//Calcular el número total de nodos en la dirección radial
void Mesher::GetTotalRadialNodes(){
int i, j;
double Sumatorio = 0.0;
double Rplus;

	if(OpcionYplus == 0){ //No activado
		NR = NR;
	}
	else if(OpcionYplus == 1){ //Activado
		Sumatorio += Radius[0]*0.50*(AxialCoord[0] + AxialCoord[1]);
		Sumatorio += Radius[NA-1]*(Xfinal - 0.50*(AxialCoord[NA-2] + AxialCoord[NA-1]));

		for(i = 1; i < NA-1; i++){
			Sumatorio += Radius[i]*0.50*(AxialCoord[i+1] - AxialCoord[i-1]);
			
		}
		
		
		Rplus = Sumatorio/Xfinal;
		DeltaYplus = 0.012186681661309455;
		
		if(OpcionR == 1){ //Regular
			NR = Rplus/(2.0*DeltaYplus);
		}
		else if(OpcionR == 2){ //Tangencial hiperbólica
			NR = StretFactorR/(StretFactorR - atanh((1.0 - (2.0*DeltaYplus)/Rplus)*tanh(StretFactorR)));
		}
		else if(OpcionR == 3){ //Senoidal
			NR = (PI/2.0)/acos(1.0 - (DeltaYplus)/Rplus);
		}
	}

}

//Alojamiento de memoria para cada matriz
void Mesher::AllocateMemory(Memory M1){

	//Matrices de coordenadas de discretización
	MP = M1.AllocateDouble(NA, NR, 2); //Coordenadas matriz de presión/temperatura
	MU = M1.AllocateDouble(NA+1, NR, 2); //Coordenadas matriz velocidades axiales
	MR = M1.AllocateDouble(NA, NR+1, 2); //Coordenadas matriz velocidades radiales

	//Matrices de distancias de volúmenes de control
	DeltasMP = M1.AllocateDouble(NA, NR, 2); //Deltas X y R de la matriz de Presión/Temperatura
	DeltasMU = M1.AllocateDouble(NA+1, NR, 2); //Deltas X y R de la matriz de velocidades axiales (U)
	DeltasMR = M1.AllocateDouble(NA, NR+1, 2); //Deltas X y R de la matriz de velocidades radiales (V)
	
	//Matrices de superficies de volúmenes de control
	SupMP = M1.AllocateDouble(NA, NR, 4); //Superficies de los volúmenes de la matriz de Presión/Temperatura
	if(Rank == 0 || Rank == Procesos-1){
		SupMU = M1.AllocateDouble(1, NR, 4); //Superficies de los volúmenes de la matriz de velocidades axiales (U)
	}
	SupMR = M1.AllocateDouble(Fx - Ix, 2, 4); //Superficies de los volúmenes de la matriz de velocidades radiales (V)

	//Matrices de volúmenes de los volúmenes de control
	VolMP = M1.AllocateDouble(NA, NR, 1); //Volúmenes de los volúmenes de la matriz de Presión/Temperatura
	if(Rank == 0 || Rank == Procesos-1){
		VolMU = M1.AllocateDouble(1, NR, 1); //Superficies de los volúmenes de la matriz de velocidades axiales (U)
	}
	VolMR = M1.AllocateDouble(Fx - Ix, 2, 1); //Volúmenes de los volúmenes de la matriz de velocidades radiales (V)

	//Matrices de los ángulos entre las velocidades y las superficies de los volúmenes de control
	AngleMR = M1.AllocateDouble(NA, NR+1, 1); //Ángulos entre las superficies de los v.c y las velocidades radiales (V)
	AngleMU = M1.AllocateDouble(NA+1, NR, 1); //Ángulos entre las velocidades axiales (U)....
	
	//Distancia mínima de cada nodo a la pared
	minDist = M1.AllocateDouble(NA, NR, 1); //Distancias mínimas de cada nodo a la pared mas cercana

	//Código:
	//0 -> West
	//1 -> East
	//2 -> South
	//3 -> North
							
}

//Creación de todas las mallas de tipo Staggered
void Mesher::Get_Mesh(){
int i, j;
double I, J;
double nr = NR;

	//Coordenadas Axiales

	//Coordenadas axiales matriz de velocidades axiales (U)
	for(i = 0; i < NA+1; i++){
		for(j = 0; j < NR; j++){
			MU[GU(i,j,0)] = AxialCoord[i];	
		}
	}

	//Coordenadas axiales matriz de velocidades radial (V)
	for(i = 0; i < NA; i++){
		for(j = 0; j < NR+1; j++){
			MR[GR(i,j,0)] = 0.5*(AxialCoord[i] + AxialCoord[i+1]);
		}
	}
	//Coordenadas axiales matriz de Presión/Temperatura
	for(i = 0; i < NA; i++){
		for(j = 0; j < NR; j++){
			MP[G(i,j,0)] = 0.5*(AxialCoord[i] + AxialCoord[i+1]);
		}
	}

	//Coordenadas Radiales
	if(Problema == "Tobera" || Problema == "Tuberia"){
		//Coordenadas radiales de la matriz de velocidades radial (V)
		if(OpcionR == 1){ //Regular
			for(i = 0; i < NA; i++){
				for(j = 0; j < NR+1; j++){
					J = j;
					MR[GR(i,j,1)] = J*(Radius[i]/nr);
				}
			}
		}
		else if(OpcionR == 2){ //Tangencial hiperbólica
			for(i = 0; i < NA; i++){
				for(j = 0; j < NR+1; j++){
					J = j;
					MR[GR(i,j,1)] = Radius[i]*((tanh(StretFactorR*((J)/nr)))/tanh(StretFactorR));
				}
			}
		}
		else if(OpcionR == 3){ //Senoidal
			for(i = 0; i < NA; i++){
				for(j = 0; j < NR+1; j++){
					J = j;
					MR[GR(i,j,1)] = Radius[i]*sind(90.0*(J/nr));
				}
			}
		}
	}
	else if(Problema == "Canal"){
		double nrtotal = NR;
		if(OpcionR == 1){ //Regular
			for(i = 0; i < NA; i++){
				for(j = 0; j < NR+1; j++){
					J = j;
					MR[GR(i,j,1)] = J*(Radius[i]/nrtotal);
				}
			}
		}
		else if(OpcionR == 2){ //Tangencial hiperbólica
		
			for(i = 0; i < NA; i++){
				for(j = 0; j < NR+1; j++){
					J = j;
					MR[GR(i,j,1)] = (Radius[i]/2.0)*(1.0 + (tanh(StretFactorR*((2.0*J - nrtotal)/nrtotal)) + tanh(StretFactorR))/tanh(StretFactorR) - 1.0);
				}
			}
		}
	}
	//Coordenadas radiales de la matriz de Presión/Temperatura
	if(Problema == "Tobera" || Problema == "Tuberia"){
		for(i = 0; i < NA; i++){
			MP[G(i,0,1)] = MR[GR(i,0,1)];
			for(j = 1; j < NR; j++){
				MP[G(i,j,1)] = 0.5*(MR[GR(i,j,1)] + MR[GR(i,j+1,1)]);
			}
		}

		//Coordenadas radiales de la matriz de velocidades axial (U)
		for(j = 0; j < NR; j++){ //Parte izda
			MU[GU(0,j,1)] = MP[G(0,j,1)];
		}

		for(i = 1; i < NA; i++){ //Centro
			MU[GU(i,0,1)] = MP[G(i,0,1)];
			for(j = 1; j < NR; j++){
				MU[GU(i,j,1)] = 0.5*(MP[G(i-1,j,1)] + MP[G(i,j,1)]);
			}
		}

		MU[GU(NA,0,1)] = MP[G(NA-1,0,1)];
		for(j = 1; j < NR; j++){ //Parte dra
			if(OpcionR == 1){ //Regular
				J = j;
				MU[GU(NA,j,1)] = 0.50*(MP[G(NA-1,j,1)] + 0.50*(J*((2.0*Rfinal - MR[GR(NA-1,NR,1)])/nr) + (J+1)*((2.0*Rfinal - MR[GR(NA-1,NR,1)])/nr)));
				
			}
			else if(OpcionR == 2){ //Tangencial hiperbólica
				J = j;
				MU[GU(NA,j,1)] = 0.50*(MP[G(NA-1,j,1)] + 0.50*((2.0*Rfinal - MR[GR(NA-1,NR,1)])*(tanh(StretFactorR*(J/nr)))/tanh(StretFactorR) + (2.0*Rfinal - MR[GR(NA-1,NR,1)])*(tanh(StretFactorR*((J+1)/nr)))/tanh(StretFactorR)));
			}
			else if(OpcionR == 3){ //Senoidal
				J = j;
				MU[GU(NA,j,1)] = 0.50*(MP[G(NA-1,j,1)] + 0.50*((2.0*Rfinal - MR[GR(NA-1,NR,1)])*sind(90.0*(J/nr)) + (2.0*Rfinal - MR[GR(NA-1,NR,1)])*sind(90.0*((J+1)/nr))));
			}
		}
	}
	else if(Problema == "Canal"){
		//Coordenadas matriz central
		for(i = 0; i < NA; i++){
			for(j = 0; j < NR; j++){
				MP[G(i,j,1)] = 0.5*(MR[GR(i,j,1)] + MR[GR(i,j+1,1)]);
			}
		}

		//Coordenadas matriz U
		for(j = 0; j < NR; j++){
			MU[GU(NA,j,1)] = MP[G(NA-1,j,1)];
			for(i = 0; i < NA; i++){
				MU[GU(i,j,1)] = MP[G(i,j,1)];
			}
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

	for(i = 0; i < NA; i++){
        	for(j = 0; j < NR; j++){
			file<<MP[G(i,j,0)]<<"\t"<<MP[G(i,j,1)]<<"\t"<<endl;
		}
        	file<<endl;;
    	}
	file.close();

stringstream InitialNameMU;
string FinalNameMU;

	FileName = "MallaMU.txt";

	InitialNameMU<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMU = InitialNameMU.str();
        file.open(FinalNameMU.c_str());

	for(i = 0; i < NA+1; i++){
        	for(j = 0; j < NR; j++){
			file<<MU[GU(i,j,0)]<<"\t"<<MU[GU(i,j,1)]<<"\t"<<endl;
		}
        	file<<endl;;
    	}
	file.close();

stringstream InitialNameMR;
string FinalNameMR;

	FileName = "MallaMR.txt";

	InitialNameMR<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMR = InitialNameMR.str();
        file.open(FinalNameMR.c_str());

	for(i = 0; i < NA; i++){
        	for(j = 0; j < NR+1; j++){
			file<<MR[GR(i,j,0)]<<"\t"<<MR[GR(i,j,1)]<<"\t"<<endl;
		}
        	file<<endl;;
    	}
	file.close();

}

//Cálculo de las distancias entre nodos en cada una de las matrices para el caso del Canal
void Mesher::Get_DeltasChannel(){
int i, j;

	//Deltas X y R de la matriz de Presión/Temperatura
	for(i = 0; i < NA; i++){
		for(j = 0; j < NR; j++){
			DeltasMP[G(i,j,0)] = MU[GU(i+1,j,0)] - MU[GU(i,j,0)]; //Deltas X
			DeltasMP[G(i,j,1)] = MR[GR(i,j+1,1)] - MR[GR(i,j,1)]; //Deltas R
		}
	}

	//Deltas X y R de la matriz U
	for(j = 0; j < NR; j++){
		//Parte izquierda
		DeltasMU[GU(0,j,0)] = MP[G(0,j,0)] - 0.0; //Deltas X
		DeltasMU[GU(0,j,1)] = MR[GR(0,j+1,1)] - MR[GR(0,j,1)]; //Deltas R

		//Parte derecha
		DeltasMU[GU(NA,j,0)] = MU[GU(NA,j,0)] - MP[G(NA-1,j,0)]; //Deltas X
		DeltasMU[GU(NA,j,1)] = MR[GR(NA-1,j+1,1)] - MR[GR(NA-1,j,1)]; //Deltas R

		for(i = 1; i < NA; i++){
			DeltasMU[GU(i,j,0)] = MP[G(i,j,0)] - MP[G(i-1,j,0)]; //Deltas X
			DeltasMU[GU(i,j,1)] = MR[GR(i,j+1,1)] - MR[GR(i,j,1)]; //Deltas R
		}
	}

	//Deltas X y R de la matriz R
	for(i = 0; i < NA; i++){	
		//Parte abajo
		DeltasMR[GR(i,0,0)] = MU[GU(i+1,0,0)] - MU[GU(i,0,0)]; //Deltas X
		DeltasMR[GR(i,0,1)] = MP[G(i,0,1)] - 0.0; //Deltas R

		//Parte arriba
		DeltasMR[GR(i,NR,0)] = MU[GU(i+1,NR-1,0)] - MU[GU(i,NR-1,0)]; //Deltas X
		DeltasMR[GR(i,NR,1)] = MR[GR(i,NR,1)] - MP[G(i,NR-1,1)]; //Deltas R

		for(j = 1; j < NR; j++){
			DeltasMR[GR(i,j,0)] = MU[GU(i+1,j,0)] - MU[GU(i,j,0)]; //Deltas X
			DeltasMR[GR(i,j,1)] = MP[G(i,j,1)] - MP[G(i,j-1,1)]; //Deltas R
		}
	}
}

//Cálculo de las distancias entre nodos en cada una de las matrices
void Mesher::Get_Deltas(){
int i, j;

	//Deltas X y R de la matriz de Presión/Temperatura
	for(i = 0; i < NA; i++){
		for(j = 0; j < NR; j++){
			DeltasMP[G(i,j,0)] = Hyp(MU[GU(i,j,0)], MU[GU(i+1,j,0)], MU[GU(i,j,1)], MU[GU(i+1,j,1)]); //Deltas X
			DeltasMP[G(i,j,1)] = MR[GR(i,j+1,1)] - MR[GR(i,j,1)]; //Deltas R
		}
	}

	//Deltas X y R de la matriz de velocidades axiales (U)	

	DeltasMU[GU(NA,0,1)] = 0.50*(0.50*(MR[GR(NA-1,0,1)] + MR[GR(NA-1,1,1)]) + MU[GU(NA,1,1)]); //Delta R (Abajo dra) 
	DeltasMU[GU(NA,NR-1,1)] = Rfinal - 0.50*(MU[GU(NA,NR-1,1)] + MU[GU(NA,NR-2,1)]); //Delta R (Arriba dra)

	for(j = 0; j < NR; j++){
		DeltasMU[GU(0,j,0)] = Hyp(MU[GU(0,j,0)], MP[G(0,j,0)], MU[GU(0,j,1)], MP[G(0,j,1)]); //Deltas X (Parte izda)
		DeltasMU[GU(NA,j,0)] = Hyp(MP[G(NA-1,j,0)], MU[GU(NA,j,0)], 0.50*(MR[GR(NA-1,j,1)] + MR[GR(NA-1,j+1,1)]), MU[GR(NA,j,1)]); //Deltas X (Parte dra)

		DeltasMU[GU(0,j,1)] = MR[GR(0,j+1,1)] - MR[GR(0,j,1)]; //Deltas R (Parte izda)	
	}
	DeltasMU[GU(NA,0,0)] = MU[GU(NA,0,0)] - MP[G(NA-1,0,0)]; //Deltas X (Parte dra)

	DeltasMU[GU(NA,1,1)] = 0.50*(MU[GU(NA,2,1)] - 0.50*(MR[GR(NA-1,0,1)] + MR[GR(NA-1,1,1)])); //Deltas R (Parte dra)
	for(j = 2; j < NR-1; j++){
		DeltasMU[GU(NA,j,1)] = 0.50*(MU[GU(NA,j+1,1)] - MU[GU(NA,j-1,1)]); //Deltas R (Parte dra)
	}
	for(i = 1; i < NA; i++){		
		DeltasMU[GU(i,0,0)] = Hyp(MP[G(i-1,0,0)], MP[G(i,0,0)], MP[G(i-1,0,1)], MP[G(i,0,1)]); //Deltas X (Parte abajo)
		DeltasMU[GU(i,NR-1,0)] = Hyp(MP[G(i-1,NR-1,0)], MP[G(i,NR-1,0)], MP[G(i-1,NR-1,1)], MP[G(i,NR-1,1)]); //Deltas X (Parte arriba)

		DeltasMU[GU(i,0,1)] = 0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]) - 0.50*(MR[GR(i-1,0,1)] + MR[GR(i,0,1)]); //Deltas R (Parte abajo)
		DeltasMU[GU(i,NR-1,1)] = 0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]) - 0.50*(MR[GR(i-1,NR-1,1)] + MR[GR(i,NR-1,1)]); //Deltas R (Parte arriba)
		for(j = 1; j < NR-1; j++){
			DeltasMU[GU(i,j,0)] = Hyp(MP[G(i-1,j,0)], MP[G(i,j,0)], MP[G(i-1,j,1)], MP[G(i,j,1)]); //Deltas X (Centro)
			DeltasMU[GU(i,j,1)] =  0.50*(MR[GR(i-1,j+1,1)] + MR[GR(i,j+1,1)]) - 0.50*(MR[GR(i-1,j,1)] + MR[GR(i,j,1)]); //Deltas R (Centro)
		}
	}
	
	//Deltas X y R de la matriz de velocidades radiales (V)
	for(j = 1; j < NR; j++){
		DeltasMR[GR(0,j,0)] = Hyp(MU[GU(0,j,0)], MU[GU(1,j,0)], 0.50*(MU[GU(0,j-1,1)] + MU[GU(0,j,1)]), 0.50*(MU[GU(1,j-1,1)] + MU[GU(1,j,1)])); //Deltas X (Parte izda)
		DeltasMR[GR(NA-1,j,0)] = Hyp(MU[GU(NA-1,j,0)], MU[GU(NA,j,0)], 0.50*(MU[GU(NA-1,j-1,1)] + MU[GU(NA-1,j,1)]), 0.50*(MU[GU(NA,j-1,1)] + MU[GU(NA,j,1)])); //Deltas X (Parte dra)
	}
	DeltasMR[GR(0,NR,0)] = Hyp(MU[GU(0,NR-1,0)], MU[GU(1,NR-1,0)], Rinicio, 0.50*(MR[GR(0,NR,1)] + MR[GR(1,NR,1)]));//Deltas X (Arriba izda)
	DeltasMR[GR(0,0,0)] = MU[GU(1,0,0)] - MU[GU(0,0,0)];//Deltas X (Abajo izda)

	DeltasMR[GR(NA-1,NR,0)] = Hyp(MU[GU(NA-1,NR-1,0)], MU[GU(NA,NR-1,0)], 0.50*(MR[GR(NA-2,NR,1)] + MR[GR(NA-1,NR,1)]), Rfinal);//Deltas X (Arriba dra)
	DeltasMR[GR(NA-1,0,0)] = MU[GU(NA,0,0)] - MU[GU(NA-1,0,0)];//Deltas X (Abajo dra)

	for(i = 1; i < NA-1; i++){ 
		DeltasMR[GR(i,0,0)] = MU[GU(i+1,0,0)] - MU[GU(i,0,0)]; //Deltas X (Parte abajo)
		DeltasMR[GR(i,NR,0)] = Hyp(MU[GU(i,NR-1,0)], MU[GU(i+1,NR-1,0)], 0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]), 0.50*(MR[GR(i,NR,1)] + MR[GR(i+1,NR,1)])); //Deltas X (Parte arriba)
		for(j = 1; j < NR; j++){
			DeltasMR[GR(i,j,0)] = Hyp(MU[GU(i,j,0)], MU[GU(i+1,j,0)], 0.50*(MU[GU(i,j-1,1)] + MU[GU(i,j,1)]), 0.50*(MU[GU(i+1,j-1,1)] + MU[GU(i+1,j,1)])); //Deltas X (Centro)
		}
	}
	for(i = 0; i < NA; i++){
		DeltasMR[GR(i,0,1)] = 0.50*(MR[GR(i,0,1)] + MR[(i,1,1)]) - MR[GR(i,0,1)]; //DeltasR (Parte abajo)
		DeltasMR[GR(i,NR,1)] = MR[GR(i,NR,1)] - MP[G(i,NR-1,1)]; //DeltasR (Parte arriba)
		for(j = 1; j < NR; j++){
			DeltasMR[GR(i,j,1)] = 0.50*(MR[GR(i,j+1,1)] + MR[GR(i,j,1)]) - 0.50*(MR[GR(i,j,1)] + MR[GR(i,j-1,1)]); //DeltasR (Centro)	
		}
	}	

}

//Cálculo de las superficies de cada uno de los volúmenes de control para el caso del canal
void Mesher::Get_SurfacesChannel(){

//Código:
//0 -> West
//1 -> East
//2 -> South
//3 -> North

int i, j;

	//Superficies matriz central
	for(i = 0; i < NA; i++){
		for(j = 0; j < NR; j++){
			SupMP[G(i,j,0)] = ChannelDepth*DeltasMU[GU(i,j,1)]; //West
			SupMP[G(i,j,1)] = ChannelDepth*DeltasMU[GU(i+1,j,1)]; //East
	
			SupMP[G(i,j,2)] = ChannelDepth*DeltasMR[GU(i,j,0)]; //South
			SupMP[G(i,j,3)] = ChannelDepth*DeltasMR[GU(i,j+1,0)]; //North
		}
	}

}

//Cálculo de las superficies de cada uno de los volúmenes de control
void Mesher::Get_Surfaces(){

//Código:
//0 -> West
//1 -> East
//2 -> South
//3 -> North

int i, j;
	
	//Superficies de los volúmenes de control de Presión/Temperatura

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NR-1; j++){
			SupMP[G(i,j,0)] = PI*(pow(0.50*(MR[GR(i,j+1,1)] + MR[GR(i-1,j+1,1)]),2.0) - pow(0.50*(MR[GR(i,j,1)] + MR[GR(i-1,j,1)]),2.0)); //West
			SupMP[G(i,j,1)] = PI*(pow(0.50*(MR[GR(i,j+1,1)] + MR[GR(i+1,j+1,1)]),2.0) - pow(0.50*(MR[GR(i,j,1)] + MR[GR(i+1,j,1)]),2.0)); //East
	
			SupMP[G(i,j,2)] = PI*DeltasMR[GR(i,j,0)]*(0.50*(MR[GR(i,j,1)] + MR[GR(i-1,j,1)]) + 0.50*(MR[GR(i,j,1)] + MR[GR(i+1,j,1)])); //South
			SupMP[G(i,j,3)] = PI*DeltasMR[GR(i,j+1,0)]*(0.50*(MR[GR(i,j+1,1)] + MR[GR(i-1,j+1,1)]) + 0.50*(MR[GR(i,j+1,1)] + MR[GR(i+1,j+1,1)])); //North
		}
	}

	for(j = 1; j < NR-1; j++){
		//Parte izquierda
		SupMP[G(0,j,0)] = PI*(pow(MR[GR(0,j+1,1)],2.0) - pow(MR[GR(0,j,1)],2.0)); //West
		SupMP[G(0,j,1)] = PI*(pow(0.50*(MR[GR(0,j+1,1)] + MR[GR(1,j+1,1)]),2.0) - pow(0.50*(MR[GR(0,j,1)] + MR[GR(1,j,1)]),2.0)); //East
	
		SupMP[G(0,j,2)] = PI*DeltasMR[GR(0,j,0)]*(MR[GR(0,j,1)] + 0.50*(MR[GR(0,j,1)] + MR[GR(1,j,1)])); //South
		SupMP[G(0,j,3)] = PI*DeltasMR[GR(0,j+1,0)]*(0.50*(MU[GU(0,j,1)] + MU[GU(0,j+1,1)]) + 0.50*(MU[GU(1,j,1)] + MU[GU(1,j+1,1)])); //North

		//Parte derecha
		SupMP[G(NA-1,j,0)] = PI*(pow(0.50*(MR[GR(NA-1,j+1,1)] + MR[GR(NA-2,j+1,1)]),2.0) - pow(0.50*(MR[GR(NA-1,j,1)] + MR[GR(NA-2,j,1)]),2.0)); //West
		SupMP[G(NA-1,j,1)] = PI*(pow(MR[GR(NA-1,j+1,1)],2.0) - pow(MR[GR(NA-1,j,1)],2.0)); //East
	
		SupMP[G(NA-1,j,2)] = PI*DeltasMR[GR(NA-1,j,0)]*(0.50*(MR[GR(NA-1,j,1)] + MR[GR(NA-2,j,1)]) + MR[GR(NA-1,j,1)]); //South
		SupMP[G(NA-1,j,3)] = PI*DeltasMR[GR(NA-1,j+1,0)]*(0.50*(MU[GU(NA-1,j,1)] + MU[GU(NA-1,j+1,1)]) + 0.50*(MU[GU(NA,j,1)] + MU[GU(NA,j+1,1)])); //North
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SupMP[G(i,0,0)] = PI*(pow(0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]),2.0) - pow(0.0,2.0)); //West
		SupMP[G(i,0,1)] = PI*(pow(0.50*(MR[GR(i,1,1)] + MR[GR(i+1,1,1)]),2.0) - pow(0.0,2.0)); //East
		
		SupMP[G(i,0,2)] = 0.0; //South
		SupMP[G(i,0,3)] = PI*DeltasMR[GR(i,1,0)]*(0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]) + 0.50*(MR[GR(i,1,1)] + MR[GR(i+1,1,1)])); //North

		//Parte arriba
		SupMP[G(i,NR-1,0)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i-1,NR,1)]),2.0) - pow(0.50*(MU[GU(i,NR-1,1)] + MU[GU(i,NR-2,1)]),2.0)); //West
		SupMP[G(i,NR-1,1)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i+1,NR,1)]),2.0) - pow(0.50*(MU[GU(i+1,NR-1,1)] + MU[GU(i+1,NR-2,1)]),2.0)); //East
	
		SupMP[G(i,NR-1,2)] = PI*DeltasMR[GR(i,NR-1,0)]*(0.50*(MU[GU(i,NR-1,1)] + MU[GU(i,NR-2,1)]) + 0.50*(MU[GU(i+1,NR-1,1)] + MU[GU(i+1,NR-2,1)])); //South
		SupMP[G(i,NR-1,3)] = PI*DeltasMR[GR(i,NR,0)]*(0.50*(MR[GR(i,NR,1)] + MR[GR(i+1,NR,1)]) + 0.50*(MR[GR(i,NR,1)] + MR[GR(i-1,NR,1)])); //North
	}

	//Esquina abajo izquierda
	SupMP[G(0,0,0)] = PI*(pow(MR[GR(0,1,1)],2.0) - pow(0.0,2.0)); //West
	SupMP[G(0,0,1)] = PI*(pow(0.50*(MR[GR(0,1,1)] + MR[GR(1,1,1)]),2.0) - pow(0.0,2.0)); //East
	
	SupMP[G(0,0,2)] = 0.0; //South
	SupMP[G(0,0,3)] = PI*DeltasMR[GR(0,1,0)]*(MR[GR(0,1,1)] + 0.50*(MR[GR(0,1,1)] + MR[GR(1,1,1)]));  //North

	//Esquina arriba izquierda
	SupMP[G(0,NR-1,0)] = PI*(pow(Rinicio,2.0) - pow(0.50*(MU[GU(0,NR-1,1)] + MU[GU(0,NR-2,1)]),2.0)); //West
	SupMP[G(0,NR-1,1)] = PI*(pow(0.50*(MR[GR(0,NR,1)] + MR[GR(1,NR,1)]),2.0) - pow(0.50*(MU[GU(1,NR-1,1)] + MU[GU(1,NR-2,1)]),2.0)); //East
	
	SupMP[G(0,NR-1,2)] = PI*DeltasMR[GR(0,NR-1,0)]*(0.50*(MU[GU(0,NR-1,1)] + MU[GU(0,NR-2,1)]) + 0.50*(MU[GU(1,NR-1,1)] + MU[GU(1,NR-2,1)])); //South
	SupMP[G(0,NR-1,3)] = PI*DeltasMR[GR(0,NR,0)]*(Rinicio + 0.50*(MR[GR(0,NR,1)] + MR[GR(1,NR,1)])); //North

	//Esquina abajo derecha
	SupMP[G(NA-1,0,0)] = PI*(pow(0.50*(MR[GR(NA-1,1,1)] + MR[GR(NA-2,1,1)]),2.0) - pow(0.0,2.0)); //West
	SupMP[G(NA-1,0,1)] = PI*(pow(MR[GR(NA-1,1,1)] + (MU[GU(NA,0,0)] - MR[GR(NA-1,1,0)])*tan(0.50*(AngleMU[GU(NA,0,0)] + AngleMU[GU(NA,1,0)])),2.0) - pow(0.0,2.0)); //East
	
	SupMP[G(NA-1,0,2)] = 0.0; //South
	SupMP[G(NA-1,0,3)] = PI*DeltasMR[GR(NA-1,1,0)]*(0.50*(MR[GR(NA-1,1,1)] + MR[GR(NA-2,1,1)]) + MR[GR(NA-1,1,1)]); //North

	//Esquina arriba derecha
	SupMP[G(NA-1,NR-1,0)] = PI*(pow(0.50*(MR[GR(NA-2,NR,1)] + MR[GR(NA-1,NR,1)]),2.0) - pow(0.50*(MU[GU(NA-1,NR-1,1)] + MU[GU(NA-1,NR-2,1)]),2.0)); //West
	SupMP[G(NA-1,NR-1,1)] = PI*(pow(Rfinal,2.0) - pow(0.50*(MU[GU(NA,NR-1,1)] + MU[GU(NA,NR-2,1)]),2.0)); //East
	
	SupMP[G(NA-1,NR-1,2)] = PI*DeltasMR[GR(NA-1,NR-1,0)]*(0.50*(MU[GU(NA-1,NR-1,1)] + MU[GU(NA-1,NR-2,1)]) + 0.50*(MU[GU(NA,NR-1,1)] + MU[GU(NA,NR-2,1)])); //South
	SupMP[G(NA-1,NR-1,3)] = PI*DeltasMR[GR(NA-1,NR,0)]*(Rfinal + 0.50*(MR[GR(NA-2,NR,1)] + MR[GR(NA-1,NR,1)])); //North

	//Superficies Nodos U
	if(Rank == 0){
		double Ri1 = MR[GR(0, 1, 1)] + ((MR[GR(1, 1, 1)] - MR[GR(0, 1, 1)])/(MR[GR(1, 1, 0)] - MR[GR(0, 1, 0)]))*(Xinicio - MR[GR(0, 1, 0)]);
		for(j = 2; j < NR - 1; j++){
			SupMU[MU(0, j, 0)] = PI*(pow(0.50*(MU[GU(0,j,1)] + MU[GU(0,j+1,1)]),2.0) - pow(0.50*(MU[GU(0,j-1,1)] + MU[GU(0,j,1)]),2.0));
			SupMU[MU(0, j, 1)] = PI*(pow(MR[GR(0,j+1,1)],2.0) - pow(MR[GR(0,j,1)],2.0)); 
			SupMU[MU(0, j, 2)] = 2*PI*(0.50*DeltasMR[GR(0,j,0)])*(0.50*(MR[GR(0,j,1)] + 0.50*(MU[GU(0,j-1,1)] + MU[GU(0,j,1)])));
			SupMU[MU(0, j, 3)] = 2*PI*(0.50*DeltasMR[GR(0,j+1,0)])*(0.50*(MR[GR(0,j+1,1)] + 0.50*(MU[GU(0,j,1)] + MU[GU(0,j+1,1)])));
		}

		//Esquina abajo izquierda
		SupMU[MU(0, 1, 0)] = PI*(pow(0.50*(MU[GU(0,1,1)] + MU[GU(0,2,1)]),2.0) - pow(Ri1,2.0));
		SupMU[MU(0, 1, 1)] = PI*(pow(MR[GR(0,2,1)],2.0) - pow(MR[GR(0,1,1)],2.0)); 
		SupMU[MU(0, 1, 2)] = 2*PI*(0.50*DeltasMR[GR(0,1,0)])*(0.50*(MR[GR(0,1,1)] + Ri1));
		SupMU[MU(0, 1, 3)] = 2*PI*(0.50*DeltasMR[GR(0,2,0)])*(0.50*(MR[GR(0,2,1)] + 0.50*(MU[GU(0,1,1)] + MU[GU(0,2,1)])));
		
		SupMU[MU(0, 0, 0)] = PI*(pow(Ri1,2.0) - 0.0);
		SupMU[MU(0, 0, 1)] = PI*(pow(MR[GR(0,1,1)],2.0) - pow(MR[GR(0,0,1)],2.0)); 
		SupMU[MU(0, 0, 2)] = 2*PI*(0.50*DeltasMR[GR(0,0,0)])*(0.50*(MR[GR(0,0,1)] + 0.0));
		SupMU[MU(0, 0, 3)] = 2*PI*(0.50*DeltasMR[GR(0,1,0)])*(0.50*(MR[GR(0,1,1)] + Ri1));
		
		//Esquina arriba izquierda
		SupMU[MU(0, NR - 1, 0)] = PI*(pow(Rinicio,2.0) - pow(0.50*(MU[GU(0,NR - 2,1)] + MU[GU(0,NR - 1,1)]),2.0));
		SupMU[MU(0, NR - 1, 1)] = PI*(pow(MR[GR(0,NR,1)],2.0) - pow(MR[GR(0,NR - 1,1)],2.0)); 
		SupMU[MU(0, NR - 1, 2)] = 2*PI*(0.50*DeltasMR[GR(0,NR - 1,0)])*(0.50*(MR[GR(0,NR - 1,1)] + 0.50*(MU[GU(0,NR - 2,1)] + MU[GU(0,NR - 1,1)])));
		SupMU[MU(0, NR - 1, 3)] = 2*PI*(0.50*DeltasMR[GR(0,NR,0)])*(0.50*(MR[GR(0,NR,1)] + Rinicio));
		
	}
	//LAS SUPERFICIES DEL CORE 1 ESTÁN VERIFICADAS, ASI COMO TODOS LOS VOLÚMENES
	else if(Rank == Procesos-1){
		double Rna = MR[GR(NA-1, 1, 1)] + ((MR[GR(NA-1, 1, 1)] - MR[GR(NA - 2, 1, 1)])/(MR[GR(NA-1, 1, 0)] - MR[GR(NA - 2, 1, 0)]))*(Xfinal - MR[GR(NA-1, 1, 0)]);
		for(j = 2; j < NR - 1; j++){
			SupMU[MU(0, j, 0)] = PI*(pow(MR[GR(NA - 1,j+1,1)],2.0) - pow(MR[GR(NA - 1,j,1)],2.0)); 
			SupMU[MU(0, j, 1)] = PI*(pow(0.50*(MU[GU(NA,j,1)] + MU[GU(NA,j+1,1)]),2.0) - pow(0.50*(MU[GU(NA,j-1,1)] + MU[GU(NA,j,1)]),2.0));
			SupMU[MU(0, j, 2)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,j,0)])*(0.50*(MR[GR(NA - 1,j,1)] + 0.50*(MU[GU(NA,j-1,1)] + MU[GU(NA,j,1)])));
			SupMU[MU(0, j, 3)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,j+1,0)])*(0.50*(MR[GR(NA - 1,j+1,1)] + 0.50*(MU[GU(NA,j,1)] + MU[GU(NA,j+1,1)])));
		}

		//Esquina abajo derecha
		SupMU[MU(0, 1, 0)] = PI*(pow(MR[GR(NA - 1,2,1)],2.0) - pow(MR[GR(NA - 1,1,1)],2.0)); 
		SupMU[MU(0, 1, 1)] = PI*(pow(0.50*(MU[GU(NA,1,1)] + MU[GU(NA,2,1)]),2.0) - pow(Rna,2.0));
		SupMU[MU(0, 1, 2)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,1,0)])*(0.50*(MR[GR(NA - 1,1,1)] + Rna));
		SupMU[MU(0, 1, 3)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,2,0)])*(0.50*(MR[GR(NA - 1,2,1)] + 0.50*(MU[GU(NA,1,1)] + MU[GU(NA,2,1)])));
		
		SupMU[MU(0, 0, 0)] = PI*(pow(MR[GR(NA - 1,1,1)],2.0) - pow(MR[GR(NA - 1,0,1)],2.0)); 
		SupMU[MU(0, 0, 1)] = PI*(pow(Rna,2.0) - 0.0);
		SupMU[MU(0, 0, 2)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,0,0)])*(0.50*(MR[GR(NA - 1,0,1)] + 0.0));
		SupMU[MU(0, 0, 3)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,1,0)])*(0.50*(MR[GR(NA - 1,1,1)] + Rna));
	
		//Esquina arriba derecha
		SupMU[MU(0, NR - 1, 0)] = PI*(pow(MR[GR(NA - 1,NR,1)],2.0) - pow(MR[GR(NA - 1,NR - 1,1)],2.0)); 
		SupMU[MU(0, NR - 1, 1)] = PI*(pow(Rfinal,2.0) - pow(0.50*(MU[GU(NA,NR - 2,1)] + MU[GU(NA,NR - 1,1)]),2.0));
		SupMU[MU(0, NR - 1, 2)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,NR - 1,0)])*(0.50*(MR[GR(NA - 1,NR - 1,1)] + 0.50*(MU[GU(NA,NR - 2,1)] + MU[GU(NA,NR - 1,1)])));
		SupMU[MU(0, NR - 1, 3)] = 2*PI*(0.50*DeltasMR[GR(NA - 1,NR,0)])*(0.50*(MR[GR(NA - 1,NR,1)] + Rfinal));
		
	}
	
	//Superficies Nodos R
	if(Rank == 0){
		double Ri1 = MR[GR(0, 1, 1)] + ((MR[GR(1, 1, 1)] - MR[GR(0, 1, 1)])/(MR[GR(1, 1, 0)] - MR[GR(0, 1, 0)]))*(Xinicio - MR[GR(0, 1, 0)]);
		for(i = Ix + 1; i < Fx; i++){
			//Parte abajo
			SupMR[MR(i - Ix, 0, 0)] = PI*(pow(0.50*0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]),2.0) - 0.0);
			SupMR[MR(i - Ix, 0, 1)] = PI*(pow(0.50*0.50*(MR[GR(i,1,1)] + MR[GR(i+1,1,1)]),2.0) - 0.0);
			SupMR[MR(i - Ix, 0, 2)] = 0.0;
			SupMR[MR(i - Ix, 0, 3)] = 2.0*PI*DeltasMP[G(i,0,0)]*0.50*(MR[GR(i,0,1)] + MR[GR(i,1,1)]);

			//Parte arriba
			SupMR[MR(i - Ix, 1, 0)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i - 1, NR, 1)]),2.0) - pow(MU[GU(i,NR - 1,1)],2.0));
			SupMR[MR(i - Ix, 1, 1)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i + 1, NR, 1)]),2.0) - pow(MU[GU(i + 1,NR - 1,1)],2.0));
			SupMR[MR(i - Ix, 1, 2)] = 2.0*PI*DeltasMP[G(i,NR - 1,0)]*MP[G(i,NR - 1,1)];
			SupMR[MR(i - Ix, 1, 3)] = 2.0*PI*DeltasMR[GR(i,NR,0)]*MR[GR(i,NR,1)];
		}

		//Esquina abajo izquierda
		SupMR[MR(0, 0, 0)] = PI*(pow(0.50*Ri1,2.0) - 0.0);
		SupMR[MR(0, 0, 1)] = PI*(pow(0.50*0.50*(MR[GR(0,1,1)] + MR[GR(1,1,1)]),2.0) - 0.0);
		SupMR[MR(0, 0, 2)] = 0.0;
		SupMR[MR(0, 0, 3)] = 2.0*PI*DeltasMP[G(0,0,0)]*0.50*(MR[GR(0,0,1)] + MR[GR(0,1,1)]);

		//Esquina arriba izquierda
		SupMR[MR(0, 1, 0)] = PI*(pow(Rinicio,2.0) - pow(MU[GU(0,NR - 1,1)],2.0));
		SupMR[MR(0, 1, 1)] = PI*(pow(0.50*(MR[GR(0,NR,1)] + MR[GR(1, NR, 1)]),2.0) - pow(MU[GU(1,NR - 1,1)],2.0));
		SupMR[MR(0, 1, 2)] = 2.0*PI*DeltasMP[G(0,NR - 1,0)]*MP[G(0,NR - 1,1)];
		SupMR[MR(0, 1, 3)] = 2.0*PI*DeltasMR[GR(0,NR,0)]*MR[GR(0,NR,1)];

	}
	else if(Rank == Procesos - 1){
		double Rna = MR[GR(NA-1, 1, 1)] + ((MR[GR(NA-1, 1, 1)] - MR[GR(NA - 2, 1, 1)])/(MR[GR(NA-1, 1, 0)] - MR[GR(NA - 2, 1, 0)]))*(Xfinal - MR[GR(NA-1, 1, 0)]);
		for(i = Ix; i < Fx - 1; i++){
			//Parte abajo
			SupMR[MR(i - Ix, 0, 0)] = PI*(pow(0.50*0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]),2.0) - 0.0);
			SupMR[MR(i - Ix, 0, 1)] = PI*(pow(0.50*0.50*(MR[GR(i,1,1)] + MR[GR(i+1,1,1)]),2.0) - 0.0);
			SupMR[MR(i - Ix, 0, 2)] = 0.0;
			SupMR[MR(i - Ix, 0, 3)] = 2.0*PI*DeltasMP[G(i,0,0)]*0.50*(MR[GR(i,0,1)] + MR[GR(i,1,1)]);

			//Parte arriba
			SupMR[MR(i - Ix, 1, 0)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i - 1, NR, 1)]),2.0) - pow(MU[GU(i,NR - 1,1)],2.0));
			SupMR[MR(i - Ix, 1, 1)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i + 1, NR, 1)]),2.0) - pow(MU[GU(i + 1,NR - 1,1)],2.0));
			SupMR[MR(i - Ix, 1, 2)] = 2.0*PI*DeltasMP[G(i,NR - 1,0)]*MP[G(i,NR - 1,1)];
			SupMR[MR(i - Ix, 1, 3)] = 2.0*PI*DeltasMR[GR(i,NR,0)]*MR[GR(i,NR,1)];
		}

		//Esquina abajo derecha
		SupMR[MR(NA - 1 - Ix, 0, 0)] = PI*(pow(0.50*0.50*(MR[GR(NA-1,1,1)] + MR[GR(NA-2,1,1)]),2.0) - 0.0);
		SupMR[MR(NA - 1 - Ix, 0, 1)] = PI*(pow(0.50*Rna,2.0) - 0.0);
		SupMR[MR(NA - 1 - Ix, 0, 2)] = 0.0;
		SupMR[MR(NA - 1 - Ix, 0, 3)] = 2.0*PI*DeltasMP[G(NA - 1,0,0)]*0.50*(MR[GR(NA - 1,0,1)] + MR[GR(NA - 1,1,1)]);

		//Esquina arriba derecha
		SupMR[MR(NA - 1 - Ix, 1, 0)] = PI*(pow(0.50*(MR[GR(NA - 1,NR,1)] + MR[GR(NA - 2, NR, 1)]),2.0) - pow(MU[GU(NA - 1,NR - 1,1)],2.0));
		SupMR[MR(NA - 1 - Ix, 1, 1)] = PI*(pow(Rfinal,2.0) - pow(MU[GU(NA,NR - 1,1)],2.0));
		SupMR[MR(NA - 1 - Ix, 1, 2)] = 2.0*PI*DeltasMP[G(NA - 1,NR - 1,0)]*MP[G(NA - 1,NR - 1,1)];
		SupMR[MR(NA - 1 - Ix, 1, 3)] = 2.0*PI*DeltasMR[GR(NA - 1,NR,0)]*MR[GR(NA - 1,NR,1)];
	}
	else{
		for(i = Ix; i < Fx; i++){
			//Parte abajo
			SupMR[MR(i - Ix, 0, 0)] = PI*(pow(0.50*0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]),2.0) - 0.0);
			SupMR[MR(i - Ix, 0, 1)] = PI*(pow(0.50*0.50*(MR[GR(i,1,1)] + MR[GR(i+1,1,1)]),2.0) - 0.0);
			SupMR[MR(i - Ix, 0, 2)] = 0.0;
			SupMR[MR(i - Ix, 0, 3)] = 2.0*PI*DeltasMP[G(i,0,0)]*0.50*(MR[GR(i,0,1)] + MR[GR(i,1,1)]);

			//Parte arriba
			SupMR[MR(i - Ix, 1, 0)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i - 1, NR, 1)]),2.0) - pow(MU[GU(i,NR - 1,1)],2.0));
			SupMR[MR(i - Ix, 1, 1)] = PI*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i + 1, NR, 1)]),2.0) - pow(MU[GU(i + 1,NR - 1,1)],2.0));
			SupMR[MR(i - Ix, 1, 2)] = 2.0*PI*DeltasMP[G(i,NR - 1,0)]*MP[G(i,NR - 1,1)];
			SupMR[MR(i - Ix, 1, 3)] = 2.0*PI*DeltasMR[GR(i,NR,0)]*MR[GR(i,NR,1)];
		}
	}

}

//Cálculo de los volúmenes de control de cada volúmen para el caso del canal
void Mesher::Get_VolumesChannel(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NR; j++){
			VolMP[G(i,j,0)] = DeltasMP[G(i,j,0)]*DeltasMP[G(i,j,1)]*ChannelDepth;
		}
	}
}

//Cálculo de los volúmenes de control de cada volúmen
void Mesher::Get_Volumes(){
int i, j;

	
	//Cálculo de los volúmenes de la matriz de Presión/Temperatura

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NR-1; j++){
			VolMP[G(i,j,0)] = (MU[GU(i+1,j,0)] - MU[GU(i,j,0)])*(PI/3.0)*(pow(0.50*(MR[GR(i+1,j+1,1)] + MR[GR(i,j+1,1)]),2.0) + pow(0.50*(MR[GR(i,j+1,1)] + MR[GR(i-1,j+1,1)]),2.0) + 0.50*(MR[GR(i+1,j+1,1)] + MR[GR(i,j+1,1)])*0.50*(MR[GR(i,j+1,1)] + MR[GR(i-1,j+1,1)]) - pow(0.50*(MR[GR(i+1,j,1)] + MR[GR(i,j,1)]),2.0) - pow(0.50*(MR[GR(i,j,1)] + MR[GR(i-1,j,1)]),2.0) - 0.50*(MR[GR(i+1,j,1)] + MR[GR(i,j,1)])*0.50*(MR[GR(i,j,1)] + MR[GR(i-1,j,1)]));
		}
	}

	//Parte abajo
	for(i = 1; i < NA-1; i++){
		VolMP[G(i,0,0)] = (MU[GU(i+1,0,0)] - MU[GU(i,0,0)])*(PI/3.0)*(pow(0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)]),2.0) + pow(0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]),2.0) + 0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)])*0.50*(MR[GR(i,1,1)] + MR[GR(i-1,1,1)]));
	}

	//Parte arriba
	for(i = 1; i < NA-1; i++){
		VolMP[G(i,NR-1,0)] = (MU[GU(i+1,NR-1,0)] - MU[GU(i,NR-1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i+1,NR,1)]),2.0) + pow(0.50*(MR[GR(i,NR,1)] + MR[GR(i-1,NR,1)]),2.0) + 0.50*(MR[GR(i,NR,1)] + MR[GR(i+1,NR,1)])*0.50*(MR[GR(i,NR,1)] + MR[GR(i-1,NR,1)]) - pow(0.50*(MU[GU(i+1,NR-1,1)] + MU[GU(i+1,NR-2,1)]),2.0) - pow(0.50*(MU[GU(i,NR-1,1)] + MU[GU(i,NR-2,1)]),2.0) - 0.50*(MU[GU(i+1,NR-1,1)] + MU[GU(i+1,NR-2,1)])*0.50*(MU[GU(i,NR-1,1)] + MU[GU(i,NR-2,1)]));
	}

	for(j = 1; j < NR-1; j++){
		//Parte izquierda
		VolMP[G(0,j,0)] = (MU[GU(1,j,0)] - MU[GU(0,j,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,j+1,1)] + MR[GR(0,j+1,1)]),2.0) + pow(MR[GR(0,j+1,1)],2.0) + 0.50*(MR[GR(1,j+1,1)] + MR[GR(0,j+1,1)])*MR[GR(0,j+1,1)] - pow(0.50*(MR[GR(1,j,1)] + MR[GR(0,j,1)]),2.0) - pow(MR[GR(0,j,1)],2.0) - 0.50*(MR[GR(1,j,1)] + MR[GR(0,j,1)])*MR[GR(i,j,1)]);

		//Parte derecha
		VolMP[G(NA-1,j,0)] = (MU[GU(NA,j,0)] - MU[GU(NA-1,j,0)])*(PI/3.0)*(pow(MR[GR(NA-1,j+1,1)],2.0) + pow(0.50*(MR[GR(NA-1,j+1,1)] + MR[GR(NA-2,j+1,1)]),2.0) + MR[GR(NA-1,j+1,1)]*0.50*(MR[GR(NA-1,j+1,1)] + MR[GR(NA-2,j+1,1)]) - pow(MR[GR(NA-1,j,1)],2.0) - pow(0.50*(MR[GR(NA-1,j,1)] + MR[GR(NA-2,j,1)]),2.0) - MR[GR(NA-1,j,1)]*0.50*(MR[GR(NA-1,j,1)] + MR[GR(NA-2,j,1)]));
	}

	//Esquina arriba izquierda
	VolMP[G(0,NR-1,0)] = (MU[GU(1,NR-1,0)] - MU[GU(0,NR-1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(0,NR,1)] + MR[GR(1,NR,1)]),2.0) + pow(Rinicio,2.0) + 0.50*(MR[GR(0,NR,1)] + MR[GR(1,NR,1)])*Rinicio - pow(0.50*(MU[GU(1,NR-1,1)] + MU[GU(1,NR-2,1)]),2.0) - pow(0.50*(MU[GU(0,NR-1,1)] + MU[GU(0,NR-2,1)]),2.0) - 0.50*(MU[GU(1,NR-1,1)] + MU[GU(1,NR-2,1)])*0.50*(MU[GU(0,NR-1,1)] + MU[GU(0,NR-2,1)]));

	//Esquina arriba derecha
	VolMP[G(NA-1,NR-1,0)] = (MU[GU(NA,NR-1,0)] - MU[GU(NA-1,NR-1,0)])*(PI/3.0)*(pow(Rfinal,2.0) + pow(0.50*(MR[GR(NA-1,NR,1)] + MR[GR(NA-2,NR,1)]),2.0) + Rfinal*0.50*(MR[GR(NA-1,NR,1)] + MR[GR(NA-2,NR,1)]) - pow(0.50*(MU[GU(NA,NR-1,1)] + MU[GU(NA,NR-2,1)]),2.0) - pow(0.50*(MU[GU(NA-1,NR-1,1)] + MU[GU(NA-1,NR-2,1)]),2.0) - 0.50*(MU[GU(NA,NR-1,1)] + MU[GU(NA,NR-2,1)])*0.50*(MU[GU(NA,NR-1,1)] + MU[GU(NA,NR-2,1)]));

	//Esquina abajo izquierda
	VolMP[G(0,0,0)] = (MU[GU(1,0,0)] - MU[GU(0,0,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)]),2.0) + pow(MR[GR(0,1,1)],2.0) + 0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)])*MR[GR(0,1,1)] - pow(0.50*(MR[GR(1,0,1)] + MR[GR(0,0,1)]),2.0) - pow(MR[GR(0,0,1)],2.0) - 0.50*(MR[GR(1,0,1)] + MR[GR(0,0,1)])*MR[GR(0,0,1)]);

	//Esquina abajo derecha
	VolMP[G(NA-1,0,0)] = (MU[GU(NA,0,0)] - MU[GU(NA-1,0,0)])*(PI/3.0)*(pow(MR[GR(NA-1,1,1)],2.0) + pow(0.50*(MR[GR(NA-1,1,1)] + MR[GR(NA-2,1,1)]),2.0) + MR[GR(NA-1,1,1)]*0.50*(MR[GR(NA-1,1,1)] + MR[GR(NA-2,1,1)]) - pow(MR[GR(NA-1,0,1)],2.0) - pow(0.50*(MR[GR(NA-1,0,1)] + MR[GR(NA-2,0,1)]),2.0) - MR[GR(NA-1,0,1)]*0.50*(MR[GR(NA-1,0,1)] + MR[GR(NA-2,0,1)]));

	//Volúmenes de control Nodos U
	if(Rank == 0){
		double Ri1 = MR[GR(0, 1, 1)] + ((MR[GR(1, 1, 1)] - MR[GR(0, 1, 1)])/(MR[GR(1, 1, 0)] - MR[GR(0, 1, 0)]))*(Xinicio - MR[GR(0, 1, 0)]);
		for(j = 2; j < NR - 1; j++){
			//Parte centro
			VolMU[j] = 
				 	 + (MU[GU(1,j,0)] - MU[GU(0,j,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,j+1,1)] + MR[GR(0,j+1,1)]),2.0) + pow(0.50*(MU[GU(0,j+1,1)] + MU[GU(0,j,1)]),2.0) + 0.50*(MR[GR(1,j+1,1)] + MR[GR(0,j+1,1)])*0.50*(MU[GU(0,j+1,1)] + MU[GU(0,j,1)]))
				 	 - (MU[GU(1,j,0)] - MU[GU(0,j,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,j,1)] + MR[GR(0,j,1)]),2.0) + pow(0.50*(MU[GU(0,j,1)] + MU[GU(0,j-1,1)]),2.0) + 0.50*(MR[GR(1,j,1)] + MR[GR(0,j,1)])*0.50*(MU[GU(0,j,1)] + MU[GU(0,j-1,1)]))
				 	 ;
		}

		//Esquina abajo izquierda
		VolMU[1] = 
				 + (MU[GU(1,1,0)] - MU[GU(0,1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,2,1)] + MR[GR(0,2,1)]),2.0) + pow(0.50*(MU[GU(0,2,1)] + MU[GU(0,1,1)]),2.0) + 0.50*(MR[GR(1,2,1)] + MR[GR(0,2,1)])*0.50*(MU[GU(0,2,1)] + MU[GU(0,1,1)]))
				 - (MU[GU(1,1,0)] - MU[GU(0,1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)]),2.0) + pow(Ri1,2.0) + 0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)])*Ri1)
				 ;

		VolMU[0] = 
				 + (MU[GU(1,0,0)] - MU[GU(0,0,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)]),2.0) + pow(Ri1,2.0) + 0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)])*Ri1)
				 - 0.0
				 ;

		//Esquina arriba izquierda
		VolMU[NR - 1] = 
				 + (MU[GU(1,NR - 1,0)] - MU[GU(0,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,NR,1)] + MR[GR(0,NR,1)]),2.0) + pow(Rinicio,2.0) + 0.50*(MR[GR(1,NR,1)] + MR[GR(0,NR,1)])*Rinicio)
				 - (MU[GU(1,NR - 1,0)] - MU[GU(0,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,NR - 1,1)] + MR[GR(0,NR - 1,1)]),2.0) + pow(0.50*(MU[GU(0,NR - 1,1)] + MU[GU(0,NR - 2,1)]),2.0) + 0.50*(MR[GR(1,NR - 1,1)] + MR[GR(0,NR - 1,1)])*0.50*(MU[GU(0,NR - 1,1)] + MU[GU(0,NR - 2,1)]))
				 ;

	}
	else if(Rank == Procesos - 1){
		double Rna = MR[GR(NA-1, 1, 1)] + ((MR[GR(NA-1, 1, 1)] - MR[GR(NA - 2, 1, 1)])/(MR[GR(NA-1, 1, 0)] - MR[GR(NA - 2, 1, 0)]))*(Xfinal - MR[GR(NA-1, 1, 0)]);
	
		for(j = 2; j < NR - 1; j++){
			//Parte centro
			VolMU[j] = 
					 + (MU[GU(NA,j,0)] - MU[GU(NA - 1,j,0)])*(PI/3.0)*(pow(0.50*(MU[GU(NA,j+1,1)] + MU[GU(NA,j,1)]),2.0) + pow(0.50*(MR[GR(NA - 1,j+1,1)] + MR[GR(NA - 2,j+1,1)]),2.0) + 0.50*(MU[GU(NA,j+1,1)] + MU[GU(NA,j,1)])*0.50*(MR[GR(NA - 1,j+1,1)] + MR[GR(NA - 2,j+1,1)]))
					 - (MU[GU(NA,j,0)] - MU[GU(NA - 1,j,0)])*(PI/3.0)*(pow(0.50*(MU[GU(NA,j,1)] + MU[GU(NA,j-1,1)]),2.0) + pow(0.50*(MR[GR(NA - 1,j,1)] + MR[GR(NA - 2,j,1)]),2.0) + 0.50*(MU[GU(NA,j,1)] + MU[GU(NA,j-1,1)])*0.50*(MR[GR(NA - 1,j,1)] + MR[GR(NA - 2,j,1)]))
					 ;
		}

		//Esquina abajo derecha
		VolMU[1] = 
				 + (MU[GU(NA,1,0)] - MU[GU(NA - 1,1,0)])*(PI/3.0)*(pow(0.50*(MU[GU(NA,2,1)] + MU[GU(NA,1,1)]),2.0) + pow(0.50*(MR[GR(NA - 1,2,1)] + MR[GR(NA - 2,2,1)]),2.0) + 0.50*(MU[GU(NA,2,1)] + MU[GU(NA,1,1)])*0.50*(MR[GR(NA - 1,2,1)] + MR[GR(NA - 2,2,1)]))
				 - (MU[GU(NA,1,0)] - MU[GU(NA - 1,1,0)])*(PI/3.0)*(pow(Rna,2.0) + pow(0.50*(MR[GR(NA - 1,1,1)] + MR[GR(NA - 2,1,1)]),2.0) + Rna*0.50*(MR[GR(NA - 1,1,1)] + MR[GR(NA - 2,1,1)]))
				 ;

		VolMU[0] = 
				 + (MU[GU(NA,0,0)] - MU[GU(NA - 1,0,0)])*(PI/3.0)*(pow(Rna,2.0) + pow(0.50*(MR[GR(NA - 1,1,1)] + MR[GR(NA - 2,1,1)]),2.0) + Rna*0.50*(MR[GR(NA - 1,1,1)] + MR[GR(NA - 2,1,1)]))
				 - 0.0
				 ;

		//Esquina arriba derecha
		VolMU[NR - 1] = 
				 + (MU[GU(NA,NR - 1,0)] - MU[GU(NA - 1,NR - 1,0)])*(PI/3.0)*(pow(Rfinal,2.0) + pow(0.50*(MR[GR(NA - 1,NR,1)] + MR[GR(NA - 2,NR,1)]),2.0) + Rfinal*0.50*(MR[GR(NA - 1,NR,1)] + MR[GR(NA - 2,NR,1)]))
				 - (MU[GU(NA,NR - 1,0)] - MU[GU(NA - 1,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MU[GU(NA,NR - 1,1)] + MU[GU(NA,NR - 2,1)]),2.0) + pow(0.50*(MR[GR(NA - 1,NR - 1,1)] + MR[GR(NA - 2,NR - 1,1)]),2.0) + 0.50*(MU[GU(NA,NR - 1,1)] + MU[GU(NA,NR - 2,1)])*0.50*(MR[GR(NA - 1,NR - 1,1)] + MR[GR(NA - 2,NR - 1,1)]))
				 ;
	}

	//Volúmenes de control Nodos R
	if(Rank == 0){
		double Ri1 = MR[GR(0, 1, 1)] + ((MR[GR(1, 1, 1)] - MR[GR(0, 1, 1)])/(MR[GR(1, 1, 0)] - MR[GR(0, 1, 0)]))*(Xinicio - MR[GR(0, 1, 0)]);
		for(i = Ix + 1; i < Fx; i++){
			//Parte abajo
			VolMR[VR(i - Ix, 0, 0)] = 
									+ (MU[GU(i + 1,0,0)] - MU[GU(i,0,0)])*(PI/3.0)*(pow(0.50*0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)]),2.0) + pow(0.50*0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]),2.0) + 0.50*0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)])*0.50*0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]))
						 			- 0.0
						 			;

			//Parte arriba
			VolMR[VR(i - Ix, 1, 0)] = 
						  			+ (MU[GU(i + 1,NR - 1,0)] - MU[GU(i,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(i+1,NR,1)] + MR[GR(i,NR,1)]),2.0) + pow(0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]),2.0) + 0.50*(MR[GR(i+1,NR,1)] + MR[GR(i,NR,1)])*0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]))
						 			- (MU[GU(i + 1,NR - 1,0)] - MU[GU(i,NR - 1,0)])*(PI/3.0)*(pow(MU[GU(i+1, NR - 1, 1)],2.0) + pow(MU[GU(i, NR - 1, 1)],2.0) + MU[GU(i, NR - 1, 1)]*MU[GU(i+1, NR - 1, 1)])
						 			;
		}

		//Esquina abajo izquierda
		VolMR[VR(0, 0, 0)] = 
						   + (MU[GU(1,0,0)] - MU[GU(0,0,0)])*(PI/3.0)*(pow(0.50*0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)]),2.0) + pow(0.50*Ri1,2.0) + 0.50*0.50*(MR[GR(1,1,1)] + MR[GR(0,1,1)])*0.50*Ri1)
						   - 0.0
						   ;

		//Esquina arriba izquierda
		VolMR[VR(0, 1, 0)] = 
						   + (MU[GU(1,NR - 1,0)] - MU[GU(0,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(1,NR,1)] + MR[GR(0,NR,1)]),2.0) + pow(Rinicio,2.0) + 0.50*(MR[GR(1,NR,1)] + MR[GR(0,NR,1)])*Rinicio)
						   - (MU[GU(1,NR - 1,0)] - MU[GU(0,NR - 1,0)])*(PI/3.0)*(pow(MU[GU(1, NR - 1, 1)],2.0) + pow(MU[GU(0, NR - 1, 1)],2.0) + MU[GU(0, NR - 1, 1)]*MU[GU(1, NR - 1, 1)])
						   ;

	}
	else if(Rank == Procesos - 1){
		double Rna = MR[GR(NA-1, 1, 1)] + ((MR[GR(NA-1, 1, 1)] - MR[GR(NA - 2, 1, 1)])/(MR[GR(NA-1, 1, 0)] - MR[GR(NA - 2, 1, 0)]))*(Xfinal - MR[GR(NA-1, 1, 0)]);
		for(i = Ix; i < Fx - 1; i++){
			//Parte abajo
			VolMR[VR(i - Ix, 0, 0)] = 
									+ (MU[GU(i + 1,0,0)] - MU[GU(i,0,0)])*(PI/3.0)*(pow(0.50*0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)]),2.0) + pow(0.50*0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]),2.0) + 0.50*0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)])*0.50*0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]))
						 			- 0.0
						 			;

			//Parte arriba
			VolMR[VR(i - Ix, 1, 0)] = 
						  			+ (MU[GU(i + 1,NR - 1,0)] - MU[GU(i,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(i+1,NR,1)] + MR[GR(i,NR,1)]),2.0) + pow(0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]),2.0) + 0.50*(MR[GR(i+1,NR,1)] + MR[GR(i,NR,1)])*0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]))
						 			- (MU[GU(i + 1,NR - 1,0)] - MU[GU(i,NR - 1,0)])*(PI/3.0)*(pow(MU[GU(i+1, NR - 1, 1)],2.0) + pow(MU[GU(i, NR - 1, 1)],2.0) + MU[GU(i, NR - 1, 1)]*MU[GU(i+1, NR - 1, 1)])
						 			;

		}

		//Esquina abajo derecha
		VolMR[VR(NA - 1 - Ix, 0, 0)] = 
									 + (MU[GU(NA,0,0)] - MU[GU(NA - 1,0,0)])*(PI/3.0)*(pow(0.50*Rna,2.0) + pow(0.50*0.50*(MR[GR(NA - 2,1,1)] + MR[GR(NA - 1,1,1)]),2.0) + 0.50*Rna*0.50*0.50*(MR[GR(NA - 2,1,1)] + MR[GR(NA - 1,1,1)]))
						 			 - 0.0
						 			 ;

		//Esquina arriba derecha
		VolMR[VR(NA - 1 - Ix, 1, 0)] = 
						  			+ (MU[GU(NA,NR - 1,0)] - MU[GU(NA - 1,NR - 1,0)])*(PI/3.0)*(pow(Rfinal,2.0) + pow(0.50*(MR[GR(NA - 2,NR,1)] + MR[GR(NA - 1,NR,1)]),2.0) + Rfinal*0.50*(MR[GR(NA - 2,NR,1)] + MR[GR(NA - 1,NR,1)]))
						 			- (MU[GU(NA,NR - 1,0)] - MU[GU(NA - 1,NR - 1,0)])*(PI/3.0)*(pow(MU[GU(NA, NR - 1, 1)],2.0) + pow(MU[GU(NA - 1, NR - 1, 1)],2.0) + MU[GU(NA - 1, NR - 1, 1)]*MU[GU(NA, NR - 1, 1)])
						 			;
	}
	else{

		for(i = Ix; i < Fx; i++){
			//Parte abajo
			VolMR[VR(i - Ix, 0, 0)] = 
									+ (MU[GU(i + 1,0,0)] - MU[GU(i,0,0)])*(PI/3.0)*(pow(0.50*0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)]),2.0) + pow(0.50*0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]),2.0) + 0.50*0.50*(MR[GR(i+1,1,1)] + MR[GR(i,1,1)])*0.50*0.50*(MR[GR(i-1,1,1)] + MR[GR(i,1,1)]))
						 			- 0.0
						 			;

			//Parte arriba
			VolMR[VR(i - Ix, 1, 0)] = 0.0
						  			+ (MU[GU(i + 1,NR - 1,0)] - MU[GU(i,NR - 1,0)])*(PI/3.0)*(pow(0.50*(MR[GR(i+1,NR,1)] + MR[GR(i,NR,1)]),2.0) + pow(0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]),2.0) + 0.50*(MR[GR(i+1,NR,1)] + MR[GR(i,NR,1)])*0.50*(MR[GR(i-1,NR,1)] + MR[GR(i,NR,1)]))
						 			- (MU[GU(i + 1,NR - 1,0)] - MU[GU(i,NR - 1,0)])*(PI/3.0)*(pow(MU[GU(i+1, NR - 1, 1)],2.0) + pow(MU[GU(i, NR - 1, 1)],2.0) + MU[GU(i, NR - 1, 1)]*MU[GU(i+1, NR - 1, 1)])
						 			;

		}

	}
	
}

//Cálculo de los ángulos entre las velocidades y las superficies de los volúmenes de control
void Mesher::Get_Angles(){
int i, j;

	//Ángulos entre las superficies y las velocidades radiales (V)	
	for(i = 0; i < NA; i++){
		AngleMR[G(i,0,0)] = 0.0; //Parte abajo
		AngleMR[G(i,NR,0)] = atan((Radius[i+1] - Radius[i])/(AxialCoord[i+1] - AxialCoord[i])); //Parte arriba
		for(j = 1; j < NR; j++){
			AngleMR[GR(i, j, 0)] = atan((0.50 * (MU[GU(i + 1, j, 1)] + MU[GU(i + 1, j + 1, 1)]) - 0.50 * (MU[GU(i, j, 1)] + MU[GU(i, j + 1, 1)])) / (0.50 * (MU[GU(i + 1, j, 0)] + MU[GU(i + 1, j + 1, 0)]) - 0.50 * (MU[GU(i, j, 0)] + MU[GU(i, j + 1, 0)]))); //Centro
		}
	}

	
	//Ángulos entre las velocidades axiales (U)....	
	for(j = 0; j < NR; j++){
		for(i = 0; i < NA+1; i++){
			AngleMU[GU(i,j,0)] = 0.0; //Centro
		}
	}
	
}

//Cálculo de la longitud característica de la geometría
void Mesher::Get_CharLength(){
int i, j;
double TotalVolume = 0.0;

	if(Problema == "Tobera"){
		for(i = 0; i < NA; i++){ 
			for(j = 0; j < NR; j++){
				TotalVolume += VolMP[G(i,j,0)];
			}
		}

		Lcar = TotalVolume/ThroatArea; //Segun Sutton (Rocket Propulsion elements)
	}
	else if(Problema == "Tuberia"){ Lcar = PipeDiameter; }
		
}

//Cálculo de la distancia mínima de cada nodo a la pared más cercana (Para el Spalart-Allmaras)
void Mesher::ClosestDist(){
int i, j, k, m;
double Minimo;
	
	for(i = 0; i < NA; i++){
		for(j = 0; j < NR; j++){
			Minimo = 2.0*Xfinal;
			if(Problema == "Tobera"){
				for(k = 0; k < NR; k++){
					if(Hyp(MP[G(i,j,0)], MU[GU(0,k,0)], MP[G(i,j,1)], MU[GU(0,k,1)]) <= Minimo){
						minDist[G(i,j,0)] = Minimo;
					}
				}
				for(m = 0; m < NA; m++){
					if(Hyp(MP[G(i,j,0)], MR[GR(m,NR,0)], MP[G(i,j,1)], MR[GR(m,NR,1)]) < Minimo){
						Minimo = Hyp(MP[G(i,j,0)], MR[GR(m,NR,0)], MP[G(i,j,1)], MR[GR(m,NR,1)]);	
					}
				}
			}
			else if(Problema == "Tuberia"){
				for(m = 0; m < NA; m++){
					if(Hyp(MP[G(i,j,0)], MR[GR(m,NR,0)], MP[G(i,j,1)], MR[GR(m,NR,1)]) < Minimo){
						Minimo = Hyp(MP[G(i,j,0)], MR[GR(m,NR,0)], MP[G(i,j,1)], MR[GR(m,NR,1)]);	
					}
				}
			}

			minDist[G(i,j,0)] = Minimo;
		}
	}



}

//Pasar los resultados a un archivo VTK en 2D
void Mesher::MallaVTK2D(string Carpeta, string Variable, string NombreFile, double *MC, int Na, int Nr){
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
        file<<"DIMENSIONS"<<"   "<<Na<<"   "<<Nr<<"   "<<1<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<Na*Nr<<"   "<<"double"<<endl;
	
	for(j = 0; j < Nr; j++){
		for(i = 0; i < Na; i++){
			file<<MC[G(i,j,0)]<<"   "<<MC[G(i,j,1)]<<"   "<<0.0<<endl;
			
		}
		
	}
        
        file<<endl;
		file<<"POINT_DATA"<<"   "<<NA*NR<<endl;
        file<<"SCALARS "<<Variable<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
        file<<endl;
        for(j = 0; j < Nr; j++){
		for(i = 0; i < Na; i++){
				file<<0.0<<" ";
			
		}
	}
        file.close();
}

//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1, Geometry G1, ParPro MPI1){
int i;

	if(Problema == "Tobera"){	
		GetTotalAxialNodes(); //Cálculo del número de nodos totales en función de las opciones de discretización	
	}
	
	//Matrices de coordenadas de discretización
	AxialCoord = M1.AllocateDouble(NA+1, 1, 1); //Matriz coordenada axial de la distribución
	Radius = M1.AllocateDouble(NA+1, 1, 1); //Matriz radios en cada coordenada axial

	MPI1.Initial_WorkSplit(NA, Ix, Fx);

	Get_AxialCoordinates(); //Cálculo del la posición axial de los nodos de velocidad axial

	Get_Radius(M1, G1); //Cálculo del radio de cada coordenada axial discretizada

	GetTotalRadialNodes(); //Calcular el número total de nodos en la dirección radial
	
	AllocateMemory(M1); //Alojamiento de memoria para cada matriz 
	
	Get_Mesh(); //Creación de todas las mallas
	
	if(Problema == "Tuberia" || Problema == "Tobera"){
		Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
	
		Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control

		Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen
	}
	else if(Problema == "Canal"){
		Get_DeltasChannel(); //Cálculo de las distancias entre nodos en cada una de las matrices para el caso del canal
	
		Get_SurfacesChannel(); //Cálculo de las superficies de cada uno de los volúmenes de control para el caso del canal

		Get_VolumesChannel(); //Cálculo de los volúmenes de control de cada volúmen para el caso del canal
	}
	
	
	Get_Angles(); //Cálculo de los ángulos entre las velocidades y las superficies de los volúmenes de control

	
	
	ClosestDist(); //Cálculo de la distancia mínima de cada nodo a la pared más cercana
	
	Get_CharLength(); //Cálculo de la longitud característica de la geometría


	if(Rank == 1){	
	/*	for(int j = NR -1; j >= 0; j--){
			for(int i = 0; i < NA; i++){
				cout<<SupMP[G(i,j,3)]<<", ";
			}
			cout<<endl;
		}*/
			
		
		PrintTxt(); //Pasar los resultados de las mallas a un txt
		MallaVTK2D("ParaviewResults/MeshResults/", "MallaPresion", "MallaMP", MP, NA, NR);
		cout<<"Mesh created."<<endl;
	}

}

