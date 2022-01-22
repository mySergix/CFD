#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <cassert>
#include <time.h>

using namespace std;

#include "/home_nobck/sergiogus/TFG/Codes/HeaderCodes/Memory.h"
#include "/home_nobck/sergiogus/TFG/Codes/HeaderCodes/ReadData.h"
#include "/home_nobck/sergiogus/TFG/Codes/HeaderCodes/Geometry.h"
#include "/home_nobck/sergiogus/TFG/Codes/HeaderCodes/Mesher.h"

#define PI 3.141592653589793

#define sind(x) sin(x * (PI/180.0)) //Cálculo seno en grados
#define cosd(x) cos(x * (PI/180.0)) //Cálculo coseno en grados
#define tand(x) tan(x * (PI/180.0)) //Cálculo tangente en grados

#define Hyp(x1, x2, y1, y2) sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0)) //Cálculo hipotenusa


//Constructor del mallador
Mesher::Mesher(Memory M1, ReadData R1, Geometry G1){
		
	//Datos sobre el problema
		Problema = R1.Get_ProblemType(); //Problema (Tobera/Tubería)
		Tobera = R1.Get_NozzleType();	//Tipo de tobera (NASA/RAO)
	
	//Datos sobre la geometría del problema

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
		NRad = R1.NumericalData[5]; //Dirección radial
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
	
}

//Calcular el número total de nodos en la dirección radial
void Mesher::GetTotalRadialNodes(){
int i, j;
double Sumatorio = 0.0;
double Rplus;

	if(OpcionYplus == 0){ //No activado
		NRad = NRad;
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
			NRad = Rplus/(2.0*DeltaYplus);
		}
		else if(OpcionR == 2){ //Tangencial hiperbólica
			NRad = StretFactorR/(StretFactorR - atanh((1.0 - (2.0*DeltaYplus)/Rplus)*tanh(StretFactorR)));
		}
		else if(OpcionR == 3){ //Senoidal
			NRad = (PI/2.0)/acos(1.0 - (DeltaYplus)/Rplus);
		}
	}

}

//Alojamiento de memoria para cada matriz
void Mesher::AllocateMemory(Memory M1){

	//Matrices de coordenadas de discretización
	MP = M1.AllocateDouble3D(NA, NRad, 2, "P"); //Coordenadas matriz de presión/temperatura
	MU = M1.AllocateDouble3D(NA, NRad, 2, "U"); //Coordenadas matriz velocidades axiales
	MR = M1.AllocateDouble3D(NA, NRad, 2, "V"); //Coordenadas matriz velocidades radiales

	//Matrices de distancias de volúmenes de control
	DeltasMP = M1.AllocateDouble3D(NA, NRad, 2, "P"); //Deltas X y R de la matriz de Presión/Temperatura
	DeltasMU = M1.AllocateDouble3D(NA, NRad, 2, "U"); //Deltas X y R de la matriz de velocidades axiales (U)
	DeltasMR = M1.AllocateDouble3D(NA, NRad, 2, "V"); //Deltas X y R de la matriz de velocidades radiales (V)
	
	//Matrices de superficies de volúmenes de control
	SupMP = M1.AllocateDouble3D(NA, NRad, 4, "P"); //Superficies de los volúmenes de la matriz de Presión/Temperatura
	SupMU = M1.AllocateDouble3D(NA, NRad, 4, "U"); //Superficies de los volúmenes de la matriz de velocidades axiales (U)
	SupMR = M1.AllocateDouble2D(NA, 4); //Superficies de los volúmenes de la matriz de velocidades radiales (V)
	
	//Matrices de volúmenes de los volúmenes de control
	VolMP = M1.AllocateDouble2D(NA, NRad); //Volúmenes de los volúmenes de la matriz de Presión/Temperatura
	VolMU = M1.AllocateDouble2D(NA+1, NRad); //Volúmenes de los volúmenes de la matriz de velocidades axiales (U)
	VolMR = M1.AllocateDouble1D(NA); //Volúmenes de los volúmenes de la matriz de velocidades radiales (V)

	//Matrices de los ángulos entre las velocidades y las superficies de los volúmenes de control
	AngleMR = M1.AllocateDouble2D(NA, NRad+1); //Ángulos entre las superficies de los v.c y las velocidades radiales (V)
	AngleMU = M1.AllocateDouble2D(NA+1, NRad); //Ángulos entre las velocidades axiales (U)....
	
	//Distancia mínima de cada nodo a la pared
	minDist = M1.AllocateDouble2D(NA, NRad); //Distancias mínimas de cada nodo a la pared mas cercana

	//Matrices de las superficies efectivas de los volúmenes de control
	eSupMP = M1.AllocateDouble3D(NA, NRad, 4, "P"); //Superficies efectiva de los volúmenes de control de Presión

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
double nr = NRad;

	//Coordenadas Axiales

	//Coordenadas axiales matriz de velocidades axiales (U)
	for(i = 0; i < NA+1; i++){
		for(j = 0; j < NRad; j++){
			MU[i][j][0] = AxialCoord[i];	
		}
	}

	//Coordenadas axiales matriz de velocidades radial (V)
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad+1; j++){
			MR[i][j][0] = 0.5*(AxialCoord[i] + AxialCoord[i+1]);
		}
	}
	//Coordenadas axiales matriz de Presión/Temperatura
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MP[i][j][0] = 0.5*(AxialCoord[i] + AxialCoord[i+1]);
		}
	}

	//Coordenadas Radiales

	//Coordenadas radiales de la matriz de velocidades radial (V)
	if(OpcionR == 1){ //Regular
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad+1; j++){
				J = j;
				MR[i][j][1] = J*(Radius[i]/nr);
			}
		}
	}
	else if(OpcionR == 2){ //Tangencial hiperbólica
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad+1; j++){
				J = j;
				MR[i][j][1] = Radius[i]*((tanh(StretFactorR*((J)/nr)))/tanh(StretFactorR));
			}
		}
	}
	else if(OpcionR == 3){ //Senoidal
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad+1; j++){
				J = j;
				MR[i][j][1] = Radius[i]*sind(90.0*(J/nr));
			}
		}
	}
	
	//Coordenadas radiales de la matriz de Presión/Temperatura
	for(i = 0; i < NA; i++){
		MP[i][0][1] = MR[i][0][1];
		for(j = 1; j < NRad; j++){
			MP[i][j][1] = 0.5*(MR[i][j][1] + MR[i][j+1][1]);
		}
	}

	//Coordenadas radiales de la matriz de velocidades axial (U)
	for(j = 0; j < NRad; j++){ //Parte izda
		MU[0][j][1] = MP[0][j][1];
	}

	for(i = 1; i < NA; i++){ //Centro
		MU[i][0][1] = MP[i][0][1];
		for(j = 0; j < NRad; j++){
			MU[i][j][1] = 0.5*(MP[i-1][j][1] + MP[i][j][1]);
		}
	}

	MU[NA][0][1] = MP[NA-1][0][1];
	for(j = 1; j < NRad; j++){ //Parte dra
		if(OpcionR == 1){ //Regular
			J = j;
			MU[NA][j][1] = 0.50*(MP[NA-1][j][1] + 0.50*(J*((2.0*Rfinal - MR[NA-1][NRad][1])/nr) + (J+1)*((2.0*Rfinal - MR[NA-1][NRad][1])/nr)));
				
		}
		else if(OpcionR == 2){ //Tangencial hiperbólica
			J = j;
			MU[NA][j][1] = 0.50*(MP[NA-1][j][1] + 0.50*((2.0*Rfinal - MR[NA-1][NRad][1])*(tanh(StretFactorR*(J/nr)))/tanh(StretFactorR) + (2.0*Rfinal - MR[NA-1][NRad][1])*(tanh(StretFactorR*((J+1)/nr)))/tanh(StretFactorR)));
		}
		else if(OpcionR == 3){ //Senoidal
			J = j;
			MU[NA][j][1] = 0.50*(MP[NA-1][j][1] + 0.50*((2.0*Rfinal - MR[NA-1][NRad][1])*sind(90.0*(J/nr)) + (2.0*Rfinal - MR[NA-1][NRad][1])*sind(90.0*((J+1)/nr))));
		}
	}
	
	

}

//Pasar los resultados de las mallas a un txt
void Mesher::PrintTxt(){
	
int i, j;

	FILE *fp1;
	fp1 = fopen("/home_nobck/sergiogus/TFG/MeshResults/GnuPlot/MallaPresion.txt","w");
	for(i = 0; i < NA; i++){
        	for(j = 0; j < NRad; j++){
			fprintf(fp1,"%f %f \n",MP[i][j][0],MP[i][j][1]);
		}
        	fprintf(fp1, "\n");
    	}
	fclose(fp1);

}

//Cálculo de las distancias entre nodos en cada una de las matrices
void Mesher::Get_Deltas(){
int i, j;

	//Deltas X y R de la matriz de Presión/Temperatura
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			DeltasMP[i][j][0] = Hyp(MU[i][j][0], MU[i+1][j][0], MU[i][j][1], MU[i+1][j][1]); //Deltas X
			DeltasMP[i][j][1] = MR[i][j+1][1] - MR[i][j][1]; //Deltas R
		}
	}

	//Deltas X y R de la matriz de velocidades axiales (U)	

	DeltasMU[NA][0][1] = 0.50*(0.50*(MR[NA-1][0][1] + MR[NA-1][1][1]) + MU[NA][1][1]); //Delta R (Abajo dra) 
	DeltasMU[NA][NRad-1][1] = Rfinal - 0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1]); //Delta R (Arriba dra)

	for(j = 0; j < NRad; j++){
		DeltasMU[0][j][0] = Hyp(MU[0][j][0], MP[0][j][0], MU[0][j][1], MP[0][j][1]); //Deltas X (Parte izda)
		DeltasMU[NA][j][0] = Hyp(MP[NA-1][j][0], MU[NA][j][0], 0.50*(MR[NA-1][j][1] + MR[NA-1][j+1][1]), MU[NA][j][1]); //Deltas X (Parte dra)

		DeltasMU[0][j][1] = MR[0][j+1][1] - MR[0][j][1]; //Deltas R (Parte izda)	
	}
	DeltasMU[NA][0][0] = MU[NA][0][0] - MP[NA-1][0][0]; //Deltas X (Parte dra)

	DeltasMU[NA][1][1] = 0.50*(MU[NA][2][1] - 0.50*(MR[NA-1][0][1] + MR[NA-1][1][1])); //Deltas R (Parte dra)
	for(j = 2; j < NRad-1; j++){
		DeltasMU[NA][j][1] = 0.50*(MU[NA][j+1][1] - MU[NA][j-1][1]); //Deltas R (Parte dra)
	}
	for(i = 1; i < NA; i++){		
		DeltasMU[i][0][0] = Hyp(MP[i-1][0][0], MP[i][0][0], MP[i-1][0][1], MP[i][0][1]); //Deltas X (Parte abajo)
		DeltasMU[i][NRad-1][0] = Hyp(MP[i-1][NRad-1][0], MP[i][NRad-1][0], MP[i-1][NRad-1][1], MP[i][NRad-1][1]); //Deltas X (Parte arriba)

		DeltasMU[i][0][1] = 0.50*(MR[i-1][1][1] + MR[i][1][1]) - 0.50*(MR[i-1][0][1] + MR[i][0][1]); //Deltas R (Parte abajo)
		DeltasMU[i][NRad-1][1] = 0.50*(MR[i-1][NRad][1] + MR[i][NRad][1]) - 0.50*(MR[i-1][NRad-1][1] + MR[i][NRad-1][1]); //Deltas R (Parte arriba)
		for(j = 1; j < NRad-1; j++){
			DeltasMU[i][j][0] = Hyp(MP[i-1][j][0], MP[i][j][0], MP[i-1][j][1], MP[i][j][1]); //Deltas X (Centro)
			DeltasMU[i][j][1] =  0.50*(MR[i-1][j+1][1] + MR[i][j+1][1]) - 0.50*(MR[i-1][j][1] + MR[i][j][1]); //Deltas R (Centro)
		}
	}
	
	//Deltas X y R de la matriz de velocidades radiales (V)
	for(j = 1; j < NRad; j++){
		DeltasMR[0][j][0] = Hyp(MU[0][j][0], MU[1][j][0], 0.50*(MU[0][j-1][1] + MU[0][j][1]), 0.50*(MU[1][j-1][1] + MU[1][j][1])); //Deltas X (Parte izda)
		DeltasMR[NA-1][j][0] = Hyp(MU[NA-1][j][0], MU[NA][j][0], 0.50*(MU[NA-1][j-1][1] + MU[NA-1][j][1]), 0.50*(MU[NA][j-1][1] + MU[NA][j][1])); //Deltas X (Parte dra)
	}
	DeltasMR[0][NRad][0] = Hyp(MU[0][NRad-1][0], MU[1][NRad-1][0], Rinicio, 0.50*(MR[0][NRad][1] + MR[1][NRad][1]));//Deltas X (Arriba izda)
	DeltasMR[0][0][0] = MU[1][0][0] - MU[0][0][0];//Deltas X (Abajo izda)

	DeltasMR[NA-1][NRad][0] = Hyp(MU[NA-1][NRad-1][0], MU[NA][NRad-1][0], 0.50*(MR[NA-2][NRad][1] + MR[NA-1][NRad][1]), Rfinal);//Deltas X (Arriba dra)
	DeltasMR[NA-1][0][0] = MU[NA][0][0] - MU[NA-1][0][0];//Deltas X (Abajo dra)

	for(i = 1; i < NA-1; i++){ 
		DeltasMR[i][0][0] = MU[i+1][0][0] - MU[i][0][0]; //Deltas X (Parte abajo)
		DeltasMR[i][NRad][0] = Hyp(MU[i][NRad-1][0], MU[i+1][NRad-1][0], 0.50*(MR[i-1][NRad][1] + MR[i][NRad][1]), 0.50*(MR[i][NRad][1] + MR[i+1][NRad][1])); //Deltas X (Parte arriba)
		for(j = 1; j < NRad; j++){
			DeltasMR[i][j][0] = Hyp(MU[i][j][0], MU[i+1][j][0], 0.50*(MU[i][j-1][1] + MU[i][j][1]), 0.50*(MU[i+1][j-1][1] + MU[i+1][j][1])); //Deltas X (Centro)
		}
	}
	for(i = 0; i < NA; i++){
		DeltasMR[i][0][1] = 0.50*(MR[i][0][1] + MR[i][1][1]) - MR[i][0][1]; //DeltasR (Parte abajo)
		DeltasMR[i][NRad][1] = MR[i][NRad][1] - MP[i][NRad-1][1]; //DeltasR (Parte arriba)
		for(j = 1; j < NRad; j++){
			DeltasMR[i][j][1] = 0.50*(MR[i][j+1][1] + MR[i][j][1]) - 0.50*(MR[i][j][1] + MR[i][j-1][1]); //DeltasR (Centro)	
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
		for(j = 1; j < NRad-1; j++){
			SupMP[i][j][0] = PI*(pow(0.50*(MR[i][j+1][1] + MR[i-1][j+1][1]),2.0) - pow(0.50*(MR[i][j][1] + MR[i-1][j][1]),2.0)); //West
			SupMP[i][j][1] = PI*(pow(0.50*(MR[i][j+1][1] + MR[i+1][j+1][1]),2.0) - pow(0.50*(MR[i][j][1] + MR[i+1][j][1]),2.0)); //East
	
			SupMP[i][j][2] = PI*DeltasMR[i][j][0]*(0.50*(MR[i][j][1] + MR[i-1][j][1]) + 0.50*(MR[i][j][1] + MR[i+1][j][1])); //South
			SupMP[i][j][3] = PI*DeltasMR[i][j+1][0]*(0.50*(MR[i][j+1][1] + MR[i-1][j+1][1]) + 0.50*(MR[i][j+1][1] + MR[i+1][j+1][1])); //North
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SupMP[0][j][0] = PI*(pow(MR[0][j+1][1],2.0) - pow(MR[0][j][1],2.0)); //West
		SupMP[0][j][1] = PI*(pow(0.50*(MR[0][j+1][1] + MR[1][j+1][1]),2.0) - pow(0.50*(MR[0][j][1] + MR[1][j][1]),2.0)); //East
	
		SupMP[0][j][2] = PI*DeltasMR[0][j][0]*(MR[0][j][1] + 0.50*(MR[0][j][1] + MR[1][j][1])); //South
		SupMP[0][j][3] = PI*DeltasMR[0][j+1][0]*(0.50*(MU[0][j][1] + MU[0][j+1][1]) + 0.50*(MU[1][j][1] + MU[1][j+1][1])); //North

		//Parte derecha
		SupMP[NA-1][j][0] = PI*(pow(0.50*(MR[NA-1][j+1][1] + MR[NA-2][j+1][1]),2.0) - pow(0.50*(MR[NA-1][j][1] + MR[NA-2][j][1]),2.0)); //West
		SupMP[NA-1][j][1] = PI*(pow(MR[NA-1][j+1][1],2.0) - pow(MR[NA-1][j][1],2.0)); //East
	
		SupMP[NA-1][j][2] = PI*DeltasMR[NA-1][j][0]*(0.50*(MR[NA-1][j][1] + MR[NA-2][j][1]) + MR[NA-1][j][1]); //South
		SupMP[NA-1][j][3] = PI*DeltasMR[NA-1][j+1][0]*(0.50*(MU[NA-1][j][1] + MU[NA-1][j+1][1]) + 0.50*(MU[NA][j][1] + MU[NA][j+1][1])); //North
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SupMP[i][0][0] = PI*(pow(0.50*(MR[i][1][1] + MR[i-1][1][1]),2.0) - pow(0.0,2.0)); //West
		SupMP[i][0][1] = PI*(pow(0.50*(MR[i][1][1] + MR[i+1][1][1]),2.0) - pow(0.0,2.0)); //East
		
		SupMP[i][0][2] = 0.0; //South
		SupMP[i][0][3] = PI*DeltasMR[i][1][0]*(0.50*(MR[i][1][1] + MR[i-1][1][1]) + 0.50*(MR[i][1][1] + MR[i+1][1][1])); //North

		//Parte arriba
		SupMP[i][NRad-1][0] = PI*(pow(0.50*(MR[i][NRad][1] + MR[i-1][NRad][1]),2.0) - pow(0.50*(MU[i][NRad-1][1] + MU[i][NRad-2][1]),2.0)); //West
		SupMP[i][NRad-1][1] = PI*(pow(0.50*(MR[i][NRad][1] + MR[i+1][NRad][1]),2.0) - pow(0.50*(MU[i+1][NRad-1][1] + MU[i+1][NRad-2][1]),2.0)); //East
	
		SupMP[i][NRad-1][2] = PI*DeltasMR[i][NRad-1][0]*(0.50*(MU[i][NRad-1][1] + MU[i][NRad-2][1]) + 0.50*(MU[i+1][NRad-1][1] + MU[i+1][NRad-2][1])); //South
		SupMP[i][NRad-1][3] = PI*DeltasMR[i][NRad][0]*(0.50*(MR[i][NRad][1] + MR[i+1][NRad][1]) + 0.50*(MR[i][NRad][1] + MR[i-1][NRad][1])); //North
	}

	//Esquina abajo izquierda
	SupMP[0][0][0] = PI*(pow(MR[0][1][1],2.0) - pow(0.0,2.0)); //West
	SupMP[0][0][1] = PI*(pow(0.50*(MR[0][1][1] + MR[1][1][1]),2.0) - pow(0.0,2.0)); //East
	
	SupMP[0][0][2] = 0.0; //South
	SupMP[0][0][3] = PI*DeltasMR[0][1][0]*(MR[0][1][1] + 0.50*(MR[0][1][1] + MR[1][1][1]));  //North

	//Esquina arriba izquierda
	SupMP[0][NRad-1][0] = PI*(pow(Rinicio,2.0) - pow(0.50*(MU[0][NRad-1][1] + MU[0][NRad-2][1]),2.0)); //West
	SupMP[0][NRad-1][1] = PI*(pow(0.50*(MR[0][NRad][1] + MR[1][NRad][1]),2.0) - pow(0.50*(MU[1][NRad-1][1] + MU[1][NRad-2][1]),2.0)); //East
	
	SupMP[0][NRad-1][2] = PI*DeltasMR[0][NRad-1][0]*(0.50*(MU[0][NRad-1][1] + MU[0][NRad-2][1]) + 0.50*(MU[1][NRad-1][1] + MU[1][NRad-2][1])); //South
	SupMP[0][NRad-1][3] = PI*DeltasMR[0][NRad][0]*(Rinicio + 0.50*(MR[0][NRad][1] + MR[1][NRad][1])); //North

	//Esquina abajo derecha
	SupMP[NA-1][0][0] = PI*(pow(0.50*(MR[NA-1][1][1] + MR[NA-2][1][1]),2.0) - pow(0.0,2.0)); //West
	SupMP[NA-1][0][1] = PI*(pow(MR[NA-1][1][1] + (MU[NA][0][0] - MR[NA-1][1][0])*tan(0.50*(AngleMU[NA][0] + AngleMU[NA][1])),2.0) - pow(0.0,2.0)); //East
	
	SupMP[NA-1][0][2] = 0.0; //South
	SupMP[NA-1][0][3] = PI*DeltasMR[NA-1][1][0]*(0.50*(MR[NA-1][1][1] + MR[NA-2][1][1]) + MR[NA-1][1][1]); //North

	//Esquina arriba derecha
	SupMP[NA-1][NRad-1][0] = PI*(pow(0.50*(MR[NA-2][NRad][1] + MR[NA-1][NRad][1]),2.0) - pow(0.50*(MU[NA-1][NRad-1][1] + MU[NA-1][NRad-2][1]),2.0)); //West
	SupMP[NA-1][NRad-1][1] = PI*(pow(Rfinal,2.0) - pow(0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1]),2.0)); //East
	
	SupMP[NA-1][NRad-1][2] = PI*DeltasMR[NA-1][NRad-1][0]*(0.50*(MU[NA-1][NRad-1][1] + MU[NA-1][NRad-2][1]) + 0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1])); //South
	SupMP[NA-1][NRad-1][3] = PI*DeltasMR[NA-1][NRad][0]*(Rfinal + 0.50*(MR[NA-2][NRad][1] + MR[NA-1][NRad][1])); //North


	//Superficies Nodos R
	for(i = 1; i < NA-1; i++){
		SupMR[i][0] = PI*(pow(0.50*(MR[i][NRad][1] + MR[i-1][NRad][1]),2.0) - pow(MU[i][NRad-1][1],2.0)); 
		SupMR[i][1] = PI*(pow(0.50*(MR[i+1][NRad][1] + MR[i][NRad][1]),2.0) - pow(MU[i+1][NRad-1][1],2.0)); 
		SupMR[i][2] = PI*DeltasMP[i][NRad-1][0]*(MU[i][NRad-1][1] + MU[i+1][NRad-1][1]);
		SupMR[i][3] = SupMP[i][NRad-1][3];
	}

	//Esquina arriba izquierda
		SupMR[0][0] = PI*(pow(Rinicio,2.0) - pow(MU[0][NRad-1][1],2.0)); 
		SupMR[0][1] = PI*(pow(0.50*(MR[1][NRad][1] + MR[0][NRad][1]),2.0) - pow(MU[1][NRad-1][1],2.0)); 
		SupMR[0][2] = PI*DeltasMP[0][NRad-1][0]*(MU[0][NRad-1][1] + MU[1][NRad-1][1]);
		SupMR[0][3] = SupMP[0][NRad-1][3];

	//Esquina arriba derecha
		SupMR[NA-1][0] = PI*(pow(0.50*(MR[NA-1][NRad][1] + MR[NA-2][NRad][1]),2.0) - pow(MU[NA-1][NRad-1][1],2.0)); 
		SupMR[NA-1][1] = PI*(pow(0.50*(Rfinal + MR[NA-1][NRad][1]),2.0) - pow(MU[NA][NRad-1][1],2.0)); 
		SupMR[NA-1][2] = PI*DeltasMP[NA-1][NRad-1][0]*(MU[NA-1][NRad-1][1] + MU[NA][NRad-1][1]);
		SupMR[NA-1][3] = SupMP[NA-1][NRad-1][3];

}

//Cálculo de los volúmenes de control de cada volúmen
void Mesher::Get_Volumes(){
int i, j;
	
	
	//Cálculo de los volúmenes de la matriz de Presión/Temperatura

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			VolMP[i][j] = (MU[i+1][j][0] - MU[i][j][0])*(PI/3.0)*(pow(0.50*(MR[i+1][j+1][1] + MR[i][j+1][1]),2.0) + pow(0.50*(MR[i][j+1][1] + MR[i-1][j+1][1]),2.0) + 0.50*(MR[i+1][j+1][1] + MR[i][j+1][1])*0.50*(MR[i][j+1][1] + MR[i-1][j+1][1]) - pow(0.50*(MR[i+1][j][1] + MR[i][j][1]),2.0) - pow(0.50*(MR[i][j][1] + MR[i-1][j][1]),2.0) - 0.50*(MR[i+1][j][1] + MR[i][j][1])*0.50*(MR[i][j][1] + MR[i-1][j][1]));
		}
	}

	//Parte abajo
	for(i = 1; i < NA-1; i++){
		VolMP[i][0] = (MU[i+1][0][0] - MU[i][0][0])*(PI/3.0)*(pow(0.50*(MR[i+1][1][1] + MR[i][1][1]),2.0) + pow(0.50*(MR[i][1][1] + MR[i-1][1][1]),2.0) + 0.50*(MR[i+1][1][1] + MR[i][1][1])*0.50*(MR[i][1][1] + MR[i-1][1][1]));
	}

	//Parte arriba
	for(i = 1; i < NA-1; i++){
		VolMP[i][NRad-1] = (MU[i+1][NRad-1][0] - MU[i][NRad-1][0])*(PI/3.0)*(pow(0.50*(MR[i][NRad][1] + MR[i+1][NRad][1]),2.0) + pow(0.50*(MR[i][NRad][1] + MR[i-1][NRad][1]),2.0) + 0.50*(MR[i][NRad][1] + MR[i+1][NRad][1])*0.50*(MR[i][NRad][1] + MR[i-1][NRad][1]) - pow(0.50*(MU[i+1][NRad-1][1] + MU[i+1][NRad-2][1]),2.0) - pow(0.50*(MU[i][NRad-1][1] + MU[i][NRad-2][1]),2.0) - 0.50*(MU[i+1][NRad-1][1] + MU[i+1][NRad-2][1])*0.50*(MU[i][NRad-1][1] + MU[i][NRad-2][1]));
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		VolMP[0][j] = (MU[1][j][0] - MU[0][j][0])*(PI/3.0)*(pow(0.50*(MR[1][j+1][1] + MR[0][j+1][1]),2.0) + pow(MR[0][j+1][1],2.0) + 0.50*(MR[1][j+1][1] + MR[0][j+1][1])*MR[0][j+1][1] - pow(0.50*(MR[1][j][1] + MR[0][j][1]),2.0) - pow(MR[0][j][1],2.0) - 0.50*(MR[1][j][1] + MR[0][j][1])*MR[i][j][1]);

		//Parte derecha
		VolMP[NA-1][j] = (MU[NA][j][0] - MU[NA-1][j][0])*(PI/3.0)*(pow(MR[NA-1][j+1][1],2.0) + pow(0.50*(MR[NA-1][j+1][1] + MR[NA-2][j+1][1]),2.0) + MR[NA-1][j+1][1]*0.50*(MR[NA-1][j+1][1] + MR[NA-2][j+1][1]) - pow(MR[NA-1][j][1],2.0) - pow(0.50*(MR[NA-1][j][1] + MR[NA-2][j][1]),2.0) - MR[NA-1][j][1]*0.50*(MR[NA-1][j][1] + MR[NA-2][j][1]));
	}

	//Esquina arriba izquierda
	VolMP[0][NRad-1] = (MU[1][NRad-1][0] - MU[0][NRad-1][0])*(PI/3.0)*(pow(0.50*(MR[0][NRad][1] + MR[1][NRad][1]),2.0) + pow(Rinicio,2.0) + 0.50*(MR[0][NRad][1] + MR[1][NRad][1])*Rinicio - pow(0.50*(MU[1][NRad-1][1] + MU[1][NRad-2][1]),2.0) - pow(0.50*(MU[0][NRad-1][1] + MU[0][NRad-2][1]),2.0) - 0.50*(MU[1][NRad-1][1] + MU[1][NRad-2][1])*0.50*(MU[0][NRad-1][1] + MU[0][NRad-2][1]));

	//Esquina arriba derecha
	VolMP[NA-1][NRad-1] = (MU[NA][NRad-1][0] - MU[NA-1][NRad-1][0])*(PI/3.0)*(pow(Rfinal,2.0) + pow(0.50*(MR[NA-1][NRad][1] + MR[NA-2][NRad][1]),2.0) + Rfinal*0.50*(MR[NA-1][NRad][1] + MR[NA-2][NRad][1]) - pow(0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1]),2.0) - pow(0.50*(MU[NA-1][NRad-1][1] + MU[NA-1][NRad-2][1]),2.0) - 0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1])*0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1]));

	//Esquina abajo izquierda
	VolMP[0][0] = (MU[1][0][0] - MU[0][0][0])*(PI/3.0)*(pow(0.50*(MR[1][1][1] + MR[0][1][1]),2.0) + pow(MR[0][1][1],2.0) + 0.50*(MR[1][1][1] + MR[0][1][1])*MR[0][1][1] - pow(0.50*(MR[1][0][1] + MR[0][0][1]),2.0) - pow(MR[0][0][1],2.0) - 0.50*(MR[1][0][1] + MR[0][0][1])*MR[0][0][1]);

	//Esquina abajo derecha
	VolMP[NA-1][0] = (MU[NA][0][0] - MU[NA-1][0][0])*(PI/3.0)*(pow(MR[NA-1][1][1],2.0) + pow(0.50*(MR[NA-1][1][1] + MR[NA-2][1][1]),2.0) + MR[NA-1][1][1]*0.50*(MR[NA-1][1][1] + MR[NA-2][1][1]) - pow(MR[NA-1][0][1],2.0) - pow(0.50*(MR[NA-1][0][1] + MR[NA-2][0][1]),2.0) - MR[NA-1][0][1]*0.50*(MR[NA-1][0][1] + MR[NA-2][0][1]));

	
	//Volúmenes de control Nodos R
	for(i = 1; i < NA-1; i++){
		VolMR[i] = (MU[i+1][NRad-1][0] - MU[i][NRad-1][0])*(PI/3.0)*( + (pow(0.50*(MR[i][NRad][1] + MR[i+1][NRad][1]),2.0) + pow(0.50*(MR[i][NRad][1] + MR[i-1][NRad][1]),2.0) + 0.50*(MR[i][NRad][1] + MR[i+1][NRad][1])*0.50*(MR[i][NRad][1] + MR[i-1][NRad][1]))
- (pow(0.50*(MR[i][NRad-1][1] + MR[i+1][NRad-1][1]),2.0) + pow(0.50*(MR[i][NRad-1][1] + MR[i-1][NRad-1][1]),2.0) + 0.50*(MR[i][NRad-1][1] + MR[i+1][NRad-1][1])*0.50*(MR[i][NRad-1][1] + MR[i-1][NRad-1][1])));
	}

	//Esquina arriba izquierda
	VolMR[0] = (MU[1][NRad-1][0] - MU[0][NRad-1][0])*(PI/3.0)*( + (pow(0.50*(MR[0][NRad][1] + MR[1][NRad][1]),2.0) + pow(0.50*(MR[0][NRad][1] + Rinicio),2.0) + 0.50*(MR[0][NRad][1] + MR[1][NRad][1])*0.50*(MR[0][NRad][1] + Rinicio))
- (pow(0.50*(MR[0][NRad-1][1] + MR[1][NRad-1][1]),2.0) + pow(0.50*(MR[0][NRad-1][1] + MR[0][NRad-1][1]),2.0) + 0.50*(MR[0][NRad-1][1] + MR[1][NRad-1][1])*0.50*(MR[0][NRad-1][1] + MR[0][NRad-1][1])));

	//Esquina arriba derecha
	VolMR[NA-1] = (MU[NA][NRad-1][0] - MU[NA-1][NRad-1][0])*(PI/3.0)*( + (pow(Rfinal,2.0) + pow(0.50*(MR[NA-1][NRad][1] + MR[NA-2][NRad][1]),2.0) + Rfinal*0.50*(MR[NA-1][NRad][1] + MR[NA-2][NRad][1])) - (pow(0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1]),2.0) + pow(0.50*(MU[NA-1][NRad-1][1] + MU[NA-1][NRad-2][1]),2.0) + 0.50*(MU[NA][NRad-1][1] + MU[NA][NRad-2][1])*0.50*(MU[NA-1][NRad-1][1] + MU[NA-1][NRad-2][1])));

}

//Cálculo de los ángulos entre las velocidades y las superficies de los volúmenes de control
void Mesher::Get_Angles(){
int i, j;

	//Ángulos entre las superficies y las velocidades radiales (V)
	
	for(i = 0; i < NA; i++){
		AngleMR[i][0] = 0.0; //Parte abajo
		AngleMR[i][NRad] = acos((MU[i+1][NRad-1][0] - MU[i][NRad-1][0])/DeltasMR[i][NRad][0]); //Parte arriba
		for(j = 1; j < NRad; j++){
			AngleMR[i][j] = acos((MU[i+1][j][0] - MU[i][j][0])/DeltasMR[i][j][0]);  //Centro
		}
	}

	
	//Ángulos entre las velocidades axiales (U)....
	
	for(j = 0; j < NRad; j++){
		AngleMU[0][j] = acos((MP[0][j][0] - 0.0)/DeltasMU[0][j][0]); //Parte izda
		AngleMU[NA][j] = 0.0;
		//AngleMU[NA][j] = acos((MU[NA][j][0] - MP[NA-1][j][0])/DeltasMU[NA][j][0]); //Parte izda
		for(i = 1; i < NA; i++){
			AngleMU[i][j] = acos((MP[i][j][0] - MP[i-1][j][0])/DeltasMU[i][j][0]); //Centro
		}
	}

}

//Cálculo de las superficies efectivas de los volúmenes de control
void Mesher::Get_EffectiveSurfaces(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			eSupMP[i][j][0] = SupMP[i][j][0]*cos(AngleMU[i][j]); //West
			eSupMP[i][j][1] = SupMP[i][j][1]*cos(AngleMU[i+1][j]); //East
			eSupMP[i][j][2] = SupMP[i][j][2]*cos(AngleMR[i][j]); //South
			eSupMP[i][j][3] = SupMP[i][j][3]*cos(AngleMR[i][j+1]); //North
		}
	}

}

//Cálculo de la longitud característica de la geometría
void Mesher::Get_CharLength(){
int i, j;
double Volume = 0.0;

	if(Problema == "Tobera"){
		for(i = 0; i < NA; i++){ 
			for(j = 0; j < NRad; j++){
				Volume += VolMP[i][j];
			}
		}

		Lcar = Volume/ThroatArea; //Segun Sutton (Rocket Propulsion elements)
	}
	else{ Lcar = PipeDiameter; }
		
}

//Cálculo de la distancia mínima de cad anodo a la pared más cercana
void Mesher::ClosestDist(){
int i, j, k, m;
double Minimo = 2.0*PipeLength;
	
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Minimo = 2.0*PipeLength;
			//for(k = 0; k < NRad; k++){
				//if(Hyp(MP[i][j][0], MU[0][k][0], MP[i][j][1], MU[0][k][1]) <= Minimo){
				//	minDist[i][j] = Minimo;
				//}
			//}
			for(m = 0; m < NA; m++){
				if(Hyp(MP[i][j][0], MR[m][NRad][0], MP[i][j][1], MR[m][NRad][1]) < Minimo){
					Minimo = Hyp(MP[i][j][0], MR[m][NRad][0], MP[i][j][1], MR[m][NRad][1]);	
				}
			}
			minDist[i][j] = Minimo;
		}
	}



}

//Pasar todos los resultados escalares en 2D a un archivo VTK
void Mesher::EscalarVTK2D(string name, string Documento, double **Matrix){
int i, j;

	ofstream file;
        stringstream name_t;
        string name_f;

	name_t<<"/home_nobck/sergiogus/TFG/MeshResults/Pipe/Paraview/"<<Problema<<"_"<<Documento<<".vtk";

	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<name<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET STRUCTURED_GRID"<<endl;
        file<<"DIMENSIONS"<<"   "<<NA<<"   "<<NRad<<"   "<<1<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<NA*NRad<<"   "<<"double"<<endl;
	
	for(j = 0; j < NRad; j++){
		for(i = 0; i < NA; i++){
			file<<MP[i][j][0]<<"   "<<MP[i][j][1]<<"   "<<0.0<<endl;
			
		}
		
	}
        
        file<<endl;
	file<<"POINT_DATA"<<"   "<<NA*NRad<<endl;
        file<<"SCALARS "<<name<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<name<<endl;
        file<<endl;
       for(j = 0; j < NRad; j++){
		for(i = 0; i < NA; i++){
				file<<Matrix[i][j]<<" ";
			
		}
	}
        file.close();


}

//Pasar todos los resultados escalares en 2D a un archivo VTK
void Mesher::EscalarVTK2D_Especial(string name, string Documento, double **Matrix){
int i, j;

	ofstream file;
        stringstream name_t;
        string name_f;

	name_t<<"/home_nobck/sergiogus/TFG/MeshResults/Pipe/Paraview/"<<Problema<<"_"<<Documento<<".vtk";

	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<name<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET STRUCTURED_GRID"<<endl;
        file<<"DIMENSIONS"<<"   "<<NA<<"   "<<NRad+1<<"   "<<1<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<(NA)*(NRad+1)<<"   "<<"double"<<endl;
	
	for(j = 0; j < NRad+1; j++){
		for(i = 0; i < NA; i++){
			file<<MR[i][j][0]<<"   "<<MR[i][j][1]<<"   "<<0.0<<endl;
			
		}
		
	}
        
        file<<endl;
	file<<"POINT_DATA"<<"   "<<(NA)*(NRad+1)<<endl;
        file<<"SCALARS "<<name<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<name<<endl;
        file<<endl;
       for(j = 0; j < NRad+1; j++){
		for(i = 0; i < NA; i++){
				file<<Matrix[i][j]<<" ";
			
		}
	}
        file.close();


}


//Pasar todos los resultados escalares en 2D a un archivo VTK
void Mesher::EscalarVTK2Dv2(string name, string Documento){
int i, j;

	ofstream file;
        stringstream name_t;
        string name_f;

	name_t<<"/home/sergio/Desktop/TFG/MeshResults/Paraview/"<<Problema<<"_"<<Documento<<".vtk";

	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<name<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET STRUCTURED_GRID"<<endl;
        file<<"DIMENSIONS"<<"   "<<NA<<"   "<<NRad<<"   "<<1<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<NA*NRad<<"   "<<"double"<<endl;
	
	for(j = 0; j < NRad; j++){
		for(i = 0; i < NA; i++){
			file<<MP[i][j][0]<<"   "<<MP[i][j][1]<<"   "<<0.0<<endl;
			
		}
		
	}
        
        file<<endl;
	file<<"POINT_DATA"<<"   "<<NA*NRad<<endl;
        file<<"SCALARS "<<name<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<name<<endl;
        file<<endl;
       for(j = 0; j < NRad; j++){
		for(i = 0; i < NA; i++){
				file<<0.0<<" ";
			
		}
	}
        file.close();


}

//Pasar todos los resultados escalares en 3D a un archivo VTK
void Mesher::EscalarVTK3D(Memory M1, string name, string Documento){
int i, j, k;
double J, K;

	ofstream file;
        stringstream name_t;
        string name_f;

	name_t<<"/home_nobck/sergiogus/TFG/MeshResults/Paraview/"<<Documento<<".vtk";

	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<name<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET STRUCTURED_GRID"<<endl;
        file<<"DIMENSIONS"<<"   "<<NAng+1<<"   "<<NRad+1<<"   "<<NA<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<NA*(NRad+1)*(NAng+1)<<"   "<<"double"<<endl;
	
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad+1; j++){
			for(k = 0; k < NAng+1; k++){
			K = k;
			file<<MR[i][j][0]<<"   "<<MR[i][j][1]*cosd((K*360.0)/NAng)<<"   "<<MR[i][j][1]*sind((K*360.0)/NAng)<<endl;
			}
		}
		
	}
        
        file<<endl;
	file<<"POINT_DATA"<<"   "<<(NRad+1)*(NAng+1)*NA<<endl;
        file<<"SCALARS "<<name<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<name<<endl;
        file<<endl;
       for(i = 0; i < NA; i++){
		for(j = 0; j < NRad+1; j++){
			for(k = 0; k < NAng+1; k++){
				file<<30.0*i+90.0*j + 60.0*k<<" ";
			}
		}
	}
        file.close();


}

//Printar por pantalla
void Mesher::PrintWindow(){
int i, j;
double Voltotal = 0.0;	
	
	for(j = NRad; j >= 0; j--){
		for(i = 0; i < NA; i++){
			//Voltotal += SupMP[i][j][1];
			cout<<DeltasMR[i][j][1]<<", ";
		}
		cout<<endl;
	}
	//MU[NA][1][0] - MU[NA-1][1][0]

	cout<<"Volumen: "<<Voltotal<<endl;
}

//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1, Geometry G1){
int i;
	
	if(Problema == "Tobera"){	
		GetTotalAxialNodes(); //Cálculo del número de nodos totales en función de las opciones de discretización	
	}
	
	//Matrices de coordenadas de discretización
	AxialCoord = M1.AllocateDouble1D(NA + 1); //Matriz coordenada axial de la distribución
	Radius = M1.AllocateDouble1D(NA + 1); //Matriz radios en cada coordenada axial

	Get_AxialCoordinates(); //Cálculo del la posición axial de los nodos de velocidad axial

	Get_Radius(M1, G1); //Cálculo del radio de cada coordenada axial discretizada

	GetTotalRadialNodes(); //Calcular el número total de nodos en la dirección radial
	cout<<"NA: "<<NA<<endl;
	cout<<"NRad: "<<NRad<<endl;

	AllocateMemory(M1); //Alojamiento de memoria para cada matriz 
	
	
	Get_Mesh(); //Creación de todas las mallas
	
	Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
	
	Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control
	
	Get_Angles(); //Cálculo de los ángulos entre las velocidades y las superficies de los volúmenes de control

	Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen
	
	//Get_EffectiveSurfaces(); //Cálculo de las superficies efectivas de los volúmenes de control

	//PrintWindow(); //Printar por pantalla
	
	ClosestDist(); //Cálculo de la distancia mínima de cad anodo a la pared más cercana
	
	PrintTxt(); //Pasar los resultados de las mallas a un txt

	Get_CharLength(); //Cálculo de la longitud característica de la geometría

	

	EscalarVTK2Dv2("Prueba geometry", "Prueba");

}

