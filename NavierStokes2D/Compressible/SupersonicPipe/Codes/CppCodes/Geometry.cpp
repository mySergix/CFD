#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <cassert>
using namespace std;

#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Memory.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/ReadData.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Geometry.h"

#define PI 3.141592653589793

#define sind(x) sin(x * (PI/180.0))
#define cosd(x) cos(x * (PI/180.0))
#define tand(x) tan(x * (PI/180.0))

//Constructor del mallador
Geometry::Geometry(Memory M1, ReadData R1){
	
	//Datos método iterativo
	Convergencia = 1e-10; //Convergencia para el método iterativo de los coeficientes (Tobera NASA)

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
	
	//Datos geometría sección divergente
	CoeficientesExpresion = M1.AllocateDouble1D(4); //Coeficientes de la expresión matemática
	CoeficientesExpresionSup = M1.AllocateDouble1D(4); //Coeficientes de la expresión matemática

	//Cálculo de datos necesarios
	m = -tand(Theta1); //Pendiente recta Sección 3
	rW = (Dt/2.0) +  Rad3*(1 - cosd(Theta3)); //Radio del punto de tangencia sección divergente
	Xt = L1 + Rad1*sind(Theta1) + Rad2*sind(Theta2) + ((Dt/2 + Rad2*(1 - cosd(Theta2))) - (D1/2 - Rad1*(1 - cosd(Theta1))))/m;
	//Coordenada axial de la garganta de la tobera
	xW = Rad3*sind(Theta3); //Coordenada axial del punto de tangencia sección divergente

	
	
	
}

void Geometry::GeometryNasa(Memory M1){
int i;
double MaxDiferencia = 2.0*Convergencia;
	
	//Inicialización de coeficientes
	CoeficientesExpresion[0] = 0.0; //Coeficiente a
	CoeficientesExpresion[1] = 0.0; //Coeficiente b
	CoeficientesExpresion[2] = 0.0; //Coeficiente c
	CoeficientesExpresion[3] = 0.0; //Coeficiente d

	//Inicialización de coeficientes supuestos
	CoeficientesExpresionSup[0] = CoeficientesExpresion[0]; //Coeficiente a
	CoeficientesExpresionSup[1] = CoeficientesExpresion[1]; //Coeficiente b
	CoeficientesExpresionSup[2] = CoeficientesExpresion[2]; //Coeficiente c
	CoeficientesExpresionSup[3] = CoeficientesExpresion[3]; //Coeficiente d

	while(MaxDiferencia >= Convergencia){
		//Cálculo a
		CoeficientesExpresion[0] = rW - sqrt(CoeficientesExpresion[1]); 

		//Cálculo b
		CoeficientesExpresion[1] = pow(rE - CoeficientesExpresion[0] - CoeficientesExpresion[3]*(xE-xW),2.0) - 	CoeficientesExpresion[2] *(xE-xW); 

		//Cálculo c
		CoeficientesExpresion[2] = 2.0*sqrt(CoeficientesExpresion[1])*(tand(Theta3) - CoeficientesExpresion[3]);

		//Cálculo d
		CoeficientesExpresion[3] = tand(ThetaE) - CoeficientesExpresion[2]/(2.0*sqrt(CoeficientesExpresion[1] + CoeficientesExpresion[2]*(xE-xW)));

		MaxDiferencia = 0.0;
		for(i = 0; i < 4; i++){
			if(abs(CoeficientesExpresion[i] - CoeficientesExpresionSup[i]) >= MaxDiferencia){
				MaxDiferencia = abs(CoeficientesExpresion[i] - CoeficientesExpresionSup[i]);
			} 
		} 
		CoeficientesExpresionSup[0] = CoeficientesExpresion[0]; //Coeficiente a
		CoeficientesExpresionSup[1] = CoeficientesExpresion[1]; //Coeficiente b
		CoeficientesExpresionSup[2] = CoeficientesExpresion[2]; //Coeficiente c
		CoeficientesExpresionSup[3] = CoeficientesExpresion[3]; //Coeficiente d
	}

}


void Geometry::GeometryRAO(Memory M1){

	//Resolución sistema de ecuaciones parábola
	CoeficientesExpresion[0] = (1/tand(Theta3) - 1/tand(ThetaE))/(2*rW - 2*rE); //Coeficiente a

	CoeficientesExpresion[1] = 1/tand(Theta3) - 2*rE*CoeficientesExpresion[0]; //Coeficiente b

	CoeficientesExpresion[2] = xW - pow(rW,2)*CoeficientesExpresion[0] - rW*CoeficientesExpresion[1]; //Coeficiente c

	CoeficientesExpresion[3] = 0.0; //Coeficiente d
	
}

//Método de cálculo de los límites de la Tobera geometría NASA
void Geometry::NasaLimits(Memory M1){
int i, j, k;
double I, J, K;

	//Datos método límites
	LimInfXE = 1.0; //Límite inferior de la coordenada axial de salida (m)
	LimSupXE = 4.0; //Límite superior de la coordenada axial de salida (m)
	AxialCoordIterations = 100; //Número total de coordenadas axiales de salida a estudiar dentro del rango

	LimInfRE = 0.0; //Límite inferior de la coordenada radial de salida (m)
	LimSupRE = 6.5; //Límite superior de la coordenada radial de salida (m)
	ExitRadiusIterations = 200; //Número total de coordenadas radiales de salida a estudiar dentro del rango

	LimInfThetaE = 0.0; //Límite superior del ángulo de salida (º)
	LimSupThetaE = 45.0; //Límite inferior del ángulo de salida (º)
	AngleIterations = 10; //Número total de ángulos a estudiar dentro del rango

	XEpositions = M1.AllocateDouble1D(AxialCoordIterations+1);
	UpperRE = M1.AllocateDouble1D(AxialCoordIterations+1);
	LowerRE = M1.AllocateDouble1D(AxialCoordIterations+1);

double A, B, C, D;
double Asup, Bsup, Csup, Dsup;

double MaxDiferencia;

double UpperLimitXE = LimInfXE; 
double LowerLimitXE = LimSupXE;
double UpperLimitRE = LimInfRE;
double LowerLimitRE = LimSupRE;

double XE, RE, TE;
double AxialIter = AxialCoordIterations;
double RadialIter = ExitRadiusIterations;
double AngleIter = AngleIterations;
int OK;
int Comprobacion;
int Numero;
	TE = ThetaE;
	FILE *fp1;
	fp1 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/NumericalResults/Geometry/LimiteCoordenadas.txt","w");
		fprintf(fp1,"Coordenada Axial \t Límite Inf RE \t Límite Sup RE \n");

		for(j = 0; j < AxialCoordIterations+1; j++){
			J = j;
			XE = LimInfXE + J*(LimSupXE - LimInfXE)/AxialIter;
			XEpositions[j] = XE;
			UpperLimitRE = LimInfRE;
			LowerLimitRE = LimSupRE;
			cout<<(J/AxialIter)*100.0<<" %"<<endl;
			
			for(k = 0; k < ExitRadiusIterations+1; k++){
				K = k;
				RE = LimInfRE + K*(LimSupRE - LimInfRE)/RadialIter;
				
				//Inicialización de coeficientes
				A = 0.0; //Coeficiente a
				B = 0.0; //Coeficiente b
				C = 0.0; //Coeficiente c
				D = 0.0; //Coeficiente d

				//Inicialización de coeficientes supuestos
				Asup = A; //Coeficiente a
				Bsup = B; //Coeficiente b
				Csup = C; //Coeficiente c
				Dsup = D; //Coeficiente d

				MaxDiferencia = 2.0*Convergencia;
				while(MaxDiferencia >= Convergencia){
					//Cálculo a
					A = rW - sqrt(B); 

					//Cálculo b
					B = pow(RE - A - D*(XE-xW),2.0) - 	C*(xE-xW); 

					//Cálculo c
					C = 2.0*sqrt(B)*(tand(Theta3) - D);

					//Cálculo d
					D = tand(ThetaE) - C/(2.0*sqrt(B + C*(XE-xW)));

					MaxDiferencia = 0.0;
					
					if(abs(A - Asup) >= MaxDiferencia){ MaxDiferencia = abs(A - Asup); }
					if(abs(B - Bsup) >= MaxDiferencia){ MaxDiferencia = abs(B - Bsup); }
					if(abs(C - Csup) >= MaxDiferencia){ MaxDiferencia = abs(C - Csup); }
					if(abs(D - Dsup) >= MaxDiferencia){ MaxDiferencia = abs(D - Dsup); }
							
					//Parada de emergencia
					if(MaxDiferencia >= 10.0 || MaxDiferencia < 0.01*Convergencia){ 
						break;
					}

					Asup = A; //Coeficiente a
					Bsup = B; //Coeficiente b
					Csup = C; //Coeficiente c
					Dsup = D; //Coeficiente d
				}

				if(MaxDiferencia <= Convergencia && MaxDiferencia >= 0.1*Convergencia){
					cout<<"MaxDif: "<<MaxDiferencia<<", RE: "<<RE<<", UpperLimit: "<<UpperLimitRE<<endl;
					if(RE >= UpperLimitRE){ UpperLimitRE = RE; }
					if(RE <= LowerLimitRE){ LowerLimitRE = RE; }
				}			
			}
				if(UpperLimitRE == LimInfRE && LowerLimitRE < LimSupRE){ //Mal UpperRE
					fprintf(fp1,"%f \t \t NO \t \t %f \n", XEpositions[j], LowerLimitRE);
				}
				else if(LowerLimitRE == LimSupRE && UpperLimitRE > LimInfRE){ //Mal LowerRE
					fprintf(fp1,"%f \t \t %f \t \t NO \n", XEpositions[j], UpperLimitRE);
				}
				else if(UpperLimitRE == LimInfRE && LowerLimitRE == LimSupRE){ //Mal ambos
					fprintf(fp1,"%f \t \t NO \t \t NO \n", XEpositions[j]);	
				}
				else{ //Bien ambos
					fprintf(fp1,"%f \t \t %f \t \t %f \n", XEpositions[j], LowerLimitRE, UpperLimitRE);
				}
			
			
		}	
			
		fclose(fp1);
		
}




