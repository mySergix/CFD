#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/ReadTXT.h"
#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/Mesh.h"

using namespace std;

//Constructor del mallador
Mallador::Mallador(ReadTXT R1){

	//Datos geométricos del problema
	Xtotal = R1.Get_GeometryData(0);
	Ytotal = R1.Get_GeometryData(1);
	Wtotal = R1.Get_GeometryData(2);
	Xcent = R1.Get_GeometryData(3);
	Ycent = R1.Get_GeometryData(4);
	Xsol = R1.Get_GeometryData(5);
	Ysol = R1.Get_GeometryData(6);

	//Datos numéricos del problema
	Problema = R1.Get_NumericalData(0);
	Progresion = R1.Get_NumericalData(1);
	NX = R1.Get_NumericalData(2);
	NY = R1.Get_NumericalData(3);
	OpcionX = R1.Get_NumericalData(4);
	OpcionY = R1.Get_NumericalData(5);
	
	StretFactor = R1.Get_PrecisionData(0);
	
	//Definición geometría de los problemas
	if(Problema == 1 || Problema == 2){
		_X1 = 0;
		_Y1 = 0;
		_X2 = Xtotal;
		_Y2 = Ytotal;
	}
	else{
		_X1 = Xcent - Xsol/2;
		_X2 = Xsol;
		_X3 = Xtotal - Xcent - Xsol/2;

		_Y1 = Ycent - Ysol/2;
		_Y2 = Ysol;
		_Y3 = Ytotal - Ycent - Ysol/2;
	}
	
	if(Problema == 1 || Problema == 2){
		NXTOT = NX;
		NYTOT = NY;
	}
	else{

		if(Progresion == 1){
			NY1 = NY/2;
			NY3 = NY/2;

			NX1 = 0.25*NX;
			NX3 = NX-NX1;
	
		 double nx1 = NX1;
		 double	nx3 = NX3;

		 double	ny1 = NY1;
		 double	ny3 = NY3;

			//Calculo de incrementos en X
			if(OpcionX == 1){//Regular
				Delta3 = _X1/nx1;
				Delta4 = _X3/nx3;
			}
			else if(OpcionX == 2){//Tangente Hiperbólica
				Delta3 = _X1-(_X1/2)*(1+(tanh(StretFactor*((nx1-2)/nx1)) + tanh(StretFactor))/tanh(StretFactor) - 1);
				Delta4 = (_X3/2)*(1+(tanh(StretFactor*((2-nx3)/nx3)) + tanh(StretFactor))/tanh(StretFactor) - 1);
			}
			else if(OpcionX == 3){//Cosenoidal
				Delta3 = _X1 - _X1*(1 - cos(((nx1-1)/nx1)*(PI/2)));
				Delta4 = _X3*(1 - cos((1/nx3)*(PI/2)));
				}
			else if(OpcionX == 4){//Senoidal
				Delta3 = _X1 - _X1*(sin(((nx1-1)/nx1)*(PI/2)));
				Delta4 = _X3*(sin((1/nx3)*(PI/2)));
			}
			else if(OpcionX == 5){//Cúbica
				Delta3 = _X1 - (_X1/2)*(1 + pow(2*((nx1)/(nx1+1))-1,3));
				Delta4 = (_X3/2)*(1 + pow(2*((1+1)/(nx3+1))-1,3));
			}

			//Calulo de incrementos en Y
			if(OpcionY == 1){//Regular
				Delta1 = _Y1/ny1;
				Delta2 = _Y3/ny3;
			}
			else if(OpcionY == 2){//Tangente Hiperbólica
				Delta1 = _Y1-(_Y1/2)*(1+(tanh(StretFactor*((ny1-2)/ny1)) + tanh(StretFactor))/tanh(StretFactor) - 1);
				Delta2 = (_Y3/2)*(1+(tanh(StretFactor*((2-ny3)/ny3)) + tanh(StretFactor))/tanh(StretFactor) - 1);
			}
			else if(OpcionY == 3){//Cosenoidal
				Delta1 = _Y1 - _Y1*(1 - cos(((ny1-1)/ny1)*(PI/2)));
				Delta2 = _Y3*(1 - cos((1/ny3)*(PI/2)));
			}
			else if(OpcionY == 4){//Senoidal
				Delta1 = _Y1 - _Y1*(sin(((ny1-1)/ny1)*(PI/2)));
				Delta2 = _Y3*(sin((1/ny3)*(PI/2)));
			}
			else if(OpcionY == 5){//Cúbica
				Delta1 = _Y1 - (_Y1/2)*(1 + pow(2*((ny1)/(nx1+1))-1,3));
				Delta2 = (_Y3/2)*(1 + pow(2*((1+1)/(nx3+1))-1,3));
			}

			//Calculo del número de nodos necesario
			NX2 = (2*_X2 - Delta3 + Delta4)/(Delta3 + Delta4);
			NY2 = (2*_Y2 - Delta1 + Delta2)/(Delta1 + Delta2);

		}
		else{ //No progresión lineal en las transiciones
			
			NY2 = 50;
			NY1 = (NY - NY2)/2;
			NY3 = NY - NY2 - NY1;

			NX2 = 50;
			NX1 = 0.25*(NX-NX2);
			NX3 = (NX - NX2 - NX1);
			
		}
	
		//Calculo del número de nodos necesario
		NXTOT = NX1 + NX2 + NX3;
		NYTOT = NY1 + NY2 + NY3;

	}
	
	//Alojamiento de memoria para las matrices
	M1 = Allocate3D(2, 1);
	M2 = Allocate3D(2, 2);
	M3 = Allocate3D(2, 3);

	DeltasYP = Allocate2D(NXTOT, NYTOT);
	DeltasXP = Allocate2D(NXTOT, NYTOT);

	DeltasYU = Allocate2D(NXTOT+1, NYTOT);
	DeltasXU = Allocate2D(NXTOT+1, NYTOT);

	DeltasYV = Allocate2D(NXTOT, NYTOT+1);
	DeltasXV = Allocate2D(NXTOT, NYTOT+1);

	SurfaceYU = Allocate2D(NXTOT+1, NYTOT);
	SurfaceXU = Allocate2D(NXTOT+1, NYTOT);

	SurfaceYV = Allocate2D(NXTOT, NYTOT+1);
	SurfaceXV = Allocate2D(NXTOT, NYTOT+1);

	SurfaceYP = Allocate2D(NXTOT, NYTOT);
	SurfaceXP = Allocate2D(NXTOT, NYTOT);

}

int Mallador::Get_NodosX(){ return NXTOT;}	
int Mallador::Get_NodosY(){ return NYTOT;}
int Mallador::Get_NX1(){ return NX1;}	
int Mallador::Get_NX2(){ return NX2;}
int Mallador::Get_NX3(){ return NX3;}
int Mallador::Get_NY1(){ return NY1;}	
int Mallador::Get_NY2(){ return NY2;}
int Mallador::Get_NY3(){ return NY3;}
double Mallador::Get_Xtotal(){ return Xtotal;}	
double Mallador::Get_Ytotal(){ return Ytotal;}
double Mallador::M1_1(int i, int j, int dim){ return M1[i][j][dim];}
double Mallador::M1_2(int i, int j, int dim){ return M2[i][j][dim];}
double Mallador::M1_3(int i, int j, int dim){ return M3[i][j][dim];}
double Mallador::DeltaXP(int i, int j){ return DeltasXP[i][j];}
double Mallador::DeltaYP(int i, int j){ return DeltasYP[i][j];}
double Mallador::DeltaXU(int i, int j){ return DeltasXU[i][j];}
double Mallador::DeltaYU(int i, int j){ return DeltasYU[i][j];}
double Mallador::DeltaXV(int i, int j){ return DeltasXV[i][j];}
double Mallador::DeltaYV(int i, int j){ return DeltasYV[i][j];}
double Mallador::SupXU(int i, int j){ return SurfaceXU[i][j];}
double Mallador::SupYU(int i, int j){ return SurfaceYU[i][j];}
double Mallador::SupXV(int i, int j){ return SurfaceXV[i][j];}
double Mallador::SupYV(int i, int j){ return SurfaceYV[i][j];}
double Mallador::SupXP(int i, int j){ return SurfaceXP[i][j];}
double Mallador::SupYP(int i, int j){ return SurfaceYP[i][j];}
	
//Memoria dinámica matriz 2D
double **Mallador::Allocate2D(int _NXTOT, int _NYTOT){
	int i, j;
	double **M1;

	M1 = new double *[_NXTOT];
	
		for(i=0;i < _NXTOT;i++){
			M1[i] = new double [_NYTOT];
		
		}
	return M1;
}


//Memoria dinámica matriz 3D
double ***Mallador::Allocate3D(int Dim, int NofMesh){
int i, j;	
double ***M1;

		if(NofMesh == 1 || NofMesh == 2){
			M1 = new double **[NXTOT];
		}
		else if(NofMesh == 3){
			M1 = new double **[NXTOT+1];
		}
	
		for (i = 0; i < NXTOT+1; i++){
			if(i < NXTOT){
				if(NofMesh == 1 || NofMesh == 3){
					M1[i] = new double *[NYTOT];
				}
				else{
					M1[i] = new double *[NYTOT +1];
				}	
			}
			else{
				if(NofMesh == 3){
				M1[i] = new double *[NYTOT];
				}
			
			}
			
		}
		for (i = 0; i < NXTOT+1; i++){
			for (j = 0; j < NYTOT+1; j++){
				if(i < NXTOT && j < NYTOT){
					M1[i][j] = new double [Dim];
				
				}
				else if(j == NYTOT && i < NXTOT){
					if(NofMesh == 2){
						M1[i][j] = new double [Dim];
					}
				}
				else if(j < NYTOT && i == NXTOT){
					if(NofMesh == 3){
						M1[i][j] = new double [Dim];
					}
				}		
			}
		}
	return M1;
}


//Creación de las mallas
void Mallador::Get_Mesh(){
int i, j;
double I, J;

	if(Problema == 1 || Problema == 2){
	
		double nxtot = NXTOT;
		double nytot = NYTOT;

		double DeltaX;
		double DeltaY;

		//Discretización vertical de la malla 2
			for (i = 0; i < NXTOT; i++){
				for (j = 0; j < NYTOT+1; j++){
					
					if(j==0){
						M2[i][j][1] = 0;
					}
					else if (j==NYTOT){
						M2[i][j][1] = Ytotal;
					}
					else{
						if(OpcionY == 1){//Regular
							DeltaY = Ytotal/nytot;
							M2[i][j][1] = DeltaY*j;
						}
						else if(OpcionY == 2){//Tangente Hiperbólica
							J = j;
							M2[i][j][1] = (Ytotal/2)*(1+(tanh(StretFactor*((2*J-nytot)/nytot)) + tanh(StretFactor))/tanh(StretFactor) - 1);
						}
						else if(OpcionY == 3){//Cosenoidal
							J = j;
							M2[i][j][1] = Ytotal*(1 - cos(((J/nytot)*(PI/2))));
						}
						else if(OpcionY == 4){//Senoidal
							J = j;
							M2[i][j][1] = Ytotal*(sin(((J/nytot)*(PI/2))));
						}
						else if(OpcionY == 5){//Cúbica
							J = j;
							M2[i][j][1] = (Ytotal/2)*(1 + pow(2*((J+1)/(nytot+1))-1,3));
						}
					}

				}
			}

		//Discretización horizontal de la malla 3
			for (i = 0; i < NXTOT+1; i++){
				for (j = 0; j < NYTOT; j++){
			
					if(i==0){
						M3[i][j][0] = 0;
					}
					else if(i == NXTOT){
						M3[i][j][0] = Xtotal;
					}
					else{
						if(OpcionX == 1){//Regular
							DeltaX = Xtotal/nxtot;
							M3[i][j][0] = DeltaX*i;
						}
						else if(OpcionX == 2){//Tangente Hiperbólica
							I = i;
							M3[i][j][0] = (Xtotal/2)*(1+(tanh(StretFactor*((2*I-nxtot)/nxtot)) + tanh(StretFactor))/tanh(StretFactor) - 1);
						}
						else if(OpcionX == 3){//Cosenoidal
							I = i;
							M3[i][j][0] = Xtotal*(1 - cos((I/nxtot)*(PI/2)));
						}
						else if(OpcionX == 4){//Senoidal
							I = i;
							M3[i][j][0] = Xtotal*(sin((I/nxtot)*(PI/2)));
						}
						else if(OpcionX == 5){//Cúbica
							I = i;
							M3[i][j][0] = (Xtotal/2)*(1 + pow(2*((I+1)/(nxtot+1))-1,3));
						}
					}
				}
			}



		//Discretización horizontal y vertical de la malla 1
			for (i = 0; i < NXTOT; i++){
				for (j = 0; j < NYTOT; j++){
					//Discretización horizontal
					M1[i][j][0] = 0.5*(M3[i][j][0] + M3[i+1][j][0]);

					//Discretización vertical
					M1[i][j][1] = 0.5*(M2[i][j][1] + M2[i][j+1][1]);
				}	
			}

		//Discretización horizontal de la malla 2
			for (i = 0; i < NXTOT; i++){
				for (j = 0; j < NYTOT+1; j++){
					M2[i][j][0] = 0.5*(M3[i][1][0] + M3[i+1][1][0]);
				}	
			}

		//Discretización vertical de la malla 3	
			for (i = 0; i < NXTOT+1; i++){
				for (j = 0; j < NYTOT; j++){
					M3[i][j][1] = 0.5*(M2[1][j][1] + M2[1][j+1][1]);
				}
			}
	}
	else if(Problema == 3){ //Problema 3
	
	int NXP, NYP, NXA, NYA;
	double X1, X2, Y1, Y2;
	int Opcion1, Opcion2;
	int Opcion;

		for(Opcion = 1; Opcion <= 9; Opcion++){

			switch (Opcion){
				case 1:
					NXP = NX1;
					NYP = NY1;
					X1 = 0;
					Y1 = 0;
					X2 = _X1;
					Y2 = _Y1;
					NXA = 0;
					NYA = 0;
					Opcion1 = OpcionX;
					Opcion2 = OpcionY;
					break;

				case 2:
					NXP = NX2;
					NYP = NY1;
					X1 = _X1;
					Y1 = 0;
					X2 = _X1 + _X2;
					Y2 = _Y1;
					NXA = NX1;
					NYA = 0;
					if(Progresion == 1){ Opcion1 = 6;}		
					else{ Opcion1 = OpcionX;}
					Opcion2 = OpcionY;
					break;

				case 3:
					NXP = NX3;
					NYP = NY1;
					X1 = _X1 + _X2;
					Y1 = 0;
					X2 = Xtotal;
					Y2 = _Y1;
					NXA = NX1 + NX2;
					NYA = 0;
					Opcion1 = OpcionX;
					Opcion2 = OpcionY;
					break;

				case 4:
					NXP = NX1;
					NYP = NY2;
					X1 = 0;
					Y1 = _Y1;
					X2 = _X1;
					Y2 = _Y1 + _Y2;
					NXA = 0;
					NYA = NY1;
					Opcion1 = OpcionX;
					if(Progresion == 1){ Opcion2 = 6;}
					else{ Opcion2 = OpcionY;}
					break;

				case 5:
					NXP = NX3;
					NYP = NY2;
					X1 = _X1 + _X2;
					Y1 = _Y1;
					X2 = Xtotal;
					Y2 = _Y1 + _Y2;
					NXA = NX1 + NX2;
					NYA = NY1;
					Opcion1 = OpcionX;
					if(Progresion == 1){ Opcion2 = 6;}
					else{ Opcion2 = OpcionY;}
					break;

				case 6:
					NXP = NX1;
					NYP = NY3;
					X1 = 0;
					Y1 = _Y1 + _Y2;
					X2 = _X1;
					Y2 = Ytotal;
					NXA = 0;
					NYA = NY1 + NY2;
					Opcion1 = OpcionX;
					Opcion2 = OpcionY;
					break;

				case 7:
					NXP = NX2;
					NYP = NY3;
					X1 = _X1;
					Y1 = _Y1 + _Y2;
					X2 = _X1 + _X2;
					Y2 = Ytotal;
					NXA = NX1;
					NYA = NY1 + NY2;
					if(Progresion == 1){ Opcion1 = 6;}
					else{ Opcion1 = OpcionX;}
					Opcion2 = OpcionY;
					break;

				case 8:
					NXP = NX3;
					NYP = NY3;
					X1 = _X1 + _X2;
					Y1 = _Y1 + _Y2;
					X2 = Xtotal;
					Y2 = Ytotal;
					NXA = NX1 + NX2;
					NYA = NY1 + NY2;
					Opcion1 = OpcionX;
					Opcion2 = OpcionY;
					break;
				
				case 9:
					NXP = NX2;
					NYP = NY2;
					X1 = _X1;
					Y1 = _Y1;
					X2 = _X1 + _X2;
					Y2 = _Y1 + _Y2;
					NXA = NX1;
					NYA = NY1;
					if(Progresion == 1){ 
						Opcion1 = 6;
						Opcion2 = 6;
					}
					else{
						Opcion1 = OpcionX;
						Opcion2 = OpcionY;
					}
					break;
			}

		double X = X2 - X1;
		double Y = Y2 - Y1;

		double nx = NXP;
		double ny = NYP;

		double DeltaX;
		double DeltaY;

		double nx2 = NX2;
		double ny2 = NY2;


		//Discretización vertical de la malla 2
			for (i = 0; i < NXP; i++){
				for (j = 0; j < NYP+1; j++){
			
					if(j == 0){
						M2[i+NXA][j+NYA][1] = Y1;
					}
					else if (j == NYP){
						M2[i+NXA][j+NYA][1] = Y2;
					}
					else{

						if(Opcion2 == 1){//Regular
							DeltaY = Y/ny;
							M2[i+NXA][j+NYA][1] = Y1 + DeltaY*j;
						}
						else if(Opcion2 == 2){//Tangente Hiperbólica
							J = j;
							M2[i+NXA][j+NYA][1] = Y1 + (Y/2)*(1+(tanh(StretFactor*((2*J-ny)/ny)) + tanh(StretFactor))/tanh(StretFactor) - 1);
						}
						else if(Opcion2 == 3){//Cosenoidal
							J = j;
							M2[i+NXA][j+NYA][1] = Y1 + Y*(1 - cos(((J/ny)*(PI/2))));
						}
						else if(Opcion2 == 4){//Senoidal
							J = j;
							M2[i+NXA][j+NYA][1] = Y1 + Y*(sin(((J/ny)*(PI/2))));
						}
						else if(Opcion2 == 5){//Cúbica
							J = j;
							M2[i+NXA][j+NYA][1] = Y1 + (Y/2)*(1 + pow(2*((J+1)/(ny+1))-1,3));
						}
						else if(Opcion2 == 6){
							if(j==0){
								M2[i+NXA][j+NYA][1] = Y1;
							}
							else if(j==NYP){
								M2[i+NXA][j+NYA][1] = Y2;
							}
							else{
								M2[i+NXA][j+NYA][1] = M2[i+NXA][j+NYA-1][1] + Delta1 + (j-1)*((Delta2 - Delta1)/(ny2));
							}
						}
					}

				}
			}

		//Discretización horizontal de la malla 3
			for (i = 0; i < NXP+1; i++){
				for (j = 0; j < NYP; j++){
			
					if(i == 0){
						M3[i+NXA][j+NYA][0] = X1;
					}
					else if(i == NXP){
						M3[i+NXA][j+NYA][0] = X2;
					}
					else{
						if(Opcion1 == 1){//Regular
							DeltaX = X/nx;
							M3[i+NXA][j+NYA][0] = X1 + DeltaX*i;
						}
						else if(Opcion1 == 2){//Tangente Hiperbólica
							I = i;
							M3[i+NXA][j+NYA][0] = X1 + (X/2)*(1+(tanh(StretFactor*((2*I-nx)/nx)) + tanh(StretFactor))/tanh(StretFactor) - 1);
						}
						else if(Opcion1 == 3){//Cosenoidal
							I = i;
							M3[i+NXA][j+NYA][0] = X1 + X*(1 - cos((I/nx)*(PI/2)));
						}
						else if(Opcion1 == 4){//Senoidal
							I = i;
							M3[i+NXA][j+NYA][0] = X1 + X*(sin((I/nx)*(PI/2)));
						}
						else if(Opcion1 == 5){//Cúbica
							I = i;
							M3[i+NXA][j+NYA][0] = X1 + (X/2)*(1 + pow(2*((I+1)/(nx+1))-1,3));
						}
						else if(Opcion1 == 6){
							if(i == 0){
								M3[i + NXA][j + NYA][0] = X1;
							}
							else if(i == NXP){
								M3[i + NXA][j + NYA][0] = X2;
							}
							else{
								M3[i + NXA][j + NYA][0] = M3[i+NXA-1][j+NYA][0] + Delta3 + (i-1)*((Delta4 - Delta3)/(nx2));
							}
				
						}
					}
				}
			}



		//Discretización horizontal y vertical de la malla 1
			for (i = 0; i < NXP; i++){
				for (j = 0; j < NYP; j++){
					//Discretización horizontal
					M1[i+NXA][j+NYA][0] = 0.5*(M3[i+NXA][j+NYA][0] + M3[i+NXA+1][j+NYA][0]);

					//Discretización vertical
					M1[i+NXA][j+NYA][1] = 0.5*(M2[i+NXA][j+NYA][1] + M2[i+NXA][j+NYA+1][1]);
				}	
			}

		//Discretización horizontal de la malla 2
			for (i = 0; i < NXP; i++){
				for (j = 0; j < NYP+1; j++){
					M2[i+NXA][j+NYA][0] = 0.5*(M3[i+NXA][1+NYA][0] + M3[i+NXA+1][1+NYA][0]);
				}	
			}

		//Discretización vertical de la malla 3	
			for (i = 0; i < NXP+1; i++){
				for (j = 0; j < NYP; j++){
					M3[i+NXA][j+NYA][1] = 0.5*(M2[1+NXA][j+NYA][1] + M2[1+NXA][j+NYA+1][1]);

				}
			}
		}
	}
}


//Longitud de las caras de los volúmenes de control
void Mallador::Get_Walls(){
int i, j;
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			DeltasYP[i][j] = (M2[i][j+1][1] - M2[i][j][1]);
			DeltasXP[i][j] = (M3[i+1][j][0] - M3[i][j][0]);
		}
	}

	for(i = 0; i < NXTOT+1; i++){
		for(j = 0; j < NYTOT; j++){
			DeltasYU[i][j] = (M2[1][j+1][1] - M2[1][j][1]);
			if(i == 0){
			DeltasXU[i][j] = M1[i][j][0];
			}
			else if(i == NXTOT){
			DeltasXU[i][j] = Xtotal - M1[i-1][j][0];
			}
			else{
			DeltasXU[i][j] = M1[i][j][0] - M1[i-1][j][0];
			}	
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT+1; j++){
			if (j == 0){
			DeltasYV[i][j] = M1[i][j][1];
			}
			else if(j == NYTOT){
			DeltasYV[i][j] = Ytotal - M1[i][j-1][1];
			}
			else{
			DeltasYV[i][j] = M1[i][j][1] - M1[i][j-1][1];
			}
			DeltasXV[i][j] = M3[i+1][1][0] - M3[i][1][0];
		}
	}
}


void Mallador::Get_Surfaces(){
int i, j;

	for(i = 0; i < NXTOT+1; i++){
		for(j = 0; j < NYTOT; j++){
			SurfaceXU[i][j] = Wtotal*DeltasXU[i][j];
			SurfaceYU[i][j] = Wtotal*DeltasYU[i][j];
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT+1; j++){
			SurfaceXV[i][j] = Wtotal*DeltasXV[i][j];
			SurfaceYV[i][j] = Wtotal*DeltasYV[i][j];
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			SurfaceXP[i][j] = Wtotal*DeltasXP[i][j];
			SurfaceYP[i][j] = Wtotal*DeltasYP[i][j];
		}
	}


}

//Ejecución del todo el mallador
void Mallador::EjecucionMesh(){

	Get_Mesh(); //Creación de las 3 mallas
	Get_Walls(); //Cálculo de las distancias de los volúmenes de control
	Get_Surfaces(); //Cálculo de las superficies de los volúmenes de control
	PrintTxt(); //Pasar los resultados de las 3 mallas a un TXT

}


//Pasar los resultados de las mallas a un txt
void Mallador::PrintTxt(){
	
	int i, j;

	FILE *fp1;
	fp1 = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Malla/M1.txt","w");
	for(i = 0; i < NXTOT; i++){
        	for(j = 0; j < NYTOT; j++){
			fprintf(fp1,"%f %f \n",M1[i][j][0],M1[i][j][1]);
		}
        	fprintf(fp1, "\n");
    	}
	fclose(fp1);

	FILE *fp2;
	fp2 = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Malla/M2.txt","w");
	for(i = 0; i < NXTOT; i++){
        	for(j = 0; j < NYTOT+1; j++){
			fprintf(fp1,"%f %f \n",M2[i][j][0],M2[i][j][1]);
		}
        	fprintf(fp2, "\n");
    	}
	fclose(fp2);

	FILE *fp3;
	fp3 = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Malla/M3.txt","w");
	for(i = 0; i < NXTOT+1; i++){
        	for(j = 0; j < NYTOT; j++){
			fprintf(fp1,"%f %f \n",M3[i][j][0],M3[i][j][1]);
		}
        	fprintf(fp3, "\n");
    	}
	fclose(fp3);

}

