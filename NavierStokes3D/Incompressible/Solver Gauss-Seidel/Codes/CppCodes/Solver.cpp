#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>
#include <chrono>

using namespace std;

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/ParPro.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/Solver.h"

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

//Coeficientes A
#define LAG(i,j,k,dim) ((NY*NZ)*((i))) + ((j) + (k)*NY)
#define LAL(i,j,k,dim) ((NY*NZ)*((i) - Ix)) + ((j) + (k)*NY)

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
Solver::Solver(Memory M1, ReadData R1, ParPro MPI1, Mesher MESH, PostProcessing POST1, string InputDirectorio){
		
	//Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3]; 
	NZ = R1.ProblemNumericalData[4];

	StepsPantalla = R1.ProblemNumericalData[8];
	StepsFile = R1.ProblemNumericalData[9];

	ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

	EsquemaLargo = R1.ConvectiveScheme1;
	EsquemaCorto = R1.ConvectiveScheme2;

	DIRECTORIO = InputDirectorio;

	//Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	Xcentroide = R1.GeometryData[3];
	Ycentroide = R1.GeometryData[4];

	Xcuadrado = R1.GeometryData[5];
	Ycuadrado = R1.GeometryData[6];

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

	//Variables del Runge - Kutta de 4o order (Capuano et All)
	/*
	c1 = 0.0;
	c3 = 0.25;
	c2 = (c3 - 1.0)/(4.0*c3 - 3.0);
	c4 = 1.0;

	b1 = 1.0/(12.0*(c3 - 1.0));
	b2 = pow(4.0*c3 - 3.0,2.0)/(12.0*(c3 - 1.0)*(2.0*c3 - 1.0));
	b3 = - 1.0/(12.0*(c3 - 1.0)*(2.0*c3 - 1.0));
	b4 = (4.0*c3 - 3.0)/(12.0*(c3 - 1.0));


	a_21 = (c3 - 1.0)/(4.0*c3 - 3.0);

	a_31 = c3 - ((2.0*c3 - 1.0)*(4.0*c3 - 3.0))/(2.0*(c3 - 1));
	a_32 = ((2.0*c3 - 1.0)*(4.0*c3 - 3.0))/(2.0*(c3 - 1));

	a_41 = - pow(2.0*c3 - 1.0,2.0)/(2.0*(c3 - 1.0)*(4.0*c3 - 3.0));
	a_42 = (6.0*pow(c3,2.0) - 8.0*c3 + 3.0)/(2.0*(c3 - 1.0)*(2.0*c3 - 1.0));
	a_43 = (c3 - 1.0)/((2.0*c3 - 1.0)*(4.0*c3 - 3.0));
*/

	c1 = 0.0;
	c2 = 0.5;
	c3 = 0.5;
	c4 = 1.0;

	b1 = 1.0/6.0;
	b2 = 1.0/3.0;
	b3 = 1.0/3.0;
	b4 = 1.0/6.0;


	a_21 = 0.5;

	a_31 = 0.0;
	a_32 = 0.5;

	a_41 = 0.0;
	a_42 = 0.0;
	a_43 = 1.0;

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

//Alojamiento de memoria para las matrices necesarias
void Solver::AllocateMatrix(Memory M1){

	//Matrices del Solver

	//Coeficientes de discretizacion
	aw = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);
	ae = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);

	as = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);
	an = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);

	ah = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);
	at = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);

	ap = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);

	//Divergencia Predictoras
	bp = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);

	if(Rank != 0 && Rank != Procesos - 1){

		//Matrices locales de propiedades

		//Presion
		PLFUT = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);
		PLSUP = M1.AllocateDouble(Fx - Ix + 2*Halo, NY, NZ, 1);

		//Velocidad U
		ULPRES = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		ULFUT = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		//Velocidad V
		VLPRES = M1.AllocateDouble(Fx - Ix + 2*Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		VLFUT = M1.AllocateDouble(Fx - Ix + 2*Halo, NY+1 + 2*HP, NZ + 2*HP, 1);

		//Velocidad W
		WLPRES = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		WLFUT = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ+1 + 2*HP, 1);

		//Matrices de las contribuciones de cada una de las ecuaciones

		ConvectiveU = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		ConvectiveV = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		ConvectiveW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);
		ConvectiveT = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		DiffusiveU = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		DiffusiveV = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		DiffusiveW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);
		DiffusiveT = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		K1_U = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		K1_V = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		K1_W = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);
		K1_T = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		K2_U = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		K2_V = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		K2_W = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);
		K2_T = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		K3_U = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		K3_V = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		K3_W = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);
		K3_T = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		K4_U = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		K4_V = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		K4_W = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);
		K4_T = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);

		UL_NEW = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		VL_NEW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		WL_NEW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);

		PredU = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
		PredV = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		PredW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);

	}
	else if(Rank == 0){

		//Matrices globales de propiedades

		//Presion
		PGPRES = M1.AllocateDouble(NX, NY, NZ, 1); //Presión P Global Presente

		//Velocidad U
		UGPRES = M1.AllocateDouble(NX + 1 + 2*HP,NY + 2*HP,NZ + 2*HP,1); //Velocidad U Global Presente
		UGFUT = M1.AllocateDouble(NX + 1 + 2*HP,NY + 2*HP,NZ + 2*HP,1); //Velocidad U Global Futuro

		//Velocidad V
		VGPRES = M1.AllocateDouble(NX + 2*HP,NY+1 + 2*HP,NZ + 2*HP,1); //Velocidad V Global Presente
		VGFUT = M1.AllocateDouble(NX + 2*HP,NY+1 + 2*HP,NZ + 2*HP,1); //Velocidad V Global Futuro

		//Velocidad W
		WGPRES = M1.AllocateDouble(NX + 2*HP,NY + 2*HP,NZ + 1 + 2*HP,1); //Velocidad W Global Presente
		WGFUT = M1.AllocateDouble(NX + 2*HP,NY + 2*HP,NZ + 1 + 2*HP,1); //Velocidad W Global Futuro

		//Matrices locales de propiedades

		//Presion
		PLFUT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY, NZ, 1);
		PLSUP = M1.AllocateDouble(Fx - Ix + HP + Halo, NY, NZ, 1);

		//Velocidad U
		ULPRES = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		ULFUT = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		//Velocidad V
		VLPRES = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		VLFUT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);

		//Velocidad W
		WLPRES = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		WLFUT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);

		//Matrices de las contribuciones de cada una de las ecuaciones
		ConvectiveU = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		ConvectiveV = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		ConvectiveW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		ConvectiveT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		DiffusiveU = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		DiffusiveV = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		DiffusiveW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		DiffusiveT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K1_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K1_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K1_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K1_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K2_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K2_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K2_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K2_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K3_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K3_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K3_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K3_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K4_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K4_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K4_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K4_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		UL_NEW = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		VL_NEW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		WL_NEW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);

		PredU = M1.AllocateDouble(Fx - Ix + 1 + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
		PredV = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 1 + 2*HP, NZ + 2*HP, 1);
		PredW = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 1 + 2*HP, 1);

	}
	else if(Rank == Procesos - 1){

		//Matrices locales de propiedades

		//Presion
		PLFUT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY, NZ, 1);
		PLSUP = M1.AllocateDouble(Fx - Ix + HP + Halo, NY, NZ, 1);

		//Velocidad U
		ULPRES = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		ULFUT = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		//Velocidad V
		VLPRES = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		VLFUT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);

		//Velocidad W
		WLPRES = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		WLFUT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);

		//Matrices de las contribuciones de cada una de las ecuaciones
		ConvectiveU = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		ConvectiveV = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		ConvectiveW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		ConvectiveT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		DiffusiveU = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		DiffusiveV = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		DiffusiveW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		DiffusiveT = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K1_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K1_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K1_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K1_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K2_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K2_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K2_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K2_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K3_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K3_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K3_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K3_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		K4_U = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		K4_V = M1.AllocateDouble(Fx - Ix + HP + Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
		K4_W = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ+1 + 2*HP, 1);
		K4_T = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);

		UL_NEW = M1.AllocateDouble(Fx - Ix + 1 + HP + Halo, NY + 2*HP, NZ + 2*HP, 1);
		VL_NEW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 1 + 2*HP, NZ + 2*HP, 1);
		WL_NEW = M1.AllocateDouble(Fx - Ix + HP + Halo, NY + 2*HP, NZ + 1 + 2*HP, 1);

		PredU = M1.AllocateDouble(Fx - Ix + 1 + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
		PredV = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 1 + 2*HP, NZ + 2*HP, 1);
		PredW = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 1 + 2*HP, 1);

	}
	

	//Matrices de condiciones de contorno
	if(Rank == 0){
		Uleft = M1.AllocateDouble(1, NY, NZ, 1);
		Vleft = M1.AllocateDouble(1, NY + 1, NZ, 1);
		Wleft = M1.AllocateDouble(1, NY, NZ + 1, 1);
		TLEFT = M1.AllocateDouble(1, NY, NZ, 1);
	}	

	if(Rank == Procesos - 1){
		Uright = M1.AllocateDouble(1, NY, NZ, 1);
		Vright = M1.AllocateDouble(1, NY + 1, NZ, 1);
		Wright = M1.AllocateDouble(1, NY, NZ + 1, 1);
		TRIGHT = M1.AllocateDouble(1, NY, NZ, 1);
	}

	//Parte Abajo
	Ubot = M1.AllocateDouble(Fx - Ix + 1, 1, NZ, 1);
	Vbot = M1.AllocateDouble(Fx - Ix + 2, 1, NZ, 1);
	Wbot = M1.AllocateDouble(Fx - Ix + 2, 1, NZ, 1);
	TBOT = M1.AllocateDouble(Fx - Ix + 1, 1, NZ, 1);

	//Parte Arriba
	Utop = M1.AllocateDouble(Fx - Ix + 1, 1, NZ, 1);
	Vtop = M1.AllocateDouble(Fx - Ix + 2, 1, NZ, 1);
	Wtop = M1.AllocateDouble(Fx - Ix + 2, 1, NZ, 1);
	TTOP = M1.AllocateDouble(Fx - Ix + 1, 1, NZ, 1);

	//Parte Here
	Uhere = M1.AllocateDouble(Fx - Ix + 1, NY, 1, 1);
	Vhere = M1.AllocateDouble(Fx - Ix + 2, NY + 1, 1, 1);
	Where = M1.AllocateDouble(Fx - Ix + 2, NY, 1, 1);
	There = M1.AllocateDouble(Fx - Ix, NY, 1, 1);

	//Parte There
	Uthere = M1.AllocateDouble(Fx - Ix + 1, NY, 1, 1);
	Vthere = M1.AllocateDouble(Fx - Ix + 2, NY + 1, 1, 1);
	Wthere = M1.AllocateDouble(Fx - Ix + 2, NY, 1, 1);
	Tthere = M1.AllocateDouble(Fx - Ix, NY, 1, 1);

}

//Seteo de las matrices iniciales
void Solver::Get_InitialConditions(){
int i, j, k;

	//Presión
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				PLFUT[LPC(i,j,k,0)] = 0.0;
				PLSUP[LPC(i,j,k,0)] = PLFUT[LPC(i,j,k,0)];
			}
		}
	}
	
	//Velocidad U
	for(i = Ix; i < Fx + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	
				ULPRES[LUC(i,j,k,0)] = 0.1;
				ULFUT[LUC(i,j,k,0)] = 0.0;
			}
		}
	}
	
	//Velocidad V
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){
				VLPRES[LVC(i,NY,k,0)] = 0.0;
				VLFUT[LVC(i,NY,k,0)] = 0.0;
			}
		}
		
	}

	//Velocidad W
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){			
				WLPRES[LWC(i,j,k,0)] = 0.0;
				WLFUT[LWC(i,j,k,0)] = 0.0;
			}
		}
	}
	
}

//Seteo y Actualización de las condiciones de contorno
void Solver::Get_BoundaryConditions(Mesher MESH){
int i, j, k;

	switch(Problema){

		//Problema Driven Cavity
		case 1:

		if(Rank == 0){

			//Parte Izquierda

			//Velocidad U
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					Uleft[ULEFT(0,j,k)] = 0.0;
				}
			}	

			//Velocidad V
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					Vleft[VLEFT(0,j,k)] = 0.0;
				}
			}

			//Velocidad W
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					Wleft[WLEFT(0,j,k)] = 0.0;
				}
			}	

		}
		else if(Rank == Procesos - 1){

			//Parte Derecha

			//Velocidad U
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					Uright[URIGHT(0,j,k)] = 0.0;
				}
			}	

			//Velocidad V
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					Vright[VRIGHT(0,j,k)] = 0.0;
				}
			}

			//Velocidad W
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					Wright[WRIGHT(0,j,k)] = 0.0;
				}
			}

		}

		//Parte Inferior

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Ubot[UBOT(i,0,k)] = 0.0;
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Vbot[VBOT(i,0,k)] = 0.0;
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ + 1; k++){
				Wbot[WBOT(i,0,k)] = 0.0;
			}
		}


		//Parte Superior

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Utop[UTOP(i,0,k)] = 1.0;
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Vtop[VTOP(i,0,k)] = 0.0;
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ + 1; k++){
				Wtop[WTOP(i,0,k)] = 0.0;
			}
		}
		

		//Parte Here

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Uhere[UHERE(i,j,0)] = ULPRES[LUC(i,j,NZ-1,0)];
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY + 1; j++){
				Vhere[VHERE(i,j,0)] = VLPRES[LVC(i,j,NZ-1,0)];
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Where[WHERE(i,j,0)] = WLPRES[LWC(i,j,NZ-1,0)];
			}
		}


		//Parte There

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Uthere[UTHERE(i,j,k)] = ULPRES[LUC(i,j,0,0)];
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY + 1; j++){
				Vthere[VTHERE(i,j,k)] = VLPRES[LVC(i,j,0,0)];
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Wthere[WTHERE(i,j,k)] = WLPRES[LWC(i,j,0,0)];
			}
		}

		break;


		//Problema Dfferentially Heated
		case 2:

		if(Rank == 0){

			//Parte Izquierda

			//Velocidad U
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					Uleft[ULEFT(0,j,k)] = 0.0;
				}
			}	

			//Velocidad V
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					Vleft[VLEFT(0,j,k)] = 0.0;
				}
			}

			//Velocidad W
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					Wleft[WLEFT(0,j,k)] = 0.0;
				}
			}	

			//Temperatura
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					TLEFT[PLEFT(0,j,k)] = Tleft;
				}
			}	

		}
		else if(Rank == Procesos - 1){

			//Parte Derecha

			//Velocidad U
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					Uright[URIGHT(0,j,k)] = 0.0;
				}
			}	

			//Velocidad V
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					Vright[VRIGHT(0,j,k)] = 0.0;
				}
			}

			//Velocidad W
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					Wright[WRIGHT(0,j,k)] = 0.0;
				}
			}

			//Temperatura
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					TRIGHT[PRIGHT(0,j,k)] = Tright;
				}
			}

		}

		//Parte Inferior

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Ubot[UBOT(i,0,k)] = 0.0;
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Vbot[VBOT(i,0,k)] = 0.0;
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ + 1; k++){
				Wbot[WBOT(i,0,k)] = 0.0;
			}
		}

		//Temperatura
		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				TBOT[PBOT(i,0,k)] = TLPRES[LPC(i,0,k,0)];
			}
		}

		//Parte Superior

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Utop[UTOP(i,0,k)] = 0.0;
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Vtop[VTOP(i,0,k)] = 0.0;
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ + 1; k++){
				Wtop[WTOP(i,0,k)] = 0.0;
			}
		}
		
		//Temperatura
		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				TTOP[PTOP(i,NY,k)] = TLPRES[LPC(i,NY - 1,k,0)];
			}
		}

		//Parte Here

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Uhere[UHERE(i,j,0)] = ULPRES[LUC(i,j,NZ-1,0)];
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY + 1; j++){
				Vhere[VHERE(i,j,0)] = VLPRES[LVC(i,j,NZ-1,0)];
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Where[WHERE(i,j,0)] = WLPRES[LWC(i,j,NZ-1,0)];
			}
		}

		//Temperatura
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				There[PHERE(i,j,0)] = TLPRES[LPC(i,j,NZ-1,0)];
			}
		}


		//Parte There

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Uthere[UTHERE(i,j,k)] = ULPRES[LUC(i,j,0,0)];
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY + 1; j++){
				Vthere[VTHERE(i,j,k)] = VLPRES[LVC(i,j,0,0)];
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Wthere[WTHERE(i,j,k)] = WLPRES[LWC(i,j,0,0)];
			}
		}

		//Temperatura
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				Tthere[PTHERE(i,j,k)] = TLPRES[LPC(i,j,0,0)];
			}
		}

		break;


		//Problema Differentially Heated
		case 3:

		if(Rank == 0){

			//Parte Izquierda

			//Velocidad U
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					Uleft[ULEFT(0,j,k)] = Uref*(1.0 - pow(MESH.MU[GU(0,j,k,1)] - Ydominio/2.0,2)/pow(Ydominio/2.0,2));
				}
			}	

			//Velocidad V
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					Vleft[VLEFT(0,j,k)] = 0.0;
				}
			}

			//Velocidad W
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					Wleft[WLEFT(0,j,k)] = 0.0;
				}
			}	

		}
		else if(Rank == Procesos - 1){

			//Parte Derecha

			//Velocidad U
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					Uright[URIGHT(0,j,k)] = ULPRES[LUR(NX - 1,j,k,0)];
				}
			}	

			//Velocidad V
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					Vright[VRIGHT(0,j,k)] = VLPRES[LVR(NX - 1,j,k,0)];
				}
			}

			//Velocidad W
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					Wright[WRIGHT(0,j,k)] = WLPRES[LWR(NX - 1,j,k,0)];
				}
			}

		}

		//Parte Inferior

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Ubot[UBOT(i,0,k)] = 0.0;
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Vbot[VBOT(i,0,k)] = 0.0;
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ + 1; k++){
				Wbot[WBOT(i,0,k)] = 0.0;
			}
		}

		//Parte Superior

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Utop[UTOP(i,0,k)] = 0.0;
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				Vtop[VTOP(i,0,k)] = 0.0;
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(k = 0; k < NZ + 1; k++){
				Wtop[WTOP(i,0,k)] = 0.0;
			}
		}
		
		//Parte Here

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Uhere[UHERE(i,j,0)] = ULPRES[LUC(i,j,NZ-1,0)];
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY + 1; j++){
				Vhere[VHERE(i,j,0)] = VLPRES[LVC(i,j,NZ-1,0)];
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Where[WHERE(i,j,0)] = WLPRES[LWC(i,j,NZ-1,0)];
			}
		}

		//Parte There

		//Velocidad U
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Uthere[UTHERE(i,j,k)] = ULPRES[LUC(i,j,0,0)];
			}
		}

		//Velocidad V
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY + 1; j++){
				Vthere[VTHERE(i,j,k)] = VLPRES[LVC(i,j,0,0)];
			}
		}

		//Velocidad W
		for(i = Ix - 1; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				Wthere[WTHERE(i,j,k)] = WLPRES[LWC(i,j,0,0)];
			}
		}

		break;

	}
	
}

//Actualización del Halo de Velocidades en cada Step
void Solver::Get_HaloVelocities(){
int i, j, k;

	//Velocidad U
	if(Rank == 0){

		//Halo Izquierda
		for(i = - HP; i < 0; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LUC(i,j,k,0)] = Uleft[ULEFT(0,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HP; j++){
					ULPRES[LUC(i,j,k,0)] = Utop[UTOP(0,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HP; j < 0; j++){
					ULPRES[LUC(i,j,k,0)] = Ubot[UBOT(0,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LUC(i,j,k,0)] = Uhere[UHERE(0,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HP; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LUC(i,j,k,0)] = Uthere[UTHERE(0,j,0)];
				}
			}

		}

	}
	else if(Rank == Procesos - 1){

		//Halo Derecha
		for(i = NX + 1; i < NX + HP + 1; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LUC(i,j,k,0)] = Uright[URIGHT(NX,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HP; j++){
					ULPRES[LUC(i,j,k,0)] = Utop[UTOP(NX,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HP; j < 0; j++){
					ULPRES[LUC(i,j,k,0)] = Ubot[UBOT(NX,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LUC(i,j,k,0)] = Uhere[UHERE(NX,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HP; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LUC(i,j,k,0)] = Uthere[UTHERE(NX,j,0)];
				}
			}

		}

	}

	//Resto de Halos
	for(i = Ix; i < Fx + 1; i++){

		//Halo Superior
		for(k = 0; k < NZ; k++){
			for(j = NY; j < NY + HP; j++){
				ULPRES[LUC(i,j,k,0)] = Utop[UTOP(i,NY,k)];
			}
		}

		//Halo Inferior
		for(k = 0; k < NZ; k++){
			for(j = - HP; j < 0; j++){
				ULPRES[LUC(i,j,k,0)] = Ubot[UBOT(i,0,k)];
			}
		}

		//Halo Here
		for(k = - HP; k < 0; k++){
			for(j = 0; j < NY; j++){
				ULPRES[LUC(i,j,k,0)] = ULPRES[LUC(i,j,NZ + k,0)];
			}
		}

		//Halo There
		for(k = NZ; k < NZ + HP; k++){
			for(j = 0; j < NY; j++){
				ULPRES[LUC(i,j,k,0)] = ULPRES[LUC(i,j,k - NZ,0)];
			}
		}

	}	

	//Velocidad V
	if(Rank == 0){

		//Halo Izquierda
		for(i = - HP; i < 0; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 1; j < NY; j++){
					VLPRES[LVC(i,j,k,0)] = Vleft[VLEFT(0,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HP + 1; j++){
					VLPRES[LVC(i,j,k,0)] = Vtop[VTOP(0,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HP; j <= 0; j++){
					VLPRES[LVC(i,j,k,0)] = Vbot[VBOT(0,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LVC(i,j,k,0)] = Vhere[VHERE(0,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HP; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LVC(i,j,k,0)] = Vthere[VTHERE(0,j,0)];
				}
			}

		}

	}
	else if(Rank == Procesos - 1){

		//Halo Derecha
		for(i = NX; i < NX + HP; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 1; j < NY; j++){
					VLPRES[LVC(i,j,k,0)] = Vright[VRIGHT(NX,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HP; j++){
					VLPRES[LVC(i,j,k,0)] = Vtop[VTOP(NX,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HP; j <= 0; j++){
					VLPRES[LVC(i,j,k,0)] = Vbot[VBOT(NX,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LVC(i,j,k,0)] = Vhere[VHERE(NX,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HP; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LVC(i,j,k,0)] = Vthere[VTHERE(NX,j,0)];
				}
			}

		}

	}

	//Resto de Halos
	for(i = Ix - 1; i < Fx + 1; i++){

		//Halo Superior
		for(k = 0; k < NZ; k++){
			for(j = NY + 1; j < NY + HP + 1; j++){
				VLPRES[LVC(i,j,k,0)] = Vtop[VTOP(i,NY,k)];
			}
		}

		//Halo Inferior
		for(k = 0; k < NZ; k++){
			for(j = - HP; j < 0; j++){
				VLPRES[LVC(i,j,k,0)] = Vbot[VBOT(i,0,k)];
			}
		}

		//Halo Here
		for(k = - HP; k < 0; k++){
			for(j = 0; j < NY + 1; j++){
				VLPRES[LVC(i,j,k,0)] = VLPRES[LVC(i,j,NZ + k,0)];
			}
		}

		//Halo There
		for(k = NZ; k < NZ + HP; k++){
			for(j = 0; j < NY + 1; j++){
				VLPRES[LVC(i,j,k,0)] = VLPRES[LVC(i,j,k - NZ,0)];
			}
		}

	}


	//Velocidad W
	if(Rank == 0){

		//Halo Izquierda
		for(i = - HP; i < 0; i++){

			//Halo Channel
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LWC(i,j,k,0)] = Wleft[WLEFT(0,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ + 1; k++){
				for(j = NY; j < NY + HP; j++){
					WLPRES[LWC(i,j,k,0)] = Wtop[WTOP(0,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ + 1; k++){
				for(j = - HP; j < 0; j++){
					WLPRES[LWC(i,j,k,0)] = Wbot[WBOT(0,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LWC(i,j,k,0)] = Where[WHERE(0,j,0)];
				}
			}

			//Halo There
			for(k = NZ + 1; k < NZ + HP + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LWC(i,j,k,0)] = Wthere[WTHERE(0,j,0)];
				}
			}

		}

	}
	else if(Rank == Procesos - 1){

		//Halo Derecha
		for(i = NX; i < NX + HP; i++){

			//Halo Channel
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LWC(i,j,k,0)] = Wright[WRIGHT(NX,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ + 1; k++){
				for(j = NY; j < NY + HP; j++){
					WLPRES[LWC(i,j,k,0)] = Wtop[WTOP(NX,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ + 1; k++){
				for(j = - HP; j < 0; j++){
					WLPRES[LWC(i,j,k,0)] = Wbot[WBOT(NX,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LWC(i,j,k,0)] = Where[WHERE(NX,j,0)];
				}
			}

			//Halo There
			for(k = NZ + 1; k < NZ + HP + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LWC(i,j,k,0)] = Wthere[WTHERE(NX,j,0)];
				}
			}

		}

	}

	//Resto de Halos
	for(i = Ix - 1; i < Fx + 1; i++){

		//Halo Superior
		for(k = 0; k < NZ + 1; k++){
			for(j = NY; j < NY + HP; j++){
				WLPRES[LWC(i,j,k,0)] = Wtop[WTOP(i,NY,k)];
			}
		}

		//Halo Inferior
		for(k = 0; k < NZ + 1; k++){
			for(j = - HP; j < 0; j++){
				WLPRES[LWC(i,j,k,0)] = Wbot[WBOT(i,0,k)];
			}
		}

		//Halo Here
		for(k = - HP; k < 0; k++){
			for(j = 0; j < NY; j++){
				WLPRES[LWC(i,j,k,0)] = WLPRES[LWC(i,j,NZ + k,0)];
			}
		}

		//Halo There
		for(k = NZ + 1; k < NZ + HP + 1; k++){
			for(j = 0; j < NY; j++){
				WLPRES[LWC(i,j,k,0)] = WLPRES[LWC(i,j,k - NZ,0)];
			}
		}

	}
	//Problema Differentially Heated
	if(Problema == 2){

		if(Rank == 0){

			//Parte Izquierda
			for(i = - HP; i < 0; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){
						TLPRES[LPC(i,j,k,0)] = TLEFT[PLEFT(0,j,k)];
					}
				}
			}

		}
		else if(Rank == Procesos - 1){

			//Parte Derecha
			for(i = NX; i < NX + HP; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){
						TLPRES[LPC(i,j,k,0)] = TRIGHT[PRIGHT(0,j,k)];
					}
				}
			}

		}

		//Resto de Halos
	
		for(i = Ix; i < Fx; i++){

			//Partes Superior e Inferior
			for(k = 0; k < NZ; k++){

				//Parte Superior
				for(j = NY; j < NY + HP; j++){
					TLPRES[LPC(i,j,k,0)] = TTOP[PTOP(i,NY,k)];
				}
				
				//Parte Inferior
				for(j = - HP; j < 0; j++){
					TLPRES[LPC(i,j,k,0)] = TBOT[PBOT(i,0,k)];
				}

			}

			//Partes Here y There
			for(j = 0; j < NY; j++){

				//Parte Here
				for(k = - HP; k < 0; k++){
					TLPRES[LPC(i,j,k,0)] = There[PHERE(i,j,0)];
				}

				//Parte There
				for(k = NZ; k < NZ + HP; k++){
					TLPRES[LPC(i,j,k,0)] = Tthere[PTHERE(i,j,NZ)];
				}

			}

		}

	}

}

//Cálculo del tiempo entre Steps
void Solver::Get_StepTime(Mesher MESH, ParPro MPI1){
int i, j , k;
DeltaT = 1000.0;
double Tpar = 0.60;

MPI_Status ST;

	//Comparacion 
	//a*(a>b) + b*(b>a)
	//b+=(a-b)*(a>b)

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){

				
				//CFL Velocidades (U, V, W)

				//Velocidad U
				DeltaT += (abs((Tpar*MESH.DeltasMP[GP(i,j,k,0)])/(0.50*(ULPRES[LUC(i,j,k,0)] + ULPRES[LUC(i+1,j,k,0)]) + 1e-10)) - DeltaT)*(abs((Tpar*MESH.DeltasMP[GP(i,j,k,0)])/(0.50*(ULPRES[LUC(i,j,k,0)] + ULPRES[LUC(i+1,j,k,0)]) + 1e-10)) <= DeltaT);

				//Velocidad V
				DeltaT += (abs((Tpar*MESH.DeltasMP[GP(i,j,k,1)])/(0.50*(VLPRES[LVC(i,j,k,0)] + VLPRES[LVC(i,j+1,k,0)]) + 1e-10)) - DeltaT)*(abs((Tpar*MESH.DeltasMP[GP(i,j,k,1)])/(0.50*(VLPRES[LVC(i,j,k,0)] + VLPRES[LVC(i,j+1,k,0)]) + 1e-10)) <= DeltaT);
				
				//Velocidad W
				DeltaT += (abs((Tpar*MESH.DeltasMP[GP(i,j,k,2)])/(0.50*(WLPRES[LWC(i,j,k,0)] + WLPRES[LWC(i,j,k+1,0)]) + 1e-10)) - DeltaT)*(abs((Tpar*MESH.DeltasMP[GP(i,j,k,2)])/(0.50*(WLPRES[LWC(i,j,k,0)] + WLPRES[LWC(i,j,k+1,0)]) + 1e-10)) <= DeltaT);
				
			
				//CFL Difusivo Velocidades

				//Difusivo U
				DeltaT += ((Tpar*Rho*pow(MESH.DeltasMP[GP(i,j,k,0)],2.0))/(mu + 1e-10) - DeltaT)*((Tpar*Rho*pow(MESH.DeltasMP[GP(i,j,k,0)],2.0))/(mu + 1e-10) <= DeltaT);

				//Difusivo V
				DeltaT += ((Tpar*Rho*pow(MESH.DeltasMP[GP(i,j,k,1)],2.0))/(mu + 1e-10) - DeltaT)*((Tpar*Rho*pow(MESH.DeltasMP[GP(i,j,k,1)],2.0))/(mu + 1e-10) <= DeltaT);

				//Difusivo W
				DeltaT += ((Tpar*Rho*pow(MESH.DeltasMP[GP(i,j,k,2)],2.0))/(mu + 1e-10) - DeltaT)*((Tpar*Rho*pow(MESH.DeltasMP[GP(i,j,k,2)],2.0))/(mu + 1e-10) <= DeltaT);

			}	
		}
	}
	
	if(Problema == 2){

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){

					//CFL Difusivo Temperatura

					//Difusivo U
					DeltaT += ((Tpar*Rho*Cp*pow(MESH.DeltasMP[GP(i,j,k,0)],2.0))/(K + 1e-10) - DeltaT)*((Tpar*Rho*Cp*pow(MESH.DeltasMP[GP(i,j,k,0)],2.0))/(K + 1e-10) <= DeltaT);

					//Difusivo V
					DeltaT += ((Tpar*Rho*Cp*pow(MESH.DeltasMP[GP(i,j,k,1)],2.0))/(K + 1e-10) - DeltaT)*((Tpar*Rho*Cp*pow(MESH.DeltasMP[GP(i,j,k,1)],2.0))/(K + 1e-10) <= DeltaT);

					//Difusivo W
					DeltaT += ((Tpar*Rho*Cp*pow(MESH.DeltasMP[GP(i,j,k,2)],2.0))/(K + 1e-10) - DeltaT)*((Tpar*Rho*Cp*pow(MESH.DeltasMP[GP(i,j,k,2)],2.0))/(K + 1e-10) <= DeltaT);

				}
			}
		}

	}

	MPI_Allreduce(&DeltaT, &DeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

//Función con los diferentes esquemas convectivos utilizados
inline double Solver::ConvectiveScheme(double CoordObjetivo, double Velocity, double Coord1, double Phi1, double Coord2, double Phi2, double Coord3, double Phi3, double Coord4, double Phi4, string Esquema){

double PhiObjetivo;

double CoordD;
double PhiD;

double CoordC;
double PhiC;

double CoordU;
double PhiU;

	if (Velocity <= 0.0 || (Phi1 == 0.0 && Coord1 == 0.0)){

		CoordD = Coord2;
		PhiD = Phi2;
		CoordC = Coord3;
		PhiC = Phi3;
		CoordU = Coord4;
		PhiU = Phi4;

	}
	else if(Velocity > 0.0 || (Phi4 == 0.0 && Coord4 == 0.0)){

		CoordD = Coord3;
		PhiD = Phi3;
		CoordC = Coord2;
		PhiC = Phi2;
		CoordU = Coord1;
		PhiU = Phi1;

	}

	//Adimensionalizacion
	double PhiAdimC;
	double AdimCoordC;
	double AdimCoordE;

	PhiAdimC = (PhiC - PhiU)/(PhiD - PhiU);

	AdimCoordC = (CoordC - CoordU)/(CoordD - CoordU);

	AdimCoordE = (CoordObjetivo - CoordU)/(CoordD - CoordU);

	//Evaluacion
	double PhiF;

	if (PhiD == PhiU){
		PhiObjetivo = PhiD;
	}
	else{
		PhiF = AdimCoordE + (((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0))))*(PhiAdimC - AdimCoordC);

		//Dimensionalizacion
		PhiObjetivo = PhiU + (PhiD - PhiU)*PhiF;
	}

	return PhiObjetivo;

}

void Solver::Get_PressureCoefficients(Mesher MESH){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Parte Central
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){		
					aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
					ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

					as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
					an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

					ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
					at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

				}
			}		
		}

		//Parte Inferior
		j = 0;

		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = 0.0;
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte Superior
		j = NY - 1;

		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){	
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = 0.0;

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte Here
		k = 0;

		for(i = Ix; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte There
		k = NZ - 1;

		for(i = Ix; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Esquinas
		for(i = Ix; i < Fx; i++){

			//Esquina Abajo Here
			j = 0;
			k = 0;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);


			//Esquina Arriba There
			j = NY - 1; 
			k = NZ - 1;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		}

	}
	else if(Rank == 0){

		//Parte Central
		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){		
					aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
					ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

					as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
					an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

					ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
					at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

				}
			}		
		}

		//Parte Inferior
		j = 0;

		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = 0.0;
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte Superior
		j = NY - 1;

		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){	
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = 0.0;

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte Here
		k = 0;

		for(i = Ix + 1; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte There
		k = NZ - 1;

		for(i = Ix + 1; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Esquinas
		for(i = Ix + 1; i < Fx; i++){

			//Esquina Abajo Here
			j = 0;
			k = 0;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);


			//Esquina Arriba There
			j = NY - 1; 
			k = NZ - 1;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

		}

		//Parte Izquierda
		i = 0;

		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = 0.0;
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}
		}

		//Parte Inferior
		j = 0;

		for(k = 1; k < NZ - 1; k++){		
			aw[LAL(i,j,k,0)] = 0.0;
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

		}		
		
		//Parte Superior
		j = NY - 1;

		for(k = 1; k < NZ - 1; k++){	
			aw[LAL(i,j,k,0)] = 0.0;
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
		}

		//Parte Here
		k = 0;

		for(j = 1; j < NY - 1; j++){		
			aw[LAL(i,j,k,0)] = 0.0;
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
		}

		//Parte There
		k = NZ - 1;

		for(j = 1; j < NY - 1; j++){		
			aw[LAL(i,j,k,0)] = 0.0;
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
		}

		//Esquina Abajo Here
		j = 0;
		k = 0;

		aw[LAL(i,j,k,0)] = 0.0;
		ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

		as[LAL(i,j,k,0)] = 0.0;
		an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
		//Esquina Arriba Here
		j = NY - 1;
		k = 0;

		aw[LAL(i,j,k,0)] = 0.0;
		ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

		as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
		an[LAL(i,j,k,0)] = 0.0;

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

		//Esquina Abajo There
		j = 0;
		k = NZ - 1;

		aw[LAL(i,j,k,0)] = 0.0;
		ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

		as[LAL(i,j,k,0)] = 0.0;
		an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);


		//Esquina Arriba There
		j = NY - 1; 
		k = NZ - 1;

		aw[LAL(i,j,k,0)] = 0.0;
		ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

		as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
		an[LAL(i,j,k,0)] = 0.0;

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);


	}
	else if(Rank == Procesos - 1){

		//Parte Central
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){		
					aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
					ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

					as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
					an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

					ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
					at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

				}
			}		
		}

		//Parte Inferior
		j = 0;

		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = 0.0;
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte Superior
		j = NY - 1;

		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){	
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = 0.0;

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte Here
		k = 0;

		for(i = Ix; i < Fx - 1; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Parte There
		k = NZ - 1;

		for(i = Ix; i < Fx - 1; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			}		
		}

		//Esquinas
		for(i = Ix; i < Fx - 1; i++){

			//Esquina Abajo Here
			j = 0;
			k = 0;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);


			//Esquina Arriba There
			j = NY - 1; 
			k = NZ - 1;

			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i+1,j,k,0)]);

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		}


		//Parte Derecha
		i = NX - 1;

		//Parte Central
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){		
				aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
				ae[LAL(i,j,k,0)] = 0.0;

				as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
				an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

				ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
			}		
		}

		//Parte Inferior
		j = 0;

		for(k = 1; k < NZ - 1; k++){		
			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = 0.0;

			as[LAL(i,j,k,0)] = 0.0;
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);		
		}

		//Parte Superior
		j = NY - 1;

		for(k = 1; k < NZ - 1; k++){	
			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = 0.0;

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = 0.0;

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
		}

		//Parte Here
		k = 0;

		for(j = 1; j < NY - 1; j++){		
			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = 0.0;

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);		
		}

		//Parte There
		k = NZ - 1;

		for(j = 1; j < NY - 1; j++){		
			aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
			ae[LAL(i,j,k,0)] = 0.0;

			as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
			an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

			ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		}

		//Esquinas
		
		//Esquina Abajo Here
		j = 0;
		k = 0;

		aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
		ae[LAL(i,j,k,0)] = 0.0;

		as[LAL(i,j,k,0)] = 0.0;
		an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);
		
		//Esquina Arriba Here
		j = NY - 1;
		k = 0;

		aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
		ae[LAL(i,j,k,0)] = 0.0;

		as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
		an[LAL(i,j,k,0)] = 0.0;

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

		//Esquina Abajo There
		j = 0;
		k = NZ - 1;

		aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
		ae[LAL(i,j,k,0)] = 0.0;

		as[LAL(i,j,k,0)] = 0.0;
		an[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j+1,k,1)]);

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);


		//Esquina Arriba There
		j = NY - 1; 
		k = NZ - 1;

		aw[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMU[GU(i,j,k,0)]);
		ae[LAL(i,j,k,0)] = 0.0;

		as[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMV[GV(i,j,k,1)]);
		an[LAL(i,j,k,0)] = 0.0;

		ah[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LAL(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(Rho*MESH.VolMP[GP(i,j,k,0)]*MESH.DeltasMW[GW(i,j,k+1,2)]);

	}

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				ap[LAL(i,j,k,0)] = aw[LAL(i,j,k,0)] + ae[LAL(i,j,k,0)] + as[LAL(i,j,k,0)] + an[LAL(i,j,k,0)] + ah[LAL(i,j,k,0)] + at[LAL(i,j,k,0)];
			}
		}
	}
}

//Cálculo del término difusivo de la velocidad U
void Solver::Get_DiffusiveU(Mesher MESH, double *UFIELD){
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
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
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
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
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
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
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
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
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
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
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
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
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
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
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
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
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
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LUC(i,j+1,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LUC(i,j,k+1,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - Uhere[UHERE(i,j,k)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LUC(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LUC(i+1,j,k,0)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(Uthere[UTHERE(i,j,k)] - UFIELD[LUC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LUC(i,j,k,0)] - UFIELD[LUC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}	
		
}

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

//Cálculo del término convectivo de la velocidad U
void Solver::Get_ConvectiveU(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD){
int i, j, k;
double uW, uE, uS, uN, uH, uT, vS, vN, wH, wT;
double uW_pred, uE_pred, uS_pred, uN_pred, uH_pred, uT_pred, vS_pred, vN_pred, wH_pred, wT_pred;

	//ESQUEMA CONVECTIVO:
	//Corto (Interpolación)
	//(CoordObjetivo, Coord1, Valor1, Coord2, Valor2)


	//Largo
	//(CoordObjetivo, Velocidad, Coord1, Valor1, Coord2, Valor2, Coord3, Valor3, Coord4, Valor4, EsquemaLargo)

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro 
		for(i = Ix; i < Fx + 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

					uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
					uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

					uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
					uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

					vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

					uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

					uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
					uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

					uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
					uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
					wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

					ConvectiveU[LUC(i,j,k,0)] = (1.0/MESH.VolMU[GU(i,j,k,0)])*(
											  - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											  + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											  - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											  + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											  - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											  + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											  );
				}	
			}

			//Partes Inferior y Superior	

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Ubot[UBOT(i,j,k)];
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = Ubot[UBOT(i,j,k)];
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Utop[UTOP(i,j,k)];

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = Utop[UTOP(i,j,k)];

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}


			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j <  NY - 1; j++){

				
				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Uhere[UHERE(i,j,k)];
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = Uhere[UHERE(i,j,k)];
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Uthere[UTHERE(i,j,k)];

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = Uthere[UTHERE(i,j,k)];

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Esquina Abajo Here
			j = 0;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

		}		

	}
	else if(Rank == 0){

		 
		for(i = Ix + 1; i < Fx + 1; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

					uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
					uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

					uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
					uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

					vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

					uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

					uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
					uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

					uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
					uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
					wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

					ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Ubot[UBOT(i,j,k)];
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = Ubot[UBOT(i,j,k)];
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Utop[UTOP(i,j,k)];

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = Utop[VTOP(i,j,k)];

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}


			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j <  NY - 1; j++){	

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Uhere[UHERE(i,j,k)];
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = Uhere[UHERE(i,j,k)];
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Uthere[UTHERE(i,j,k)];

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = Uthere[UTHERE(i,j,k)];

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Esquina Abajo Here
			j = 0;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );


			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );


		}	

	}
	else if(Rank == Procesos - 1){

		 
		for(i = Ix; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

					uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
					uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

					uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
					uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

					vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

					uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

					uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
					uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

					uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
					uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
					wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

					ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

				}
			}


			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Ubot[UBOT(i,j,k)];
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = Ubot[UBOT(i,j,k)];
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Utop[UTOP(i,j,k)];

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = Utop[UTOP(i,j,k)];

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Partes Here y There

			//Parte Here
			k = 0;

			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Uhere[UHERE(i,j,k)];
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = Uhere[UHERE(i,j,k)];
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte There
			k = NZ - 1;

			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uT_pred = Uthere[UTHERE(i,j,k)];

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uT = Uthere[UTHERE(i,j,k)];

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Esquina Abajo Here
			j = 0;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)]);

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LUC(i,j+2,k,0)], EsquemaLargo);

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LVC(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LVC(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LVC(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)]);

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LUC(i,j,k+2,0)], EsquemaLargo);

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LWC(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LWC(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LWC(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LUC(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LUC(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LUC(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LUC(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

		}	
	}
}

//Cálculo del término convectivo de la velocidad V
void Solver::Get_ConvectiveV(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD){
int i, j, k;
double vW, vE, vS, vN, vH, vT, uW, uE, wH, wT;
double vW_pred, vE_pred, vS_pred, vN_pred, vH_pred, vT_pred, uW_pred, uE_pred, wH_pred, wT_pred;

	//ESQUEMA CONVECTIVO:
	//Corto (Interpolación)
	//(CoordObjetivo, Coord1, Valor1, Coord2, Valor2)


	//Largo
	//(CoordObjetivo, Velocidad, Coord1, Valor1, Coord2, Valor2, Coord3, Valor3, Coord4, Valor4, EsquemaLargo)

	if(Rank != 0 && Rank != Procesos - 1){

		for(i = Ix; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
					vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

					vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
					vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

					uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

					wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

					vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
					vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

					vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
					vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				}
			}

			//Partes Here y There
			//Parte Here
			k = 0;

			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

			
				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

				
				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				
			}

		}
		
	}	
	else if(Rank == 0){

		
		for(i = Ix + 1; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
					vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

					vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
					vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

					uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

					wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

					vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
					vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

					vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
					vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);


				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY; j++){
				
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}	

		}

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){

				vW_pred = Vleft[VLEFT(i,j,k)];
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = Vleft[VLEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
			}
		}
		
		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){

				

				vW_pred = Vleft[VLEFT(i,j,k)];
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = Vleft[VLEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

				
				

				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
		}
		
		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){
					
				vW_pred = Vleft[VLEFT(i,j,k)];
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = Vleft[VLEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);


				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

		}

	}
	else if(Rank == Procesos - 1){
		
		
		for(i = Ix; i < Fx - 1; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
					vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

					vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
					vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

					uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

					wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

					vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
					vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

					vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
					vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;

			for(j = 1; j < NY; j++){
				
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,0)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
			}

			//Parte There
			k = NZ - 1;

			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LVC(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LUC(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LUC(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LUC(i+1,j+1,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,0)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}
		}

		//Parte Derecha
		i = NX - 1;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Vright[VRIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = Vright[VRIGHT(i,j,k)];

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
			}
		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){
				
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Vright[VRIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = Vright[VRIGHT(i,j,k)];

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LVC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LWC(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LWC(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LWC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(  
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
		}


		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)]);
				vE_pred = Vright[VRIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LVC(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LVC(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LVC(i+1,j,k,0)], EsquemaLargo);
				vE = Vright[VRIGHT(i,j,k)];

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LVC(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LVC(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LVC(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LUC(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LUC(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LUC(i,j+1,k,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LVC(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*( 
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				
		}

	}

}

//Cálculo del término convectivo de la velocidad W
void Solver::Get_ConvectiveW(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD){
int i, j, k;
double wW, wE, wS, wN, wH, wT, uW, uE, vS, vN;
double wW_pred, wE_pred, wS_pred, wN_pred, wH_pred, wT_pred, uW_pred, uE_pred, vS_pred, vN_pred;

	//ESQUEMA CONVECTIVO:
	//Corto (Interpolación)
	//(CoordObjetivo, Coord1, Valor1, Coord2, Valor2)


	//Largo
	//(CoordObjetivo, Velocidad, Coord1, Valor1, Coord2, Valor2, Coord3, Valor3, Coord4, Valor4, EsquemaLargo)

	if(Rank != 0 && Rank != Procesos - 1){

		
		for(i = Ix; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
					wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

					wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
					wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

					uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

					wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
					wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

					wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
					wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior 
			j = 0;
			for(k = 1; k < NZ; k++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k-1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k-1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}

		}
	}	
	else if(Rank == 0){

		
		for(i = Ix + 1; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
					wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

					wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
					wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

					uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

					wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
					wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

					wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
					wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ; k++){
				
				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
			}


		}

		//Parte Izquierda
		i = 0;

		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){

				wW_pred = Wleft[WLEFT(i,j,k)];
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = Wleft[WLEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ; k++){
				
				

				wW_pred = Wleft[WLEFT(i,j,k)];
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = Wleft[WLEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){
				
				wW_pred = Wleft[WLEFT(,j,k)];
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = Wleft[WLEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
		}

	}
	else if (Rank == Procesos - 1){

		
		for(i = Ix; i < Fx - 1; i++){

			//Centro
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
					wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

					wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
					wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

					uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

					wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
					wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

					wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
					wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}


			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
			}


		}

		//Parte Derecha
		i = NX - 1;

		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LWC(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LUC(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LUC(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LUC(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LUC(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ; k++){
				
				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Wright[WRIGHT(i,j,k)];

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = Wright[WRIGHT(i,j,k)];

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LWC(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LVC(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LVC(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LVC(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LVC(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
		}


		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)]);
				wE_pred = Wright[WRIGHT(i,j,k)];

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LWC(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LWC(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LWC(i+1,j,k,0)], EsquemaLargo);
				wE = Wright[WRIGHT(i,j,k)];

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LWC(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LWC(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LWC(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LWC(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LWC(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LWC(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LUC(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LUC(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LUC(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LUC(i,j,k+1,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LVC(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LVC(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LVC(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LVC(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LWC(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

		}

	}

}

//Cálculo del término convectivo de la temperatura T
void Solver::Get_ConvectiveT(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD, double *TFIELD){
int i, j, k;
double Tw, Te, Ts, Tn, Th, Tt;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro 
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );
				}
			}
		}

		for(i = Ix; i < Fx; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

		}


	}
	else if(Rank == 0){

		//Centro 
		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );
				}
			}
		}

		for(i = Ix + 1; i < Fx; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){
	
				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

		}

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){

					Tw = TLEFT[PLEFT(0,j,k)];
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );
			}
		}

		//Partes Superior e Inferior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ - 1; k++){

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ - 1; k++){

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY - 1; j++){

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

		}

		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY - 1; j++){

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

		}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );


	}
	else if(Rank == Procesos - 1){

		//Centro 
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );
				}
			}
		}

		for(i = Ix; i < Fx - 1; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TFIELD[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*UFIELD[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

		}


		//Parte Derecha
		i = NX - 1;

		//Centro 
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );
			}
		}

		for(i = Ix; i < Fx - 1; i++){

			//Partes Superior e Inferior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TFIELD[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TFIELD[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WFIELD[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TFIELD[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TFIELD[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TFIELD[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TFIELD[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TFIELD[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TFIELD[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TFIELD[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TFIELD[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TFIELD[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TFIELD[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*UFIELD[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

		}

	}

}

//Cálculo de término de Boussinesq de la ecuación de Cantidad de Movimiento
void Solver::Get_BoussinesqTerm(Mesher MESH){
int i, j, k;
double Tv;

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 1; j < NY; j++){

				Tv = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
					
				BoussinesqV[LVC(i,j,k,0)] = gy*(1.0 - Beta*(Tv - To));

			}
		}
	}
	
}

//Ejecución de la integración por Runge - Kutta de las velocidades
void Solver::Get_RungeKuttaVelocities(Mesher MESH, ParPro MPI1){
int i, j, k;

	//CÁLCULO DE K1

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLU(ULPRES, ULPRES, Ix, Fx); //Enviar Velocidades U
	MPI1.CommunicateDataLV(VLPRES, VLPRES, Ix, Fx); //Enviar Velocidades V
	MPI1.CommunicateDataLW(WLPRES, WLPRES, Ix, Fx); //Enviar Velocidades W

	Get_DiffusiveU(MESH, ULPRES);
	Get_DiffusiveV(MESH, VLPRES);
	Get_DiffusiveW(MESH, WLPRES);

	Get_ConvectiveU(MESH, ULPRES, VLPRES, WLPRES);
	Get_ConvectiveV(MESH, ULPRES, VLPRES, WLPRES);
	Get_ConvectiveW(MESH, ULPRES, VLPRES, WLPRES);

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K1_U[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)];
			}
		}
	}

	//Velocidades V
	switch(Problema){
		case 1:
		case 3:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K1_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)];
				}
			}
		}

		break;

		
		case 2:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K1_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)] + BoussinesqV[LVC(i,j,k,0)];
				}
			}
		}

		break;
	}
	

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				K1_W[LWC(i,j,k,0)] = - ConvectiveW[LWC(i,j,k,0)] + DiffusiveW[LWC(i,j,k,0)];
			}
		}
	}


	//CÁLCULO DE K2

	//Cálculo de las nuevas velocidades intermedias

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				UL_NEW[LUC(i,j,k,0)] = ULPRES[LUC(i,j,k,0)] + c2*DeltaT*a_21*K1_U[LUC(i,j,k,0)];
			}
		}
	}

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				VL_NEW[LVC(i,j,k,0)] = VLPRES[LVC(i,j,k,0)] + c2*DeltaT*a_21*K1_V[LVC(i,j,k,0)];
			}
		}
	}

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				WL_NEW[LWC(i,j,k,0)] = WLPRES[LWC(i,j,k,0)] + c2*DeltaT*a_21*K1_W[LWC(i,j,k,0)];
			}
		}
	}


	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLU(UL_NEW, UL_NEW, Ix, Fx); //Enviar Velocidades U
	MPI1.CommunicateDataLV(VL_NEW, VL_NEW, Ix, Fx); //Enviar Velocidades V
	MPI1.CommunicateDataLW(WL_NEW, WL_NEW, Ix, Fx); //Enviar Velocidades W

	Get_DiffusiveU(MESH, UL_NEW);
	Get_DiffusiveV(MESH, VL_NEW);
	Get_DiffusiveW(MESH, WL_NEW);

	Get_ConvectiveU(MESH, UL_NEW, VL_NEW, WL_NEW);
	Get_ConvectiveV(MESH, UL_NEW, VL_NEW, WL_NEW);
	Get_ConvectiveW(MESH, UL_NEW, VL_NEW, WL_NEW);

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K2_U[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)];
			}
		}
	}

	//Velocidades V
	switch(Problema){
		case 1:
		case 3:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K2_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)];
				}
			}
		}

		break;

		case 2:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K2_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)] + BoussinesqV[LVC(i,j,k,0)];
				}
			}
		}
	
		break;
	}
	

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				K2_W[LWC(i,j,k,0)] = - ConvectiveW[LWC(i,j,k,0)] + DiffusiveW[LWC(i,j,k,0)];
			}
		}
	}


	//CÁLCULO DE K3

	//Cálculo de las nuevas velocidades intermedias

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				UL_NEW[LUC(i,j,k,0)] = ULPRES[LUC(i,j,k,0)] + c3*DeltaT*(a_31*K1_U[LUC(i,j,k,0)] + a_32*K2_U[LUC(i,j,k,0)]);
			}
		}
	}

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				VL_NEW[LVC(i,j,k,0)] = VLPRES[LVC(i,j,k,0)] + c3*DeltaT*(a_31*K1_V[LVC(i,j,k,0)] + a_32*K2_V[LVC(i,j,k,0)]);
			}
		}
	}

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				WL_NEW[LWC(i,j,k,0)] = WLPRES[LWC(i,j,k,0)] + c3*DeltaT*(a_31*K1_W[LWC(i,j,k,0)] + a_32*K2_W[LWC(i,j,k,0)]);
			}
		}
	}

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLU(UL_NEW, UL_NEW, Ix, Fx); //Enviar Velocidades U
	MPI1.CommunicateDataLV(VL_NEW, VL_NEW, Ix, Fx); //Enviar Velocidades V
	MPI1.CommunicateDataLW(WL_NEW, WL_NEW, Ix, Fx); //Enviar Velocidades W

	Get_DiffusiveU(MESH, UL_NEW);
	Get_DiffusiveV(MESH, VL_NEW);
	Get_DiffusiveW(MESH, WL_NEW);

	Get_ConvectiveU(MESH, UL_NEW, VL_NEW, WL_NEW);
	Get_ConvectiveV(MESH, UL_NEW, VL_NEW, WL_NEW);
	Get_ConvectiveW(MESH, UL_NEW, VL_NEW, WL_NEW);

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K3_U[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)];
			}
		}
	}

	//Velocidades V
	switch(Problema){
		case 1:
		case 3:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K3_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)];
				}
			}
		}
		break;

		case 2:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K3_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)] + BoussinesqV[LVC(i,j,k,0)];
				}
			}
		}

		break;
	}
	

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				K3_W[LWC(i,j,k,0)] = - ConvectiveW[LWC(i,j,k,0)] + DiffusiveW[LWC(i,j,k,0)];
			}
		}
	}


	//CÁLCULO DE K4

	//Cálculo de las nuevas velocidades intermedias

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				UL_NEW[LUC(i,j,k,0)] = ULPRES[LUC(i,j,k,0)] + c4*DeltaT*(a_41*K1_U[LUC(i,j,k,0)] + a_42*K2_U[LUC(i,j,k,0)] + a_43*K3_U[LUC(i,j,k,0)]);
			}
		}
	}

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				VL_NEW[LVC(i,j,k,0)] = VLPRES[LVC(i,j,k,0)] + c4*DeltaT*(a_41*K1_V[LVC(i,j,k,0)] + a_42*K2_V[LVC(i,j,k,0)] + a_43*K3_V[LVC(i,j,k,0)]);
			}
		}
	}

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				WL_NEW[LWC(i,j,k,0)] = WLPRES[LWC(i,j,k,0)] + c4*DeltaT*(a_41*K1_W[LWC(i,j,k,0)] + a_42*K2_W[LWC(i,j,k,0)] + a_43*K3_W[LWC(i,j,k,0)]);
			}
		}
	}

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLU(UL_NEW, UL_NEW, Ix, Fx); //Enviar Velocidades U
	MPI1.CommunicateDataLV(VL_NEW, VL_NEW, Ix, Fx); //Enviar Velocidades V
	MPI1.CommunicateDataLW(WL_NEW, WL_NEW, Ix, Fx); //Enviar Velocidades W

	Get_DiffusiveU(MESH, UL_NEW);
	Get_DiffusiveV(MESH, VL_NEW);
	Get_DiffusiveW(MESH, WL_NEW);

	Get_ConvectiveU(MESH, UL_NEW, VL_NEW, WL_NEW);
	Get_ConvectiveV(MESH, UL_NEW, VL_NEW, WL_NEW);
	Get_ConvectiveW(MESH, UL_NEW, VL_NEW, WL_NEW);

	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K4_U[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)];
			}
		}
	}

	//Velocidades V
	switch(Problema){
		case 1:
		case 3:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K4_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)];
				}
			}
		}

		break;
		
		case 2:

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY + 1; j++){
					K4_V[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)] + BoussinesqV[LVC(i,j,k,0)];
				}
			}
		}

		break;
	}
	

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				K4_W[LWC(i,j,k,0)] = - ConvectiveW[LWC(i,j,k,0)] + DiffusiveW[LWC(i,j,k,0)];
			}
		}
	}

	//CÁLCULO DE LAS VELOCIDADES PREDICTORAS CON RUNGE KUTTA
	
	//Velocidades U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				PredU[LUC(i,j,k,0)] = ULPRES[LUC(i,j,k,0)] + DeltaT*(b1*K1_U[LUC(i,j,k,0)] + b2*K2_U[LUC(i,j,k,0)] + b3*K3_U[LUC(i,j,k,0)] + b4*K4_U[LUC(i,j,k,0)]);
			}
		}
	}

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				PredV[LVC(i,j,k,0)] = VLPRES[LVC(i,j,k,0)] + DeltaT*(b1*K1_V[LVC(i,j,k,0)] + b2*K2_V[LVC(i,j,k,0)] + b3*K3_V[LVC(i,j,k,0)] + b4*K4_V[LVC(i,j,k,0)]);
			}
		}
	}

	//Velocidades W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				PredW[LWC(i,j,k,0)] = WLPRES[LWC(i,j,k,0)] + DeltaT*(b1*K1_W[LWC(i,j,k,0)] + b2*K2_W[LWC(i,j,k,0)] + b3*K3_W[LWC(i,j,k,0)] + b4*K4_W[LWC(i,j,k,0)]);
			}
		}
	}

}

//Ejecución de la integración por Runge - Kutta de las temperaturas
void Solver::Get_RungeKuttaTemperature(Mesher MESH, ParPro MPI1){
int i, j, k;

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLP(TLPRES, TLPRES, Ix, Fx); //Comunicación Temperaturas

	Get_DiffusiveT(MESH, TLPRES);

	Get_ConvectiveT(MESH, ULFUT, VLFUT, WLFUT, TLPRES);
	
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K1_T[LPC(i,j,k,0)] = - ConvectiveT[LPC(i,j,k,0)] + DiffusiveT[LPC(i,j,k,0)];
			}
		}
	}

	//CÁLCULO DE K2

	//Cálculo de las nuevas temperaturas intermedias

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				TL_NEW[LPC(i,j,k,0)] = TLPRES[LPC(i,j,k,0)] + c2*DeltaT*a_21*K1_T[LPC(i,j,k,0)];
			}
		}
	}

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLP(TL_NEW, TL_NEW, Ix, Fx); //Comunicación Temperaturas

	Get_DiffusiveT(MESH, TL_NEW);

	Get_ConvectiveT(MESH, ULFUT, VLFUT, WLFUT, TL_NEW);
	
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K2_T[LPC(i,j,k,0)] = - ConvectiveT[LPC(i,j,k,0)] + DiffusiveT[LPC(i,j,k,0)];
			}
		}
	}

	//CÁLCULO DE K3

	//Cálculo de las nuevas temperaturas intermedias

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				TL_NEW[LPC(i,j,k,0)] = TLPRES[LPC(i,j,k,0)] + c3*DeltaT*(a_31*K1_T[LPC(i,j,k,0)] + a_32*K2_T[LPC(i,j,k,0)]);
			}
		}
	}

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLP(TL_NEW, TL_NEW, Ix, Fx); //Comunicación Temperaturas

	Get_DiffusiveT(MESH, TL_NEW);

	Get_ConvectiveT(MESH, ULFUT, VLFUT, WLFUT, TL_NEW);

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K3_T[LPC(i,j,k,0)] = - ConvectiveT[LPC(i,j,k,0)] + DiffusiveT[LPC(i,j,k,0)];
			}
		}
	}

	//CÁLCULO DE K4

	//Cálculo de las nuevas temperaturas intermedias

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				TL_NEW[LPC(i,j,k,0)] = TLPRES[LPC(i,j,k,0)] + c4*DeltaT*(a_41*K1_T[LPC(i,j,k,0)] + a_42*K2_T[LPC(i,j,k,0)] + a_43*K3_T[LPC(i,j,k,0)]);
			}
		}
	}

	//Comunicación de velocidades entre los procesos
	MPI1.CommunicateDataLP(TL_NEW, TL_NEW, Ix, Fx); //Comunicación Temperaturas

	Get_DiffusiveT(MESH, TL_NEW);

	Get_ConvectiveT(MESH, ULFUT, VLFUT, WLFUT, TL_NEW);

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				K4_T[LPC(i,j,k,0)] = - ConvectiveT[LPC(i,j,k,0)] + DiffusiveT[LPC(i,j,k,0)];
			}
		}
	}


	//Cálculo de la temperatura final
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				TLFUT[LPC(i,j,k,0)] = TLPRES[LPC(i,j,k,0)] + DeltaT*(b1*K1_T[LPC(i,j,k,0)] + b2*K2_T[LPC(i,j,k,0)] + b3*K3_T[LPC(i,j,k,0)] + b4*K4_T[LPC(i,j,k,0)]);
			}
		}
	}

}

//Cálculo de la divergencia de la velocidad predictora (Término bp)
void Solver::Get_PredictorsDivergence(Mesher MESH, ParPro MPI1){
int i, j, k;
//(Rho/(DeltaT*MESH.VolMP[GP(i,j,k,0)]))
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				bp[LAL(i,j,k,0)] = - (1.0/(DeltaT))*(
								   - MESH.SupMP[GP(i,j,k,0)]*PredU[LUC(i,j,k,0)]
								   + MESH.SupMP[GP(i,j,k,1)]*PredU[LUC(i + 1,j,k,0)]
								   - MESH.SupMP[GP(i,j,k,2)]*PredV[LVC(i,j,k,0)]
								   + MESH.SupMP[GP(i,j,k,3)]*PredV[LVC(i,j + 1,k,0)]
								   - MESH.SupMP[GP(i,j,k,4)]*PredW[LWC(i,j,k,0)]
								   + MESH.SupMP[GP(i,j,k,5)]*PredW[LWC(i,j,k + 1,0)]
								   );
			}
		}
	}

	//MPI1.SendMatrixToZeroBP(bp, bpGlobal, NX, NY, NZ, Procesos, Ix, Fx);
}

//Cálculo de la máxima diferencia en Gauss-Seidel
void Solver::Get_MaxDifGS(ParPro MPI1){
int i, j, k;
MaxDiffGS = 0.0;
MPI_Status ST;

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				if(abs(PLFUT[LPC(i,j,k,0)] - PLSUP[LPC(i,j,k,0)]) >= MaxDiffGS){
					MaxDiffGS = abs(PLFUT[LPC(i,j,k,0)] - PLSUP[LPC(i,j,k,0)]);
				}
			}
		}
	}

	MPI_Allreduce(&MaxDiffGS, &MaxDiffGS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

//Resolución del sistema con Gauss-Seidel
void Solver::Get_GaussSeidel(ParPro MPI1){
int i, j, k;
MaxDiffGS = 2.0*ConvergenciaGS;

	while(MaxDiffGS >= ConvergenciaGS){

		if(Rank != 0 && Rank != Procesos - 1){

			//Parte Central
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
					}
				}
			}
			
			//Parte Inferior
			j = 0;
			for(i = Ix; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Superior
			j = NY - 1;
			for(i = Ix; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Here
			k = 0;
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
						PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte There
			k = NZ - 1;
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Esquinas
			for(i = Ix; i < Fx; i++){

				//Esquina Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
			}

		}
        else if(Rank == 0){

        	//Parte Central
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
					}
				}
			}
			
			//Parte Inferior
			j = 0;
			for(i = Ix + 1; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Superior
			j = NY - 1;
			for(i = Ix + 1; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Here
			k = 0;
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
						PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte There
			k = NZ - 1;
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Esquinas
			for(i = Ix + 1; i < Fx; i++){

				//Esquina Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
			}
			

			//Parte Izquierda
			i = 0;

			//Centro
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Esquina Abajo Here
			j = 0;
			k = 0;
			PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;
			PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
			//Esquina Abajo There
			j = 0;
			k = NZ - 1;
			PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;
			PLFUT[LPC(i,j,k,0)] = (ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
		}
		else if(Rank == Procesos - 1){

			//Parte Central
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
					}
				}
			}
			
			//Parte Inferior
			j = 0;
			for(i = Ix; i < Fx - 1; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Superior
			j = NY - 1;
			for(i = Ix; i < Fx - 1; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Here
			k = 0;
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
						PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte There
			k = NZ - 1;
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Esquinas
			for(i = Ix; i < Fx - 1; i++){

				//Esquina Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
				//Esquina Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LAL(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				
			}

			//Parte Derecha
			i = NX - 1;

			//Centro
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
				}
			}

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,NZ - 1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LPC(i,j,k,0)] = (aw[LAL(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LAL(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LAL(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LAL(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LAL(i,j,k,0)]*PLFUT[LPC(i,j,0,0)] + bp[LAL(i,j,k,0)])/ap[LAL(i,j,k,0)];
			}

		}

		Get_MaxDifGS(MPI1);

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					PLSUP[LPC(i,j,k,0)] = PLFUT[LPC(i,j,k,0)];
				}
			}
		}
		
		MPI1.CommunicateDataLP(PLFUT, PLFUT, Ix, Fx);
	}

}

//Cálculo de las velocidades del sistema
void Solver::Get_Velocities(Mesher MESH, ParPro MPI1){
int i, j, k;

	//Comunicación de las nuevas presiones
	MPI1.CommunicateDataLP(PLFUT, PLFUT, Ix, Fx); //Comunicar Presiones P

	//Velocidades U
	if(Rank != 0 && Rank != Procesos - 1){

		for(i = Ix; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULFUT[LUC(i,j,k,0)] = PredU[LUC(i,j,k,0)] - (DeltaT/Rho)*((PLFUT[LPC(i,j,k,0)] - PLFUT[LPC(i - 1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]);
				}
			}
		}

	}
	else if(Rank == 0){

		for(i = Ix + 1; i < Fx + 1; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULFUT[LUL(i,j,k,0)] = PredU[LUL(i,j,k,0)] - (DeltaT/Rho)*((PLFUT[LPC(i,j,k,0)] - PLFUT[LPC(i - 1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]);
				}
			}
		}

		//Parte Izquierda
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				ULFUT[LUL(0,j,k,0)] = Uleft[UHERE(0,j,k)];
			}
		}

	}
	else if(Rank == Procesos - 1){

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULFUT[LUR(i,j,k,0)] = PredU[LUR(i,j,k,0)] - (DeltaT/Rho)*((PLFUT[LPC(i,j,k,0)] - PLFUT[LPC(i - 1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]);
				}
			}
		}

		//Parte Derecha
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				ULFUT[LUR(NX,j,k,0)] = Uright[URIGHT(0,j,k)];
			}
		}
			
	}

	//Velocidad V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){

			//Parte Abajo
			VLFUT[LVC(i,0,k,0)] = Vbot[VBOT(i,0,k)];

			//Parte Arriba
			VLFUT[LVC(i,NY,k,0)] = Vtop[VTOP(i,NY,k)];
			
			for(j = 1; j < NY; j++){
				VLFUT[LVC(i,j,k,0)] = PredV[LVC(i,j,k,0)] - (DeltaT/Rho)*((PLFUT[LPC(i,j,k,0)] - PLFUT[LPC(i,j - 1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]);
			}
		}
	}

	//Velocidad W
	for(i = Ix; i < Fx; i++){
		for(k = 1; k < NZ; k++){
			for(j = 0; j < NY; j++){
				WLFUT[LWC(i,j,k,0)] = PredW[LWC(i,j,k,0)] - (DeltaT/Rho)*((PLFUT[LPC(i,j,k,0)] - PLFUT[LPC(i,j,k - 1,0)])/MESH.DeltasMW[GW(i,j,k,2)]);
			}
		}
	}

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){

			//Parte Here
			WLFUT[LWC(i,j,0,0)] = Where[WHERE(i,j,0)];

			//Parte There
			WLFUT[LWC(i,j,NZ,0)] = Wthere[WTHERE(i,j,NZ)];

		}
	}
		
}

//Calcular la diferencia entre steps para comprobar la convergencia
void Solver::Get_Stop(ParPro MPI1){
int i, j, k;
MaxDiffGlobal = 0.0;

	//Comprobación Velocidad U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				MaxDiffGlobal += (abs((ULFUT[LUC(i,j,k,0)] - ULPRES[LUC(i,j,k,0)])/(ULPRES[LUC(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((ULFUT[LUC(i,j,k,0)] - ULPRES[LUC(i,j,k,0)])/(ULPRES[LUC(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	//Comprobación Velocidad V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				MaxDiffGlobal += (abs((VLFUT[LVC(i,j,k,0)] - VLPRES[LVC(i,j,k,0)])/(VLPRES[LVC(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((VLFUT[LVC(i,j,k,0)] - VLPRES[LVC(i,j,k,0)])/(VLPRES[LVC(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	//Comprobación Velocidad W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				MaxDiffGlobal += (abs((WLFUT[LWC(i,j,k,0)] - WLPRES[LWC(i,j,k,0)])/(WLPRES[LWC(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((WLFUT[LWC(i,j,k,0)] - WLPRES[LWC(i,j,k,0)])/(WLPRES[LWC(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	MPI_Allreduce(&MaxDiffGlobal, &MaxDiffGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

//Actualizar los valores de todas las matrices de propiedades
void Solver::Get_Update(){
int i, j, k;

	//Actualización Velocidad U
	for(i = Ix; i < Fx + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				ULPRES[LUC(i,j,k,0)] = ULFUT[LUC(i,j,k,0)];
			}
		}
	}

	//Actualización Velocidad V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				VLPRES[LVC(i,j,k,0)] = VLFUT[LVC(i,j,k,0)];
			}
		}
	}

	//Actualización Velocidad W
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ + 1; k++){
			for(j = 0; j < NY; j++){
				WLPRES[LWC(i,j,k,0)] = WLFUT[LWC(i,j,k,0)];
			}
		}
	}

	switch(Problema){
		case 2:

		//Actualización Temperaturas
		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					TLPRES[LPC(i,j,k,0)] = TLFUT[LPC(i,j,k,0)];
				}
			}
		}

	}
}

//Ejecución del solver del problema
void Solver::ExecuteSolver(Memory M1, ParPro MPI1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	auto Inicio_Total = std::chrono::high_resolution_clock::now();
	double Time = 0.0;
	int Step = 0;

	char FileName_1[300];
	char FileName_2[300];
	char FileName_3[300];

	MaxDiffGlobal = 2.0*ConvergenciaGlobal;

	if(Rank == 0){
		cout<<"Generando memoria para el solver... "<<endl;
	}

	auto Inicio_Preprograma = std::chrono::high_resolution_clock::now();
	
	AllocateMatrix(M1); //Alojamiento de Memoria

	if(Rank == 0){
		cout<<"Memoria generada "<<endl;
	}

	Get_InitialConditions(); //Condiciones iniciales de las matrices
	
	if(Rank == 0){
		//Generar los coeficientes de la matriz laplaciana
		std::cout<<"Generando los coeficientes de la matriz laplaciana...\n";
	}
		
	Get_PressureCoefficients(MESH); //Cálculo de los coeficientes de presión
	
	if(Rank == 0){
		//Generar los coeficientes de la matriz laplaciana
		std::cout<<"Coeficientes generados \n";	
	}

	auto Fin_Preprograma = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> Tiempo_Preprograma = Fin_Preprograma - Inicio_Preprograma;

	std::chrono::duration<double> Duracion_InStep, Tiempo_InStep, Duracion_RungeKutta, Tiempo_RungeKutta;
	std::chrono::duration<double> Duracion_Predictors, Tiempo_Predictors, Duracion_GS, Tiempo_GS;
	std::chrono::duration<double> Duracion_Stop, Tiempo_Stop, Duracion_Post, Tiempo_Post;
	
	while(MaxDiffGlobal >= ConvergenciaGlobal){
		
		auto Inicio_InStep = std::chrono::high_resolution_clock::now();
		Step++;
		Get_BoundaryConditions(MESH); //Actualizar condiciones de contorno
		Get_HaloVelocities(); //Actualizar el Halo de las velocidades
		Get_StepTime(MESH, MPI1); //Get Time Step
		Time += DeltaT;
		auto Fin_InStep = std::chrono::high_resolution_clock::now();

		Duracion_InStep = (Fin_InStep - Inicio_InStep);
		Tiempo_InStep += Duracion_InStep;

		auto Inicio_RungeKutta = std::chrono::high_resolution_clock::now();
		switch(Problema){
			case 1: 
			case 3:

			Get_RungeKuttaVelocities(MESH, MPI1);
			break;

			case 2:
			Get_RungeKuttaVelocities(MESH, MPI1);
			Get_BoussinesqTerm(MESH);
			Get_RungeKuttaVelocities(MESH, MPI1);
			break;
		}
		auto Fin_RungeKutta = std::chrono::high_resolution_clock::now();
		Duracion_RungeKutta = (Fin_RungeKutta - Inicio_RungeKutta);
		Tiempo_RungeKutta += Duracion_RungeKutta;

		auto Inicio_Predictors = std::chrono::high_resolution_clock::now();
		Get_PredictorsDivergence(MESH, MPI1); //Calculo bp
		auto Fin_Predictors = std::chrono::high_resolution_clock::now();
		Duracion_Predictors = (Fin_Predictors - Inicio_Predictors);
		Tiempo_Predictors += Duracion_Predictors;

		auto Inicio_GS = std::chrono::high_resolution_clock::now();
		Get_GaussSeidel(MPI1); //Resolucion por GS
		auto Fin_GS = std::chrono::high_resolution_clock::now();
		Duracion_GS = (Fin_GS - Inicio_GS);
		Tiempo_GS += Duracion_GS;


		Get_Velocities(MESH, MPI1); //Calculo de las velocidades

		switch(Problema){
			case 2:
			Get_RungeKuttaTemperature(MESH, MPI1);
			break;
		}
	
		auto Inicio_Stop = std::chrono::high_resolution_clock::now();
		if(Step%100 == 0){
			Get_Stop(MPI1);		
			if(Rank == 0){
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<endl;
			}
		}
		auto Fin_Stop = std::chrono::high_resolution_clock::now();
		Duracion_Stop = (Fin_Stop - Inicio_Stop);
		Tiempo_Stop += Duracion_Stop;

		auto Inicio_Post = std::chrono::high_resolution_clock::now();
		if(Step%500 == 0){

			//Envío Matrices Locales a Globales Presente
			MPI1.SendMatrixToZeroMU(ULFUT, UGFUT, NX, NY, NZ, Procesos, Ix, Fx);
			MPI1.SendMatrixToZeroMV(VLFUT, VGFUT, NX, NY, NZ, Procesos, Ix, Fx);
			MPI1.SendMatrixToZeroMW(WLFUT, WGFUT, NX, NY, NZ, Procesos, Ix, Fx);

			if(Rank == 0){

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				sprintf(FileName_2, "MapaVelocidades_Step_%d", Step);

				switch(Problema){

					case 1: 
					POST1.EscalarVTK3D("DrivenCavity/", "Presion", FileName_1, PGPRES, MESH.MP, NX, NY, NZ);
					POST1.VectorialVTK3D(MESH, "DrivenCavity/", "Velocidad", FileName_2, UGFUT, VGFUT, WGFUT, MESH.MP, NX, NY, NZ);
					//POST1.Get_DrivenResults(MESH, Time, Procesos, Step, UGFUT, VGFUT);
					break;

					case 2:
					POST1.EscalarVTK3D("DifferentiallyHeated/", "Presion", FileName_1, PGPRES, MESH.MP, NX, NY, NZ);
					POST1.VectorialVTK3D(MESH, "DifferentiallyHeated/", "Velocidad", FileName_2, UGFUT, VGFUT, WGFUT, MESH.MP, NX, NY, NZ);
					sprintf(FileName_3, "MapaTemperaturas_Step_%d", Step);
					POST1.EscalarVTK3D("DifferentiallyHeated/", "Temperatura", FileName_3, TGFUT, MESH.MP, NX, NY, NZ);
					POST1.Get_DifferentiallyResults(MESH, Time, Procesos, Step, UGFUT, VGFUT, TGFUT);
					break;

					case 3:
					POST1.EscalarVTK3D("SquareCylinder/", "Presion", FileName_1, PGPRES, MESH.MP, NX, NY, NZ);
					POST1.VectorialVTK3D(MESH, "SquareCylinder/", "Velocidad", FileName_2, UGFUT, VGFUT, WGFUT, MESH.MP, NX, NY, NZ);
					
					break;

				}
				
			}
		}
		
		auto Fin_Post = std::chrono::high_resolution_clock::now();
		Duracion_Post = (Fin_Post - Inicio_Post);
		Tiempo_Post += Duracion_Post;

		Get_Update();
	
	}
	
	auto Fin_Total = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> Tiempo_Total = Fin_Total - Inicio_Total;

	if(Rank == 0){
		std::cout << "Tiempo_Total: " << Tiempo_Total.count() << " s\n";
		std::cout << "Tiempo_Preprograma: " << Tiempo_Preprograma.count() << " s, "<<100*(Tiempo_Preprograma.count()/Tiempo_Total.count())<<" %\n";
		std::cout << "Tiempo_InStep: " << Tiempo_InStep.count() << " s, "<<100*(Tiempo_InStep.count()/Tiempo_Total.count())<<" %\n";
		std::cout << "Tiempo_RungeKutta: " << Tiempo_RungeKutta.count() << " s, "<<100*(Tiempo_RungeKutta.count()/Tiempo_Total.count())<<" %\n";
		std::cout << "Tiempo_Predictors: " << Tiempo_Predictors.count() << " s, "<<100*(Tiempo_Predictors.count()/Tiempo_Total.count())<<" %\n";
		std::cout << "Tiempo_GS: " << Tiempo_GS.count() << " s, "<<100*(Tiempo_GS.count()/Tiempo_Total.count())<<" %\n";
		std::cout << "Tiempo_Stop: " << Tiempo_Stop.count() << " s, "<<100*(Tiempo_Stop.count()/Tiempo_Total.count())<<" %\n";
		std::cout << "Tiempo_Post: " << Tiempo_Post.count() << " s, "<<100*(Tiempo_Post.count()/Tiempo_Total.count())<<" %\n";
	}
	
	if(Rank == 0){
		cout<<"Solver Completed"<<endl;
	}

}


			