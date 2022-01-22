#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ReadData.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ParPro.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Geometry.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Mesher.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Solver.h"

#define PI 3.141592653589793

#define G(i, j, dim) (((j) + (i)*NR) + NA * NR * (dim))				   //Global Index
#define GU(i, j, dim) (((j) + (i)*NR) + (NA + 1) * NR * (dim))		   //Global Index Matriz U
#define GR(i, j, dim) (((j) + (i) * (NR + 1)) + NA * (NR + 1) * (dim)) //Global Index Matriz R

#define MU(i,j,dim) ((j) + NR*(dim)) 
#define MR(i,j,dim) ((i) + (Fx - Ix)*(dim) + (Fx - Ix)*(4)*(j))

#define VR(i,j,dim) ((i) + (Fx - Ix)*(j))

#define LU(i, j, dim) (((j) + ((i)-Ix) * NR)) //Local Index Axial Velocity Nodes
#define LR(i, j, dim) (((j) + ((i)-Ix) * (NR + 1))) //Local Index Radial Velocity Nodes

#define LNH(i, j, dim) (((j) + ((i)-Ix + Halo) * NR))								 //Local Index No Halo
#define LSH(i, j, dim) (((j) + ((i)-Ix) * NY) + NY * (Fx - Ix + 2 * Halo) * dim) //Local Index Si Halo

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/"

//Constructor del mallador
Solver::Solver(Memory M1, ReadData R1, ParPro MPI1, Mesher MESH){
		
	//Datos del problema
	NA = MESH.ReturnNA();
	NR = MESH.ReturnNR();

	TypeMesh = R1.Get_MeshingType(); //Tipo de mallado seleccionado (Collocated/Staggered)
	EsquemaAmplio = "CDS";
	Problema = R1.Get_ProblemType(); //Problema (Tobera/Tubería)

	//Datos y variables de computación paralela
	Rank = MESH.Rank;
	Procesos = MPI1.Get_Processes();
	Ix = MESH.Ix;
	Fx = MESH.Fx;
	Halo = MPI1.Get_Halo();

	BulkRe = R1.ProblemPhysicalData[0];
	BulkM = R1.ProblemPhysicalData[1];
	Pr = R1.ProblemPhysicalData[2];
	Twall = R1.ProblemPhysicalData[3];
	Gamma = R1.ProblemPhysicalData[4];
	Rideal = R1.ProblemPhysicalData[5];
	CpBase = R1.ProblemPhysicalData[6];

	PipeLength = R1.GeometryData[0];   //Longitud de la tubería
	PipeDiameter = R1.GeometryData[1]; //Diametro de la tubería

	H = 0.50*R1.GeometryData[15];

	//Datos sobre el modelo de viscosidad
//	B = 1.0;
//	n = 0.7;
//	muW = B * pow(Twall, n);

	//muW(276 K) = 1.754e-5;
	//muW(500 K) = 2.760e-5;

	muW = 1.754e-5;
	S = 110.4;

	SoundSpeed = sqrt(Gamma*Rideal*Twall);
	RhoM = (muW*BulkRe)/(BulkM*H*SoundSpeed);
	
	TimeBetta = 0.5; //Adams-Bashforth
	Convergencia = 1e-6;

}

//Alojamiento de memoria para las matrices necesarias
void Solver::AllocateMatrix(Memory M1){

	//Matrices globales de propiedades
	if(Rank == 0){

		//Step Futuro
		UglobalFut = M1.AllocateDouble(NA, NR, 1);
		VglobalFut = M1.AllocateDouble(NA, NR, 1);
		TglobalFut = M1.AllocateDouble(NA, NR, 1);
		RhoGlobalFut = M1.AllocateDouble(NA, NR, 1);
		PresionGlobal = M1.AllocateDouble(NA, NR, 1);

		TauRZglobal = M1.AllocateDouble(NA, NR, 1);
		UwallsMUglobal = M1.AllocateDouble(NA+1, NR, 1);

		//Step Presente
		UglobalPres = M1.AllocateDouble(NA, NR, 1);
		VglobalPres = M1.AllocateDouble(NA, NR, 1);
		TglobalPres = M1.AllocateDouble(NA, NR, 1);
		RhoGlobalPres = M1.AllocateDouble(NA, NR, 1);

		PDT = M1.AllocateDouble(Procesos, 1, 1);
	}

	//Matrices de propiedades locales
	//Densidad
	RhoLocalPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	RhoLocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	RhoLocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Presión
	Pressure = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Velocidad Axial
	UlocalPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	UlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	UlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Velocidad Radial
	VlocalPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	VlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	VlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Energía Interna
	ElocalPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	ElocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	ElocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Temperatura
	TlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	TlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Matrices de propiedades en las paredes de los volúmenes de control
	//Densidad
	RhoWallsMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	RhoWallsMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Presión
	PwallsMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	PwallsMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Velocidad Axial
	UwallsMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	UwallsMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Velocidad Radial
	VwallsMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	VwallsMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Conductividad Térmica
	KwallsMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	KwallsMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Energía Interna
	EwallsMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	EwallsMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Viscosidad Dinámica
	muTotalMU = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	muTotalMR = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Matrices de cálculo de contribuciones a las propiedades
	//Densidad
	FRhoPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	FRhoPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Velocidad Axial
	FUPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	FUPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Velocidad Radial
	FVPrev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	FVPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Energía Interna
	FEprev = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	FEpres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	
	//Términos de las ecuaciones de Navier-Stokes
	//Momentum Axial
	MomentumConvectiveAxial = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	MomentumDiffusiveAxial = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	PressureGradientAxial = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Momentum Radial
	MomentumConvectiveRadial = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	MomentumDiffusiveRadial = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	PressureGradientRadial = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	//Energía
	EnergyConvective = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1); //Término convectivo de la ecuación de energía
	EnergyDiffusive = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);	 //Término difusivo de la ecuación de conservación de la energía
	EnergyPressureTerm = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1); //Término de presión de la ecuación de la energía
	EnergyViscous = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);	   //Término de disipación viscosa en la ecuación de energía

	//Matrices de Esfuerzos Viscosos
	Divergence = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	DivergenceMR = M1.AllocateDouble(Fx - Ix, 2, 1);
	if(Rank == 0 || Rank == Procesos - 1){
		DivergenceMU = M1.AllocateDouble(1, NR, 1);
	}
	TauRR = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	TauRZ = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	TauZZ = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	TauThetaTheta = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);

	TauRR_mu = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	TauRR_mr = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	TauRZ_mu = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	TauRZ_mr = M1.AllocateDouble(Fx - Ix, NR + 1, 1);
		
	TauZZ_mu = M1.AllocateDouble(Fx - Ix + 1, NR, 1);
	TauZZ_mr = M1.AllocateDouble(Fx - Ix, NR + 1, 1);

	//Matrices locales de propiedades físicas del fluido
	//Propiedades térmicas
	Cp = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1); //Calor específico de los gases en cada nodo
	K = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);  //Conductividad térmica en los nodos

	//Viscosidad
	muBase = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	muTurb = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1);
	muTotal = M1.AllocateDouble(Fx - Ix + 2 * Halo, NR, 1); //Viscosidad dinámica total en cada uno de los nodos

	//Matrices y arrays de condiciones de contorno	

	if(Rank == 0){
		RhoLeft = M1.AllocateDouble(NR, 1, 1); //Densidad pared izquierda
		Uleft = M1.AllocateDouble(NR, 1, 1);   //Velocidad axial pared izquierda
		Vleft = M1.AllocateDouble(NR, 1, 1);   //Velocidad radial pared izquierda
		Tleft = M1.AllocateDouble(NR, 1, 1);   //Temperaturas pared izquierda
		Pleft = M1.AllocateDouble(NR, 1, 1);   //Presión pared izquierda
		muLeft = M1.AllocateDouble(NR, 1, 1);  //Viscosidad dinámica pared izquierda
		KLeft = M1.AllocateDouble(NR, 1, 1);   //Conductividad térmica pared izquierda
		Eleft = M1.AllocateDouble(NR, 1, 1);   //Energía Interna pared izquierda
	}

	if(Rank == Procesos - 1){
		RhoRight = M1.AllocateDouble(NR, 1, 1); //Densidad pared derecha
		Uright = M1.AllocateDouble(NR, 1, 1);	//Velocidad axial pared derecha
		Vright = M1.AllocateDouble(NR, 1, 1);	//Velocidad radial pared derecha
		Tright = M1.AllocateDouble(NR, 1, 1);	//Temperaturas pared derecha
		Pright = M1.AllocateDouble(NR, 1, 1);	//Presión axial pared derecha
		muRight = M1.AllocateDouble(NR, 1, 1);	//Viscosidad dinámica pared derecha
		KRight = M1.AllocateDouble(NR, 1, 1);	//Conductividad térmica pared derecha
		Eright = M1.AllocateDouble(NR, 1, 1);		//Energía Interna pared derecha
	}

	//Variables parte de arriba
	RhoUp = M1.AllocateDouble(Fx - Ix, 1, 1);   //Densidad pared superior
	Uup = M1.AllocateDouble(Fx - Ix, 1, 1);		//Velocidad axial pared superior
	Vup = M1.AllocateDouble(Fx - Ix, 1, 1);		//Velocidad radial pared superior
	Tup = M1.AllocateDouble(Fx - Ix, 1, 1);		//Temperaturas pared superior
	Pup = M1.AllocateDouble(Fx - Ix, 1, 1);		//Presión axial pared superior
	muUp = M1.AllocateDouble(Fx - Ix, 1, 1);	//Viscosidad dinámica pared superior
	KUp = M1.AllocateDouble(Fx - Ix, 1, 1);		//Conductividad térmica pared superior
	Eup = M1.AllocateDouble(Fx - Ix, 1, 1);		//Energía Interna pared superior

	//Variables parte de abajo
	RhoDown = M1.AllocateDouble(Fx - Ix, 1, 1);   //Densidad pared inferior
	Udown = M1.AllocateDouble(Fx - Ix, 1, 1);		//Velocidad axial pared inferior
	Vdown = M1.AllocateDouble(Fx - Ix, 1, 1);		//Velocidad radial pared inferior
	Tdown = M1.AllocateDouble(Fx - Ix, 1, 1);		//Temperaturas pared inferior
	Pdown = M1.AllocateDouble(Fx - Ix, 1, 1);		//Presión axial pared inferior
	muDown = M1.AllocateDouble(Fx - Ix, 1, 1);	//Viscosidad dinámica pared inferior
	Kdown = M1.AllocateDouble(Fx - Ix, 1, 1);		//Conductividad térmica pared inferior
	Edown = M1.AllocateDouble(Fx - Ix, 1, 1);		//Energía Interna pared inferior

	//Código;
	//1 -> West (aW)
	//2 -> East (aE)
	//3 -> South (aS)
	//4 -> North (aN)
	//5 -> Central (aP)
	
}

//Inicialización de los campos de temperaturas
void Solver::InitializeFields(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){

			//Temperatura
			TlocalPres[LNH(i, j, 0)] = 0.80*Twall;

			//Calor específico
			Cp[LNH(i, j, 0)] = CpBase;

			//Densidad
			RhoLocalPrev[LNH(i, j, 0)] = 0.5*RhoM;
			RhoLocalPres[LNH(i, j, 0)] = 0.5*RhoM;

			//Presión
			Pressure[LNH(i, j, 0)] = RhoLocalPres[LNH(i, j, 0)] * TlocalPres[LNH(i, j, 0)] * Rideal;

			//Velocidad Axial
			UlocalPrev[LNH(i, j, 0)] = ((3.0*muW*BulkRe)/(2.0*RhoM*pow(H,3.0)))*(H - (MESH.MP[G(0,j,1)] - H))*(H + (MESH.MP[G(0,j,1)] - H));
			//*(H - (MESH.MP[G(0,j,1)] - H))*(H + (MESH.MP[G(0,j,1)] - H));
			UlocalPres[LNH(i, j, 0)] = ((3.0*muW*BulkRe)/(2.0*RhoM*pow(H,3.0)))*(H - (MESH.MP[G(0,j,1)] - H))*(H + (MESH.MP[G(0,j,1)] - H));
			//;

			//Velocidad Radial
			VlocalPrev[LNH(i, j, 0)] = 0.0;
			VlocalPres[LNH(i, j, 0)] = 0.0;

			//Energía Interna
			ElocalPrev[LNH(i, j, 0)] = (Cp[LNH(i, j, 0)] / Gamma) * TlocalPres[LNH(i, j, 0)];
			ElocalPres[LNH(i, j, 0)] = (Cp[LNH(i, j, 0)] / Gamma) * TlocalPres[LNH(i, j, 0)];
		}
	}		
}


//Asignación de temperaturas a las condiciones de contorno
void Solver::UpdateBoundaryConditions(Mesher MESH){
int i, j;

	if (Rank == 0){
		for(j = 0; j < NR; j++){
			RhoLeft[j] = RhoM; 									//Densidad pared izquierda
			Uleft[j] = ((3.0*muW*BulkRe)/(2.0*RhoM*pow(H,3.0)))*(H - (MESH.MP[G(0,j,1)] - H))*(H + (MESH.MP[G(0,j,1)] - H));	 									//Velocidad axial pared izquierda
			Vleft[j] = 0.0;							    		//Velocidad radial pared izquierda
			Tleft[j] = 0.50*Twall;									//Temperaturas pared izquierda
			Pleft[j] = Pressure[LNH(0, j, 0)];				//Presión pared izquierda
			muLeft[j] = muTotal[LNH(Ix, j, 0)];					//Viscosidad dinámica pared izquierda
			KLeft[j] = K[LNH(Ix, j, 0)];						//Conductividad térmica pared izquierda
			Eleft[j] = (Cp[LNH(0, j, 0)] / Gamma) * Tleft[j];  //Energía Interna pared izquierda
		}
	}

	if (Rank == Procesos - 1){
		for(j = 0; j < NR; j++){
			RhoRight[j] = RhoLocalPres[LNH(Fx - 1, j, 0)]; //Densidad pared derecha
			Uright[j] =  UlocalPres[LNH(Fx - 1, j, 0)];		   //Velocidad axial pared derecha
			Vright[j] = VlocalPres[LNH(Fx - 1, j, 0)];	   //Velocidad radial pared derecha
			Tright[j] = TlocalPres[LNH(Fx - 1, j, 0)];	   //Temperaturas pared derecha
			Pright[j] = Pressure[LNH(Fx - 1, j, 0)];	   //Presión axial pared derecha
			muRight[j] = muTotal[LNH(Fx - 1, j, 0)];	   //Viscosidad dinámica pared derecha
			KRight[j] = K[LNH(Fx - 1, j, 0)];			   //Conductividad térmica pared derecha
			Eright[j] = ElocalPres[LNH(Fx - 1, j, 0)];	   //Energía Interna pared derecha
		}
	}

	for(i = Ix; i < Fx; i++){

			//Variables Parte de arriba
			RhoUp[i - Ix] = RhoLocalPres[LNH(i, NR-1, 0)];
			Uup[i - Ix] = 0.0;
			Vup[i - Ix] = 0.0;
			Tup[i - Ix] = Twall;
			Pup[i - Ix] = Pressure[LNH(i, NR - 1, 0)];
			muUp[i - Ix] = muTotal[LNH(i, NR - 1, 0)];
			KUp[i - Ix] = K[LNH(i, NR - 1, 0)];
			Eup[i - Ix] = (Cp[LNH(i, NR - 1, 0)] / Gamma) * Tup[i];

			//Variables Parte de abajo
			RhoDown[i - Ix] = RhoLocalPres[LNH(i, 0, 0)];
			Udown[i - Ix] = 0.0;
			Vdown[i - Ix] = 0.0;
			Tdown[i - Ix] = Twall;
			Pdown[i - Ix] = Pressure[LNH(i, 0, 0)];
			muDown[i - Ix] = muTotal[LNH(i, 0, 0)];
			Kdown[i - Ix] = K[LNH(i, 0, 0)];
			Edown[i - Ix] = (Cp[LNH(i, 0, 0)] / Gamma) * Tdown[i];

	}

}

//Seteo inicial de las Contribuciones de las ecuaciones
void Solver::InitialF(){
int i, j;

	for (i = Ix; i < Fx; i++){
		for (j = 0; j < NR; j++){
			//Densidad
			FRhoPrev[LNH(i, j, 0)] = 0.0;
			FRhoPres[LNH(i, j, 0)] = 0.0;

			//Velocidad Axial
			FUPrev[LNH(i, j, 0)] = 0.0;
			FUPres[LNH(i, j, 0)] = 0.0;

			//Velocidad Radial
			FVPrev[LNH(i, j, 0)] = 0.0;
			FVPres[LNH(i, j, 0)] = 0.0;

			//Energía Interna
			FEprev[LNH(i, j, 0)] = 0.0;
			FEpres[LNH(i, j, 0)] = 0.0;
		}
	}
}

void Solver::Get_FluidViscosity(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			muBase[LNH(i, j, 0)] = muW*pow(TlocalPres[LNH(i,j,0)]/Twall,1.5)*((Twall + S)/(TlocalPres[LNH(i,j,0)] + S));
		}
	}

}

//Cálculo del DeltaT para el siguiente Step
void Solver::Get_StepTime(Mesher MESH, ParPro MPI1){
int i, j;
double DeltasT = 1000.0;
double Tpar = 0.05;
MPI_Status ST;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){

			//CFL Temperatura
			if (abs((Tpar * RhoLocalPres[LNH(i,j,0)] * Cp[LNH(i,j,0)] * pow(MESH.DeltasMP[G(i,j,0)], 2.0)) / (K[LNH(i,j,0)] + 1e-10)) < DeltasT){
				DeltasT = abs((Tpar * RhoLocalPres[LNH(i, j, 0)] * Cp[LNH(i, j, 0)] * pow(MESH.DeltasMP[G(i, j, 0)], 2.0)) / (K[LNH(i, j, 0)] + 1e-10));
			}
			if (abs((Tpar * RhoLocalPres[LNH(i, j, 0)] * Cp[LNH(i, j, 0)] * pow(MESH.DeltasMP[G(i, j, 1)], 2.0)) / (K[LNH(i, j, 0)] + 1e-10)) < DeltasT){
				DeltasT = abs((Tpar * RhoLocalPres[LNH(i, j, 0)] * Cp[LNH(i, j, 0)] * pow(MESH.DeltasMP[G(i, j, 1)], 2.0)) / (K[LNH(i, j, 0)] + 1e-10));
			}

			//CFL Velocidades (U + V)
			if(abs((Tpar*MESH.DeltasMP[G(i,j,0)])/(SoundSpeed + abs(UlocalPres[LNH(i,j,0)])) < DeltasT)){
				DeltasT = (Tpar * MESH.DeltasMP[G(i, j, 0)]) / (SoundSpeed + abs(UlocalPres[LNH(i, j, 0)]));
			}
			if(abs((Tpar*MESH.DeltasMP[G(i,j,1)])/(SoundSpeed + abs(VlocalPres[LNH(i,j,0)])) < DeltasT)){
				DeltasT = (Tpar * MESH.DeltasMP[G(i, j, 1)]) / (SoundSpeed + abs(VlocalPres[LNH(i, j, 0)]));
			}

			//CFL Difusivo
			if((Tpar*RhoLocalPres[LNH(i,j,0)]*pow(MESH.DeltasMP[G(i,j,0)],2.0))/(muTotal[LNH(i,j,0)] + 1e-10) < DeltasT){
				DeltasT = (Tpar * RhoLocalPres[LNH(i, j, 0)] * pow(MESH.DeltasMP[G(i, j, 0)], 2.0)) / (muTotal[LNH(i, j, 0)] + 1e-10);
			}
			if((Tpar*RhoLocalPres[LNH(i,j,0)]*pow(MESH.DeltasMP[G(i,j,1)],2.0))/(muTotal[LNH(i,j,0)] + 1e-10) < DeltasT){
				DeltasT = (Tpar * RhoLocalPres[LNH(i, j, 0)] * pow(MESH.DeltasMP[G(i, j, 1)], 2.0)) / (muTotal[LNH(i, j, 0)] + 1e-10);
			}
		}
	}
	
	MPI1.SendDataToZero(DeltasT, PDT);

	double DT;
	if(Rank == 0){
		DT = PDT[0];
		for(i = 1; i < Procesos; i++){
			if(PDT[i] <= DT){ 
				DT = PDT[i];
			}
		}
	}

	MPI1.SendDataToAll(DT, DeltaT);

}


//Cálculo de las velocidades en las paredes de los volúmenes de control
void Solver::Get_VelocityWalls(Mesher MESH){
int i, j;

	//Nodos U
	if(Rank == 0){
		for(j = 0; j < NR; j++){
			UwallsMU[LU(0,j,0)] = Uleft[j];
			VwallsMU[LU(0,j,0)] = Vleft[j];

			for(i = Ix + 1; i < Fx + 1; i++){
				UwallsMU[LU(i, j, 0)] = 0.50 * (UlocalPres[LNH(i, j, 0)] + UlocalPres[LNH(i - 1, j, 0)]);
				VwallsMU[LU(i, j, 0)] = 0.50 * (VlocalPres[LNH(i, j, 0)] + VlocalPres[LNH(i - 1, j, 0)]);
			}
		}
	}
	else if(Rank == Procesos-1){
		for (j = 0; j < NR; j++){
			UwallsMU[LU(NA, j, 0)] = Uright[j];
			VwallsMU[LU(NA, j, 0)] = Vright[j];

			for (i = Ix; i < NA; i++){
				UwallsMU[LU(i, j, 0)] = 0.50 * (UlocalPres[LNH(i, j, 0)] + UlocalPres[LNH(i - 1, j, 0)]);
				VwallsMU[LU(i, j, 0)] = 0.50 * (VlocalPres[LNH(i, j, 0)] + VlocalPres[LNH(i - 1, j, 0)]);
			}
		}
	}
	else{
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NR; j++){
				UwallsMU[LU(i, j, 0)] = 0.50 * (UlocalPres[LNH(i, j, 0)] + UlocalPres[LNH(i - 1, j, 0)]);
				VwallsMU[LU(i, j, 0)] = 0.50 * (VlocalPres[LNH(i, j, 0)] + VlocalPres[LNH(i - 1, j, 0)]);
			}
		}
	}
	
	//Nodos R
	for(i = Ix; i < Fx; i++){
		UwallsMR[LR(i, 0, 0)] = Udown[i - Ix];
		UwallsMR[LR(i, NR, 0)] = Uup[i - Ix];

		VwallsMR[LR(i, 0, 0)] = Vdown[i - Ix];
		VwallsMR[LR(i, NR, 0)] = Vup[i - Ix];

		for(j = 1; j < NR; j++){
			UwallsMR[LR(i, j, 0)] = 0.50 * (UlocalPres[LNH(i, j, 0)] + UlocalPres[LNH(i, j - 1, 0)]);
			VwallsMR[LR(i, j, 0)] = 0.50 * (VlocalPres[LNH(i, j, 0)] + VlocalPres[LNH(i, j - 1, 0)]);
		}
	}
	
}

void Solver::Get_DensityWalls(Mesher MESH){
int i, j;

	//Nodos U
	if(Rank == 0){
		for (j = 0; j < NR; j++){
			RhoWallsMU[LU(0, j, 0)] = RhoLeft[j];

			for (i = Ix + 1; i < Fx + 1; i++){
				RhoWallsMU[LU(i, j, 0)] = sqrt(RhoLocalPres[LNH(i, j, 0)]) * sqrt(RhoLocalPres[LNH(i - 1, j, 0)]);
			}
		}
	}
	else if(Rank == Procesos - 1){
		for (j = 0; j < NR; j++){
			RhoWallsMU[LU(NA, j, 0)] = RhoRight[j];

			for (i = Ix; i < NA; i++){
				RhoWallsMU[LU(i, j, 0)] = sqrt(RhoLocalPres[LNH(i, j, 0)]) * sqrt(RhoLocalPres[LNH(i - 1, j, 0)]);
			}
		}
	}
	else{
		for (i = Ix; i < Fx + 1; i++){
			for (j = 0; j < NR; j++){
				RhoWallsMU[LU(i, j, 0)] = sqrt(RhoLocalPres[LNH(i, j, 0)]) * sqrt(RhoLocalPres[LNH(i - 1, j, 0)]);
			}
		}
	}
	
	//Nodos R
	for (i = Ix; i < Fx + 1; i++){
		RhoWallsMR[LR(i, 0, 0)] = RhoDown[i - Ix];
		RhoWallsMR[LR(i, NR, 0)] = RhoUp[i - Ix];

		for (j = 1; j < NR; j++){
			RhoWallsMR[LR(i, j, 0)] = sqrt(RhoLocalPres[LNH(i, j, 0)]) * sqrt(RhoLocalPres[LNH(i, j - 1, 0)]);
		}
	}
}


//Calculo de las FRho (Densidad, Conservación de masa)
void Solver::Get_DensityConvective(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){ 
			FRhoPres[LNH(i,j,0)] = 
								 - (1.0/MESH.VolMP[G(i,j,0)])*(
								 + MESH.SupMP[G(i,j,0)]*RhoWallsMU[LU(i,j,0)]*(UwallsMU[LU(i,j,0)]*cos(PI + MESH.AngleMU[GU(i,j,0)]) + VwallsMU[LU(i,j,0)]*sin(PI + MESH.AngleMU[GU(i,j,0)])) 
								 + MESH.SupMP[G(i,j,1)]*RhoWallsMU[LU(i + 1,j,0)]*(UwallsMU[LU(i + 1,j,0)]*cos(MESH.AngleMU[GU(i + 1,j,0)]) + VwallsMU[LU(i + 1,j,0)]*sin(MESH.AngleMU[GU(i + 1,j,0)])) 
								 + MESH.SupMP[G(i,j,2)]*RhoWallsMR[LR(i,j,0)]*(UwallsMR[LR(i,j,0)]*cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + VwallsMR[LR(i,j,0)]*sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])) 
								 + MESH.SupMP[G(i,j,3)]*RhoWallsMR[LR(i,j + 1,0)]*(UwallsMR[LR(i,j + 1,0)]*cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]) + VwallsMR[LR(i,j + 1,0)]*sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]))
								 );			
		}
	}

}

//Cálculo de las presiones en las paredes de los volúmenes de control
void Solver::Get_PressureWalls(Mesher MESH){
int i, j;

	//Nodos U 
	if(Rank == 0){
		for(j = 0; j < NR; j++){
			PwallsMU[LU(0,j,0)] = Pleft[j];
			for(i = Ix + 1; i < Fx+1; i++){
				PwallsMU[LU(i, j, 0)] = Pressure[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((Pressure[LNH(i,j,0)] - Pressure[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}
	else if(Rank == Procesos - 1){
		for (j = 0; j < NR; j++){
			PwallsMU[LU(NA, j, 0)] = Pright[j];
			for (i = Ix; i < NA; i++){
				PwallsMU[LU(i, j, 0)] = Pressure[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((Pressure[LNH(i,j,0)] - Pressure[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}
	else{
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NR; j++){
				PwallsMU[LU(i,j,0)] = Pressure[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((Pressure[LNH(i,j,0)] - Pressure[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}

	//Nodos V
	for(i = Ix; i < Fx; i++){
		PwallsMR[LR(i,0,0)] = Pdown[i - Ix];
		PwallsMR[LR(i,NR,0)] = Pup[i - Ix];

		for(j = 1; j < NR; j++){
			PwallsMR[LR(i,j,0)] = Pressure[LNH(i,j-1,0)] + (MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)])*((Pressure[LNH(i,j,0)] - Pressure[LNH(i,j-1,0)])/(MESH.MP[G(i,j,1)] - MESH.MP[G(i,j-1,1)]));
		}
	}

}

//Cálculo de los gradientes de presión en cada nodo en ambas direciones (U y V)
void Solver::Get_PressureGradients(Mesher MESH){
int i, j;

	//Gradiente Axial
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			PressureGradientAxial[LNH(i,j,0)] = 
											  - (1.0/MESH.VolMP[G(i,j,0)])*(
											  + MESH.SupMP[G(i,j,0)]*PwallsMU[LU(i,j,0)]*cos(PI + MESH.AngleMU[GU(i,j,0)])
											  + MESH.SupMP[G(i,j,1)]*PwallsMU[LU(i + 1,j,0)]*cos(MESH.AngleMU[GU(i + 1,j,0)])
											  + MESH.SupMP[G(i,j,2)]*PwallsMR[LR(i,j,0)]*cos(1.50*PI + MESH.AngleMR[GR(i,j,0)])
											  + MESH.SupMP[G(i,j,3)]*PwallsMR[LR(i,j + 1,0)]*cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])
											  );
		}
	}

	//Gradiente Radial
	for(i = Ix; i < Fx; i++){
		for (j = 1; j < NR-1; j++){
			PressureGradientRadial[LNH(i, j, 0)] =
				                                 - (1.0 / MESH.VolMP[G(i, j, 0)]) * (
												 + MESH.SupMP[G(i, j, 0)] * PwallsMU[LU(i, j, 0)] * sin(PI + MESH.AngleMU[GU(i, j, 0)]) 
												 + MESH.SupMP[G(i, j, 1)] * PwallsMU[LU(i + 1, j, 0)] * sin(MESH.AngleMU[GU(i + 1, j, 0)]) 
												 + MESH.SupMP[G(i, j, 2)] * PwallsMR[LR(i, j, 0)] * sin(1.50 * PI + MESH.AngleMR[GR(i, j, 0)]) 
												 + MESH.SupMP[G(i, j, 3)] * PwallsMR[LR(i, j + 1, 0)] * sin(0.50 * PI + MESH.AngleMR[GR(i, j + 1, 0)])
												 )
											//	 + Pressure[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
												 ;
		}
	}
}

//Cálculo del mapa de presiones con la Ley de los Gases Ideales
void Solver::Get_Pressure(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			Pressure[LNH(i, j, 0)] = RhoLocalPres[LNH(i, j, 0)] * TlocalPres[LNH(i, j, 0)] * Rideal;
		}
	}

}

//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
void Solver::Get_MomentumConvectiveAxial(Mesher MESH){
int i, j;

	for (i = Ix; i < Fx; i++){
		for (j = 0; j < NR; j++){
			MomentumConvectiveAxial[LNH(i,j,0)] = 
									 		    - (1.0 / MESH.VolMP[G(i,j,0)]) * (
									            + MESH.SupMP[G(i,j,0)] * (UwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)]) + VwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])) * RhoWallsMU[LU(i,j,0)] * UwallsMU[LU(i,j,0)] 
									            + MESH.SupMP[G(i,j,1)] * (UwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)]) + VwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)])) * RhoWallsMU[LU(i + 1,j,0)] * UwallsMU[LU(i + 1,j,0)] 
									            + MESH.SupMP[G(i,j,2)] * (UwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + VwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])) * RhoWallsMR[LR(i,j,0)] * UwallsMR[LR(i,j,0)] 
									            + MESH.SupMP[G(i,j,3)] * (UwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]) + VwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])) * RhoWallsMR[LR(i,j + 1,0)] * UwallsMR[LR(i,j + 1,0)] 
												);
		}
	}
}

//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomentumConvectiveRadial(Mesher MESH){
int i, j;

	for (i = Ix; i < Fx; i++){
		for (j = 1; j < NR; j++){
			MomentumConvectiveRadial[LNH(i,j,0)] = 
										  		 - (1.0 / MESH.VolMP[G(i,j,0)]) * (
										   		 + MESH.SupMP[G(i,j,0)] * (UwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)]) + VwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])) * RhoWallsMU[LU(i,j,0)] * VwallsMU[LU(i,j,0)] 
										   		 + MESH.SupMP[G(i,j,1)] * (UwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)]) + VwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)])) * RhoWallsMU[LU(i + 1,j,0)] * VwallsMU[LU(i + 1,j,0)] 
										   		 + MESH.SupMP[G(i,j,2)] * (UwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + VwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])) * RhoWallsMR[LR(i,j,0)] * VwallsMR[LR(i,j,0)] 
										   		 + MESH.SupMP[G(i,j,3)] * (UwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]) + VwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])) * RhoWallsMR[LR(i,j + 1,0)] * VwallsMR[LR(i,j + 1,0)] 
										   		 );
		}
	}
}

void Solver::Get_Stresses(Mesher MESH, ParPro MPI1){
int i, j;

	//Cálculo de la Divergencia en cada VC
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			Divergence[LNH(i,j,0)] = 
								   + (1.0 / MESH.VolMP[G(i,j,0)]) * ( 
								   + MESH.SupMP[G(i,j,0)] * (UwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)]) + VwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])) 
								   + MESH.SupMP[G(i,j,1)] * (UwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)]) + VwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)]))
								   + MESH.SupMP[G(i,j,2)] * (UwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + VwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])) 
								   + MESH.SupMP[G(i,j,3)] * (UwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]) + VwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]))
								   );
		}
	}

	//Cálculo de la divergencia en los volúmenes del inlet y outlet
	if(Rank == 0){
		//Parte central
		for(j = 0; j < NR; j++){
			DivergenceMU[j] = 
							+ (2.0 / MESH.VolMP[G(0,j,0)]) * (
							+ MESH.SupMP[G(0,j,0)] * (UwallsMU[LU(0,j,0)] * cos(PI + MESH.AngleMU[GU(0,j,0)]) + VwallsMU[LU(0,j,0)] * sin(PI + MESH.AngleMU[GU(0,j,0)]))  
							+ MESH.SupMP[G(0,j,1)] * (UlocalPres[LNH(0,j,0)] * cos(MESH.AngleMU[GU(1,j,0)]) + VlocalPres[LNH(0,j,0)] * sin(MESH.AngleMU[GU(1,j,0)]))
							+ 0.50*MESH.SupMP[G(0,j,2)] * (UwallsMR[LR(0,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(0,j,0)]) + VwallsMR[LR(0,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(0,j,0)])) 
							+ 0.50*MESH.SupMP[G(0,j,3)] * (UwallsMR[LR(0,j+1,0)] * cos(0.50*PI + MESH.AngleMR[GR(0,j+1,0)]) + VwallsMR[LR(0,j+1,0)] * sin(0.50*PI + MESH.AngleMR[GR(0,j+1,0)])) 
							);
		}	
	}
	else if(Rank == Procesos - 1){
		//Parte central
		for(j = 0; j < NR; j++){
			DivergenceMU[j] = 
							+ (2.0 / MESH.VolMP[G(NA-1,j,0)]) * (
							+ MESH.SupMP[G(NA-1,j,0)] * (UlocalPres[LNH(NA - 1,j,0)] * cos(PI + MESH.AngleMU[GU(NA - 1,j,0)]) + VlocalPres[LNH(NA - 1,j,0)] * sin(PI + MESH.AngleMU[GU(NA - 1,j,0)]))  
							+ MESH.SupMP[G(NA-1,j,1)] * (UwallsMU[LU(NA,j,0)] * cos(MESH.AngleMU[GU(NA,j,0)]) + VwallsMU[LU(NA,j,0)] * sin(MESH.AngleMU[GU(NA,j,0)]))
							+ 0.50*MESH.SupMP[G(NA-1,j,2)] * (UwallsMR[LR(NA - 1,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(NA - 1,j,0)]) + VwallsMR[LR(NA - 1,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(NA - 1,j,0)])) 
							+ 0.50*MESH.SupMP[G(NA-1,j,3)] * (UwallsMR[LR(NA - 1,j+1,0)] * cos(0.50*PI + MESH.AngleMR[GR(NA - 1,j+1,0)]) + VwallsMR[LR(NA - 1,j+1,0)] * sin(0.50*PI + MESH.AngleMR[GR(NA - 1,j+1,0)])) 
							);
		}	
	}

	//Cálculo de la divergencia en los volúmenes del eje y el radio máximo 
	for(i = Ix; i < Fx; i++){
		//Parte abajo
		DivergenceMR[VR(i - Ix, 0, 0)] = 
								   	   + (2.0 / MESH.VolMP[G(i,0,0)]) * ( 
								       + 0.50*MESH.SupMP[G(i,0,0)] * (UwallsMU[LU(i,0,0)] * cos(PI + MESH.AngleMU[GU(i,0,0)]) + VwallsMU[LU(i,0,0)] * sin(PI + MESH.AngleMU[GU(i,0,0)])) 
								       + 0.50*MESH.SupMP[G(i,0,1)] * (UwallsMU[LU(i + 1,0,0)] * cos(MESH.AngleMU[GU(i + 1,0,0)]) + VwallsMU[LU(i + 1,0,0)] * sin(MESH.AngleMU[GU(i + 1,0,0)]))
								       + MESH.SupMP[G(i,0,2)] * (UwallsMR[LR(i,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,0,0)]) + VwallsMR[LR(i,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])) 
								       + MESH.SupMP[G(i,0,3)] * (UlocalPres[LNH(i,0,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,1,0)]) + VlocalPres[LNH(i,0,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)]))
								       );

		//Parte arriba
		DivergenceMR[VR(i - Ix, 1, 0)] =
								       + (2.0 / MESH.VolMP[G(i,NR-1,0)]) * ( 
								   	   + 0.50*MESH.SupMP[G(i,NR-1,0)] * (UwallsMU[LU(i,NR - 1,0)] * cos(PI + MESH.AngleMU[GU(i,NR - 1,0)]) + VwallsMU[LU(i,NR - 1,0)] * sin(PI + MESH.AngleMU[GU(i,NR - 1,0)])) 
								       + 0.50*MESH.SupMP[G(i,NR-1,1)] * (UwallsMU[LU(i + 1,NR - 1,0)] * cos(MESH.AngleMU[GU(i + 1,NR - 1,0)]) + VwallsMU[LU(i + 1,NR - 1,0)] * sin(MESH.AngleMU[GU(i + 1,NR - 1,0)]))
								       + MESH.SupMP[G(i,NR-1,2)] * (UlocalPres[LNH(i,NR - 1,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,NR - 1,0)]) + VlocalPres[LR(i,NR - 1,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,NR - 1,0)]))   
								       + MESH.SupMP[G(i,NR-1,3)] * (UwallsMR[LR(i,NR,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,NR,0)]) + VwallsMR[LR(i,NR,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,NR,0)]))   
								       );

	}

	//Cálculo Esfuerzos Viscosos

	//Cálculo Esfuerzo Viscoso Radial
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			TauRR[LNH(i,j,0)] = 
							  + muTotal[LNH(i,j,0)]*(
							  + 2.0 * ( (1.0 / MESH.VolMP[G(i,j,0)]) * (
							  + MESH.SupMP[G(i,j,0)] * VwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])
							  + MESH.SupMP[G(i,j,1)] * VwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)])
							  + MESH.SupMP[G(i,j,2)] * VwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  + MESH.SupMP[G(i,j,3)] * VwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  )
							  //- VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
							  )
							  - (2.0/3.0)*Divergence[LNH(i,j,0)]
							  );
		}
	}

	//Cálculo Esfuerzo Viscoso Axial
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
				TauZZ[LNH(i,j,0)] = 
							      + muTotal[LNH(i,j,0)]*(
							      + 2.0 * ( (1.0 / MESH.VolMP[G(i,j,0)]) * (
							      + MESH.SupMP[G(i,j,0)] * UwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)])
							      + MESH.SupMP[G(i,j,1)] * UwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)])
							      + MESH.SupMP[G(i,j,2)] * UwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							      + MESH.SupMP[G(i,j,3)] * UwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							      )
							  
							      )
							      - (2.0/3.0)*Divergence[LNH(i,j,0)]
							      );
		}
	}

	//Cálculo Esfuerzo Viscoso Axial-Radial
	for(i = Ix; i < Fx; i++){
		TauRZ[LNH(i,0,0)] = 
						  + muTotal[LNH(i,0,0)]*(
						  + ((1.0 / MESH.VolMP[G(i,0,0)]) * (
						  + MESH.SupMP[G(i,0,0)] * UwallsMU[LU(i,0,0)] * sin(PI + MESH.AngleMU[GU(i,0,0)])
						  + MESH.SupMP[G(i,0,1)] * UwallsMU[LU(i + 1,0,0)] * sin(MESH.AngleMU[GU(i + 1,0,0)])
						  + MESH.SupMP[G(i,0,2)] * UwallsMR[LR(i,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])
						  + MESH.SupMP[G(i,0,3)] * UwallsMR[LR(i,1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
						  )
						//  - (0.50*(UwallsMR[LR(i,0,0)] + UwallsMR[LR(i,1,0)]))/(0.50*(MESH.MR[GR(i,0,1)] + MESH.MR[GR(i,1,1)]))
						  )

						  + (1.0 / MESH.VolMP[G(i,0,0)]) * (
						  + MESH.SupMP[G(i,0,0)] * VwallsMU[LU(i,0,0)] * cos(PI + MESH.AngleMU[GU(i,0,0)])
						  + MESH.SupMP[G(i,0,1)] * VwallsMU[LU(i + 1,0,0)] * cos(MESH.AngleMU[GU(i + 1,0,0)])
						  + MESH.SupMP[G(i,0,2)] * VwallsMR[LR(i,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,0,0)])
						  + MESH.SupMP[G(i,0,3)] * VwallsMR[LR(i,1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
						  )
						  );

		for(j = 1; j < NR; j++){
			TauRZ[LNH(i,j,0)] = 
							  + muTotal[LNH(i,j,0)]*(
							  + ((1.0 / MESH.VolMP[G(i,j,0)]) * (
						      + MESH.SupMP[G(i,j,0)] * UwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])
							  + MESH.SupMP[G(i,j,1)] * UwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)])
							  + MESH.SupMP[G(i,j,2)] * UwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  + MESH.SupMP[G(i,j,3)] * UwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  )
							//  - UlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
							  )

							  + (1.0 / MESH.VolMP[G(i,j,0)]) * (
							  + MESH.SupMP[G(i,j,0)] * VwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)])
							  + MESH.SupMP[G(i,j,1)] * VwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)])
							  + MESH.SupMP[G(i,j,2)] * VwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  + MESH.SupMP[G(i,j,3)] * VwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  )
							  );
		}
	}

	//Cálculo Esfuerzo Viscoso Angular
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			TauThetaTheta[LNH(i,j,0)] = 
									  - muTotal[LNH(i,j,0)]* (
									  + 2.0*(VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)])
									  - (2.0/3.0)*Divergence[LNH(i,j,0)]
									  );
		}
	}

	//Envío de los Halos de los esfuerzos viscosos
	MPI1.SendData(TauRR, Ix, Fx);
	MPI1.ReceiveData(TauRR, Ix, Fx);
	

	MPI1.SendData(TauRZ, Ix, Fx);
	MPI1.ReceiveData(TauRZ, Ix, Fx);
		

	MPI1.SendData(TauZZ, Ix, Fx);
	MPI1.ReceiveData(TauZZ, Ix, Fx);
	

	//Nodos U
	if(Rank == 0){
		//Esquina abajo izquierda
		//TauRR_mu[LU(0,j,0)] = 0.0;
	//	TauRZ_mu[LU(0,0,0)] = TauRZ[LNH(0,0,0)];
					/*		- muTotal[LNH(0,0,0)]*(
							+ ((1.0 / MESH.VolMU[0]) * (
							+ MESH.SupMU[MU(0,0,0)] * UwallsMU[LU(0,0,0)] * sin(PI + MESH.AngleMU[GU(0,0,0)])
							+ MESH.SupMU[MU(0,0,1)] * UlocalPres[LNH(0,0,0)] * sin(MESH.AngleMU[GU(1,0,0)])
						//	+ MESH.SupMU[MU(0,0,2)] * UwallsMR[LR(0,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(0,0,0)])
						///	+ MESH.SupMU[MU(0,0,3)] * UwallsMR[LR(0,1,0)] * sin(0.50*PI + MESH.AngleMR[GR(0,1,0)])	  
							)
						//	- (0.50*(UwallsMR[LR(0,0,0)] + UwallsMR[LR(0,1,0)]))/(0.50*(MESH.MR[GR(0,0,1)] + MESH.MR[GR(0,1,1)]))
							)

							+ (1.0 / MESH.VolMU[0]) * (
							+ MESH.SupMU[MU(0,0,0)] * VwallsMU[LU(0,0,0)] * cos(PI + MESH.AngleMU[GU(0,0,0)])
							+ MESH.SupMU[MU(0,0,1)] * VlocalPres[LNH(0,0,0)] * cos(MESH.AngleMU[GU(1,0,0)])
							+ MESH.SupMU[MU(0,0,2)] * VwallsMR[LR(0,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(0,0,0)])
							+ MESH.SupMU[MU(0,0,3)] * VwallsMR[LR(0,1,0)] * cos(0.50*PI + MESH.AngleMR[GR(0,1,0)])	  
							)
							);*/

	//	TauZZ_mu[LU(0,0,0)] = TauZZ[LNH(0,0,0)];
						/*	- muTotal[LNH(0,0,0)]*(
							+ 2.0 * ( (1.0 / MESH.VolMU[0]) * (
							+ MESH.SupMU[MU(0,0,0)] * UwallsMU[LU(0,0,0)] * cos(PI + MESH.AngleMU[GU(0,0,0)])
							+ MESH.SupMU[MU(0,0,1)] * UlocalPres[LNH(0,0,0)] * cos(MESH.AngleMU[GU(1,0,0)])
							+ MESH.SupMU[MU(0,0,2)] * UwallsMR[LR(0,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(0,0,0)])
							+ MESH.SupMU[MU(0,0,3)] * UwallsMR[LR(0,1,0)] * cos(0.50*PI + MESH.AngleMR[GR(0,1,0)])	  
							)
							  
							)
							- (2.0/3.0)*DivergenceMU[0]
							);*/

		//Parte abajo
	/*	for(i = Ix + 1; i < Fx + 1; i++){
				TauRR_mu[LU(i, 0, 0)] = 0.50 * (TauRR[LNH(i, 0, 0)] + TauRR[LNH(i - 1, 0, 0)]);
				TauRZ_mu[LU(i, 0, 0)] = 0.50 * (TauRZ[LNH(i, 0, 0)] + TauRZ[LNH(i - 1, 0, 0)]);
				TauZZ_mu[LU(i, 0, 0)] = TauZZ_mu[LU(i-1,0,0)] + (MESH.MU[GU(i,0,0)] - MESH.MU[GU(i-1,0,0)])*((TauZZ_mu[LU(i,0,0)] - TauZZ_mu[LU(i-1,0,0)])/(MESH.MU[GU(i,0,0)] - MESH.MU[GU(i-1,0,0)]));
		}*/

		for(j = 0; j < NR; j++){
			//Parte izquierda
			TauRR_mu[LU(0,j,0)] =  
								+ muTotal[LNH(0,j,0)]*(
							    + 2.0 * ( (2.0 / MESH.VolMP[G(0,j,0)]) * (
							    + MESH.SupMP[G(0,j,0)] * VwallsMU[LU(0,j,0)] * sin(PI + MESH.AngleMU[GU(0,j,0)])
							    + MESH.SupMP[G(0,j,1)] * VlocalPres[LNH(0,j,0)] * sin(MESH.AngleMU[GU(1,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,2)] * VwallsMR[LR(0,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(0,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,3)] * VwallsMR[LR(0,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(0,j + 1,0)])	  
							    )
							    //- VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
							    )
							    - (2.0/3.0)*DivergenceMU[j]
							    );

			TauRZ_mu[LU(0,j,0)] = 
								+ muTotal[LNH(0,j,0)]*(
							    + ((2.0 / MESH.VolMP[G(0,j,0)]) * (
							    + MESH.SupMP[G(0,j,0)] * UwallsMU[LU(0,j,0)] * sin(PI + MESH.AngleMU[GU(0,j,0)])
							    + MESH.SupMP[G(0,j,1)] * UlocalPres[LNH(0,j,0)] * sin(MESH.AngleMU[GU(1,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,2)] * UwallsMR[LR(0,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(0,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,3)] * UwallsMR[LR(0,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(0,j + 1,0)])	  
							    )
							    //-  UwallsMU[LU(0,j,0)]/MESH.MU[GU(0,j,1)]
								)

							    + (2.0 / MESH.VolMP[G(0,j,0)]) * (
							    + MESH.SupMP[G(0,j,0)] * VwallsMU[LU(0,j,0)] * cos(PI + MESH.AngleMU[GU(0,j,0)])
							    + MESH.SupMP[G(0,j,1)] * VlocalPres[LNH(0,j,0)] * cos(MESH.AngleMU[GU(1,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,2)] * VwallsMR[LR(0,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(0,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,3)] * VwallsMR[LR(0,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(0,j + 1,0)])	  
							    )
							    );

			TauZZ_mu[LU(0,j,0)] =
								+ muTotal[LNH(0,j,0)]*(
							    + 2.0 * ( (2.0 / MESH.VolMP[G(0,j,0)]) * (
							    + MESH.SupMP[G(0,j,0)] * UwallsMU[LU(0,j,0)] * cos(PI + MESH.AngleMU[GU(0,j,0)])
							    + MESH.SupMP[G(0,j,1)] * UlocalPres[LNH(0,j,0)] * cos(MESH.AngleMU[GU(1,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,2)] * UwallsMR[LR(0,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(0,j,0)])
							    + 0.50*MESH.SupMP[G(0,j,3)] * UwallsMR[LR(0,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(0,j + 1,0)])	  
							    )
							  
							    )
							    - (2.0/3.0)*DivergenceMU[j]
							    );
			//Parte centro
			for(i = Ix + 1; i < Fx + 1; i++){
				TauRR_mu[LU(i, j, 0)] = TauRR[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauRR[LNH(i,j,0)] - TauRR[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				TauRZ_mu[LU(i, j, 0)] = TauRZ[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauRZ[LNH(i,j,0)] - TauRZ[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				TauZZ_mu[LU(i, j, 0)] = TauZZ[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauZZ[LNH(i,j,0)] - TauZZ[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}
	else if(Rank == Procesos-1){

		//Esquina abajo derecha
		//TauRR_mu[LU(0,j,0)] = 0.0;
//		TauRZ_mu[LU(NA,0,0)] = TauRZ[LNH(NA-1,0,0)];
						/*	- muTotal[LNH(NA - 1,0,0)]*(
							+ ((1.0 / MESH.VolMU[0]) * (
							+ MESH.SupMU[MU(0,0,0)] * UlocalPres[LNH(NA - 1,0,0)] * sin(PI + MESH.AngleMU[GU(NA - 1,0,0)])				
							+ MESH.SupMU[MU(0,0,1)] * UwallsMU[LU(NA,0,0)] * sin(MESH.AngleMU[GU(NA,0,0)])
							+ MESH.SupMU[MU(0,0,2)] * UwallsMR[LR(NA - 1,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(NA - 1,0,0)])
							+ MESH.SupMU[MU(0,0,3)] * UwallsMR[LR(NA - 1,1,0)] * sin(0.50*PI + MESH.AngleMR[GR(NA - 1,1,0)])	  
							)
						//	- (0.50*(UwallsMR[LR(NA - 1,0,0)] + UwallsMR[LR(NA - 1,1,0)]))/(0.50*(MESH.MR[GR(NA - 1,0,1)] + MESH.MR[GR(NA - 1,1,1)]))
							)

							+ (1.0 / MESH.VolMU[0]) * (
							+ MESH.SupMU[MU(0,0,0)] * VlocalPres[LNH(NA - 1,0,0)] * cos(PI + MESH.AngleMU[GU(NA - 1,0,0)])
							+ MESH.SupMU[MU(0,0,1)] * VwallsMU[LU(NA,0,0)] * cos(MESH.AngleMU[GU(NA,0,0)])
							+ MESH.SupMU[MU(0,0,2)] * VwallsMR[LR(NA - 1,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(NA - 1,0,0)])
							+ MESH.SupMU[MU(0,0,3)] * VwallsMR[LR(NA - 1,1,0)] * cos(0.50*PI + MESH.AngleMR[GR(NA - 1,1,0)])	  
							)
							);*/

	//	TauZZ_mu[LU(NA,0,0)] = TauZZ[LNH(NA-1,0,0)];
						/*	 - muTotal[LNH(NA - 1,0,0)]*(
							 + 2.0 * ( (1.0 / MESH.VolMU[0]) * (
							 + MESH.SupMU[MU(0,0,0)] * UlocalPres[LNH(NA - 1,0,0)] * cos(PI + MESH.AngleMU[GU(NA - 1,0,0)])
							 + MESH.SupMU[MU(0,0,1)] * UwallsMU[LU(NA,0,0)] * cos(MESH.AngleMU[GU(NA,0,0)])
							 + MESH.SupMU[MU(0,0,2)] * UwallsMR[LR(NA - 1,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(NA - 1,0,0)])
							 + MESH.SupMU[MU(0,0,3)] * UwallsMR[LR(NA - 1,1,0)] * cos(0.50*PI + MESH.AngleMR[GR(NA - 1,1,0)])	  
							 )
							  
							 )
							 - (2.0/3.0)*DivergenceMU[0]
							 );*/

		//Parte abajo
	/*	for(i = Ix; i < Fx; i++){
			TauRR_mu[LU(i, 0, 0)] = 0.50 * (TauRR[LNH(i, 0, 0)] + TauRR[LNH(i - 1, 0, 0)]);
			TauRZ_mu[LU(i, 0, 0)] = 0.50 * (TauRZ[LNH(i, 0, 0)] + TauRZ[LNH(i - 1, 0, 0)]);
			TauZZ_mu[LU(i, 0, 0)] = TauZZ_mu[LU(i-1,0,0)] + (MESH.MU[GU(i,0,0)] - MESH.MU[GU(i-1,0,0)])*((TauZZ_mu[LU(i,0,0)] - TauZZ_mu[LU(i-1,0,0)])/(MESH.MU[GU(i,0,0)] - MESH.MU[GU(i-1,0,0)]));
		}*/

		for (j = 0; j < NR; j++){
			//Parte derecha
			TauRR_mu[LU(NA,j,0)] = 
								 + muTotal[LNH(NA-1,j,0)]*(
							     + 2.0 * ( (2.0 / MESH.VolMP[G(NA-1,j,0)]) * (
							     + MESH.SupMP[G(NA-1,j,0)] * VlocalPres[LNH(NA-1,j,0)] * sin(PI + MESH.AngleMU[GU(NA-1,j,0)])
							     + MESH.SupMP[G(NA-1,j,1)] * VwallsMU[LU(NA,j,0)] * sin(MESH.AngleMU[GU(NA,j,0)])
							     + 0.50*MESH.SupMP[G(NA-1,j,2)] * VwallsMR[LR(NA-1,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(NA-1,j,0)])
							     + 0.50*MESH.SupMP[G(NA-1,j,3)] * VwallsMR[LR(NA-1,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(NA-1,j + 1,0)])	  
							     )
							     //- VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
							     )
							     - (2.0/3.0)*DivergenceMU[j]
							     );

			TauRZ_mu[LU(NA,j,0)] = 
							     + muTotal[LNH(NA - 1,j,0)]*(
							     + ((2.0 / MESH.VolMP[G(NA-1,j,0)]) * (
							     + MESH.SupMP[G(NA-1,j,0)] * UlocalPres[LNH(NA - 1,j,0)] * sin(PI + MESH.AngleMU[GU(NA - 1,j,0)])
							     + MESH.SupMP[G(NA-1,j,1)] * UwallsMU[LU(NA,j,0)] * sin(MESH.AngleMU[GU(NA,j,0)])
							     + 0.50*MESH.SupMP[G(NA-1,j,2)] * UwallsMR[LR(NA - 1,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(NA - 1,j,0)])
							     + 0.50*MESH.SupMP[G(NA-1,j,3)] * UwallsMR[LR(NA - 1,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(NA - 1,j + 1,0)])	  
							     )
							     //- UwallsMU[LU(NA,j,0)]/MESH.MU[GU(NA,j,1)]
							     )

							     + (2.0 / MESH.VolMP[G(NA-1,j,0)]) * (
							     + MESH.SupMP[G(NA-1,j,0)] * VlocalPres[LNH(NA - 1,j,0)] * cos(PI + MESH.AngleMU[GU(NA - 1,j,0)])
							     + MESH.SupMP[G(NA-1,j,1)] * VwallsMU[LU(NA,j,0)] * cos(MESH.AngleMU[GU(NA,j,0)]) 
							     + 0.50*MESH.SupMP[G(NA-1,j,2)] * VwallsMR[LR(NA - 1,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(NA - 1,j,0)])
							     + 0.50*MESH.SupMP[G(NA-1,j,3)] * VwallsMR[LR(NA - 1,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(NA - 1,j + 1,0)])	  
							     )
							     );

			TauZZ_mu[LU(NA,j,0)] = 
								 + muTotal[LNH(NA - 1,j,0)]*(
								 + 2.0 * ( (2.0 / MESH.VolMP[G(NA-1,j,0)]) * (
								 + MESH.SupMP[G(NA-1,j,0)] * UlocalPres[LNH(NA - 1,j,0)] * cos(PI + MESH.AngleMU[GU(NA - 1,j,0)])
								 + MESH.SupMP[G(NA-1,j,1)] * UwallsMU[LU(NA,j,0)] * cos(MESH.AngleMU[GU(NA,j,0)])
								 + 0.50*MESH.SupMP[G(NA-1,j,2)] * UwallsMR[LR(NA - 1,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(NA - 1,j,0)])
								 + 0.50*MESH.SupMP[G(NA-1,j,3)] * UwallsMR[LR(NA - 1,j+1,0)] * cos(0.50*PI + MESH.AngleMR[GR(NA - 1,j+1,0)])	  
								 )
							  
								 )
								 - (2.0/3.0)*DivergenceMU[j]
								 ); 

			for (i = Ix; i < Fx; i++){
				TauRR_mu[LU(i, j, 0)] = TauRR[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauRR[LNH(i,j,0)] - TauRR[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				TauRZ_mu[LU(i, j, 0)] = TauRZ[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauRZ[LNH(i,j,0)] - TauRZ[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				TauZZ_mu[LU(i, j, 0)] = TauZZ[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauZZ[LNH(i,j,0)] - TauZZ[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}
	else{

		for(i = Ix; i < Fx + 1; i++){
			//Parte centro
			for(j = 0; j < NR; j++){
				TauRR_mu[LU(i, j, 0)] = TauRR[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauRR[LNH(i,j,0)] - TauRR[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				TauRZ_mu[LU(i, j, 0)] = TauRZ[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauRZ[LNH(i,j,0)] - TauRZ[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				TauZZ_mu[LU(i, j, 0)] = TauZZ[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((TauZZ[LNH(i,j,0)] - TauZZ[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}

	}

	//Nodos R
	for(i = Ix; i < Fx; i++){
		//Parte abajo
		TauRR_mr[LR(i, 0, 0)] = 
							  + muTotal[LNH(i,0,0)]*(
							  + 2.0 * ( (2.0 / MESH.VolMP[G(i,0,0)]) * (
							  + 0.50*MESH.SupMP[G(i,0,0)] * VwallsMU[LU(i,0,0)] * sin(PI + MESH.AngleMU[GU(i,0,0)])
							  + 0.50*MESH.SupMP[G(i,0,1)] * VwallsMU[LU(i+1,0,0)] * sin(MESH.AngleMU[GU(i+1,0,0)])
							  + MESH.SupMP[G(i,0,2)] * VwallsMR[LR(i,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])
							  + MESH.SupMP[G(i,0,3)] *  VlocalPres[LNH(i,0,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
							  )
							  //- VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
							  )
							  - (2.0/3.0)*DivergenceMR[VR(i - Ix,0,0)]
							  );

		TauRZ_mr[LR(i, 0, 0)] = 
							  + muTotal[LNH(i,0,0)]*(
						  	  + ((2.0 / MESH.VolMP[G(i,0,0)]) * (
						  	  + 0.50*MESH.SupMP[G(i,0,0)] * UwallsMU[LU(i,0,0)] * sin(PI + MESH.AngleMU[GU(i,0,0)])
						  	  + 0.50*MESH.SupMP[G(i,0,1)] * UwallsMU[LU(i + 1,0,0)] * sin(MESH.AngleMU[GU(i + 1,0,0)])
						  	  + MESH.SupMP[G(i,0,2)] * UwallsMR[LR(i,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])
						  	  + MESH.SupMP[G(i,0,3)] * UlocalPres[LNH(i,0,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
						  	  )
						  	  //- (0.50*(UwallsMR[LR(i,0,0)] + 0.50*(UwallsMR[LR(i,0,0)] + UwallsMR[LR(i,1,0)])))/(0.50*(MESH.MR[GR(i,0,1)] + 0.50*(MESH.MR[GR(i,0,1)] + MESH.MR[GR(i,1,1)])))
							  )

						  	  + (2.0 / MESH.VolMP[G(i,0,0)]) * (
						  	  + 0.50*MESH.SupMP[G(i,0,0)] * VwallsMU[LU(i,0,0)] * cos(PI + MESH.AngleMU[GU(i,0,0)])
						  	  + 0.50*MESH.SupMP[G(i,0,1)] * VwallsMU[LU(i + 1,0,0)] * cos(MESH.AngleMU[GU(i + 1,0,0)])
						  	  + MESH.SupMP[G(i,0,2)] * VwallsMR[LR(i,0,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,0,0)])
						  	  + MESH.SupMP[G(i,0,3)] * VlocalPres[LNH(i,0,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
						  	  )
						  	  );

		TauZZ_mr[LR(i, 0, 0)] = 
							  + muTotal[LNH(i,0,0)]*(
							  + 2.0 * ( (2.0 / MESH.VolMP[G(i,0,0)]) * (
							  + 0.50*MESH.SupMP[G(i,0,0)] * UwallsMU[LU(i,0,0)] * sin(PI + MESH.AngleMU[GU(i,0,0)])
							  + 0.50*MESH.SupMP[G(i,0,1)] * UwallsMU[LU(i + 1,0,0)] * sin(MESH.AngleMU[GU(i + 1,0,0)])
							  + MESH.SupMP[G(i,0,2)] * UwallsMR[LR(i,0,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])
							  + MESH.SupMP[G(i,0,3)] * UlocalPres[LNH(i,0,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
							  )
							  
							  )
							  - (2.0/3.0)*DivergenceMR[VR(i - Ix,0,0)]
							  );
	
		//Parte arriba
		TauRR_mr[LR(i, NR, 0)] =
							   + muTotal[LNH(i,NR-1,0)]*(
							   + 2.0 * ( (2.0 / MESH.VolMP[G(i,NR-1,0)]) * (
							   + 0.50*MESH.SupMP[G(i,NR-1,0)] * VwallsMU[LU(i,NR-1,0)] * sin(PI + MESH.AngleMU[GU(i,NR-1,0)])
							   + 0.50*MESH.SupMP[G(i,NR-1,1)] * VwallsMU[LU(i+1,NR-1,0)] * sin(MESH.AngleMU[GU(i+1,NR-1,0)])
							   + MESH.SupMP[G(i,NR-1,2)] * VlocalPres[LNH(i,NR-1,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,NR-1,0)])
							   + MESH.SupMP[G(i,NR-1,3)] * VwallsMR[LR(i,NR,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,NR,0)])	  
							   )
							   //- VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
							   )
							   - (2.0/3.0)*DivergenceMR[VR(i - Ix,1,0)]
							   );

		TauRZ_mr[LR(i, NR, 0)] = 
							   + muTotal[LNH(i,NR-1,0)]*(
						  	   + ((2.0 / MESH.VolMP[G(i,NR-1,0)]) * (
						  	   + 0.50*MESH.SupMP[G(i,NR-1,0)] * UwallsMU[LU(i,NR-1,0)] * sin(PI + MESH.AngleMU[GU(i,NR-1,0)])
						  	   + 0.50*MESH.SupMP[G(i,NR-1,1)] * UwallsMU[LU(i + 1,NR-1,0)] * sin(MESH.AngleMU[GU(i + 1,NR-1,0)])
						  	   + MESH.SupMP[G(i,NR-1,2)] * UlocalPres[LNH(i,NR-1,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,NR-1,0)])
						  	   + MESH.SupMP[G(i,NR-1,3)] * UwallsMR[LR(i,NR,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,NR,0)])	  
						  	   )
						  	   //- (0.50*(UwallsMR[LR(i,0,0)] + 0.50*(UwallsMR[LR(i,0,0)] + UwallsMR[LR(i,1,0)])))/(0.50*(MESH.MR[GR(i,0,1)] + 0.50*(MESH.MR[GR(i,0,1)] + MESH.MR[GR(i,1,1)])))
							   )

						  	   + (2.0 / MESH.VolMP[G(i,NR-1,0)]) * (
						  	   + 0.50*MESH.SupMP[G(i,NR-1,0)] * VwallsMU[LU(i,NR-1,0)] * cos(PI + MESH.AngleMU[GU(i,NR-1,0)])
						  	   + 0.50*MESH.SupMP[G(i,NR-1,1)] * VwallsMU[LU(i + 1,NR-1,0)] * cos(MESH.AngleMU[GU(i + 1,NR-1,0)])
						  	   + MESH.SupMP[G(i,NR-1,2)] * VlocalPres[LNH(i,NR-1,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,NR-1,0)])
						  	   + MESH.SupMP[G(i,NR-1,3)] * VwallsMR[LR(i,NR,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,NR,0)])	  
						  	   )
						  	   );

		TauZZ_mr[LR(i, NR, 0)] =  
							   + muTotal[LNH(i,NR - 1,0)]*(
							   + 2.0 * ( (2.0 / MESH.VolMP[G(i,NR-1,0)]) * (
							   + 0.50*MESH.SupMP[G(i,NR-1,0)] * UwallsMU[LU(i,NR - 1,0)] * sin(PI + MESH.AngleMU[GU(i,NR - 1,0)])
							   + 0.50*MESH.SupMP[G(i,NR-1,1)] * UwallsMU[LU(i + 1,NR - 1,0)] * sin(MESH.AngleMU[GU(i + 1,NR - 1,0)])
							   + MESH.SupMP[G(i,NR-1,2)] * UlocalPres[LNH(i,NR - 1,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,NR - 1,0)])
							   + MESH.SupMP[G(i,NR-1,3)] * UwallsMR[LU(i,NR,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	  
							   )
							  
							   )
							   - (2.0/3.0)*DivergenceMR[VR(i - Ix,1,0)]
							   );

		//Parte centro
		for(j = 1; j < NR; j++){
			TauRR_mr[LR(i, j, 0)] = TauRR[LNH(i,j-1,0)] + (MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)])*((TauRR[LNH(i,j,0)] - TauRR[LNH(i,j-1,0)])/(MESH.MP[G(i,j,1)] - MESH.MP[G(i,j-1,1)]));
			TauRZ_mr[LR(i, j, 0)] = TauRZ[LNH(i,j-1,0)] + (MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)])*((TauRZ[LNH(i,j,0)] - TauRZ[LNH(i,j-1,0)])/(MESH.MP[G(i,j,1)] - MESH.MP[G(i,j-1,1)]));
			TauZZ_mr[LR(i, j, 0)] = TauZZ[LNH(i,j-1,0)] + (MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)])*((TauZZ[LNH(i,j,0)] - TauZZ[LNH(i,j-1,0)])/(MESH.MP[G(i,j,1)] - MESH.MP[G(i,j-1,1)]));
		}
	}

}

//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveAxial(Mesher MESH){
int i, j;
	
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			MomentumDiffusiveAxial[LNH(i,j,0)] = 
									   		   + (1.0/MESH.VolMP[G(i,j,0)])*( 
							                   + MESH.SupMP[G(i,j,0)]*(TauZZ_mu[LU(i,j,0)]*cos(PI + MESH.AngleMU[GU(i,j,0)]) + TauRZ_mu[LU(i,j,0)]*sin(PI + MESH.AngleMU[GU(i,j,0)]))
							                   + MESH.SupMP[G(i,j,1)]*(TauZZ_mu[LU(i+1,j,0)]*cos(MESH.AngleMU[GU(i + 1,j,0)]) + TauRZ_mu[LU(i + 1,j,0)]*sin(MESH.AngleMU[GU(i + 1,j,0)]))  
							                   + MESH.SupMP[G(i,j,2)]*(TauZZ_mr[LR(i,j,0)]*cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + TauRZ_mr[LR(i,j,0)]*sin(1.50*PI + MESH.AngleMR[GR(i,j,0)]))
							                   + MESH.SupMP[G(i,j,3)]*(TauZZ_mr[LR(i,j+1,0)]*cos(0.50*PI + MESH.AngleMR[GR(i,j+1,0)]) + TauRZ_mr[LR(i,j+1,0)]*sin(0.50*PI + MESH.AngleMR[GR(i,j+1,0)]))
							                   );
		}
	}	
	 
}

//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveRadial(Mesher MESH){
int i, j;
	
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			MomentumDiffusiveRadial[LNH(i,j,0)] = 
										 		+ (1.0/MESH.VolMP[G(i,j,0)])*(
									     		+ MESH.SupMP[G(i,j,0)]*(TauRZ_mu[LU(i,j,0)]*cos(PI + MESH.AngleMU[GU(i,j,0)]) + TauRR_mu[LU(i,j,0)]*sin(PI + MESH.AngleMU[GU(i,j,0)]))
									     		+ MESH.SupMP[G(i,j,1)]*(TauRZ_mu[LU(i+1,j,0)]*cos(MESH.AngleMU[GU(i+1,j,0)]) + TauRR_mu[LU(i+1,j,0)]*sin(MESH.AngleMU[GU(i+1,j,0)]))  
									     		+ MESH.SupMP[G(i,j,2)]*(TauRZ_mr[LR(i,j,0)]*cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + TauRR_mr[LR(i,j,0)]*sin(1.50*PI + MESH.AngleMR[GR(i,j,0)]))
									     		+ MESH.SupMP[G(i,j,3)]*(TauRZ_mr[LR(i,j+1,0)]*cos(0.50*PI + MESH.AngleMR[GR(i,j+1,0)]) + TauRR_mr[LR(i,j+1,0)]*sin(0.50*PI + MESH.AngleMR[GR(i,j+1,0)]))
									     		)
									     	//	- TauThetaTheta[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
									     		;
		}
	}
}

//Cálculo de la conductividad térmica de los gases en cada nodo y paredes
void Solver::Get_K(Mesher MESH, ParPro MPI1){
int i, j;

	//Cálculo de K en los nodos
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			K[LNH(i,j,0)] = (Cp[LNH(i,j,0)]*muTotal[LNH(i,j,0)])/Pr;
		}
	}

	//Envío de los Halos de la conductividad térmica
	MPI1.SendData(K, Ix, Fx);
	MPI1.ReceiveData(K, Ix, Fx);

	//Cálculo de K en las paredes de los volúmenes de control	
	//Nodos U
	if(Rank == 0){
		for(j = 0; j < NR; j++){
			KwallsMU[LU(0,j,0)] = KLeft[j];
			for(i = Ix + 1; i < Fx + 1; i++){
				KwallsMU[LU(i, j, 0)] = K[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((K[LNH(i,j,0)] - K[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}
	else if(Rank == Procesos - 1){
		for (j = 0; j < NR; j++){
			KwallsMU[LU(NA, j, 0)] = KRight[j];
			for (i = Ix; i < Fx; i++){
				KwallsMU[LU(i, j, 0)] = K[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((K[LNH(i,j,0)] - K[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}
	else{
		for(j = 0; j < NR; j++){
			for(i = Ix; i < Fx + 1; i++){
				KwallsMU[LU(i, j, 0)] = K[LNH(i-1,j,0)] + (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])*((K[LNH(i,j,0)] - K[LNH(i-1,j,0)])/(MESH.MP[G(i,j,0)] - MESH.MP[G(i-1,j,0)]));
			}
		}
	}

	//Nodos V
	for(i = Ix; i < Fx; i++){
		KwallsMR[LR(i, 0, 0)] = Kdown[i - Ix];
		KwallsMR[LR(i, NR, 0)] = KUp[i - Ix];
			for (j = 1; j < NR; j++){
			KwallsMR[LR(i,j,0)] = K[LNH(i,j-1,0)] + (MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)])*((K[LNH(i,j,0)] - K[LNH(i,j-1,0)])/(MESH.MP[G(i,j,1)] - MESH.MP[G(i,j-1,1)]));
		}
	}
	
}

//Cálculo del término difusivo de la ecuación de conservación de la energía
void Solver::Get_EnergyDiffusive(Mesher MESH){
int i, j;
	
	if(Rank == 0){
		for(j = 1; j < NR - 1; j++){
			EnergyDiffusive[LNH(0,j,0)] = 
									    + (1.0/MESH.VolMP[G(0,j,0)]) * (
										+ KwallsMU[LU(0,j,0)] * (MESH.SupMP[G(0,j,0)]/MESH.DeltasMU[GU(0,j,0)]) *	(TlocalPres[LNH(0,j,0)] - Tleft[j]) * cos(PI + MESH.AngleMU[GU(0,j,0)])
										+ KwallsMU[LU(1,j,0)] * (MESH.SupMP[G(0,j,1)]/MESH.DeltasMU[GU(1,j,0)]) *	(TlocalPres[LNH(1,j,0)] - TlocalPres[LNH(0,j,0)]) * cos(MESH.AngleMU[GU(1,j,0)])
										+ KwallsMR[LR(0,j,0)] * (MESH.SupMP[G(0,j,2)]/MESH.DeltasMR[GR(0,j,1)]) *	(TlocalPres[LNH(0,j,0)] - TlocalPres[LNH(0,j - 1,0)]) * sin(1.50*PI + MESH.AngleMR[GR(0,j,0)])
										+ KwallsMR[LR(0,j + 1,0)] * (MESH.SupMP[G(0,j,3)]/MESH.DeltasMR[GR(0,j + 1,1)]) *	(TlocalPres[LNH(0,j + 1,0)] - TlocalPres[LNH(0,j,0)]) * sin(0.50*PI + MESH.AngleMR[GR(0,j + 1,0)])
										);

			for(i = Ix + 1; i < Fx; i++){

				EnergyDiffusive[LNH(i,j,0)] = 
										    + (1.0/MESH.VolMP[G(i,j,0)]) * (
										    + KwallsMU[LU(i,j,0)] * (MESH.SupMP[G(i,j,0)]/MESH.DeltasMU[GU(i,j,0)]) *	(TlocalPres[LNH(i,j,0)] - TlocalPres[LNH(i - 1,j,0)]) * cos(PI + MESH.AngleMU[GU(i,j,0)])
										    + KwallsMU[LU(i + 1,j,0)] * (MESH.SupMP[G(i,j,1)]/MESH.DeltasMU[GU(i + 1,j,0)]) *	(TlocalPres[LNH(i + 1,j,0)] - TlocalPres[LNH(i,j,0)]) * cos(MESH.AngleMU[GU(i + 1,j,0)])
										    + KwallsMR[LR(i,j,0)] * (MESH.SupMP[G(i,j,2)]/MESH.DeltasMR[GR(i,j,1)]) *	(TlocalPres[LNH(i,j,0)] - TlocalPres[LNH(i,j - 1,0)]) * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
										    + KwallsMR[LR(i,j + 1,0)] * (MESH.SupMP[G(i,j,3)]/MESH.DeltasMR[GR(i,j + 1,1)]) *	(TlocalPres[LNH(i,j + 1,0)] - TlocalPres[LNH(i,j,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])
										    );
										
			}
		}
		for(i = Ix + 1; i < Fx; i++){
			EnergyDiffusive[LNH(i,0,0)] =
										+ (1.0/MESH.VolMP[G(i,0,0)]) * (
										+ KwallsMU[LU(i,0,0)] * (MESH.SupMP[G(i,0,0)]/MESH.DeltasMU[GU(i,0,0)]) *	(TlocalPres[LNH(i,0,0)] - TlocalPres[LNH(i - 1,0,0)]) * cos(PI + MESH.AngleMU[GU(i,0,0)])
										+ KwallsMU[LU(i + 1,0,0)] * (MESH.SupMP[G(i,0,1)]/MESH.DeltasMU[GU(i + 1,0,0)]) *	(TlocalPres[LNH(i + 1,0,0)] - TlocalPres[LNH(i,0,0)]) * cos(MESH.AngleMU[GU(i + 1,0,0)])
										+ KwallsMR[LR(i,0,0)] * (MESH.SupMP[G(i,0,2)]/MESH.DeltasMR[GR(i,0,1)]) *	(TlocalPres[LNH(i,0,0)] - Tdown[i - Ix]) * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])	     
										+ KwallsMR[LR(i,1,0)] * (MESH.SupMP[G(i,0,3)]/MESH.DeltasMR[GR(i,1,1)]) *	(TlocalPres[LNH(i,1,0)] - TlocalPres[LNH(i,0,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	     
									    );
			EnergyDiffusive[LNH(i,NR - 1,0)] = 
											 + (1.0/MESH.VolMP[G(i,NR - 1,0)]) * (
										     + KwallsMU[LU(i,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,0)]/MESH.DeltasMU[GU(i,NR - 1,0)]) *	(TlocalPres[LNH(i,NR - 1,0)] - TlocalPres[LNH(i - 1,NR - 1,0)]) * cos(PI + MESH.AngleMU[GU(i,NR - 1,0)])
										     + KwallsMU[LU(i + 1,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,1)]/MESH.DeltasMU[GU(i + 1,NR - 1,0)]) *	(TlocalPres[LNH(i + 1,NR - 1,0)] - TlocalPres[LNH(i,NR - 1,0)]) * cos(MESH.AngleMU[GU(i + 1,NR - 1,0)])
										     + KwallsMR[LR(i,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,2)]/MESH.DeltasMR[GR(i,NR - 1,1)]) *	(TlocalPres[LNH(i,NR - 1,0)] - TlocalPres[LNH(i,NR - 2,0)]) * sin(1.50*PI + MESH.AngleMR[GR(i,NR - 1,0)])
										     + KwallsMR[LR(i,NR,0)] * (MESH.SupMP[G(i,NR - 1,3)]/MESH.DeltasMR[GR(i,NR,1)]) *	(Tup[i - Ix] - TlocalPres[LNH(i,NR - 1,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,NR,0)])
									  		 );
		}
		EnergyDiffusive[LNH(0,0,0)] = 
									+ (1.0/MESH.VolMP[G(0,0,0)]) * (
									+ KwallsMU[LU(0,0,0)] * (MESH.SupMP[G(0,0,0)]/MESH.DeltasMU[GU(0,0,0)]) *	(TlocalPres[LNH(0,0,0)] - Tleft[0]) * cos(PI + MESH.AngleMU[GU(0,0,0)])
									+ KwallsMU[LU(1,0,0)] * (MESH.SupMP[G(0,0,1)]/MESH.DeltasMU[GU(1,0,0)]) *	(TlocalPres[LNH(1,0,0)] - TlocalPres[LNH(0,0,0)]) * cos(MESH.AngleMU[GU(1,0,0)])
									+ KwallsMR[LR(0,0,0)] * (MESH.SupMP[G(0,0,2)]/MESH.DeltasMR[GR(0,0,1)]) *	(TlocalPres[LNH(0,0,0)] - Tdown[0]) * sin(1.50*PI + MESH.AngleMR[GR(0,0,0)])	     
									+ KwallsMR[LR(0,1,0)] * (MESH.SupMP[G(0,0,3)]/MESH.DeltasMR[GR(0,1,1)]) *	(TlocalPres[LNH(0,1,0)] - TlocalPres[LNH(0,0,0)]) * sin(0.50*PI + MESH.AngleMR[GR(0,1,0)])				   
								    );
		EnergyDiffusive[LNH(0,NR - 1,0)] = 
										 + (1.0/MESH.VolMP[G(0,NR - 1,0)]) * (
										 + KwallsMU[LU(0,NR - 1,0)] * (MESH.SupMP[G(0,NR - 1,0)]/MESH.DeltasMU[GU(0,NR - 1,0)]) *	(TlocalPres[LNH(0,NR - 1,0)] - Tleft[NR - 1]) * cos(PI + MESH.AngleMU[GU(0,NR - 1,0)])
										 + KwallsMU[LU(1,NR - 1,0)] * (MESH.SupMP[G(0,NR - 1,1)]/MESH.DeltasMU[GU(1,NR - 1,0)]) *	(TlocalPres[LNH(1,NR - 1,0)] - TlocalPres[LNH(0,NR - 1,0)]) * cos(MESH.AngleMU[GU(1,NR - 1,0)])
										 + KwallsMR[LR(0,NR - 1,0)] * (MESH.SupMP[G(0,NR - 1,2)]/MESH.DeltasMR[GR(0,NR - 1,1)]) *	(TlocalPres[LNH(0,NR - 1,0)] - TlocalPres[LNH(0,NR - 2,0)]) * sin(1.50*PI + MESH.AngleMR[GR(0,NR - 1,0)])
										 + KwallsMR[LR(0,NR,0)] * (MESH.SupMP[G(0,NR - 1,3)]/MESH.DeltasMR[GR(0,NR,1)]) *	(Tup[0] - TlocalPres[LNH(0,NR - 1,0)]) * sin(0.50*PI + MESH.AngleMR[GR(0,NR,0)])
										 );
	}
	else if(Rank == Procesos - 1){
		for(j = 1; j < NR - 1; j++){
			EnergyDiffusive[LNH(NA - 1,j,0)] = 
											 + (1.0/MESH.VolMP[G(NA - 1,j,0)]) * (
										     + KwallsMU[LU(NA - 1,j,0)] * (MESH.SupMP[G(NA - 1,j,0)]/MESH.DeltasMU[GU(NA - 1,j,0)]) *	(TlocalPres[LNH(NA - 1,j,0)] - TlocalPres[LNH(NA - 2,j,0)]) * cos(PI + MESH.AngleMU[GU(NA - 1,j,0)])
										     + KwallsMU[LU(NA,j,0)] * (MESH.SupMP[G(NA - 1,j,1)]/MESH.DeltasMU[GU(NA,j,0)]) *	(Tright[j] - TlocalPres[LNH(NA - 1,j,0)]) * cos(MESH.AngleMU[GU(NA,j,0)])
										     + KwallsMR[LR(NA - 1,j,0)] * (MESH.SupMP[G(NA - 1,j,2)]/MESH.DeltasMR[GR(NA - 1,j,1)]) *	(TlocalPres[LNH(NA - 1,j,0)] - TlocalPres[LNH(NA - 1,j - 1,0)]) * sin(1.50*PI + MESH.AngleMR[GR(NA - 1,j,0)])
										     + KwallsMR[LR(NA - 1,j + 1,0)] * (MESH.SupMP[G(NA - 1,j,3)]/MESH.DeltasMR[GR(NA - 1,j + 1,1)]) *	(TlocalPres[LNH(NA - 1,j + 1,0)] - TlocalPres[LNH(NA - 1,j,0)]) * sin(0.50*PI + MESH.AngleMR[GR(NA - 1,j + 1,0)])
										     );
			for(i = Ix; i < Fx - 1; i++){
				EnergyDiffusive[LNH(i,j,0)] = 
										    + (1.0/MESH.VolMP[G(i,j,0)]) * (
										    + KwallsMU[LU(i,j,0)] * (MESH.SupMP[G(i,j,0)]/MESH.DeltasMU[GU(i,j,0)]) *	(TlocalPres[LNH(i,j,0)] - TlocalPres[LNH(i - 1,j,0)]) * cos(PI + MESH.AngleMU[GU(i,j,0)])
										    + KwallsMU[LU(i + 1,j,0)] * (MESH.SupMP[G(i,j,1)]/MESH.DeltasMU[GU(i + 1,j,0)]) *	(TlocalPres[LNH(i + 1,j,0)] - TlocalPres[LNH(i,j,0)]) * cos(MESH.AngleMU[GU(i + 1,j,0)])
										    + KwallsMR[LR(i,j,0)] * (MESH.SupMP[G(i,j,2)]/MESH.DeltasMR[GR(i,j,1)]) *	(TlocalPres[LNH(i,j,0)] - TlocalPres[LNH(i,j - 1,0)]) * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
										    + KwallsMR[LR(i,j + 1,0)] * (MESH.SupMP[G(i,j,3)]/MESH.DeltasMR[GR(i,j + 1,1)]) *	(TlocalPres[LNH(i,j + 1,0)] - TlocalPres[LNH(i,j,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])
										    );
			}
		}	
		for(i = Ix; i < Fx - 1; i++){
			EnergyDiffusive[LNH(i,0,0)] = 
										+ (1.0/MESH.VolMP[G(i,0,0)]) * (
										+ KwallsMU[LU(i,0,0)] * (MESH.SupMP[G(i,0,0)]/MESH.DeltasMU[GU(i,0,0)]) *	(TlocalPres[LNH(i,0,0)] - TlocalPres[LNH(i - 1,0,0)]) * cos(PI + MESH.AngleMU[GU(i,0,0)])
										+ KwallsMU[LU(i + 1,0,0)] * (MESH.SupMP[G(i,0,1)]/MESH.DeltasMU[GU(i + 1,0,0)]) *	(TlocalPres[LNH(i + 1,0,0)] - TlocalPres[LNH(i,0,0)]) * cos(MESH.AngleMU[GU(i + 1,0,0)])
										+ KwallsMR[LR(i,0,0)] * (MESH.SupMP[G(i,0,2)]/MESH.DeltasMR[GR(i,0,1)]) *	(TlocalPres[LNH(i,0,0)] - Tdown[i - Ix]) * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])	     
										+ KwallsMR[LR(i,1,0)] * (MESH.SupMP[G(i,0,3)]/MESH.DeltasMR[GR(i,1,1)]) *	(TlocalPres[LNH(i,1,0)] - TlocalPres[LNH(i,0,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	   
									    );
			EnergyDiffusive[LNH(i,NR - 1,0)] = 
											 + (1.0/MESH.VolMP[G(i,NR - 1,0)]) * (
										     + KwallsMU[LU(i,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,0)]/MESH.DeltasMU[GU(i,NR - 1,0)]) *	(TlocalPres[LNH(i,NR - 1,0)] - TlocalPres[LNH(i - 1,NR - 1,0)]) * cos(PI + MESH.AngleMU[GU(i,NR - 1,0)])
										     + KwallsMU[LU(i + 1,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,1)]/MESH.DeltasMU[GU(i + 1,NR - 1,0)]) *	(TlocalPres[LNH(i + 1,NR - 1,0)] - TlocalPres[LNH(i,NR - 1,0)]) * cos(MESH.AngleMU[GU(i + 1,NR - 1,0)])
										     + KwallsMR[LR(i,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,2)]/MESH.DeltasMR[GR(i,NR - 1,1)]) *	(TlocalPres[LNH(i,NR - 1,0)] - TlocalPres[LNH(i,NR - 2,0)]) * sin(1.50*PI + MESH.AngleMR[GR(i,NR - 1,0)])
										     + KwallsMR[LR(i,NR,0)] * (MESH.SupMP[G(i,NR - 1,3)]/MESH.DeltasMR[GR(i,NR,1)]) *	(Tup[i - Ix] - TlocalPres[LNH(i,NR - 1,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,NR,0)])
										     );
		}
		EnergyDiffusive[LNH(NA - 1,0,0)] = 
										 + (1.0/MESH.VolMP[G(NA - 1,0,0)]) * (
										 + KwallsMU[LU(NA - 1,0,0)] * (MESH.SupMP[G(NA - 1,0,0)]/MESH.DeltasMU[GU(NA - 1,0,0)]) *	(TlocalPres[LNH(NA - 1,0,0)] - TlocalPres[LNH(NA - 2,0,0)]) * cos(PI + MESH.AngleMU[GU(NA - 1,0,0)])
										 + KwallsMU[LU(NA,0,0)] * (MESH.SupMP[G(NA - 1,0,1)]/MESH.DeltasMU[GU(NA,0,0)]) *	(Tright[0] - TlocalPres[LNH(NA - 1,0,0)]) * cos(MESH.AngleMU[GU(NA,0,0)])
										 + KwallsMR[LR(NA - 1,0,0)] * (MESH.SupMP[G(NA - 1,0,2)]/MESH.DeltasMR[GR(NA - 1,0,1)]) *	(TlocalPres[LNH(NA - 1,0,0)] - Tdown[NA - 1 - Ix]) * sin(1.50*PI + MESH.AngleMR[GR(NA - 1,0,0)])	     
										 + KwallsMR[LR(NA - 1,1,0)] * (MESH.SupMP[G(NA - 1,0,3)]/MESH.DeltasMR[GR(NA - 1,1,1)]) *	(TlocalPres[LNH(NA - 1,1,0)] - TlocalPres[LNH(NA - 1,0,0)]) * sin(0.50*PI + MESH.AngleMR[GR(NA - 1,1,0)])
										 );
		EnergyDiffusive[LNH(NA - 1, NR - 1,0)] = 
											   + (1.0/MESH.VolMP[G(NA - 1,NR - 1,0)]) * (
										       + KwallsMU[LU(NA - 1,NR - 1,0)] * (MESH.SupMP[G(NA - 1,NR - 1,0)]/MESH.DeltasMU[GU(NA - 1,NR - 1,0)]) *	(TlocalPres[LNH(NA - 1,NR - 1,0)] - TlocalPres[LNH(NA - 2,NR - 1,0)]) * cos(PI + MESH.AngleMU[GU(NA - 1,NR - 1,0)])
										       + KwallsMU[LU(NA,NR - 1,0)] * (MESH.SupMP[G(NA - 1,NR - 1,1)]/MESH.DeltasMU[GU(NA,NR - 1,0)]) *	(Tright[NR - 1] - TlocalPres[LNH(NA - 1,NR - 1,0)]) * cos(MESH.AngleMU[GU(NA,NR - 1,0)])
										       + KwallsMR[LR(NA - 1,NR - 1,0)] * (MESH.SupMP[G(NA - 1,NR - 1,2)]/MESH.DeltasMR[GR(NA - 1,NR - 1,1)]) *	(TlocalPres[LNH(NA - 1,NR - 1,0)] - TlocalPres[LNH(NA - 1,NR - 2,0)]) * sin(1.50*PI + MESH.AngleMR[GR(NA - 1,NR - 1,0)])
										       + KwallsMR[LR(NA - 1,NR,0)] * (MESH.SupMP[G(NA - 1,NR - 1,3)]/MESH.DeltasMR[GR(NA - 1,NR,1)]) *	(Tup[i - Ix] - TlocalPres[LNH(NA - 1,NR - 1,0)]) * sin(0.50*PI + MESH.AngleMR[GR(NA - 1,NR,0)])
										       );
	}
	else{
		for(i = Ix; i < Fx; i++){
			EnergyDiffusive[LNH(i,0,0)] = 
										+ (1.0/MESH.VolMP[G(i,0,0)]) * (
										+ KwallsMU[LU(i,0,0)] * (MESH.SupMP[G(i,0,0)]/MESH.DeltasMU[GU(i,0,0)]) *	(TlocalPres[LNH(i,0,0)] - TlocalPres[LNH(i - 1,0,0)]) * cos(PI + MESH.AngleMU[GU(i,0,0)])
										+ KwallsMU[LU(i + 1,0,0)] * (MESH.SupMP[G(i,0,1)]/MESH.DeltasMU[GU(i + 1,0,0)]) *	(TlocalPres[LNH(i + 1,0,0)] - TlocalPres[LNH(i,0,0)]) * cos(MESH.AngleMU[GU(i + 1,0,0)])
										+ KwallsMR[LR(i,0,0)] * (MESH.SupMP[G(i,0,2)]/MESH.DeltasMR[GR(i,0,1)]) *	(TlocalPres[LNH(i,0,0)] - Tdown[i - Ix]) * sin(1.50*PI + MESH.AngleMR[GR(i,0,0)])	    
										+ KwallsMR[LR(i,1,0)] * (MESH.SupMP[G(i,0,3)]/MESH.DeltasMR[GR(i,1,1)]) *	(TlocalPres[LNH(i,1,0)] - TlocalPres[LNH(i,0,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,1,0)])	    
									    );
			EnergyDiffusive[LNH(i,NR - 1,0)] = 
											 + (1.0/MESH.VolMP[G(i,NR - 1,0)]) * (
										     + KwallsMU[LU(i,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,0)]/MESH.DeltasMU[GU(i,NR - 1,0)]) *	(TlocalPres[LNH(i,NR - 1,0)] - TlocalPres[LNH(i - 1,NR - 1,0)]) * cos(PI + MESH.AngleMU[GU(i,NR - 1,0)])
										     + KwallsMU[LU(i + 1,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,1)]/MESH.DeltasMU[GU(i + 1,NR - 1,0)]) * (TlocalPres[LNH(i + 1,NR - 1,0)] - TlocalPres[LNH(i,NR - 1,0)]) * cos(MESH.AngleMU[GU(i + 1,NR - 1,0)])
										     + KwallsMR[LR(i,NR - 1,0)] * (MESH.SupMP[G(i,NR - 1,2)]/MESH.DeltasMR[GR(i,NR - 1,1)]) *	(TlocalPres[LNH(i,NR - 1,0)] - TlocalPres[LNH(i,NR - 2,0)]) * sin(1.50*PI + MESH.AngleMR[GR(i,NR - 1,0)])
										     + KwallsMR[LR(i,NR,0)] * (MESH.SupMP[G(i,NR - 1,3)]/MESH.DeltasMR[GR(i,NR,1)]) *	(Tup[i - Ix] - TlocalPres[LNH(i,NR - 1,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,NR,0)])    
										     );
			for(j = 1; j < NR - 1; j++){
				EnergyDiffusive[LNH(i,j,0)] = 
											+ (1.0/MESH.VolMP[G(i,j,0)]) * (
										    + KwallsMU[LU(i,j,0)] * (MESH.SupMP[G(i,j,0)]/MESH.DeltasMU[GU(i,j,0)]) *	(TlocalPres[LNH(i,j,0)] - TlocalPres[LNH(i - 1,j,0)]) * cos(PI + MESH.AngleMU[GU(i,j,0)])
										    + KwallsMU[LU(i + 1,j,0)] * (MESH.SupMP[G(i,j,1)]/MESH.DeltasMU[GU(i + 1,j,0)]) *	(TlocalPres[LNH(i + 1,j,0)] - TlocalPres[LNH(i,j,0)]) * cos(MESH.AngleMU[GU(i + 1,j,0)])
										    + KwallsMR[LR(i,j,0)] * (MESH.SupMP[G(i,j,2)]/MESH.DeltasMR[GR(i,j,1)]) *	(TlocalPres[LNH(i,j,0)] - TlocalPres[LNH(i,j - 1,0)]) * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
										    + KwallsMR[LR(i,j + 1,0)] * (MESH.SupMP[G(i,j,3)]/MESH.DeltasMR[GR(i,j + 1,1)]) *	(TlocalPres[LNH(i,j + 1,0)] - TlocalPres[LNH(i,j,0)]) * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])
										    );
			}
		}
	}
	
}

void Solver::Get_EnergyWalls(Mesher MESH){
int i, j;

		//Nodos U
		if(Rank == 0){
			for(j = 0; j < NR; j++){
				EwallsMU[LU(0,j,0)] = Eleft[j];
				for(i = Ix + 1; i < Fx + 1; i++){
					EwallsMU[LU(i, j, 0)] = sqrt(ElocalPres[LNH(i, j, 0)] * sqrt(ElocalPres[LNH(i - 1, j, 0)]));
				}
			}
		}
		else if(Rank == Procesos - 1){
			for(j = 0; j < NR; j++){
				EwallsMU[LU(NA, j, 0)] = Eright[j];
				for (i = Ix + 1; i < NA; i++){
					EwallsMU[LU(i, j, 0)] = sqrt(ElocalPres[LNH(i, j, 0)] * sqrt(ElocalPres[LNH(i - 1, j, 0)]));
				}
			}
		}
		else{
			for(j = 0; j < NR; j++){
				for(i = Ix; i < Fx + 1; i++){
					EwallsMU[LU(i, j, 0)] = sqrt(ElocalPres[LNH(i, j, 0)] * sqrt(ElocalPres[LNH(i - 1, j, 0)]));
				}
			}
		}
		
		//Nodos R
		for(i = Ix; i < Fx + 1; i++){
			EwallsMR[LR(i, 0, 0)] = ElocalPres[LNH(i,0,0)];
			EwallsMR[LR(i, NR, 0)] = Eup[i - Ix];
			for(j = 1; j < NR; j++){
				EwallsMR[LR(i, j, 0)] = sqrt(ElocalPres[LNH(i, j, 0)] * sqrt(ElocalPres[LNH(i, j - 1, 0)]));
			}
		}

}

//Cálculo del término convectivo de la energía
void Solver::Get_EnergyConvective(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			EnergyConvective[LNH(i,j,0)] = 
										 - (1.0/MESH.VolMP[G(i,j,0)])*(
								  		 + MESH.SupMP[G(i,j,0)]*(UwallsMU[LU(i,j,0)]*cos(PI + MESH.AngleMU[GU(i,j,0)]) + VwallsMU[LU(i,j,0)]*sin(PI + MESH.AngleMU[GU(i,j,0)]))*RhoWallsMU[LU(i,j,0)]*EwallsMU[LU(i,j,0)]
								  		 + MESH.SupMP[G(i,j,1)]*(UwallsMU[LU(i + 1,j,0)]*cos(MESH.AngleMU[GU(i + 1,j,0)]) + VwallsMU[LU(i + 1,j,0)]*sin(MESH.AngleMU[GU(i + 1,j,0)]))*RhoWallsMU[LU(i + 1,j,0)]*EwallsMU[LU(i + 1,j,0)]
								  		 + MESH.SupMP[G(i,j,2)]*(UwallsMR[LR(i,j,0)]*cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + VwallsMR[LR(i,j,0)]*sin(1.50*PI + MESH.AngleMR[GR(i,j,0)]))*RhoWallsMR[LR(i,j,0)]*EwallsMR[LR(i,j,0)]
								  		 + MESH.SupMP[G(i,j,3)]*(UwallsMR[LR(i,j + 1,0)]*cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]) + VwallsMR[LR(i,j + 1,0)]*sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]))*RhoWallsMR[LR(i,j + 1,0)]*EwallsMR[LR(i,j + 1,0)]
								 		 );
		}
	}

}

//Cálculo del término de presión en la ecuación de conservación de la energía
void Solver::Get_EnergyPressureTerm(Mesher MESH){
int i, j;
		
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NR; j++){
				EnergyPressureTerm[LNH(i, j, 0)] = 
												 - Pressure[LNH(i,j,0)]*(1.0/MESH.VolMP[G(i,j,0)])*(
												 + MESH.SupMP[G(i,j,0)]*(UwallsMU[LU(i,j,0)]*cos(PI + MESH.AngleMU[GU(i,j,0)]) + VwallsMU[LU(i,j,0)]*sin(PI + MESH.AngleMU[GU(i,j,0)]))
								  		 		 + MESH.SupMP[G(i,j,1)]*(UwallsMU[LU(i + 1,j,0)]*cos(MESH.AngleMU[GU(i + 1,j,0)]) + VwallsMU[LU(i + 1,j,0)]*sin(MESH.AngleMU[GU(i + 1,j,0)]))
								  		 		 + MESH.SupMP[G(i,j,2)]*(UwallsMR[LR(i,j,0)]*cos(1.50*PI + MESH.AngleMR[GR(i,j,0)]) + VwallsMR[LR(i,j,0)]*sin(1.50*PI + MESH.AngleMR[GR(i,j,0)]))
								  		 		 + MESH.SupMP[G(i,j,3)]*(UwallsMR[LR(i,j + 1,0)]*cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]) + VwallsMR[LR(i,j + 1,0)]*sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)]))
								 				 );							
			}
		}

	
}

//Cálculo del término viscoso de la ecuación de conservación de la energía
void Solver::Get_EnergyViscousTerm(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			EnergyViscous[LNH(i,j,0)] = 
									  + TauRR[LNH(i,j,0)]*(
									  + (1.0 / MESH.VolMP[G(i,j,0)]) * (
							  		  + MESH.SupMP[G(i,j,0)] * VwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])
							  		  + MESH.SupMP[G(i,j,1)] * VwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)])
							  		  + MESH.SupMP[G(i,j,2)] * VwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  		  + MESH.SupMP[G(i,j,3)] * VwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  		  )
							  		  - VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
									  )
									
									  + TauRZ[LNH(i,j,0)]*(
									  + (1.0 / MESH.VolMP[G(i,j,0)]) * (
							  		  + MESH.SupMP[G(i,j,0)] * VwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)])
							  		  + MESH.SupMP[G(i,j,1)] * VwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)])
							  		  + MESH.SupMP[G(i,j,2)] * VwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  		  + MESH.SupMP[G(i,j,3)] * VwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  		  )
									  )

									  + TauThetaTheta[LNH(i,j,0)]*(VlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,0)])
							    
									  + TauRZ[LNH(i,j,0)]*(
									  + (1.0 / MESH.VolMP[G(i,j,0)]) * (
							  		  + MESH.SupMP[G(i,j,0)] * UwallsMU[LU(i,j,0)] * sin(PI + MESH.AngleMU[GU(i,j,0)])
							  		  + MESH.SupMP[G(i,j,1)] * UwallsMU[LU(i + 1,j,0)] * sin(MESH.AngleMU[GU(i + 1,j,0)])
							  		  + MESH.SupMP[G(i,j,2)] * UwallsMR[LR(i,j,0)] * sin(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  		  + MESH.SupMP[G(i,j,3)] * UwallsMR[LR(i,j + 1,0)] * sin(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  		  )
							  		  - UlocalPres[LNH(i,j,0)]/MESH.MP[G(i,j,1)]
									  )

									  + TauZZ[LNH(i,j,0)]*(
									  + (1.0 / MESH.VolMP[G(i,j,0)]) * (
							  		  + MESH.SupMP[G(i,j,0)] * UwallsMU[LU(i,j,0)] * cos(PI + MESH.AngleMU[GU(i,j,0)])
							  		  + MESH.SupMP[G(i,j,1)] * UwallsMU[LU(i + 1,j,0)] * cos(MESH.AngleMU[GU(i + 1,j,0)])
							  		  + MESH.SupMP[G(i,j,2)] * UwallsMR[LR(i,j,0)] * cos(1.50*PI + MESH.AngleMR[GR(i,j,0)])
							  		  + MESH.SupMP[G(i,j,3)] * UwallsMR[LR(i,j + 1,0)] * cos(0.50*PI + MESH.AngleMR[GR(i,j + 1,0)])	  
							  		  )
									  );
		}
	}
	
}

//Cálculo de las velocidades del step siguiente
void Solver::Get_Velocities(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			FUPres[LNH(i,j,0)] = 
							   + MomentumConvectiveAxial[LNH(i,j,0)]
							   + MomentumDiffusiveAxial[LNH(i,j,0)]
							   + PressureGradientAxial[LNH(i,j,0)]  
							   ;
					
			FVPres[LNH(i,j,0)] = 
							   + MomentumConvectiveRadial[LNH(i,j,0)]
							   + MomentumDiffusiveRadial[LNH(i,j,0)]
							   + PressureGradientRadial[LNH(i,j,0)]  
							   ;
		}
	}

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){

			//Cálculo de velocidades axiales (U)
			UlocalFut[LNH(i,j,0)] = (2.0*TimeBetta*UlocalPres[LNH(i,j,0)] - (TimeBetta - 0.50)*UlocalPrev[LNH(i,j,0)])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FUPres[LNH(i,j,0)]/RhoLocalPres[LNH(i,j,0)]) - TimeBetta*(FUPrev[LNH(i,j,0)]/RhoLocalPrev[LNH(i,j,0)]));
			
			//Cálculo de velocidades radiales (V)
			VlocalFut[LNH(i,j,0)] = (2.0*TimeBetta*VlocalPres[LNH(i,j,0)] - (TimeBetta - 0.50)*VlocalPrev[LNH(i,j,0)])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FVPres[LNH(i,j,0)]/RhoLocalPres[LNH(i,j,0)]) - TimeBetta*(FVPrev[LNH(i,j,0)]/RhoLocalPrev[LNH(i,j,0)]));
			
		}
		//UlocalFut[LNH(i,0,0)] = (2.0*TimeBetta*UlocalPres[LNH(i,0,0)] - (TimeBetta - 0.50)*UlocalPrev[LNH(i,0,0)])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FUPres[LNH(i,0,0)]/RhoLocalPres[LNH(i,0,0)]) - TimeBetta*(FUPrev[LNH(i,0,0)]/RhoLocalPrev[LNH(i,0,0)]));
	//	VlocalFut[LNH(i,0,0)] = VlocalFut[LNH(i,1,0)];
	}

}

//Cálculo del mapa de densidades futuro
void Solver::Get_Densities(){
int i, j;
	
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			RhoLocalFut[LNH(i,j,0)] = (2.0*TimeBetta*RhoLocalPres[LNH(i,j,0)] - (TimeBetta - 0.50)*RhoLocalPrev[LNH(i,j,0)])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*FRhoPres[LNH(i,j,0)] - TimeBetta*FRhoPrev[LNH(i,j,0)]);
		}	
	}

}

//Cálculo de la contribución total de energía (Fe)
void Solver::Get_Temperatures(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			FEpres[LNH(i,j,0)] =
							   + EnergyConvective[LNH(i,j,0)]
						       + EnergyDiffusive[LNH(i,j,0)]
						       + EnergyPressureTerm[LNH(i,j,0)]
					           + EnergyViscous[LNH(i,j,0)]
							   ;
		}
	}

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			ElocalFut[LNH(i,j,0)] = (2.0*TimeBetta*ElocalPres[LNH(i,j,0)] - (TimeBetta - 0.50)*ElocalPrev[LNH(i,j,0)])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FEpres[LNH(i,j,0)]/RhoLocalPres[LNH(i,j,0)]) - TimeBetta*(FEprev[LNH(i,j,0)]/RhoLocalPrev[LNH(i,j,0)]));
		}
	}	

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			TlocalFut[LNH(i,j,0)] = TlocalPres[LNH(i,j,0)] + (ElocalFut[LNH(i,j,0)] - ElocalPres[LNH(i,j,0)])/(Cp[LNH(i,j,0)]/Gamma + 1e-10);
		}
	}

}

//Cálculo de la diferencia de resultados entre Steps
void Solver::Get_Stop(ParPro MPI1){
int i, j;
MaxDifference = 0.0;

	if(Rank == 0){
		for(i = 0; i < NA; i++){
			for(j = 0; j < NR; j++){

				//Mapa de densidades
				if(abs(RhoGlobalFut[G(i,j,0)] - RhoGlobalPres[G(i,j,0)])/(RhoGlobalPres[G(i,j,0)] + 1e-10) >= MaxDifference){ 
					MaxDifference = abs(RhoGlobalFut[G(i,j,0)] - RhoGlobalPres[G(i,j,0)])/(RhoGlobalPres[G(i,j,0)] + 1e-10); 
				}
				
				//Mapa de temperaturas
				if(abs(TglobalFut[G(i,j,0)] - TglobalPres[G(i,j,0)])/(TglobalPres[G(i,j,0)] + 1e-10) >= MaxDifference){ 
					MaxDifference = abs(TglobalFut[G(i,j,0)] - TglobalPres[G(i,j,0)])/(TglobalPres[G(i,j,0)] + 1e-10); 
				}	
	
				//Mapa de velocidades axiales (U)
				if(abs(UglobalFut[G(i,j,0)] - UglobalPres[G(i,j,0)])/(UglobalPres[G(i,j,0)] + 1e-10) >= MaxDifference){ 
					MaxDifference = abs(UglobalFut[G(i,j,0)] - UglobalPres[G(i,j,0)])/(UglobalPres[G(i,j,0)] + 1e-10); 
				}	
	
				//Mapa de velocidades radiales (V)
				if(abs(VglobalFut[G(i,j,0)] - VglobalPres[G(i,j,0)])/(VglobalPres[G(i,j,0)] + 1e-10) >= MaxDifference){ 
					MaxDifference = abs(VglobalFut[G(i,j,0)] - VglobalPres[G(i,j,0)])/(VglobalPres[G(i,j,0)] + 1e-10); 
				}
	
			}
		}
	}

	MPI1.SendDataToAll(MaxDifference, MaxDifference);
	
}

void Solver::UpdateParameters(){
int i, j;
	
	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){

			//Densidades
			RhoLocalPrev[LNH(i,j,0)] = RhoLocalPres[LNH(i,j,0)];
			RhoLocalPres[LNH(i,j,0)] = RhoLocalFut[LNH(i,j,0)];
			FRhoPrev[LNH(i,j,0)] = FRhoPres[LNH(i,j,0)];

			//Velocidades axiales (U)
			UlocalPrev[LNH(i,j,0)] = UlocalPres[LNH(i,j,0)];
			UlocalPres[LNH(i,j,0)] = UlocalFut[LNH(i,j,0)];
			FUPrev[LNH(i,j,0)] = FUPres[LNH(i,j,0)];

			//Velocidades radiales (V)
			VlocalPrev[LNH(i,j,0)] = VlocalPres[LNH(i,j,0)];
			VlocalPres[LNH(i,j,0)] = VlocalFut[LNH(i,j,0)];
			FVPrev[LNH(i,j,0)] = FVPres[LNH(i,j,0)];

			//Temperaturas
			TlocalPres[LNH(i,j,0)] = TlocalFut[LNH(i,j,0)];

			//Energías Internas
			ElocalPrev[LNH(i,j,0)] = ElocalPres[LNH(i,j,0)];
			ElocalPres[LNH(i,j,0)] = ElocalFut[LNH(i,j,0)];
			FEprev[LNH(i,j,0)] = FEpres[LNH(i,j,0)];
			
		}
	}
	
}

//Sumar las contribuciones de todas las viscosidades
void Solver::Get_TotalViscosity(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NR; j++){
			muTotal[LNH(i,j,0)] = muBase[LNH(i,j,0)];
		}
	}

}

//Pasar los resultados a un archivo VTK en 2D
void Solver::EscalarVTK2D(Mesher MESH, string Carpeta, string Variable, string NombreFile, double *MC, double *GlobalField, int Na, int Nr){
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
		file<<"POINT_DATA"<<"   "<<Na*Nr<<endl;
        file<<"SCALARS "<<Variable<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
        file<<endl;

       for(j = 0; j < Nr; j++){
			for(i = 0; i < Na; i++){
				file<<GlobalField[G(i,j,0)]<<" ";
		    }
	    }

    file.close();

}
//Pasar a un .vtk los resultados de campos vectoriales
void Solver::VectorialVTK2D(Mesher MESH, string Carpeta, string Variable, string NombreFile, double *MC, double *GlobalField1, double *GlobalField2, int Na, int Nr){
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
        file<<"POINT_DATA"<<"   "<<Na*Nr<<endl;
        file<<"VECTORS "<<Variable<<" double"<<endl;
        file<<endl;

     	for(j = 0; j < Nr; j++){
			for(i = 0; i < Na; i++){
		        file<<GlobalField1[G(i,j,0)]<<"   "<<GlobalField2[G(i,j,0)]<<"   "<<"0.0"<<endl;
            }
        }

    file.close();
}

//Ejecución de todos los procesos del solver
void Solver::ExecuteSolver(Memory M1, ParPro MPI1, Mesher MESH){
int i, j;
int Step = 0;
double Time = 0.0;
char FileName[300];
MaxDifference = 2.0*Convergencia;

	//Alojamiento de la memoria de todas las matrices utilizadas
	AllocateMatrix(M1);
	
	//Seteos Iniciales
	InitializeFields(MESH); 
	InitialF();
	
	//Mandar todas las matrices al ZERO
	//Matrices Presentes
	MPI1.SendMatrixToZero(RhoLocalPres, RhoGlobalPres, NA, NR, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(UlocalPres, UglobalPres, NA, NR, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(VlocalPres, VglobalPres, NA, NR, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(TlocalPres, TglobalPres, NA, NR, Procesos, Ix, Fx);

	if(Rank == 0){

		//Pasar los resultados a archivos .VTK Step Inicial
		sprintf(FileName, "MapaDensidades_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Densidades", FileName, MESH.MP, RhoGlobalPres, NA, NR);

		sprintf(FileName, "MapaPresiones_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Presión", FileName, MESH.MP, PresionGlobal, NA, NR);

		sprintf(FileName, "MapaTemperaturas_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Temperaturas", FileName, MESH.MP, TglobalPres, NA, NR);

		sprintf(FileName, "MapaVelocidades_Step_%d", Step);
		VectorialVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Velocidades", FileName, MESH.MP, UglobalPres, VglobalPres, NA, NR);

	}

	while(MaxDifference >= Convergencia){

		UpdateBoundaryConditions(MESH); //Actualización Condiciones de Contorno

		Get_FluidViscosity(); //Cálculo de la viscosidad base del fluido
		Get_TotalViscosity(); //Cálculo de la viscosidad total del fluido para la simulación

		Get_StepTime(MESH, MPI1); //Cálculo del valor del Time Step para todos los procesos
		Step += 1;
		Time += DeltaT;

		//Comunicación Velocidades Axiales
		MPI1.SendData(UlocalPres, Ix, Fx);
		MPI1.ReceiveData(UlocalPres, Ix, Fx);
	//	MPI_Barrier(MPI_COMM_WORLD);		

		//Comunicación Velocidades Radiales
		MPI1.SendData(VlocalPres, Ix, Fx);
		MPI1.ReceiveData(VlocalPres, Ix, Fx);
	//	MPI_Barrier(MPI_COMM_WORLD);	

		//Comunicación Densidades
		MPI1.SendData(RhoLocalPres, Ix, Fx);
		MPI1.ReceiveData(RhoLocalPres, Ix, Fx);
	//	MPI_Barrier(MPI_COMM_WORLD);	

		Get_VelocityWalls(MESH); //Cálculo de las velocidades en las paredes de los volúmenes de control
		Get_DensityWalls(MESH); //Cálculo de las densidades en las paredes de los volúmenes de control

		//Cálculo del campo de densidades
		Get_DensityConvective(MESH); //Cálculo de las contribuciones de la ecuación de continuidad
		Get_Densities(); //Cálculo del mapa de densidades futuro

		//Cálculo del campo de presiones
		//Get_Pressure();	//Cálculo del mapa de presiones con la Ley de los Gases Ideales

		//Comunicación presiones
		//MPI1.SendData(Pressure, Ix, Fx);
		//MPI1.ReceiveData(Pressure, Ix, Fx);

		//Get_PressureWalls(MESH); //Cálculo de las presiones en las paredes de los volúmenes de control
		
		//Cálculo del campo de velocidades
		Get_MomentumConvectiveAxial(MESH); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
		Get_MomentumConvectiveRadial(MESH); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)

		Get_Stresses(MESH, MPI1);
		Get_MomemtumDifusiveAxial(MESH); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Axial U)
		Get_MomemtumDifusiveRadial(MESH); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)

		
		/*if(Rank == 2 && Step == 150){
			for(i = Ix; i < Fx; i++){
				cout<<MESH.DeltasMR[GR(i,NR-1,1)]<<", ";
			}
			cout<<endl;
			
		}*/
		//Get_PressureGradients(MESH); //Cálculo de los gradientes de presión en cada nodo en ambas direciones (U y V)

		Get_Velocities(); //Cálculo del mapa de velocidades futuro

		//Cálculo del campo de temperaturas

		//Comunicación Energías Internas
		MPI1.SendData(ElocalPres, Ix, Fx);
		MPI1.ReceiveData(ElocalPres, Ix, Fx);

		Get_EnergyWalls(MESH); //Cálculo de la energía interna del fluido en las paredes de los volúmenes de control
		Get_EnergyConvective(MESH); //Cálculo del término convectivo de la ecuación de energía

		//Comunicación Temperaturas
		MPI1.SendData(TlocalPres, Ix, Fx);
		MPI1.ReceiveData(TlocalPres, Ix, Fx);

		Get_K(MESH, MPI1);  //Cálculo de las conductividades térmicas en las paredes de los volúmenes de control
		Get_EnergyDiffusive(MESH); //Cálculo del término difusivo de la ecuación de energía
		
		//Get_EnergyPressureTerm(MESH); //Cálculo del término de presión de la ecuación de energía
		//Get_EnergyViscousTerm(MESH); //Cálculo del término viscoso de la ecuación de energía

		Get_Temperatures(); //Cálculo del mapa de temperaturas futuro
	 
		if(Step%1000 == 0){

				//Mandar todas las matrices al ZERO
			//	MPI1.SendMatrixToZero(TauRZ, TauRZglobal, NA, NR, Procesos, Ix, Fx);
				//Matrices Presentes
				MPI1.SendMatrixToZero(RhoLocalPres, RhoGlobalPres, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(UlocalPres, UglobalPres, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(VlocalPres, VglobalPres, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(TlocalPres, TglobalPres, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(Pressure, PresionGlobal, NA, NR, Procesos, Ix, Fx);

				//Matrices Futuras
				MPI1.SendMatrixToZero(RhoLocalFut, RhoGlobalFut, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(UlocalFut, UglobalFut, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(VlocalFut, VglobalFut, NA, NR, Procesos, Ix, Fx);
				MPI1.SendMatrixToZero(TlocalFut, TglobalFut, NA, NR, Procesos, Ix, Fx);

				Get_Stop(MPI1);

				if(Rank == 0){
					cout<<"Step: "<<Step<<", DeltaT: "<<DeltaT<<", Tiempo: "<<Time<<", MaxDif: "<<MaxDifference<<endl;

					//Pasar los resultados a archivos .VTK
				//	sprintf(FileName, "MapaEsfuerzos_Step_%d", Step);
					//EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Esfuerzos", FileName, MESH.MP, TauRZglobal, NA, NR);

					sprintf(FileName, "MapaDensidades_Step_%d", Step);
					EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Densidades", FileName, MESH.MP, RhoGlobalFut, NA, NR);

					sprintf(FileName, "MapaPresiones_Step_%d", Step);
					EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Presión", FileName, MESH.MP, PresionGlobal, NA, NR);

					sprintf(FileName, "MapaTemperaturas_Step_%d", Step);
					EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Temperaturas", FileName, MESH.MP, TglobalFut, NA, NR);

					sprintf(FileName, "MapaVelocidades_Step_%d", Step);
					VectorialVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Velocidades", FileName, MESH.MP, UglobalFut, VglobalFut, NA, NR);

				}

		}

		UpdateParameters();
		MPI_Barrier(MPI_COMM_WORLD);	
	}

}
	