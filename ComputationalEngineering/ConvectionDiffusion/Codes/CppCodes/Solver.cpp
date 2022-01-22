#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <chrono>
#include <mpi.h>

#include "/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/Codes/HeaderCodes/Mesher.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/Codes/HeaderCodes/Solver.h"

using namespace std;

#define PI 3.141592653589793

#define G(i,j,dim) (((j) + (i)*NY) + NX*NY*(dim)) //Global Index
#define GU(i,j,dim) (((j) + (i)*NY) + (NX+1)*NY*(dim)) //Global Index Matriz U
#define GR(i,j,dim) (((j) + (i)*(NY+1)) + NX*(NY+1)*(dim)) //Global Index Matriz R

#define MU(i,j,dim) ((j) + NY*(dim)) 
#define MR(i,j,dim) ((i) + (Fx - Ix)*(dim) + (Fx - Ix)*(4)*(j))

#define VR(i,j,dim) ((i) + (Fx - Ix)*(j))

#define LU(i, j, dim) (((j) + ((i)-Ix) * NY)) //Local Index Axial Velocity Nodes
#define LR(i, j, dim) (((j) + ((i)-Ix) * (NY + 1))) //Local Index Radial Velocity Nodes

#define LNH(i,j,dim) (((j) + ((i) - Ix + Halo)*NY)) //Local Index No Halo 
#define LSH(i,j,dim) (((j) + ((i) - Ix)*NY) + NY*(Fx-Ix + 2*Halo)*dim) //Local Index Si Halo

#define DIRECTORIO "/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/"

//Constructor del mallador
Solver::Solver(Memory M1, ReadData R1, ParPro MPI1, Mesher MESH, int i){
		
	//Datos del problema
	NX = MESH.NX;
	NY = MESH.NY;

	//Datos para la computacion en paralelo
	Procesos = MPI1.Procesos;
	Rank = MPI1.Rank;
	Ix = MESH.Ix;
	Fx = MESH.Fx;
	Halo = MPI1.Get_Halo();

	ConvergenciaGlobal = 1e-8;
	ConvergenciaGS = 1e-7;
	
	//Geometría del canal
	ChannelLength = R1.GeometryData[0]; //Longitud del canal
	ChannelHeight = R1.GeometryData[1]; //Altura del canal
	ChannelDepth = R1.GeometryData[2]; //Profundidad del canal

	Problema = R1.NumericalData[4];

	//Datos Físicos del problema
	Peclet = R1.ProblemPhysicalData[4];
	RhoRef = R1.ProblemPhysicalData[2];

	Uref = Peclet/(RhoRef*ChannelLength);
	Vref = Peclet/(RhoRef*ChannelLength);
	
	Alpha =  R1.ProblemPhysicalData[3];
	AlphaRad = Alpha*PI/180;

	PhiIz = R1.ProblemPhysicalData[5];
	PhiDer = R1.ProblemPhysicalData[6];
		
	EsquemaLargo = R1.Esquema;

	FR = R1.FactorRelax[i];
}

//Alojamiento de memoria para las matrices necesarias
void Solver::AllocateMatrix(Memory M1){
	
	//Matrices de los mapas de propiedades del problema
	
		//Matrices globales de propiedades

		if(Rank == 0){

			//Step Presente
			PhiGlobalPres = M1.AllocateDouble(NX, NY, 1);
			UglobalPres = M1.AllocateDouble(NX, NY, 1);
			VglobalPres = M1.AllocateDouble(NX, NY, 1);

			//Step Futuro
			PhiGlobalFut = M1.AllocateDouble(NX, NY, 1);
			UglobalFut = M1.AllocateDouble(NX, NY, 1);
			VglobalFut = M1.AllocateDouble(NX, NY, 1);
			
			PDiff = M1.AllocateDouble(Procesos, 1, 1);
			PDT = M1.AllocateDouble(Procesos, 1, 1);

			CoordenadasAnalitico = M1.AllocateDouble(11, 1, 1);
			PhiAnalitico = M1.AllocateDouble(11, 1, 1);

			PhiReal = M1.AllocateDouble(11, 1, 1);

			MinAbajo = M1.AllocateDouble(11, 1, 1);
			MinArriba = M1.AllocateDouble(11, 1, 1);

			PhiAbajo = M1.AllocateDouble(11, 1, 1);
			PhiArriba = M1.AllocateDouble(11, 1, 1);

		}
		

		//Matrices locales de propiedades

		//Streamline
		PhiLocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		PhiLocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		PhiLocalSup = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		//Velocidad Horizontal
		UlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		UlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		//Velocidad Vertical
		VlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		VlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		
		//Matrices de valores de las propiedades en las paredes de los volúmenes de control
		UwallsMU = M1.AllocateDouble(Fx - Ix + 1, NY, 1);
		VwallsMR = M1.AllocateDouble(Fx - Ix, NY + 1, 1);

		PhiWallsMU = M1.AllocateDouble(Fx - Ix + 1, NY, 1);
		PhiWallsMR = M1.AllocateDouble(Fx - Ix, NY + 1, 1);

		//Matrices de cálculo de contribuciones a las propiedades

		Dw = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		De = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		Ds = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		Dn = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		Cw = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		Ce = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		Cs = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		Cn = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		aw = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		ae = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		as = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		an = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		ap = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		
		bp = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		//Matrices y arrays de condiciones de contorno

		if(Rank == 0){
			
			PhiLeft = M1.AllocateDouble(1, NY, 1); //Presión pared izquierda
			Uleft = M1.AllocateDouble(1, NY, 1); //Velocidad axial pared izquierda
			Vleft = M1.AllocateDouble(1, NY, 1); //Velocidad radial pared izquierda

		}
		else if(Rank == Procesos - 1){

			PhiRight = M1.AllocateDouble(1, NY, 1); //Presión pared derecha
			Uright = M1.AllocateDouble(1, NY, 1); //Velocidad axial pared derecha
			Vright = M1.AllocateDouble(1, NY, 1); //Velocidad radial pared derecha
			
		}
		
		PhiUp = M1.AllocateDouble(Fx - Ix, 1, 1); 
		PhiDown = M1.AllocateDouble(Fx - Ix, 1, 1); 

		Uup = M1.AllocateDouble(Fx - Ix, 1, 1); 
		Udown = M1.AllocateDouble(Fx - Ix, 1, 1); 
 
		Vup = M1.AllocateDouble(Fx - Ix, 1, 1);
		Vdown = M1.AllocateDouble(Fx - Ix, 1, 1); 

}

//Inicialización de los campos de temperaturas
void Solver::InitializeFields(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){

			PhiLocalPres[LNH(i,j,0)] = 0.0;
			PhiLocalFut[LNH(i,j,0)] =  0.0;

			PhiLocalSup[LNH(i,j,0)] =  0.0;

		}
	}
		
}


//Asignación de temperaturas a las condiciones de contorno
void Solver::UpdateBoundaryConditions(Mesher MESH){
int i, j;

	if(Rank == 0){
		if(Problema == 1){
			for(j = 0; j < NY; j++){	
				PhiLeft[j] = PhiIz;
			}
		}
		else if(Problema == 2){
			for(j = 0; j < NY; j++){	
				PhiLeft[j] = PhiIz;
			}
		}
		else if(Problema == 3){
			for(j = 0; j < NY; j++){	
				PhiLeft[j] = PhiIz;
			}
		}
		else if(Problema == 4){
			for(j = 0; j < NY; j++){	
				PhiLeft[j] = 1.0 - tanh(10.0);
			}
		}	
	}
	else if(Rank == Procesos - 1){
		if(Problema == 1){
			for(j = 0; j < NY; j++){	
				PhiRight[j] = PhiDer;
			}
		}
		else if(Problema == 2){
			for(j = 0; j < NY; j++){	
				PhiRight[j] = PhiDer;
			}
		}
		else if(Problema == 3){
			for(j = 0; j < NY; j++){	
				PhiRight[j] = PhiDer;
			}
		}
		else if(Problema == 4){
			for(j = 0; j < NY; j++){	
				PhiRight[j] = 1.0 - tanh(10.0);
			}
		}	
	}

	if(Problema == 1){
		for(i = Ix; i < Fx; i++){
			PhiUp[i - Ix] = PhiLocalFut[LNH(i,NY-1,0)];
			PhiDown[i - Ix] = PhiLocalFut[LNH(i,0,0)];
		}	
	}
	else if(Problema == 2){
		for(i = Ix; i < Fx; i++){
			PhiUp[i - Ix] = PhiLocalFut[LNH(i,NY-1,0)];
			PhiDown[i - Ix] = PhiLocalFut[LNH(i,0,0)];
		}		
	}
	else if(Problema == 3){
		for(i = Ix; i < Fx; i++){
			PhiUp[i - Ix] = PhiIz;
			PhiDown[i - Ix] = PhiDer;
		}		
	}
	else if(Problema == 4){
		for(i = Ix; i < Fx; i++){
			PhiUp[i - Ix] = 1.0 - tanh(10.0);
			if((MESH.MR[GR(i,0,0)] - 0.50*ChannelLength) <= 0.0){
				PhiDown[i - Ix] = 1.0 + tanh((2.0*(MESH.MR[GR(i,0,0)] - 0.50*ChannelLength) + 1.0)*10.0);
			}
			else if((MESH.MR[GR(i,0,0)] - 0.50*ChannelLength) > 0.0){
				PhiDown[i - Ix] = PhiLocalFut[LNH(i,0,0)];
			}
			
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
		for(j = 0; j < NY; j++){

			//CFL Velocidades (U + V)
			if(abs((Tpar*MESH.DeltasMP[G(i,j,0)])/(abs(UlocalFut[LNH(i,j,0)])) < DeltasT)){
				DeltasT = (Tpar * MESH.DeltasMP[G(i, j, 0)]) / (abs(UlocalFut[LNH(i, j, 0)]));
			}
			if(abs((Tpar*MESH.DeltasMP[G(i,j,1)])/(abs(VlocalFut[LNH(i,j,0)])) < DeltasT)){
				DeltasT = (Tpar * MESH.DeltasMP[G(i, j, 1)]) / (abs(VlocalFut[LNH(i, j, 0)]));
			}

			//CFL Difusivo
			if((Tpar*RhoRef*pow(MESH.DeltasMP[G(i,j,0)],2.0))/(1.0 + 1e-10) < DeltasT){
				DeltasT = (Tpar * RhoRef * pow(MESH.DeltasMP[G(i, j, 0)], 2.0)) / (1.0 + 1e-10);
			}
			if((Tpar*RhoRef*pow(MESH.DeltasMP[G(i,j,1)],2.0))/(1.0 + 1e-10) < DeltasT){
				DeltasT = (Tpar * RhoRef * pow(MESH.DeltasMP[G(i, j, 1)], 2.0)) / (1.0 + 1e-10);
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

//Función con los diferentes esquemas convectivos utilizados
double Solver::ConvectiveScheme(double CoordObjetivo, double Velocity, double Coord1, double Phi1, double Coord2, double Phi2, double Coord3, double Phi3, double Coord4, double Phi4, string Esquema){

double PhiObjetivo;

double CoordD;
double PhiD;

double CoordC;
double PhiC;

double CoordU;
double PhiU;

	if (Velocity < 0.0 || (Phi1 == 0.0 && Coord1 == 0.0)){

		CoordD = Coord2;
		PhiD = Phi2;
		CoordC = Coord3;
		PhiC = Phi3;
		CoordU = Coord4;
		PhiU = Phi4;

	}
	else if(Velocity >= 0.0 || (Phi4 == 0.0 && Coord4 == 0.0)){

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
	//double Xade;

	PhiAdimC = (PhiC - PhiU)/(PhiD - PhiU);

	AdimCoordC = (CoordC - CoordU)/(CoordD - CoordU);

	AdimCoordE = (CoordObjetivo - CoordU)/(CoordD - CoordU);

	//Evaluacion
	double PhiF;

	if (PhiD == PhiU){
		PhiObjetivo = PhiD;
	}
	else{
		if(Esquema == "CDS"){
			PhiF = ((AdimCoordE - AdimCoordC)/(1.0 - AdimCoordC)) + ((AdimCoordE - 1.0)/(AdimCoordC - 1.0))*PhiAdimC;	
		}
		else if(Esquema == "UDS"){
			PhiF = PhiAdimC;	
		}
		else if(Esquema == "SUDS"){
			PhiF = (AdimCoordE/AdimCoordC)*PhiAdimC;
		}
		else if(Esquema == "QUICK"){
			PhiF = AdimCoordE + (((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0))))*(PhiAdimC - AdimCoordC);
		}
		else if(Esquema == "SMART"){
			if(PhiAdimC > 0 && PhiAdimC < AdimCoordC/3.0){
				PhiF = -((AdimCoordE*(1.0 - 3.0*AdimCoordC + 2.0*AdimCoordE))/(AdimCoordC*(AdimCoordC - 1.0)))*PhiAdimC;
			}
			else if(PhiAdimC > AdimCoordC/6.0 && PhiAdimC < (AdimCoordC/AdimCoordE)*(1.0 + AdimCoordE - AdimCoordC)){
				PhiF = ((AdimCoordE*(AdimCoordE - AdimCoordC))/(1.0 - AdimCoordC)) + ((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0)))*PhiAdimC;
			}
			else if(PhiAdimC > (AdimCoordC/AdimCoordE)*(1.0 + AdimCoordE - AdimCoordC) && PhiAdimC < 1.0){
				PhiF = 1.0;
			}
			else{
				PhiF = PhiAdimC;
			}
		}

		//Dimensionalizacion
		PhiObjetivo = PhiU + (PhiD - PhiU)*PhiF;
	}

	return PhiObjetivo;

}



//Cálculo de las densidades en las paredes de lo volúmenes de control
void Solver::Get_PhiWalls(Mesher MESH, ParPro MPI1){
int i, j;

	//Comunicación de densidades entre los procesos
	MPI1.SendData(PhiLocalFut, Ix, Fx);
	MPI1.ReceiveData(PhiLocalFut, Ix, Fx);

	//Nodos R
	for(i = Ix; i < Fx; i++){
		//Parte abajo
		PhiWallsMR[LR(i,0,0)] = PhiDown[i - Ix];
		PhiWallsMR[LR(i,1,0)] = ConvectiveScheme(MESH.MR[GR(i,1,1)], VwallsMR[LR(i,1,0)], MESH.MR[GR(i,0,1)], PhiDown[i - Ix], MESH.MP[G(i,0,1)], PhiLocalFut[LNH(i,0,0)], MESH.MP[G(i,1,1)], PhiLocalFut[LNH(i,1,0)], MESH.MP[G(i,2,1)], PhiLocalFut[LNH(i,2,0)], EsquemaLargo);

		//Parte arriba
		PhiWallsMR[LR(i,NY,0)] = PhiUp[i - Ix];
		PhiWallsMR[LR(i,NY - 1,0)] = ConvectiveScheme(MESH.MR[GR(i,NY - 1,1)], VwallsMR[LR(i,NY - 1,0)], MESH.MP[G(i,NY-3,1)], PhiLocalFut[LNH(i,NY-3,0)], MESH.MP[G(i,NY-2,1)], PhiLocalFut[LNH(i,NY-2,0)], MESH.MP[G(i,NY-1,1)], PhiLocalFut[LNH(i,NY-1,0)], MESH.MR[GR(i,NY,1)], PhiUp[i - Ix], EsquemaLargo);

		for(j = 2; j < NY - 1; j++){
			PhiWallsMR[LR(i,j,0)] = ConvectiveScheme(MESH.MR[GR(i,j,1)], VwallsMR[LR(i,j,0)], MESH.MP[G(i,j-2,1)], PhiLocalFut[LNH(i,j-2,0)], MESH.MP[G(i,j-1,1)], PhiLocalFut[LNH(i,j-1,0)], MESH.MP[G(i,j,1)], PhiLocalFut[LNH(i,j,0)], MESH.MP[G(i,j+1,1)], PhiLocalFut[LNH(i,j+1,0)], EsquemaLargo);
		}
	}

	//Nodos U
	if(Rank != 0 && Rank != Procesos - 1){

			for(i = Ix; i < Fx + 1; i++){
				for(j = 0; j < NY; j++){
					PhiWallsMU[LU(i,j,0)] = ConvectiveScheme(MESH.MU[GU(i,j,0)], UwallsMU[LU(i,j,0)], MESH.MP[G(i-2,j,0)], PhiLocalFut[LNH(i-2,j,0)], MESH.MP[G(i-1,j,0)], PhiLocalFut[LNH(i-1,j,0)], MESH.MP[G(i,j,0)], PhiLocalFut[LNH(i,j,0)], MESH.MP[G(i+1,j,0)], PhiLocalFut[LNH(i+1,j,0)], EsquemaLargo);
				}
			}

	}
	else if(Rank == 0){

			for(j = 0; j < NY; j++){
				//Parte izquierda
				PhiWallsMU[LU(0,j,0)] = PhiLeft[j];
				PhiWallsMU[LU(1,j,0)] = ConvectiveScheme(MESH.MU[GU(1,j,0)], UwallsMU[LU(1,j,0)], MESH.MU[GU(0,j,0)], PhiLeft[j], MESH.MP[G(0,j,0)], PhiLocalFut[LNH(0,j,0)], MESH.MP[G(1,j,0)], PhiLocalFut[LNH(1,j,0)], MESH.MP[G(2,j,0)], PhiLocalFut[LNH(2,j,0)], EsquemaLargo);

				for(i = Ix + 2; i < Fx + 1; i++){
					PhiWallsMU[LU(i,j,0)] = ConvectiveScheme(MESH.MU[GU(i,j,0)], UwallsMU[LU(i,j,0)], MESH.MP[G(i-2,j,0)], PhiLocalFut[LNH(i-2,j,0)], MESH.MP[G(i-1,j,0)], PhiLocalFut[LNH(i-1,j,0)], MESH.MP[G(i,j,0)], PhiLocalFut[LNH(i,j,0)], MESH.MP[G(i+1,j,0)], PhiLocalFut[LNH(i+1,j,0)], EsquemaLargo);
				}
			}

	}
	else if(Rank == Procesos - 1){
			
			for(j = 0; j < NY; j++){
				//Parte derecha
				PhiWallsMU[LU(NX,j,0)] = PhiRight[j];
				PhiWallsMU[LU(NX - 1,j,0)] = ConvectiveScheme(MESH.MU[GU(NX-1,j,0)], UwallsMU[LU(NX-1,j,0)], MESH.MP[G(NX-3,j,0)], PhiLocalFut[LNH(NX-3,j,0)], MESH.MP[G(NX-2,j,0)], PhiLocalFut[LNH(NX-2,j,0)], MESH.MP[G(NX-1,j,0)], PhiLocalFut[LNH(NX-1,j,0)], MESH.MU[GU(NX,j,0)], PhiRight[j], EsquemaLargo);

				for(i = Ix; i < Fx - 1; i++){
					PhiWallsMU[LU(i,j,0)] = ConvectiveScheme(MESH.MU[GU(i,j,0)], UwallsMU[LU(i,j,0)], MESH.MP[G(i-2,j,0)], PhiLocalFut[LNH(i-2,j,0)], MESH.MP[G(i-1,j,0)], PhiLocalFut[LNH(i-1,j,0)], MESH.MP[G(i,j,0)], PhiLocalFut[LNH(i,j,0)], MESH.MP[G(i+1,j,0)], PhiLocalFut[LNH(i+1,j,0)], EsquemaLargo);
				}
			}
	}

}

//Cálculo de las densidades en las paredes de lo volúmenes de control
void Solver::Get_PhiWalls2(Mesher MESH, ParPro MPI1){
int i, j;

	//Comunicación de densidades entre los procesos
	MPI1.SendData(PhiLocalFut, Ix, Fx);
	MPI1.ReceiveData(PhiLocalFut, Ix, Fx);

	//Nodos R
	for(i = Ix; i < Fx; i++){
		//Parte abajo
		PhiWallsMR[LR(i,0,0)] = PhiDown[i - Ix];
		
		//Parte arriba
		PhiWallsMR[LR(i,NY,0)] = PhiUp[i - Ix];
		
		for(j = 1; j < NY; j++){
			PhiWallsMR[LR(i,j,0)] = 0.50*(PhiLocalFut[LNH(i,j,0)] + PhiLocalFut[LNH(i,j - 1,0)]);
				}
	}

	//Nodos U
	if(Rank != 0 && Rank != Procesos - 1){

			for(i = Ix; i < Fx + 1; i++){
				for(j = 0; j < NY; j++){
					PhiWallsMU[LU(i,j,0)] = 0.50*(PhiLocalFut[LNH(i,j,0)] + PhiLocalFut[LNH(i -  1,j,0)]);
					}
			}

	}
	else if(Rank == 0){

			for(j = 0; j < NY; j++){
				//Parte izquierda
				PhiWallsMU[LU(0,j,0)] = PhiLeft[j];

				for(i = Ix + 1; i < Fx + 1; i++){
					PhiWallsMU[LU(i,j,0)] = 0.50*(PhiLocalFut[LNH(i,j,0)] + PhiLocalFut[LNH(i - 1,j,0)]);
				}
			}

	}
	else if(Rank == Procesos - 1){
			
			for(j = 0; j < NY; j++){
				//Parte derecha
				PhiWallsMU[LU(NX,j,0)] = PhiRight[j];
				
				for(i = Ix; i < Fx; i++){
					PhiWallsMU[LU(i,j,0)] = 0.50*(PhiLocalFut[LNH(i,j,0)] + PhiLocalFut[LNH(i - 1,j,0)]);
				}
			}
	}

}

//Cálculo del término difusivo de la ecuación de Convección-Difusión
void Solver::Get_Diffusive(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			Dw[LNH(i,j,0)] = MESH.SupMP[G(i,j,0)]/MESH.DeltasMU[GU(i,j,0)];
			De[LNH(i,j,0)] = MESH.SupMP[G(i,j,1)]/MESH.DeltasMU[GU(i + 1,j,0)];
			Ds[LNH(i,j,0)] = MESH.SupMP[G(i,j,2)]/MESH.DeltasMR[GR(i,j,1)];
			Dn[LNH(i,j,0)] = MESH.SupMP[G(i,j,3)]/MESH.DeltasMR[GR(i,j+1,1)];
		}
	}
}

//Cálculo del término convectivo de la ecuación de Convección-Difusión
void Solver::Get_Convective(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			Cw[LNH(i,j,0)] = RhoRef*UwallsMU[LU(i,j,0)]*MESH.SupMP[G(i,j,0)];
			Ce[LNH(i,j,0)] = - RhoRef*UwallsMU[LU(i+1,j,0)]*MESH.SupMP[G(i,j,1)];
			Cs[LNH(i,j,0)] = RhoRef*VwallsMR[LR(i,j,0)]*MESH.SupMP[G(i,j,2)];
			Cn[LNH(i,j,0)] = - RhoRef*VwallsMR[LR(i,j+1,0)]*MESH.SupMP[G(i,j,3)];
		}
	}

}

//Cálculo de los coeficientes de discretización A
void Solver::Get_Coeficients(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			aw[LNH(i,j,0)] = Dw[LNH(i,j,0)];
			ae[LNH(i,j,0)] = De[LNH(i,j,0)];
			as[LNH(i,j,0)] = Ds[LNH(i,j,0)];
			an[LNH(i,j,0)] = Dn[LNH(i,j,0)];

			ap[LNH(i,j,0)] = aw[LNH(i,j,0)] + ae[LNH(i,j,0)] + as[LNH(i,j,0)] + an[LNH(i,j,0)] + (RhoRef*MESH.VolMP[G(i,j,0)])/DeltaT;

		}		
	}

}

//Cálculo de los coeficientes de discretización Betas
void Solver::Get_BetaCoefficient(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			bp[LNH(i,j,0)] = (RhoRef*MESH.VolMP[G(i,j,0)]*PhiLocalFut[LNH(i,j,0)])/DeltaT + Cw[LNH(i,j,0)]*PhiWallsMU[LU(i,j,0)] + Ce[LNH(i,j,0)]*PhiWallsMU[LU(i+1,j,0)] + Cs[LNH(i,j,0)]*PhiWallsMR[LR(i,j,0)] + Cn[LNH(i,j,0)]*PhiWallsMR[LR(i,j+1,0)];
		}
	}

}

void Solver::Get_MaxDifGS(ParPro MPI1){
int i, j;
MaxDiffGS = 0.0;
MPI_Status ST;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			if(abs(PhiLocalFut[LNH(i,j,0)] - PhiLocalSup[LNH(i,j,0)]) >= MaxDiffGS){
				MaxDiffGS = abs(PhiLocalFut[LNH(i,j,0)] - PhiLocalSup[LNH(i,j,0)]);
			//	PhiLocalSup[LNH(i,j,0)] = PhiLocalFut[LNH(i,j,0)];
			}
		}
	}

	MPI1.SendDataToZero(MaxDiffGS, PDiff);

	double Diff;
	if(Rank == 0){
		Diff = PDiff[0];
		for(i = 1; i < Procesos; i++){
			if(PDiff[i] >= Diff){
				Diff = PDiff[i];
			}
		}
	}

	MPI1.SendDataToAll(Diff, MaxDiffGS);

}

//Resolución de las ecuaciones con Gauss-Seidel
void Solver::Get_StreamlinesFR(ParPro MPI1){
int i, j;
MaxDiffGS = 2.0*ConvergenciaGS;

	while(MaxDiffGS >= ConvergenciaGS){

		if(Rank != 0 && Rank != Procesos - 1){

			for(i = Ix; i < Fx; i++){
				//Parte abajo
				PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + as[LNH(i,0,0)]*PhiDown[i - Ix] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)] + bp[LNH(i,0,0)])/ap[LNH(i,0,0)];
				
				//Parte arriba
				PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix] + bp[LNH(i,NY-1,0)])/ap[LNH(i,NY-1,0)];
			
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)] + bp[LNH(i,j,0)])/ap[LNH(i,j,0)];
				}
			}

		}
		else if(Rank == 0){

				for(i = Ix + 1; i < Fx; i++){
					//Parte abajo
					PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + as[LNH(i,0,0)]*PhiDown[i - Ix] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)] + bp[LNH(i,0,0)])/ap[LNH(i,0,0)];
			
					//Parte arriba
					PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix] + bp[LNH(i,NY-1,0)])/ap[LNH(i,NY-1,0)];
			
					for(j = 1; j < NY - 1; j++){
						PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)] + bp[LNH(i,j,0)])/ap[LNH(i,j,0)];
					}
				}

				//Parte izquierda
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(0,j,0)] = (aw[LNH(0,j,0)]*PhiLeft[j] + ae[LNH(0,j,0)]*PhiLocalFut[LNH(1,j,0)] + as[LNH(0,j,0)]*PhiLocalFut[LNH(0,j-1,0)] + an[LNH(0,j,0)]*PhiLocalFut[LNH(0,j+1,0)] + bp[LNH(0,j,0)])/ap[LNH(0,j,0)];
				}

				//Esquina abajo izquierda
				PhiLocalFut[LNH(0,0,0)] = (aw[LNH(0,0,0)]*PhiLeft[0] + ae[LNH(0,0,0)]*PhiLocalFut[LNH(1,0,0)] +as[LNH(0,0,0)]*PhiDown[0 - Ix] + an[LNH(0,0,0)]*PhiLocalFut[LNH(0,1,0)] + bp[LNH(0,0,0)])/ap[LNH(0,0,0)];
				
				//Esquina arriba izquierda
				PhiLocalFut[LNH(0,NY-1,0)] = (aw[LNH(0,NY-1,0)]*PhiLeft[NY-1] + ae[LNH(0,NY-1,0)]*PhiLocalFut[LNH(1,NY-1,0)] + as[LNH(0,NY-1,0)]*PhiLocalFut[LNH(0,NY-2,0)] + an[LNH(0,NY-1,0)]*PhiUp[0 - Ix] + bp[LNH(0,NY-1,0)])/ap[LNH(0,NY-1,0)];
				
		}
		else if(Rank == Procesos - 1){

			for(i = Ix; i < Fx - 1; i++){
				//Parte abajo
				PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + as[LNH(i,0,0)]*PhiDown[i - Ix] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)] + bp[LNH(i,0,0)])/ap[LNH(i,0,0)];
			
				//Parte arriba
				PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix] + bp[LNH(i,NY-1,0)])/ap[LNH(i,NY-1,0)];
			
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)] + bp[LNH(i,j,0)])/ap[LNH(i,j,0)];
				}
			}

			//Parte derecha
			for(j = 1; j < NY - 1; j++){
				PhiLocalFut[LNH(NX-1,j,0)] = (aw[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-2,j,0)] + ae[LNH(NX-1,j,0)]*PhiRight[j] + as[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-1,j-1,0)] + an[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-1,j+1,0)] + bp[LNH(NX-1,j,0)])/ap[LNH(NX-1,j,0)];
			}

			//Esquina abajo derecha
			PhiLocalFut[LNH(NX-1,0,0)] = (aw[LNH(NX-1,0,0)]*PhiLocalFut[LNH(NX-2,0,0)] + ae[LNH(NX-1,0,0)]*PhiRight[0] + an[LNH(NX-1,0,0)]*PhiLocalFut[LNH(NX-1,1,0)] + as[LNH(NX-1,0,0)]*PhiDown[NX-1 - Ix] + bp[LNH(NX-1,j,0)])/ap[LNH(NX-1,0,0)];
		
			//Esquina arriba derecha
			PhiLocalFut[LNH(NX-1,NY-1,0)] = (aw[LNH(NX-1,NY-1,0)]*PhiLocalFut[LNH(NX-2,NY-1,0)] + ae[LNH(NX-1,NY-1,0)]*PhiRight[NY-1] + as[LNH(NX-1,NY-1,0)]*PhiLocalFut[LNH(NX-1,NY-2,0)] + an[LNH(NX-1,NY-1,0)]*PhiUp[NX-1 - Ix] + bp[LNH(NX-1,NY-1,0)])/ap[LNH(NX-1,NY-1,0)];
		
		}

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				PhiLocalFut[LNH(i,j,0)] = PhiLocalSup[LNH(i,j,0)] + FR*(PhiLocalFut[LNH(i,j,0)] - PhiLocalSup[LNH(i,j,0)]);
			}
		}
		Get_MaxDifGS(MPI1);

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				PhiLocalSup[LNH(i,j,0)] = PhiLocalFut[LNH(i,j,0)];
			}
		}

		MPI1.SendData(PhiLocalFut, Ix, Fx);
		MPI1.ReceiveData(PhiLocalFut, Ix, Fx);

	}

}

//Resolución de las ecuaciones con Gauss-Seidel
void Solver::Get_Streamlines(ParPro MPI1){
int i, j;
MaxDiffGS = 2.0*ConvergenciaGS;
FR = 1.0;
	while(MaxDiffGS >= ConvergenciaGS){

		if(Rank != 0 && Rank != Procesos - 1){

			for(i = Ix; i < Fx; i++){
				//Parte abajo
				PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + as[LNH(i,0,0)]*PhiDown[i - Ix] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)] + bp[LNH(i,0,0)])/ap[LNH(i,0,0)];
				
				//Parte arriba
				PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix] + bp[LNH(i,NY-1,0)])/ap[LNH(i,NY-1,0)];
			
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)] + bp[LNH(i,j,0)])/ap[LNH(i,j,0)];
				}
			}

		}
		else if(Rank == 0){

				for(i = Ix + 1; i < Fx; i++){
					//Parte abajo
					PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + as[LNH(i,0,0)]*PhiDown[i - Ix] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)] + bp[LNH(i,0,0)])/ap[LNH(i,0,0)];
			
					//Parte arriba
					PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix] + bp[LNH(i,NY-1,0)])/ap[LNH(i,NY-1,0)];
			
					for(j = 1; j < NY - 1; j++){
						PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)] + bp[LNH(i,j,0)])/ap[LNH(i,j,0)];
					}
				}

				//Parte izquierda
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(0,j,0)] = (aw[LNH(0,j,0)]*PhiLeft[j] + ae[LNH(0,j,0)]*PhiLocalFut[LNH(1,j,0)] + as[LNH(0,j,0)]*PhiLocalFut[LNH(0,j-1,0)] + an[LNH(0,j,0)]*PhiLocalFut[LNH(0,j+1,0)] + bp[LNH(0,j,0)])/ap[LNH(0,j,0)];
				}

				//Esquina abajo izquierda
				PhiLocalFut[LNH(0,0,0)] = (aw[LNH(0,0,0)]*PhiLeft[0] + ae[LNH(0,0,0)]*PhiLocalFut[LNH(1,0,0)] +as[LNH(0,0,0)]*PhiDown[0 - Ix] + an[LNH(0,0,0)]*PhiLocalFut[LNH(0,1,0)] + bp[LNH(0,0,0)])/ap[LNH(0,0,0)];
				
				//Esquina arriba izquierda
				PhiLocalFut[LNH(0,NY-1,0)] = (aw[LNH(0,NY-1,0)]*PhiLeft[NY-1] + ae[LNH(0,NY-1,0)]*PhiLocalFut[LNH(1,NY-1,0)] + as[LNH(0,NY-1,0)]*PhiLocalFut[LNH(0,NY-2,0)] + an[LNH(0,NY-1,0)]*PhiUp[0 - Ix] + bp[LNH(0,NY-1,0)])/ap[LNH(0,NY-1,0)];
				
		}
		else if(Rank == Procesos - 1){

			for(i = Ix; i < Fx - 1; i++){
				//Parte abajo
				PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + as[LNH(i,0,0)]*PhiDown[i - Ix] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)] + bp[LNH(i,0,0)])/ap[LNH(i,0,0)];
			
				//Parte arriba
				PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix] + bp[LNH(i,NY-1,0)])/ap[LNH(i,NY-1,0)];
			
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)] + bp[LNH(i,j,0)])/ap[LNH(i,j,0)];
				}
			}

			//Parte derecha
			for(j = 1; j < NY - 1; j++){
				PhiLocalFut[LNH(NX-1,j,0)] = (aw[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-2,j,0)] + ae[LNH(NX-1,j,0)]*PhiRight[j] + as[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-1,j-1,0)] + an[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-1,j+1,0)] + bp[LNH(NX-1,j,0)])/ap[LNH(NX-1,j,0)];
			}

			//Esquina abajo derecha
			PhiLocalFut[LNH(NX-1,0,0)] = (aw[LNH(NX-1,0,0)]*PhiLocalFut[LNH(NX-2,0,0)] + ae[LNH(NX-1,0,0)]*PhiRight[0] + an[LNH(NX-1,0,0)]*PhiLocalFut[LNH(NX-1,1,0)] + as[LNH(NX-1,0,0)]*PhiDown[NX-1 - Ix] + bp[LNH(NX-1,j,0)])/ap[LNH(NX-1,0,0)];
		
			//Esquina arriba derecha
			PhiLocalFut[LNH(NX-1,NY-1,0)] = (aw[LNH(NX-1,NY-1,0)]*PhiLocalFut[LNH(NX-2,NY-1,0)] + ae[LNH(NX-1,NY-1,0)]*PhiRight[NY-1] + as[LNH(NX-1,NY-1,0)]*PhiLocalFut[LNH(NX-1,NY-2,0)] + an[LNH(NX-1,NY-1,0)]*PhiUp[NX-1 - Ix] + bp[LNH(NX-1,NY-1,0)])/ap[LNH(NX-1,NY-1,0)];
		
		}

		Get_MaxDifGS(MPI1);

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				PhiLocalSup[LNH(i,j,0)] = PhiLocalFut[LNH(i,j,0)];
			}
		}

		MPI1.SendData(PhiLocalFut, Ix, Fx);
		MPI1.ReceiveData(PhiLocalFut, Ix, Fx);

	}

}

//Cálculo de los campos de velocidades
void Solver::Get_Velocities(Mesher MESH){
int i, j;

	if(Problema == 1){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				UlocalFut[LNH(i,j,0)] = Uref;
				VlocalFut[LNH(i,j,0)] = 0.0;
			}
		}

	}
	else if(Problema == 2){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				UlocalFut[LNH(i,j,0)] = 0.0;
				VlocalFut[LNH(i,j,0)] = -Vref;
			}
		}

	}
	else if(Problema == 3){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				UlocalFut[LNH(i,j,0)] = Uref*cos(AlphaRad);
				VlocalFut[LNH(i,j,0)] = Vref*sin(AlphaRad);
			}
		}

	}
	else if(Problema == 4){
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				UlocalFut[LNH(i,j,0)] = 2.0*MESH.MP[G(i,j,1)]*(1.0  - pow(MESH.MP[G(i,j,0)] - 0.50*ChannelLength,2.0));
				VlocalFut[LNH(i,j,0)] = -2.0*(MESH.MP[G(i,j,0)] - 0.50*ChannelLength)*(1.0  - pow(MESH.MP[G(i,j,1)],2.0));
			}
		}
	}
	
	//Nodos R
	if(Problema == 1){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY + 1; j++){
				VwallsMR[LR(i,j,0)] = 0.0;
			}
		}

	}
	else if(Problema == 2){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY + 1; j++){
				VwallsMR[LR(i,j,0)] = -Vref;
			}
		}

	}
	else if(Problema == 3){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY + 1; j++){
				VwallsMR[LR(i,j,0)] = Vref*sin(AlphaRad);
			}
		}

	}
	else if(Problema == 4){
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY + 1; j++){
				VwallsMR[LR(i,j,0)] = -2.0*(MESH.MR[GR(i,j,0)] - 0.50*ChannelLength)*(1.0 - pow(MESH.MR[GR(i,j,1)],2.0));
			}
		}
	}

	//Nodos U
	if(Problema == 1){

		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				UwallsMU[LU(i,j,0)] = Uref;
			}
		}

	}
	else if(Problema == 2){

		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				UwallsMU[LU(i,j,0)] = 0.0;
			}
		}

	}
	else if(Problema == 3){

		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				UwallsMU[LU(i,j,0)] = Uref*cos(AlphaRad);
			}
		}

	}
	else if(Problema == 4){
		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				UwallsMU[LU(i,j,0)] = 2.0*MESH.MU[GU(i,j,1)]*(1.0 - pow(MESH.MU[GU(i,j,0)] - 0.50*ChannelLength,2.0));
			}
		}
	}

}

//Cálculo de la diferencia de resultados entre Steps
void Solver::Get_Stop(ParPro MPI1){
int i, j;
MaxDiffGlobal = 0.0;

	if(Rank == 0){
		for(i = 0; i < NX; i++){
			for(j = 0; j < NY; j++){
				if(abs((PhiGlobalFut[G(i,j,0)] - PhiGlobalPres[G(i,j,0)])/(PhiGlobalPres[G(i,j,0)] + 1e-10)) >= MaxDiffGlobal){
					MaxDiffGlobal = abs((PhiGlobalFut[G(i,j,0)] - PhiGlobalPres[G(i,j,0)])/(PhiGlobalPres[G(i,j,0)] + 1e-10));
				}
			}
		}
	}

	MPI1.SendDataToAll(MaxDiffGlobal, MaxDiffGlobal);

}

void Solver::UpdatePropertiesFields(){
int i, j;	

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			PhiLocalPres[LNH(i,j,0)] = PhiLocalFut[LNH(i,j,0)];
		}
	}
}

void Solver::Get_NumericalResults(Mesher MESH, string Carpeta, string Problema){
int i,j;

ofstream file;
string FileName;
stringstream InitialNameMP;
string FinalNameMP;
int RHO = RhoRef;
string txt = ".txt";

	InitialNameMP<<DIRECTORIO<<Carpeta<<Problema<<"_"<<EsquemaLargo<<"_"<<RHO<<txt;

	FinalNameMP = InitialNameMP.str();
        file.open(FinalNameMP.c_str());

		file<<0.0<<"\t"<<2.000<<"\t"<<endl;

		for(i = 0; i < NX; i++){
			if(MESH.MR[GR(i,0,0)] - 0.50*ChannelLength >= 0.0){
				file<<(MESH.MR[GR(i,0,0)] - 0.50*ChannelLength)<<"\t"<<PhiGlobalFut[G(i,0,0)]<<"\t"<<endl;
			}			
    	}
		file<<1.0<<"\t"<<0.0<<"\t"<<endl;

	file.close();

}

void Solver::Get_AnalyticalResults(){
int i;

	for(i = 0; i < 11; i++){
		CoordenadasAnalitico[i] = 0.1*i;
	}

	if(RhoRef == 10){
		PhiAnalitico[0] = 1.989;
		PhiAnalitico[1] = 1.402;
		PhiAnalitico[2] = 1.146;
		PhiAnalitico[3] = 0.946;
		PhiAnalitico[4] = 0.775;
		PhiAnalitico[5] = 0.621;
		PhiAnalitico[6] = 0.480;
		PhiAnalitico[7] = 0.349;
		PhiAnalitico[8] = 0.227;
		PhiAnalitico[9] = 0.111;
		PhiAnalitico[10] = 0.00;
	}
	else if(RhoRef == 1000){
		PhiAnalitico[0] = 2.0000;
		PhiAnalitico[1] = 1.9997;
		PhiAnalitico[2] = 1.9990;
		PhiAnalitico[3] = 1.9850;
		PhiAnalitico[4] = 1.8410;
		PhiAnalitico[5] = 0.9510;
		PhiAnalitico[6] = 0.1540;
		PhiAnalitico[7] = 0.0000;
		PhiAnalitico[8] = 0.0000;
		PhiAnalitico[10] = 0.0000;
	}
	else if(RhoRef == 1000000){
		PhiAnalitico[0] = 2.000;
		PhiAnalitico[1] = 2.000;
		PhiAnalitico[2] = 2.000;
		PhiAnalitico[3] = 1.999;
		PhiAnalitico[4] = 1.964;
		PhiAnalitico[5] = 1.000;
		PhiAnalitico[6] = 0.036;
		PhiAnalitico[7] = 0.000;
		PhiAnalitico[8] = 0.000;
		PhiAnalitico[10] = 0.000;
	}

}

void Solver::RelativeError(Mesher MESH, string Carpeta, string Problema){
int i,j;
double ErrorRelativoTotal = 0.0;
double MediaErrorRelativo;

	
	//Búsqueda de los puntos más cercanos a los analíticos
	for(j = 0; j < 9; j++){
		MinAbajo[j] = ChannelLength;
		MinArriba[j] = ChannelLength;
		for(i = 0.40*NX; i < NX; i++){
			if((MESH.MP[G(i,0,0)] - 0.50*ChannelLength) - CoordenadasAnalitico[j] < 0.0){
				if(abs((MESH.MP[G(i,0,0)] - 0.50*ChannelLength) - CoordenadasAnalitico[j]) < abs((MinAbajo[j] - 0.50*ChannelLength) - CoordenadasAnalitico[j])){
					MinAbajo[j] = MESH.MP[G(i,0,0)];
					PhiAbajo[j] = PhiGlobalFut[G(i,0,0)];
				}
			}
			if((MESH.MP[G(i,0,0)] - 0.50*ChannelLength) - CoordenadasAnalitico[j] > 0.0){
				MinArriba[j] = MESH.MP[G(i,0,0)];
				PhiArriba[j] = PhiGlobalFut[G(i,0,0)];
				break;
			
			}
			
		}
	
	}


	//Cálculo de la propiedad en esos puntos
	for(j = 1; j < 9; j++){
		PhiReal[j] = PhiAbajo[j] + ((PhiArriba[j] - PhiAbajo[j])/(MinArriba[j] - MinAbajo[j]))*((CoordenadasAnalitico[j] + 0.50*ChannelLength) - MinAbajo[j]);
	}

	cout<<"MinAbajo \t MinArriba"<<endl;
	for(i = 1; i < 9; i++){
		cout<<MinAbajo[i]<<"\t"<<MinArriba[i]<<endl;
	}
	cout<<"PhiAbajo \t PhiArriba"<<endl;
	for(i = 1; i < 9; i++){
		cout<<PhiAbajo[i]<<"\t"<<PhiArriba[i]<<endl;
	}
	cout<<"PhiAnalitico \t PhiReal"<<endl;
	for(i = 1; i < 9; i++){
		cout<<PhiAnalitico[i]<<"\t"<<PhiReal[i]<<endl;
	}
	for(i = 1; i < 9; i++){
		if(PhiAnalitico[i] == 0.0){
			ErrorRelativoTotal += abs(PhiAnalitico[i] - PhiReal[i])/(PhiAnalitico[i] + 1.0);
		}
		else{
			ErrorRelativoTotal += abs(PhiAnalitico[i] - PhiReal[i])/(PhiAnalitico[i] + 1e-10);
		}
		
	}

	MediaErrorRelativo = ErrorRelativoTotal/8.0;


int RHO = RhoRef;

	char Directorio[200];
	sprintf(Directorio,"/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/NumericalResults/SmithHutton_MeshStudy_%d",RHO);

	FILE *fp1;
		fp1 = fopen(Directorio,"a");
			fprintf(fp1,"%d \t %d \t %d \t %f \n", NX, NY, NX*NY, MediaErrorRelativo);
					
	fclose(fp1);

}

//Calcular el outlet como la simetría
void Solver::Get_SymetricOutlet(Mesher MESH, string Carpeta, string Problema){
int i;

ofstream file;
string FileName;
stringstream InitialNameMP;
string FinalNameMP;
int RHO = RhoRef;
string txt = ".txt";

	InitialNameMP<<DIRECTORIO<<Carpeta<<Problema<<"_"<<EsquemaLargo<<"_"<<RHO<<"Symmetry"<<txt;

	FinalNameMP = InitialNameMP.str();
        file.open(FinalNameMP.c_str());

		file<<0.0<<"\t"<<2.000<<"\t"<<endl;
		for(i = NX/2 + 1; i < NX; i++){
				file<<(MESH.MR[GR(i,0,0)] - 0.50*ChannelLength)<<"\t"<<(1.0 + tanh((2.0*(MESH.MR[GR(NX-i,0,0)] - 0.50*ChannelLength) + 1.0)*10.0))<<"\t"<<endl;			
    	}
		file<<1.0<<"\t"<<0.0<<"\t"<<endl;

	file.close();

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
void Solver::ExecuteSolver(Memory M1, ReadData R1, ParPro MPI1, Mesher MESH, int i){
int Step = 0;
char FileName[300];
double Time = 0.0;

MaxDiffGlobal = 2.0*ConvergenciaGlobal; 

 auto start = std::chrono::high_resolution_clock::now();

// Portion of code to be timed



	AllocateMatrix(M1);
	InitializeFields(MESH);
	Get_Velocities(MESH);
	if(Rank == 0){ Get_AnalyticalResults(); }
	
	//Pasar todas las matrices al ZERO
/*	//Matrices
	MPI1.SendMatrixToZero(UlocalFut, UglobalFut, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(PhiLocalFut, PhiGlobalFut, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(VlocalFut, VglobalFut, NX, NY, Procesos, Ix, Fx);

	if(Rank == 0){
		sprintf(FileName, "MapaStreamFunction_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "StreamFunctions", FileName, MESH.MP, PhiGlobalFut, NX, NY);

		sprintf(FileName, "MapaVelocidades_Step_%d", Step);
		VectorialVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Velocidades", FileName, MESH.MP, UglobalFut, VglobalFut, NX, NY);
	}
*/
	
Get_SymetricOutlet(MESH, "GnuPlotResults/NumericalResults/", "SmithHutton");
	while(MaxDiffGlobal >= ConvergenciaGlobal){

		Step++;

		UpdateBoundaryConditions(MESH);
		
		Get_StepTime(MESH, MPI1);
		Time += DeltaT;

		Get_PhiWalls(MESH, MPI1);

	/*	if(Step == 2 && Rank == 0){
			for(j = NY-1; j>= 0; j--){
				for(i = Ix; i < Fx; i++){
					cout<<1<<", ";
				}
				cout<<endl;
			}
		}*/

		Get_Diffusive(MESH);
		Get_Convective(MESH);

		Get_Coeficients(MESH); //Cálculo de los coeficientes de discretización
		Get_BetaCoefficient(MESH);

		Get_Streamlines(MPI1);
		//Get_StreamlinesFR(MPI1);

		
		if(Step%1000 == 0){
			
			//Pasar todas las matrices al ZERO

			//Matrices Step Presente
			MPI1.SendMatrixToZero(PhiLocalPres, PhiGlobalPres, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(UlocalPres, UglobalPres, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(VlocalPres, VglobalPres, NX, NY, Procesos, Ix, Fx);

			//Matrices Step Futuro
			MPI1.SendMatrixToZero(UlocalFut, UglobalFut, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(PhiLocalFut, PhiGlobalFut, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(VlocalFut, VglobalFut, NX, NY, Procesos, Ix, Fx);

			Get_Stop(MPI1);

			if(Rank == 0){
				cout<<"Simulation: "<<i<<", Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<endl;

				if(Problema == 4){
					Get_NumericalResults(MESH, "GnuPlotResults/NumericalResults/", "SmithHutton");
					
				}	
				
				sprintf(FileName, "MapaStreamFunction_Step_%d", Step);
				EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "StreamFunctions", FileName, MESH.MP, PhiGlobalFut, NX, NY);

				sprintf(FileName, "MapaVelocidades_Step_%d", Step);
				VectorialVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Velocidades", FileName, MESH.MP, UglobalFut, VglobalFut, NX, NY);
			}
		}

		UpdatePropertiesFields();

	}

	if(Rank == 0){
		RelativeError(MESH, "NumericalResults/", "SmithHutton");
	}
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	if(Rank == 0){
		int RHO = RhoRef;
		char Directorio[200];
		sprintf(Directorio,"/home/sergiogus/Desktop/ComputationalEngineering/ConvectionDiffusion/NumericalResults/SmithHutton_FRStudy_%d",RHO);

		FILE *fp1;
			fp1 = fopen(Directorio,"a");
				fprintf(fp1,"%f \t %f \n", FR, elapsed.count());
					
		fclose(fp1);
	}


	
}

