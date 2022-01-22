#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <chrono>
#include <mpi.h>

#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Mesher.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Solver.h"

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

#define DIRECTORIO "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/"

#define sind(x) sin((PI/180)*x)
#define cosd(x) cos((PI/180)*x)

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

	Problema = R1.Problema;  //Problema (Canal/Canal con Cilindro/Perfil)
	ConvergenciaGlobal = 1e-6;
	ConvergenciaGS = 1e-6;
	
	//Geometría del canal
	ChannelLength = R1.GeometryData[0]; //Longitud del canal
	ChannelHeight = R1.GeometryData[1]; //Altura del canal
	ChannelDepth = R1.GeometryData[2]; //Profundidad del canal

	CylinderRadius = 0.25; //Radio del cilindro (m)

	//Datos físicos del problema
	Uref = R1.ProblemPhysicalData[0];
	Vref = R1.ProblemPhysicalData[1];
	RhoRef = R1.ProblemPhysicalData[2]; 
	Tref = R1.ProblemPhysicalData[3];
	
	
	Cp = R1.ProblemPhysicalData[5];
	Gamma = R1.ProblemPhysicalData[6];
	Rideal = R1.ProblemPhysicalData[7];
	Pref = R1.ProblemPhysicalData[4];
	
	NumericalCirculation = 0.0;

	AnguloAtaque = R1.AttackAngle[i];
	PhiCilindro = 0.10*Uref*ChannelHeight;
	
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
			TglobalPres = M1.AllocateDouble(NX, NY, 1);
			PglobalPres = M1.AllocateDouble(NX, NY, 1);

			//Step Futuro
			PhiGlobalFut = M1.AllocateDouble(NX, NY, 1);
			UglobalFut = M1.AllocateDouble(NX, NY, 1);
			VglobalFut = M1.AllocateDouble(NX, NY, 1);
			TglobalFut = M1.AllocateDouble(NX, NY, 1);
			RhoGlobalFut = M1.AllocateDouble(NX, NY, 1);
			PglobalFut = M1.AllocateDouble(NX, NY, 1);

			UwallsMR_Global = M1.AllocateDouble(NX, NY + 1, 1);
			VwallsMU_Global = M1.AllocateDouble(NX + 1, NY, 1);

			PDiff = M1.AllocateDouble(Procesos, 1, 1);
			PCirc = M1.AllocateDouble(Procesos, 1, 1);

			GlobalDensityIndicator = M1.AllocateInt(NX, NY, 1);

			awGlobal = M1.AllocateDouble(NX, NY, 1);
			aeGlobal = M1.AllocateDouble(NX, NY, 1);
			asGlobal = M1.AllocateDouble(NX, NY, 1);
			anGlobal = M1.AllocateDouble(NX, NY, 1);
			apGlobal  = M1.AllocateDouble(NX, NY, 1);
		}
		

		//Matrices locales de propiedades

		DensityIndicator = M1.AllocateInt(Fx - Ix + 2 * Halo, NY, 1);

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

		//Temperatura
		TlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		TlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		//Densidad
		RhoLocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);

		//Presión
		PlocalPres = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		PlocalFut = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		
		//Matrices de valores de las propiedades en las paredes de los volúmenes de control
		UwallsMR = M1.AllocateDouble(Fx - Ix, NY + 1, 1);
		VwallsMU = M1.AllocateDouble(Fx - Ix + 1, NY, 1);

		RhoWallsMU = M1.AllocateDouble(Fx - Ix + 1, NY, 1);
		RhoWallsMR = M1.AllocateDouble(Fx - Ix, NY + 1, 1);

		//Matrices de cálculo de contribuciones a las propiedades

		aw = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		ae = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		as = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		an = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		ap = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY, 1);
		
		//Matrices y arrays de condiciones de contorno

		if(Rank == 0){
			PhiLeft = M1.AllocateDouble(1, NY, 1); //Presión pared izquierda
			RhoLeft = M1.AllocateDouble(1, NY, 1); //Presión pared izquierda
			Uleft = M1.AllocateDouble(1, NY, 1); //Velocidad axial pared izquierda
			Vleft = M1.AllocateDouble(1, NY, 1); //Velocidad radial pared izquierda
			Tleft = M1.AllocateDouble(1, NY, 1); //Temperaturas pared izquierda
			Pleft = M1.AllocateDouble(1, NY, 1); //Presión pared izquierda
		}
		else if(Rank == Procesos - 1){
			PhiRight = M1.AllocateDouble(1, NY, 1); //Presión pared derecha
			RhoRight = M1.AllocateDouble(1, NY, 1); //Presión pared derecha
			Uright = M1.AllocateDouble(1, NY, 1); //Velocidad axial pared derecha
			Vright = M1.AllocateDouble(1, NY, 1); //Velocidad radial pared derecha
			Tright = M1.AllocateDouble(1, NY, 1); //Temperaturas pared derecha
			Pright = M1.AllocateDouble(1, NY, 1); //Presión pared derecha
		}
		
		PhiUp = M1.AllocateDouble(Fx - Ix, 1, 1); 
		PhiDown = M1.AllocateDouble(Fx - Ix, 1, 1); 

		RhoUp = M1.AllocateDouble(Fx - Ix, 1, 1);
		RhoDown = M1.AllocateDouble(Fx - Ix, 1, 1); 

		Uup = M1.AllocateDouble(Fx - Ix, 1, 1); 
		Udown = M1.AllocateDouble(Fx - Ix, 1, 1); 
 
		Vup = M1.AllocateDouble(Fx - Ix, 1, 1);
		Vdown = M1.AllocateDouble(Fx - Ix, 1, 1); 

		Tup = M1.AllocateDouble(Fx - Ix, 1, 1); 
		Tdown = M1.AllocateDouble(Fx - Ix, 1, 1); 

		Pup = M1.AllocateDouble(Fx - Ix, 1, 1);
		Pdown = M1.AllocateDouble(Fx - Ix, 1, 1); 

}

//Inicialización de los campos de temperaturas
void Solver::InitializeFields(Mesher MESH){
int i, j;
double J;
double ny = NY;
double Xcentro = 0.50*ChannelLength;
double Ycentro = 0.50*ChannelHeight;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			J = j;
			PhiLocalPres[LNH(i,j,0)] = 0.5*Uref*ChannelHeight;
			PhiLocalFut[LNH(i,j,0)] =  0.5*Uref*ChannelHeight;

			PhiLocalSup[LNH(i,j,0)] =  0.5*Uref*ChannelHeight;

			UlocalPres[LNH(i,j,0)] = Uref;
			UlocalFut[LNH(i,j,0)] = Uref;

			VlocalPres[LNH(i,j,0)] = 0.0;
			VlocalFut[LNH(i,j,0)] = 0.0;

			TlocalPres[LNH(i,j,0)] = Tref;
			TlocalFut[LNH(i,j,0)] = Tref;

			PlocalFut[LNH(i,j,0)] = Pref;

			
				if(DensityIndicator[LNH(i,j,0)] == 1){
					RhoLocalFut[LNH(i,j,0)] = 0.0;
				}
				else{
					RhoLocalFut[LNH(i,j,0)] = RhoRef;
				}
		
		}
	}
		
}


//Asignación de temperaturas a las condiciones de contorno
void Solver::UpdateBoundaryConditions(Mesher MESH){
int i, j;

	if(Rank == 0){

		for(j = 0; j < NY; j++){	
			RhoLeft[j] = RhoRef;
			Pleft[j] = Pref;
			Vleft[j] = Vref;
			Uleft[j] = Uref;
			Tleft[j] = Tref;
			PhiLeft[j] = Uleft[j]*MESH.MP[G(0,j,1)];
		}
		
	}
	else if(Rank == Procesos - 1){

		for(j = 0; j < NY; j++){
			RhoRight[j] = RhoLocalFut[LNH(NX-1,j,0)];
			Pright[j] = PlocalFut[LNH(NX - 1,j,0)];
			Vright[j] = 0.0;
			Uright[j] = UlocalFut[LNH(NX - 1,j,0)];
			Tright[j] = TlocalFut[LNH(NX - 1,j,0)];
			PhiRight[j] = PhiLocalFut[LNH(NX - 1, j, 0)];
		}
	
	}

	for(i = Ix; i < Fx; i++){
		
		RhoUp[i - Ix] = RhoLocalFut[LNH(i,NY-1,0)];
		Pup[i - Ix] = 0.0;
		Vup[i - Ix] = Vref;
		Uup[i - Ix] = Uref;
		Tup[i - Ix] = 0.0;
		PhiUp[i - Ix] = MESH.MR[GR(i, NY, 1)]*Uref;
	
		RhoDown[i - Ix] = RhoLocalFut[LNH(i,0,0)];
		Pdown[i - Ix] = 0.0;
		Vdown[i - Ix] = Vref;
		Udown[i - Ix] = Uref;
		Tdown[i - Ix] = 0.0;
		PhiDown[i - Ix] = 0.0;
	}
	
}

//Seteo matriz de 0 y 1 para el cálculo de la densidad en los sólidos
void Solver::Get_DensityIndicatorCylinder(Mesher MESH){
int i,j;
double Xcentro = 0.50*ChannelLength;
double Ycentro = 0.50*ChannelHeight;

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){			
				if(pow(MESH.MP[G(i,j,0)] - Xcentro,2.0) + pow(MESH.MP[G(i,j,1)] - Ycentro,2.0) <= pow(CylinderRadius,2.0)){
					DensityIndicator[LNH(i,j,0)] = 1;
				}
				else{
					DensityIndicator[LNH(i,j,0)] = 0;
				}
			}
		}		
}

//Seteo de matriz de 0 para el caso del perfil aerodinámico
void Solver::Get_DensityIndicatorAirfoil(Mesher MESH){
int i, j;
double Xcentro = 0.50*ChannelLength;
double Ycentro = 0.50*ChannelHeight;
double m = 0.04; //Curvatura máxima en (%)
double p = 0.4; //Posición de la curvatura máxima (%)
double t = 0.12; //Espesor relativo máximo
double Cuerda = 1.0; //Cuerda (c) (m)
double Alpha = AnguloAtaque; //Ángulo de ataque
double UpLimit;
double DownLimit;


		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){	
				if(MESH.MP[G(i,j,0)] >= Xcentro && MESH.MP[G(i,j,0)] <= Xcentro + Cuerda*cosd(Alpha)){
						UpLimit = (MESH.MP[G(i,j,0)] - Xcentro)*sind(Alpha) + Ycentro + (5*t*Cuerda*(0.2969*sqrt((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.1260*((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.3516*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,2.0) + 0.2843*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,3.0) - 0.1015*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,4.0)))*cosd(Alpha);
						DownLimit = (MESH.MP[G(i,j,0)] - Xcentro)*sind(Alpha) + Ycentro - (5*t*Cuerda*(0.2969*sqrt((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.1260*((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.3516*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,2.0) + 0.2843*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,3.0) - 0.1015*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,4.0)))*cosd(Alpha);
						if(MESH.MP[G(i,j,1)] <= UpLimit && MESH.MP[G(i,j,1)] >= DownLimit){
							DensityIndicator[LNH(i,j,0)] = 1;
						}
						else{
							DensityIndicator[LNH(i,j,0)] = 0;
						}
				
				}
				else{
					DensityIndicator[LNH(i,j,0)] = 0;
				}

			}
		}		
}

//Cálculo de las densidades en las paredes de lo volúmenes de control
void Solver::Get_DensityWalls(Mesher MESH, ParPro MPI1){
int i, j;

	//Comunicación de densidades entre los procesos
	MPI1.SendData(RhoLocalFut, Ix, Fx);
	MPI1.ReceiveData(RhoLocalFut, Ix, Fx);

	MPI1.SendDataInt(DensityIndicator, Ix, Fx);
	MPI1.ReceiveDataInt(DensityIndicator, Ix, Fx);

	//Nodos U
	if(Rank != 0 && Rank != Procesos - 1){

		for(i = Ix; i < Fx + 1; i++){
			for(j = 0; j < NY; j++){
				if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i-1,j,0)] == 1){
					RhoWallsMU[LU(i,j,0)] = 0.0;
				}
				if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i-1,j,0)] == 0){
					RhoWallsMU[LU(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i-1,j,0)])*(MESH.DeltasMU[GU(i,j,0)]/(MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				}
				if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i-1,j,0)] == 1){
					RhoWallsMU[LU(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i,j,0)])*(MESH.DeltasMU[GU(i,j,0)]/(MESH.MP[G(i,j,0)] - MESH.MU[GU(i,j,0)]));
				}
				if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i-1,j,0)] == 0){
					RhoWallsMU[LU(i,j,0)] =  MESH.DeltasMU[GU(i,j,0)]/(  (MESH.MP[G(i,j,0)] - MESH.MU[GU(i,j,0)])/(RhoRef/RhoLocalFut[LNH(i,j,0)])   +    (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])/(RhoRef/RhoLocalFut[LNH(i-1,j,0)])   );
				}		 
			}
		}

	}
	else if(Rank == 0){

		for(j = 0; j < NY; j++){

			//Parte izquierda
			RhoWallsMU[LU(0,j,0)] = RhoRef/RhoLeft[j];

			for(i = Ix + 1; i < Fx + 1; i++){
				if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i-1,j,0)] == 1){
						RhoWallsMU[LU(i,j,0)] = 0.0;
				}
				if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i-1,j,0)] == 0){
					RhoWallsMU[LU(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i-1,j,0)])*(MESH.DeltasMU[GU(i,j,0)]/(MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				}
				if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i-1,j,0)] == 1){
					RhoWallsMU[LU(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i,j,0)])*(MESH.DeltasMU[GU(i,j,0)]/(MESH.MP[G(i,j,0)] - MESH.MU[GU(i,j,0)]));
				}
				if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i-1,j,0)] == 0){
					RhoWallsMU[LU(i,j,0)] =  MESH.DeltasMU[GU(i,j,0)]/(  (MESH.MP[G(i,j,0)] - MESH.MU[GU(i,j,0)])/(RhoRef/RhoLocalFut[LNH(i,j,0)])   +    (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])/(RhoRef/RhoLocalFut[LNH(i-1,j,0)])   );
				}	
			}
		}
	}
	else if(Rank == Procesos - 1){

		for(j = 0; j < NY; j++){

			//Parte derecha
			RhoWallsMU[LU(NX,j,0)] = RhoRef/RhoLocalFut[LNH(NX-1,j,0)];

			for(i = Ix-1; i < NX; i++){
				if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i-1,j,0)] == 1){
						RhoWallsMU[LU(i,j,0)] = 0.0;
				}
				if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i-1,j,0)] == 0){
					RhoWallsMU[LU(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i-1,j,0)])*(MESH.DeltasMU[GU(i,j,0)]/(MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)]));
				}
				if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i-1,j,0)] == 1){
					RhoWallsMU[LU(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i,j,0)])*(MESH.DeltasMU[GU(i,j,0)]/(MESH.MP[G(i,j,0)] - MESH.MU[GU(i,j,0)]));
				}
				if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i-1,j,0)] == 0){
					RhoWallsMU[LU(i,j,0)] =  MESH.DeltasMU[GU(i,j,0)]/(  (MESH.MP[G(i,j,0)] - MESH.MU[GU(i,j,0)])/(RhoRef/RhoLocalFut[LNH(i,j,0)])   +    (MESH.MU[GU(i,j,0)] - MESH.MP[G(i-1,j,0)])/(RhoRef/RhoLocalFut[LNH(i-1,j,0)])   );
				}				
			}
		}

	}

	//Nodos R

	for(i = Ix; i < Fx; i++){

		//Parte abajo
		RhoWallsMR[LR(i,0,0)] = RhoRef/RhoDown[i - Ix];

		//Parte arriba
		RhoWallsMR[LR(i,NY,0)] = RhoRef/RhoUp[i - Ix];

		for(j = 1; j < NY; j++){  
			if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i,j-1,0)] == 1){
				RhoWallsMR[LR(i,j,0)] = 0.0;
			}
			if(DensityIndicator[LNH(i,j,0)] == 1 && DensityIndicator[LNH(i,j-1,0)] == 0){
				RhoWallsMR[LR(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i,j-1,0)])*(MESH.DeltasMR[GR(i,j,1)]/(MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)]));
			}
			if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i,j-1,0)] == 1){
				RhoWallsMR[LR(i,j,0)] = (RhoRef/RhoLocalFut[LNH(i,j,0)])*(MESH.DeltasMR[GR(i,j,1)]/(MESH.MP[G(i,j,1)] - MESH.MR[GR(i,j,1)]));
			}
			if(DensityIndicator[LNH(i,j,0)] == 0 && DensityIndicator[LNH(i,j-1,0)] == 0){
				RhoWallsMR[LR(i,j,0)] =  MESH.DeltasMR[GR(i,j,1)]/(  (MESH.MP[G(i,j,1)] - MESH.MR[GR(i,j,1)])/(RhoRef/RhoLocalFut[LNH(i,j,0)])   +    (MESH.MR[GR(i,j,1)] - MESH.MP[G(i,j-1,1)])/(RhoRef/RhoLocalFut[LNH(i,j-1,0)])   );
			}
				
		}          
	}
}

//Cálculo de los coeficientes de discretización
void Solver::Get_Coeficients(Mesher MESH){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){

			aw[LNH(i,j,0)] = RhoWallsMU[LU(i,j,0)]*(MESH.DeltasMU[GU(i,j,1)]/MESH.DeltasMU[GU(i,j,0)]);
			ae[LNH(i,j,0)] = RhoWallsMU[LU(i+1,j,0)]*(MESH.DeltasMU[GU(i+1,j,1)]/MESH.DeltasMU[GU(i+1,j,0)]);
			as[LNH(i,j,0)] = RhoWallsMR[LR(i,j,0)]*(MESH.DeltasMR[GR(i,j,0)]/MESH.DeltasMR[GR(i,j,1)]);
			an[LNH(i,j,0)] = RhoWallsMR[LR(i,j+1,0)]*(MESH.DeltasMR[GR(i,j+1,0)]/MESH.DeltasMR[GR(i,j+1,1)]);

		}		
	}

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			ap[LNH(i,j,0)] = aw[LNH(i,j,0)] + ae[LNH(i,j,0)] + as[LNH(i,j,0)] + an[LNH(i,j,0)] + 1e-10;
		}		
	}



}

void Solver::Send_Coefficients(Mesher MESH, ParPro MPI1){
int i, j;
char FileName[300];

	MPI1.SendMatrixToZero(aw, awGlobal, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(ae, aeGlobal, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(as, asGlobal, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(an, anGlobal, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(ap, apGlobal, NX, NY, Procesos, Ix, Fx);

	if(Rank == 0){

		
		sprintf(FileName, "MapaAw_Step_%d", 10);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "aw", FileName, MESH.MP, awGlobal, NX, NY);

			sprintf(FileName, "MapaAe_Step_%d", 10);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "ae", FileName, MESH.MP, aeGlobal, NX, NY);

		sprintf(FileName, "MapaAs_Step_%d", 10);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "as", FileName, MESH.MP, asGlobal, NX, NY);

		sprintf(FileName, "MapaAn_Step_%d", 10);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "an", FileName, MESH.MP, anGlobal, NX, NY);

			sprintf(FileName, "MapaAp_Step_%d", 10);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "ap", FileName, MESH.MP, apGlobal, NX, NY);
	}
	
}
void Solver::Get_MaxDifGS(ParPro MPI1){
int i, j;
MaxDiffGS = 0.0;
MPI_Status ST;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			if(abs((PhiLocalFut[LNH(i,j,0)] - PhiLocalSup[LNH(i,j,0)])) >= MaxDiffGS){
				MaxDiffGS = abs((PhiLocalFut[LNH(i,j,0)] - PhiLocalSup[LNH(i,j,0)]));
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
void Solver::Get_Streamlines(ParPro MPI1, Mesher MESH){
int i, j;
MaxDiffGS = 2.0*ConvergenciaGS;

	while(MaxDiffGS >= ConvergenciaGS){

		if(Rank != 0 && Rank != Procesos - 1){

			for(i = Ix; i < Fx; i++){
				//Parte abajo
				PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)])/ap[LNH(i,0,0)];
				//Parte arriba
				PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix])/ap[LNH(i,NY-1,0)];

				for(j = 1; j < NY - 1; j++){
				//	if(DensityIndicator[LNH(i,j,0)] == 1){
					//	PhiLocalFut[LNH(i,j,0)] = PhiCilindro;
				//	}
				//	else{
						PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)])/ap[LNH(i,j,0)];
						
				//	}
					
				}
			}

		}
		else if(Rank == 0){

				for(i = Ix + 1; i < Fx; i++){
					//Parte abajo
					PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)])/ap[LNH(i,0,0)];
					
					//Parte arriba
					PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix])/ap[LNH(i,NY-1,0)];
					
					for(j = 1; j < NY - 1; j++){
						//if(DensityIndicator[LNH(i,j,0)] == 1){
					//		PhiLocalFut[LNH(i,j,0)] = PhiCilindro;
					//	}
					//	else{
							PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)])/ap[LNH(i,j,0)];
					//	}
					}
				}

				//Parte izquierda
				for(j = 1; j < NY - 1; j++){
					PhiLocalFut[LNH(0,j,0)] = (aw[LNH(0,j,0)]*PhiLeft[j] + ae[LNH(0,j,0)]*PhiLocalFut[LNH(1,j,0)] + as[LNH(0,j,0)]*PhiLocalFut[LNH(0,j-1,0)] + an[LNH(0,j,0)]*PhiLocalFut[LNH(0,j+1,0)])/ap[LNH(0,j,0)];
				}

				//Esquina abajo izquierda
				PhiLocalFut[LNH(0,0,0)] = (aw[LNH(0,0,0)]*PhiLeft[0] + ae[LNH(0,0,0)]*PhiLocalFut[LNH(1,0,0)] + an[LNH(0,0,0)]*PhiLocalFut[LNH(0,1,0)])/ap[LNH(0,0,0)];
				
				//Esquina arriba izquierda
				PhiLocalFut[LNH(0,NY-1,0)] = (aw[LNH(0,NY-1,0)]*PhiLeft[NY-1] + ae[LNH(0,NY-1,0)]*PhiLocalFut[LNH(1,NY-1,0)] + as[LNH(0,NY-1,0)]*PhiLocalFut[LNH(0,NY-2,0)] + an[LNH(0,NY-1,0)]*PhiUp[0])/ap[LNH(0,NY-1,0)];
				
		}
		else if(Rank == Procesos - 1){

			for(i = Ix; i < Fx - 1; i++){
				//Parte abajo
				PhiLocalFut[LNH(i,0,0)] = (aw[LNH(i,0,0)]*PhiLocalFut[LNH(i-1,0,0)] + ae[LNH(i,0,0)]*PhiLocalFut[LNH(i+1,0,0)] + an[LNH(i,0,0)]*PhiLocalFut[LNH(i,1,0)])/ap[LNH(i,0,0)];
			
				//Parte arriba
				PhiLocalFut[LNH(i,NY-1,0)] = (aw[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i-1,NY-1,0)] + ae[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i+1,NY-1,0)] + as[LNH(i,NY-1,0)]*PhiLocalFut[LNH(i,NY-2,0)] + an[LNH(i,NY-1,0)]*PhiUp[i - Ix])/ap[LNH(i,NY-1,0)];
			
				for(j = 1; j < NY - 1; j++){
				//	if(DensityIndicator[LNH(i,j,0)] == 1){
					//	PhiLocalFut[LNH(i,j,0)] = PhiCilindro;
				//	}
					//else{
						PhiLocalFut[LNH(i,j,0)] = (aw[LNH(i,j,0)]*PhiLocalFut[LNH(i-1,j,0)] + ae[LNH(i,j,0)]*PhiLocalFut[LNH(i+1,j,0)] + as[LNH(i,j,0)]*PhiLocalFut[LNH(i,j-1,0)] + an[LNH(i,j,0)]*PhiLocalFut[LNH(i,j+1,0)])/ap[LNH(i,j,0)];
				//	}
				}
			}

			//Parte derecha
			for(j = 1; j < NY - 1; j++){
				PhiLocalFut[LNH(NX-1,j,0)] = (aw[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-2,j,0)] + ae[LNH(NX-1,j,0)]*PhiRight[j] + as[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-1,j-1,0)] + an[LNH(NX-1,j,0)]*PhiLocalFut[LNH(NX-1,j+1,0)])/ap[LNH(NX-1,j,0)];
			}

			//Esquina abajo derecha
			PhiLocalFut[LNH(NX-1,0,0)] = (aw[LNH(NX-1,0,0)]*PhiLocalFut[LNH(NX-2,0,0)] + ae[LNH(NX-1,0,0)]*PhiRight[0] + an[LNH(NX-1,0,0)]*PhiLocalFut[LNH(NX-1,1,0)])/ap[LNH(NX-1,0,0)];
		
			//Esquina arriba derecha
			PhiLocalFut[LNH(NX-1,NY-1,0)] = (aw[LNH(NX-1,NY-1,0)]*PhiLocalFut[LNH(NX-2,NY-1,0)] + ae[LNH(NX-1,NY-1,0)]*PhiRight[NY-1] + as[LNH(NX-1,NY-1,0)]*PhiLocalFut[LNH(NX-1,NY-2,0)] + an[LNH(NX-1,NY-1,0)]*PhiUp[NX-1 - Ix])/ap[LNH(NX-1,NY-1,0)];
		
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
void Solver::Get_Velocities(Mesher MESH, ParPro MPI1){
int i, j;

	MPI1.SendData(PhiLocalFut, Ix, Fx);
	MPI1.ReceiveData(PhiLocalFut, Ix, Fx);

	//Velocidades Horizontales (Nodos R)
	for(i = Ix; i < Fx; i++){

		//Parte arriba
		UwallsMR[LR(i,NY,0)] = Uref;

		//Parte abajo
		UwallsMR[LR(i,0,0)] = Uref;

		for(j = 1; j < NY; j++){
			UwallsMR[LR(i,j,0)] = (RhoWallsMR[LR(i,j,0)]*(PhiLocalFut[LNH(i,j,0)] - PhiLocalFut[LNH(i,j-1,0)]))/MESH.DeltasMR[GR(i,j,1)];
		}
	}


	//Velocidades Verticales (Nodos U)

	if(Rank != 0 && Rank != Procesos - 1){
		for(j = 0; j < NY; j++){
			for(i = Ix; i < Fx + 1; i++){
				VwallsMU[LU(i,j,0)] = - (RhoWallsMU[LU(i,j,0)]*(PhiLocalFut[LNH(i,j,0)] - PhiLocalFut[LNH(i-1,j,0)]))/MESH.DeltasMU[GU(i,j,0)];
			}
		}
	}
	else if(Rank == 0){
		for(j = 0; j < NY; j++){

			//Parte izquierda
			VwallsMU[LU(0,j,0)] = Vleft[j];
			
			for(i = Ix + 1; i < Fx + 1; i++){
				VwallsMU[LU(i,j,0)] = - (RhoWallsMU[LU(i,j,0)]*(PhiLocalFut[LNH(i,j,0)] - PhiLocalFut[LNH(i-1,j,0)]))/MESH.DeltasMU[GU(i,j,0)];
			}
		}
	}
	else if(Rank == Procesos - 1){
		for(j = 0; j < NY; j++){
			
			//Parte derecha
			VwallsMU[LU(NX,j,0)] = Vright[j];
			
			for(i = Ix; i < Fx; i++){
				VwallsMU[LU(i,j,0)] = - (RhoWallsMU[LU(i,j,0)]*(PhiLocalFut[LNH(i,j,0)] - PhiLocalFut[LNH(i-1,j,0)]))/MESH.DeltasMU[GU(i,j,0)];	
			}
		}
	}

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){

				UlocalFut[LNH(i,j,0)] = 0.50*(UwallsMR[LR(i,j,0)] + UwallsMR[LR(i,j+1,0)]);			
				VlocalFut[LNH(i,j,0)] = 0.50*(VwallsMU[LU(i,j,0)] + VwallsMU[LU(i+1,j,0)]);
			
		}
	}
}

//Cálculo del mapa de temperaturas futuro
void Solver::Get_Temperatures(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			TlocalFut[LNH(i,j,0)] = Tref + ((pow(Uref,2.0) + pow(Vref,2.0)) - (pow(UlocalFut[LNH(i,j,0)],2.0) + pow(VlocalFut[LNH(i,j,0)],2.0)))/(2.0*Cp);
		}
	}

}

//Cálculo del mapa de presiones futuro
void Solver::Get_Pressure(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			PlocalFut[LNH(i,j,0)] = Pref*pow(TlocalFut[LNH(i,j,0)]/Tref,Gamma/(Gamma-1));
			//Pref + (1 - (pow(UlocalFut[LNH(i,j,0)],2.0) + pow(VlocalFut[LNH(i,j,0)],2.0))/pow(Uref,2.0))*0.50*RhoRef*pow(Uref,2.0);		
		}
	}

}

//Cálculo del mapa de densidades futuro
void Solver::Get_Density(){
int i, j;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){	
			if(DensityIndicator[LNH(i,j,0)] == 1){
				RhoLocalFut[LNH(i,j,0)] = 0.0;
			}
			else{
				RhoLocalFut[LNH(i,j,0)] = PlocalFut[LNH(i,j,0)]/(Rideal*TlocalFut[LNH(i,j,0)]);
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
				if(abs((TglobalFut[G(i,j,0)] - TglobalPres[G(i,j,0)])) >= MaxDiffGlobal){
					MaxDiffGlobal = abs((TglobalFut[G(i,j,0)] - TglobalPres[G(i,j,0)]));
				}
				if(abs((PglobalFut[G(i,j,0)] - PglobalPres[G(i,j,0)])) >= MaxDiffGlobal){
					MaxDiffGlobal = abs((PglobalFut[G(i,j,0)] - PglobalPres[G(i,j,0)]));
				}
		
				if(abs((PhiGlobalFut[G(i,j,0)] - PhiGlobalPres[G(i,j,0)])) >= MaxDiffGlobal){
					MaxDiffGlobal = abs((PhiGlobalFut[G(i,j,0)] - PhiGlobalPres[G(i,j,0)]));
				}
				if(abs((UglobalFut[G(i,j,0)] - UglobalPres[G(i,j,0)])) >= MaxDiffGlobal){
					MaxDiffGlobal = abs((UglobalFut[G(i,j,0)] - UglobalPres[G(i,j,0)]));
				}
				if(abs((VglobalFut[G(i,j,0)] - VglobalPres[G(i,j,0)])) >= MaxDiffGlobal){
					MaxDiffGlobal = abs((VglobalFut[G(i,j,0)] - VglobalPres[G(i,j,0)]));
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

			TlocalPres[LNH(i,j,0)] = TlocalFut[LNH(i,j,0)];
			PlocalPres[LNH(i,j,0)] = PlocalFut[LNH(i,j,0)];
			PhiLocalPres[LNH(i,j,0)] = PhiLocalFut[LNH(i,j,0)];
			UlocalPres[LNH(i,j,0)] = UlocalFut[LNH(i,j,0)];
			VlocalPres[LNH(i,j,0)] = VlocalFut[LNH(i,j,0)];

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

//Pasar los resultados a un archivo VTK en 2D
void Solver::EscalarVTK2DInt(Mesher MESH, string Carpeta, string Variable, string NombreFile, double *MC, int *GlobalField, int Na, int Nr){
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

//Cálculo de la circulación total en el cilindro
void Solver::Get_Circulation(Mesher MESH, ParPro MPI1){
int i,j;
double Circu = 0.0;

	MPI1.SendData(RhoLocalFut, Ix, Fx);
	MPI1.ReceiveData(RhoLocalFut, Ix, Fx);

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			if(DensityIndicator[LNH(i,j,0)] == 1){
				Circu += MESH.DeltasMU[GU(i+1,j,1)]*VwallsMU[LU(i+1,j,0)] - MESH.DeltasMU[GU(i,j,1)]*VwallsMU[LU(i,j,0)] - MESH.DeltasMR[GR(i,j+1,0)]*UwallsMR[LR(i,j+1,0)] + MESH.DeltasMR[GR(i,j,0)]*UwallsMR[LR(i,j,0)];
			}
		}
	}

	MPI1.SendDataToZero(Circu, PCirc);

	
	if(Rank == 0){
		for(i = 0; i < Procesos; i++){
			NumericalCirculation += PCirc[i];
		}

		//cout<<"Numerical Circulation: "<<CirculationNumerical<<endl;
	}

}

//Cálculo del Cd y del Cl del cilindro
void Solver::Get_AeroCoefficients(Mesher MESH){
int i, j;
double DragForce = 0.0;
double LiftForce = 0.0;
double Xcentro = 0.50*ChannelLength;
double Ycentro = 0.50*ChannelHeight;

double m = 0.04; //Curvatura máxima en (%)
double p = 0.4; //Posición de la curvatura máxima (%)
double t = 0.12; //Espesor relativo máximo
double Cuerda = 1.0; //Cuerda (c) (m)
double Alpha = AnguloAtaque; //Ángulo de ataque
double UpLimit;
double DownLimit;


		for(i = 0; i < NX; i++){
			for(j = 0; j < NY; j++){	
				if(MESH.MP[G(i,j,0)] >= Xcentro && MESH.MP[G(i,j,0)] <= Xcentro + Cuerda*cosd(Alpha)){
						UpLimit = (MESH.MP[G(i,j,0)] - Xcentro)*sind(Alpha) + Ycentro + (5*t*Cuerda*(0.2969*sqrt((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.1260*((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.3516*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,2.0) + 0.2843*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,3.0) - 0.1015*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,4.0)))*cosd(Alpha);
						DownLimit = (MESH.MP[G(i,j,0)] - Xcentro)*sind(Alpha) + Ycentro - (5*t*Cuerda*(0.2969*sqrt((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.1260*((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda) - 0.3516*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,2.0) + 0.2843*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,3.0) - 0.1015*pow((MESH.MP[G(i,j,0)] - Xcentro)/Cuerda,4.0)))*cosd(Alpha);
						if(MESH.MP[G(i,j,1)] <= UpLimit && MESH.MP[G(i,j,1)] >= DownLimit){
							GlobalDensityIndicator[G(i,j,0)] = 1;
						}
						else{
							GlobalDensityIndicator[G(i,j,0)] = 0;
						}
				
				}
				else{
					GlobalDensityIndicator[G(i,j,0)] = 0;
				}

			}
		}		
/*
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
				if(pow(MESH.MP[G(i,j,0)] - Xcentro,2.0) + pow(MESH.MP[G(i,j,1)] - Ycentro,2.0) <= pow(CylinderRadius,2.0)){
					GlobalDensityIndicator[G(i,j,0)] = 1;
					PglobalFut[G(i,j,0)] = 0.0;
				}
				else{
					GlobalDensityIndicator[G(i,j,0)] = 0;
				}		
		}
	}
*/
	FILE *fp5;
	fp5 = fopen("/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/NumericalResults/CpDown.txt","a");

	FILE *fp9;
	fp9 = fopen("/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/NumericalResults/CpUp.txt","a");
	for(i = 1; i < NX-1; i++){
		for(j = 1; j < NY-1; j++){
			if(GlobalDensityIndicator[G(i,j,0)] == 1){
				if(GlobalDensityIndicator[G(i+1,j,0)] == 0){
					
					DragForce += - MESH.DeltasMP[G(i,j,1)]*PglobalFut[G(i+1,j,0)];
				}
				if(GlobalDensityIndicator[G(i-1,j,0)] == 0){
					DragForce += MESH.DeltasMP[G(i,j,1)]*PglobalFut[G(i-1,j,0)];
				}

				if(GlobalDensityIndicator[G(i,j-1,0)] == 0){	
					fprintf(fp5,"%f \t %f \n",MESH.MP[G(i,j,0)] - 0.50*ChannelLength, (PglobalFut[G(i,j-1,0)] - Pref)/(0.50*RhoRef*pow(Uref,2)));		
					LiftForce += MESH.DeltasMP[G(i,j,0)]*PglobalFut[G(i,j-1,0)];
				}
				if(GlobalDensityIndicator[G(i,j+1,0)] == 0){
					fprintf(fp9,"%f \t %f \n",MESH.MP[G(i,j,0)] - 0.50*ChannelLength, (PglobalFut[G(i,j+1,0)] - Pref)/(0.50*RhoRef*pow(Uref,2)));	
					LiftForce += - MESH.DeltasMP[G(i,j,0)]*PglobalFut[G(i,j+1,0)];
				}
			}
		}
	}

	fclose(fp5);
	fclose(fp9);

	Cd = DragForce/(0.50*RhoRef*pow(Uref,2)*1.0);
	Cl = LiftForce/(0.50*RhoRef*pow(Uref,2)*1.0);

	FILE *fp1;
		fp1 = fopen("/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/NumericalResults/AeroCoefficients.txt","a");
			fprintf(fp1,"%f \t %f \t %f \t %f \t %f \t %f \n", DragForce, LiftForce, Cd, Cl, NumericalCirculation, -RhoRef*Uref*NumericalCirculation);
					
	fclose(fp1);

	FILE *fp2;
		fp2 = fopen("/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/NumericalResults/AttackAngleStudy.txt","a");
			fprintf(fp2,"%f \t %f \t %f \t %f \n",-AnguloAtaque, Cl, LiftForce, -RhoRef*Uref*NumericalCirculation);
					
	fclose(fp2);

}

//Ejecución de todos los procesos del solver
void Solver::ExecuteSolver(Memory M1, ReadData R1, ParPro MPI1, Mesher MESH){
int i, j;
int Step = 0;
char FileName[300];
MaxDiffGlobal = 2.0*ConvergenciaGlobal; 

	AllocateMatrix(M1);
	Get_DensityIndicatorAirfoil(MESH);
	//Get_DensityIndicatorCylinder(MESH);
	InitializeFields(MESH);
	
	//Pasar todas las matrices al ZERO
	//Matrices
	MPI1.SendMatrixToZero(TlocalFut, TglobalFut, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(RhoLocalFut, RhoGlobalFut, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(PlocalFut, PglobalFut, NX, NY, Procesos, Ix, Fx);	
	MPI1.SendMatrixToZero(UlocalFut, UglobalFut, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(PhiLocalFut, PhiGlobalFut, NX, NY, Procesos, Ix, Fx);
	MPI1.SendMatrixToZero(VlocalFut, VglobalFut, NX, NY, Procesos, Ix, Fx);

	if(Rank == 0){

		
		sprintf(FileName, "MapaStreamFunction_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "StreamFunctions", FileName, MESH.MP, PhiGlobalFut, NX, NY);

		sprintf(FileName, "MapaPresiones_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Presiones", FileName, MESH.MP, PglobalFut, NX, NY);

		sprintf(FileName, "MapaTemperaturas_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Temperaturas", FileName, MESH.MP, TglobalFut, NX, NY);

		sprintf(FileName, "MapaDensidades_Step_%d", Step);
		EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Densidades", FileName, MESH.MP, RhoGlobalFut, NX, NY);

		sprintf(FileName, "MapaVelocidades_Step_%d", Step);
		VectorialVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Velocidades", FileName, MESH.MP, UglobalFut, VglobalFut, NX, NY);
	}
	
	auto start = std::chrono::high_resolution_clock::now();
	while(MaxDiffGlobal >= ConvergenciaGlobal){

		Step++;
		UpdateBoundaryConditions(MESH);
		
		Get_DensityWalls(MESH, MPI1);
		Get_Coeficients(MESH);
		
		Get_Streamlines(MPI1, MESH);

		Get_Velocities(MESH, MPI1);
		Get_Temperatures(); 
		Get_Pressure();
		
		Get_Density();
		
		if(Step%500 == 0){
			
			//Pasar todas las matrices al ZERO

			//Matrices Step Presente
			MPI1.SendMatrixToZero(TlocalPres, TglobalPres, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(PlocalPres, PglobalPres, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(PhiLocalPres, PhiGlobalPres, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(UlocalPres, UglobalPres, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(VlocalPres, VglobalPres, NX, NY, Procesos, Ix, Fx);

			//Matrices Step Futuro
			MPI1.SendMatrixToZero(TlocalFut, TglobalFut, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(RhoLocalFut, RhoGlobalFut, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(PlocalFut, PglobalFut, NX, NY, Procesos, Ix, Fx);
		
			MPI1.SendMatrixToZero(UlocalFut, UglobalFut, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(PhiLocalFut, PhiGlobalFut, NX, NY, Procesos, Ix, Fx);
			MPI1.SendMatrixToZero(VlocalFut, VglobalFut, NX, NY, Procesos, Ix, Fx);

		
			Get_Stop(MPI1);

			if(Rank == 0){
	
				cout<<"Step: "<<Step<<", MaxDif: "<<MaxDiffGlobal<<endl;

			
			}
		}
		UpdatePropertiesFields();

	}
	
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";

	if(Rank == 0){

				sprintf(FileName, "MapaStreamFunction_Step_%d", 2);
				EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "StreamFunctions", FileName, MESH.MP, PhiGlobalFut, NX, NY);

				sprintf(FileName, "MapaPresiones_Step_%d", 2);
				EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Presiones", FileName, MESH.MP, PglobalFut, NX, NY);

				sprintf(FileName, "MapaTemperaturas_Step_%d", 2);
				EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Temperaturas", FileName, MESH.MP, TglobalFut, NX, NY);

				sprintf(FileName, "MapaDensidades_Step_%d", 2);
				EscalarVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Densidades", FileName, MESH.MP, RhoGlobalFut, NX, NY);

				sprintf(FileName, "MapaVelocidades_Step_%d", 2);
				VectorialVTK2D(MESH, "ParaviewResults/PropertiesResults/", "Velocidades", FileName, MESH.MP, UglobalFut, VglobalFut, NX, NY);
				
			}

	Get_Circulation(MESH, MPI1);
	if(Rank == 0){
		Get_AeroCoefficients(MESH);
	}


}

