#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
using namespace std;

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/"

SpalartAllmaras::SpalartAllmaras(){
	
}

//Cálculo del término 1 de la ecuación de Spalart-Allmaras
void Solver::SA_Term1(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino1[i][j] = -(1.0/MESH.VolMP[i][j])*(
								  MESH.SupMP[i][j][0]*0.50*(vpresent[i][j] + vpresent[i-1][j])*0.50*(RhoPres[i][j] + RhoPres[i-1][j])*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))
								+ MESH.SupMP[i][j][1]*0.50*(vpresent[i][j] + vpresent[i+1][j])*0.50*(RhoPres[i][j] + RhoPres[i+1][j])*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))
								+ MESH.SupMP[i][j][2]*0.50*(vpresent[i][j] + vpresent[i][j-1])*0.50*(RhoPres[i][j] + RhoPres[i][j-1])*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
								+ MESH.SupMP[i][j][3]*0.50*(vpresent[i][j] + vpresent[i][j+1])*0.50*(RhoPres[i][j] + RhoPres[i][j+1])*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino1[0][j] = -(1.0/MESH.VolMP[0][j])*(
								MESH.SupMP[0][j][0]*vleft[j]*RhoLeft[j]*(UwallsMU[0][j]*cos(PI + MESH.AngleMU[0][j]) + VwallsMU[0][j]*sin(PI + MESH.AngleMU[0][j]))
								+ MESH.SupMP[0][j][1]*0.50*(vpresent[0][j] + vpresent[1][j])*0.50*(RhoPres[0][j] + RhoPres[1][j])*(UwallsMU[1][j]*cos(MESH.AngleMU[1][j]) + VwallsMU[1][j]*sin(MESH.AngleMU[1][j]))
								+ MESH.SupMP[0][j][2]*0.50*(vpresent[0][j] + vpresent[0][j-1])*0.50*(RhoPres[0][j] + RhoPres[0][j-1])*(UwallsMR[0][j]*cos(1.50*PI + MESH.AngleMR[0][j]) + VwallsMR[0][j]*sin(1.50*PI + MESH.AngleMR[0][j]))
								+ MESH.SupMP[0][j][3]*0.50*(vpresent[0][j] + vpresent[0][j+1])*0.50*(RhoPres[0][j] + RhoPres[0][j+1])*(UwallsMR[0][j+1]*cos(0.50*PI + MESH.AngleMR[0][j+1]) + VwallsMR[0][j+1]*sin(0.50*PI + MESH.AngleMR[0][j+1]))
							  );

		//Parte derecha
		SA_Termino1[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(
								MESH.SupMP[NA-1][j][0]*0.50*(vpresent[NA-1][j] + vpresent[NA-2][j])*0.50*(RhoPres[NA-1][j] + RhoPres[NA-2][j])*(UwallsMU[NA-1][j]*cos(PI + MESH.AngleMU[NA-1][j]) + VwallsMU[NA-1][j]*sin(PI + MESH.AngleMU[NA-1][j]))
								+ MESH.SupMP[NA-1][j][1]*vright[j]*RhoRight[j]*(UwallsMU[NA][j]*cos(MESH.AngleMU[NA][j]) + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][j]))
								+ MESH.SupMP[NA-1][j][2]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1])*0.50*(RhoPres[NA-1][j] + RhoPres[NA-1][j-1])*(UwallsMR[NA-1][j]*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + VwallsMR[NA-1][j]*sin(1.50*PI + MESH.AngleMR[NA-1][j]))
								+ MESH.SupMP[NA-1][j][3]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1])*0.50*(RhoPres[NA-1][j] + RhoPres[NA-1][j+1])*(UwallsMR[NA-1][j+1]*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + VwallsMR[NA-1][j+1]*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))
							  );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino1[i][0] = -(1.0/MESH.VolMP[i][0])*(
								MESH.SupMP[i][0][0]*0.50*(vpresent[i][0] + vpresent[i-1][0])*0.50*(RhoPres[i][0] + RhoPres[i-1][0])*(UwallsMU[i][0]*cos(PI + MESH.AngleMU[i][0]) + VwallsMU[i][0]*sin(PI + MESH.AngleMU[i][0]))
								+ MESH.SupMP[i][0][1]*0.50*(vpresent[i][0] + vpresent[i+1][0])*0.50*(RhoPres[i][0] + RhoPres[i+1][0])*(UwallsMU[i+1][0]*cos(MESH.AngleMU[i+1][0]) + VwallsMU[i+1][0]*sin(MESH.AngleMU[i+1][0]))
								+ MESH.SupMP[i][0][2]*vdown[i]*RhoDown[i]*(UwallsMR[i][0]*cos(1.50*PI + MESH.AngleMR[i][0]) + VwallsMR[i][0]*sin(1.50*PI + MESH.AngleMR[i][0]))
								+ MESH.SupMP[i][0][3]*0.50*(vpresent[i][0] + vpresent[i][1])*0.50*(RhoPres[i][0] + RhoPres[i][1])*(UwallsMR[i][1]*cos(0.50*PI + MESH.AngleMR[i][1]) + VwallsMR[i][1]*sin(0.50*PI + MESH.AngleMR[i][1]))
							  );

		//Parte arriba
		SA_Termino1[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(
								MESH.SupMP[i][NRad-1][0]*0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1])*0.50*(RhoPres[i][NRad-1] + RhoPres[i-1][NRad-1])*(UwallsMU[i][NRad-1]*cos(PI + MESH.AngleMU[i][NRad-1]) + VwallsMU[i][NRad-1]*sin(PI + MESH.AngleMU[i][NRad-1]))
								+ MESH.SupMP[i][NRad-1][1]*0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1])*0.50*(RhoPres[i][NRad-1] + RhoPres[i+1][NRad-1])*(UwallsMU[i+1][NRad-1]*cos(MESH.AngleMU[i+1][NRad-1]) + VwallsMU[i+1][NRad-1]*sin(MESH.AngleMU[i+1][NRad-1]))
								+ MESH.SupMP[i][NRad-1][2]*0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2])*0.50*(RhoPres[i][NRad-1] + RhoPres[i][NRad-2])*(UwallsMR[i][NRad-1]*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + VwallsMR[i][NRad-1]*sin(1.50*PI + MESH.AngleMR[i][NRad-1]))
								+ MESH.SupMP[i][NRad-1][3]*vup[i]*RhoUp[i]*(UwallsMR[i][NRad]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + VwallsMR[i][NRad]*sin(0.50*PI + MESH.AngleMR[i][NRad]))
							  );
	}

	//Esquina abajo izquierda
	SA_Termino1[0][0] = -(1.0/MESH.VolMP[0][0])*(
								MESH.SupMP[0][0][0]*vleft[0]*RhoLeft[0]*(UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0]) + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0]))
								+ MESH.SupMP[0][0][1]*0.50*(vpresent[0][0] + vpresent[1][0])*0.50*(RhoPres[0][0] + RhoPres[1][0])*(UwallsMU[1][0]*cos(MESH.AngleMU[1][0]) + VwallsMU[1][0]*sin(MESH.AngleMU[1][0]))
								+ MESH.SupMP[0][0][2]*vdown[0]*RhoDown[0]*(UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))
								+ MESH.SupMP[0][0][3]*0.50*(vpresent[0][0] + vpresent[0][1])*0.50*(RhoPres[0][0] + RhoPres[0][1])*(UwallsMR[0][1]*cos(0.50*PI + MESH.AngleMR[0][1]) + VwallsMR[0][1]*sin(0.50*PI + MESH.AngleMR[0][1]))
							  );

	//Esquina arriba izquierda
	SA_Termino1[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
								MESH.SupMP[0][NRad-1][0]*vleft[NRad-1]*RhoLeft[NRad-1]*(UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))
								+ MESH.SupMP[0][NRad-1][1]*0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1])*0.50*(RhoPres[0][NRad-1] + RhoPres[1][NRad-1])*(UwallsMU[1][NRad-1]*cos(MESH.AngleMU[1][NRad-1]) + VwallsMU[1][NRad-1]*sin(MESH.AngleMU[1][NRad-1]))
								+ MESH.SupMP[0][NRad-1][2]*0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2])*0.50*(RhoPres[0][NRad-1] + RhoPres[0][NRad-2])*(UwallsMR[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + VwallsMR[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))
								+ MESH.SupMP[0][NRad-1][3]*vup[0]*RhoUp[0]*(UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad]))
							  );

	//Esquina abajo derecha
	SA_Termino1[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								MESH.SupMP[NA-1][0][0]*0.50*(vpresent[NA-1][0] + vpresent[NA-2][0])*0.50*(RhoPres[NA-1][0] + RhoPres[NA-2][0])*(UwallsMU[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0]) + VwallsMU[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0]))
								+ MESH.SupMP[NA-1][0][1]*vright[0]*RhoRight[0]*(UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0]) + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0]))
								+ MESH.SupMP[NA-1][0][2]*vdown[NA-1]*RhoDown[NA-1]*(UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))
								+ MESH.SupMP[NA-1][0][3]*0.50*(vpresent[NA-1][0] + vpresent[NA-1][1])*0.50*(RhoPres[NA-1][0] + RhoPres[NA-1][1])*(UwallsMR[NA-1][1]*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + VwallsMR[NA-1][1]*sin(0.50*PI + MESH.AngleMR[NA-1][1]))
							  );

	//Esquina arriba derecha
	SA_Termino1[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
								MESH.SupMP[NA-1][NRad-1][0]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1])*0.50*(RhoPres[NA-1][NRad-1] + RhoPres[NA-2][NRad-1])*(UwallsMU[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + VwallsMU[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1]))
								+ MESH.SupMP[NA-1][NRad-1][1]*vright[NRad-1]*RhoRight[NRad-1]*(UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + VwallsMU[NA][NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))
								+ MESH.SupMP[NA-1][NRad-1][2]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2])*0.50*(RhoPres[NA-1][NRad-1] + RhoPres[NA-1][NRad-2])*(UwallsMR[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + VwallsMR[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))
								+ MESH.SupMP[NA-1][NRad-1][3]*vup[NA-1]*RhoUp[NA-1]*(UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))
							  );
}

//Cálculo del término 3 de la ecuación de Spalart-Allmaras
void Solver::SA_Term2(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			SA_Termino2[i][j] = Cb1*Smodel[i][j]*RhoPres[i][j]*vpresent[i][j];
		}
	}
}

//Cálculo del término 2 de la ecuación de Spalart-Allmaras
void Solver::SA_Term3(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino3[i][j] = (1.0/(Sigma*MESH.VolMP[i][j]))*(
								(MESH.SupMP[i][j][0]/MESH.DeltasMU[i][j][0])*muWallsMU[i][j]*(vpresent[i][j] - vpresent[i-1][j])*cos(PI + MESH.AngleMU[i][j])
							  + (MESH.SupMP[i][j][1]/MESH.DeltasMU[i+1][j][0])*muWallsMU[i+1][j]*(vpresent[i+1][j] - vpresent[i][j])*cos(MESH.AngleMU[i+1][j])
							  + (MESH.SupMP[i][j][2]/MESH.DeltasMR[i][j][1])*muWallsMR[i][j]*(vpresent[i][j] - vpresent[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j])
							  + (MESH.SupMP[i][j][3]/MESH.DeltasMR[i][j+1][1])*muWallsMR[i][j+1]*(vpresent[i][j+1] - vpresent[i][j])*sin(0.50*PI + MESH.AngleMR[i][j+1])	
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino3[0][j] = (1.0/(Sigma*MESH.VolMP[0][j]))*(
								(MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*muWallsMU[0][j]*(vpresent[0][j] - vleft[j])*cos(PI + MESH.AngleMU[0][j])
							  + (MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*muWallsMU[1][j]*(vpresent[1][j] - vpresent[0][j])*cos(MESH.AngleMU[1][j])
							  + (MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*muWallsMR[0][j]*(vpresent[0][j] - vpresent[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j])
							  + (MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*muWallsMR[0][j+1]*(vpresent[0][j+1] - vpresent[0][j])*sin(0.50*PI + MESH.AngleMR[0][j+1])	
							  );

		//Parte derecha
		SA_Termino3[NA-1][j] = (1.0/(Sigma*MESH.VolMP[NA-1][j]))*(
								(MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*muWallsMU[NA-1][j]*(vpresent[NA-1][j] - vpresent[NA-2][j])*cos(PI + MESH.AngleMU[i][j])
							  + (MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*muWallsMU[NA][j]*(vright[j] - vpresent[NA-1][j])*cos(MESH.AngleMU[NA][j])
							  + (MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*muWallsMR[NA-1][j]*(vpresent[NA-1][j] - vpresent[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])
							  + (MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*muWallsMR[NA-1][j+1]*(vpresent[NA-1][j+1] - vpresent[NA-1][j])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])	
							  );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino3[i][0] = (1.0/(Sigma*MESH.VolMP[i][0]))*(
								(MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*muWallsMU[i][0]*(vpresent[i][0] - vpresent[i-1][0])*cos(PI + MESH.AngleMU[i][0])
							  + (MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*muWallsMU[i+1][0]*(vpresent[i+1][0] - vpresent[i][0])*cos(MESH.AngleMU[i+1][0])
							  + (MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*muWallsMR[i][0]*(vpresent[i][0] - vdown[i])*sin(1.50*PI + MESH.AngleMR[i][0])
							  + (MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*muWallsMR[i][1]*(vpresent[i][1] - vpresent[i][0])*sin(0.50*PI + MESH.AngleMR[i][1])	
							  );

		//Parte arriba
		SA_Termino3[i][NRad-1] = (1.0/(Sigma*MESH.VolMP[i][NRad-1]))*(
								(MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*muWallsMU[i][NRad-1]*(vpresent[i][NRad-1] - vpresent[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1])
							  + (MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*muWallsMU[i+1][NRad-1]*(vpresent[i+1][NRad-1] - vpresent[i][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1])
							  + (MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*muWallsMR[i][NRad-1]*(vpresent[i][NRad-1] - vpresent[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])
							  + (MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*muWallsMR[i][NRad]*(vup[i] - vpresent[i][NRad-1])*sin(0.50*PI + MESH.AngleMR[i][NRad])	
							  );
	}

	//Esquina abajo izquierda
	SA_Termino3[0][0] = (1.0/(Sigma*MESH.VolMP[0][0]))*(
								(MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*muWallsMU[0][0]*(vpresent[0][0] - vleft[0])*cos(PI + MESH.AngleMU[0][0])
							  + (MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*muWallsMU[1][0]*(vpresent[1][0] - vpresent[0][0])*cos(MESH.AngleMU[1][0])
							  + (MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*muWallsMR[0][0]*(vpresent[0][0] - vdown[0])*sin(1.50*PI + MESH.AngleMR[0][0])
							  + (MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*muWallsMR[0][1]*(vpresent[0][1] - vpresent[0][0])*sin(0.50*PI + MESH.AngleMR[0][1])	
							  );

	//Esquina arriba izquierda
	SA_Termino3[0][NRad-1] = (1.0/(Sigma*MESH.VolMP[0][NRad-1]))*(
								(MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*muWallsMU[0][NRad-1]*(vpresent[0][NRad-1] - vleft[NRad-1])*cos(PI + MESH.AngleMU[0][NRad-1])
							  + (MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*muWallsMU[1][NRad-1]*(vpresent[1][NRad-1] - vpresent[0][NRad-1])*cos(MESH.AngleMU[1][NRad-1])
							  + (MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*muWallsMR[0][NRad-1]*(vpresent[0][NRad-1] - vpresent[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])
							  + (MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*muWallsMR[0][NRad]*(vup[0] - vpresent[0][NRad-1])*sin(0.50*PI + MESH.AngleMR[0][NRad])	
							  );

	//Esquina abajo derecha
	SA_Termino3[NA-1][0] = (1.0/(Sigma*MESH.VolMP[NA-1][0]))*(
								(MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*muWallsMU[NA-1][0]*(vpresent[NA-1][0] - vpresent[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0])
							  + (MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*muWallsMU[NA][0]*(vright[0] - vpresent[NA-1][0])*cos(MESH.AngleMU[NA][0])
							  + (MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*muWallsMR[NA-1][0]*(vpresent[NA-1][0] - vdown[NA-1])*sin(1.50*PI + MESH.AngleMR[NA-1][0])
							  + (MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*muWallsMR[NA-1][1]*(vpresent[NA-1][1] - vpresent[NA-1][0])*sin(0.50*PI + MESH.AngleMR[NA-1][1])	
							  );

	//Esquina arriba derecha
	SA_Termino3[NA-1][NRad-1] = (1.0/(Sigma*MESH.VolMP[NA-1][NRad-1]))*(
								(MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*muWallsMU[NA-1][NRad-1]*(vpresent[NA-1][NRad-1] - vpresent[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1])
							  + (MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*muWallsMU[NA][NRad-1]*(vright[NRad-1] - vpresent[NA-1][NRad-1])*cos(MESH.AngleMU[NA][NRad-1])
							  + (MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*muWallsMR[NA-1][NRad-1]*(vpresent[NA-1][NRad-1] - vpresent[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])
							  + (MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*muWallsMR[NA-1][NRad]*(vup[NA-1] - vpresent[NA-1][NRad-1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])	
							  );

}

//Cálculo del término 4 de la ecuación de Spalart-Allmaras
void Solver::SA_Term4(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino4[i][j] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[i][j] + vpresent[i+1][j]) - 0.50*(vpresent[i][j] + vpresent[i-1][j]))/MESH.DeltasMP[i][j][0],2.0)
							  + pow((0.50*(vpresent[i][j] + vpresent[i][j+1]) - 0.50*(vpresent[i][j] + vpresent[i][j-1]))/MESH.DeltasMP[i][j][1],2.0)		
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino4[0][j] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[0][j] + vpresent[1][j]) - vleft[j])/MESH.DeltasMP[0][j][0],2.0)
							  + pow((0.50*(vpresent[0][j] + vpresent[0][j+1]) - 0.50*(vpresent[0][j] + vpresent[0][j-1]))/MESH.DeltasMP[0][j][1],2.0)		
							  );

		//Parte derecha
		SA_Termino4[NA-1][j] = ((Cb2)/Sigma)*(
								pow((vright[j] - 0.50*(vpresent[NA-1][j] + vpresent[NA-2][j]))/MESH.DeltasMP[NA-1][j][0],2.0)
							  + pow((0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1]) - 0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1]))/MESH.DeltasMP[i][j][1],2.0)		
							  );

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino4[i][0] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[i][0] + vpresent[i+1][0]) - 0.50*(vpresent[i][0] + vpresent[i-1][0]))/MESH.DeltasMP[i][0][0],2.0)
							  + pow((0.50*(vpresent[i][0] + vpresent[i][1]) - vdown[i])/MESH.DeltasMP[i][0][1],2.0)		
							  );

		//Parte arriba
		SA_Termino4[i][NRad-1] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1]) - 0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1]))/MESH.DeltasMP[i][NRad-1][0],2.0)
							  + pow((vup[i] - 0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2]))/MESH.DeltasMP[i][NRad-1][1],2.0)		
							  );
	
	}

	//Esquina abajo izquierda
	SA_Termino4[0][0] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[0][0] + vpresent[1][0]) - vleft[0])/MESH.DeltasMP[0][0][0],2.0)
							  + pow((0.50*(vpresent[0][0] + vpresent[0][1]) - vdown[0])/MESH.DeltasMP[0][0][1],2.0)		
							  );

	//Esquina arriba izquierda
	SA_Termino4[0][NRad-1] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1]) - vleft[NRad-1])/MESH.DeltasMP[0][NRad-1][0],2.0)
							  + pow((vup[0] - 0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2]))/MESH.DeltasMP[0][NRad-1][1],2.0)		
							  );

	//Esquina abajo derecha
	SA_Termino4[NA-1][0] = ((Cb2)/Sigma)*(
								pow((vright[0] - 0.50*(vpresent[NA-1][0] + vpresent[NA-2][0]))/MESH.DeltasMP[NA-1][0][0],2.0)
							  + pow((0.50*(vpresent[NA-1][0] + vpresent[NA-1][1]) - vdown[NA-1])/MESH.DeltasMP[NA-1][0][1],2.0)		
							  );

	//Esquina arriba derecha
	SA_Termino4[NA-1][NRad-1] = ((Cb2)/Sigma)*(
								pow((vright[NRad-1] - 0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1]))/MESH.DeltasMP[NA-1][NRad-1][0],2.0)
							  + pow((vup[NA-1] - 0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2]))/MESH.DeltasMP[NA-1][NRad-1][1],2.0)		
							  );

}

//Cálculo del término 5 de la ecuación de Spalart-Allmaras
void Solver::SA_Term5(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			SA_Termino5[i][j] = -Cw1*RhoPres[i][j]*fw[i][j]*pow((vpresent[i][j])/MESH.minDist[i][j],2.0);
		}
	}

}

//Cálculo de todas las contribuciones del modelo Spalart-Allmaras
void Solver::FmuSA(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Fmupresent[i][j] = 
							 + SA_Termino1[i][j] 
							 + SA_Termino2[i][j] 
							 + SA_Termino3[i][j] 
							 + SA_Termino4[i][j] 
							 + SA_Termino5[i][j]					
							 ;
		}
	}

}

//Inicialización matrices del modelo Spalart Allmaras
void Solver::InitialSpalartAllmaras(){
int i,j;	
	
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad; j++){
				X[i][j] = vpresent[i][j]/(muBase[i][j]/RhoPres[i][j] + 1e-10);
				fv1[i][j] = pow(X[i][j],3.0)/(pow(X[i][j],3.0) + pow(Cv1,3.0) + 1e-10);

			}
		}
}

//Cálculo de variables y matrices necesarias para aplicar el modelo Spalart-Allmaras
void Solver::SpalartAllmarasPreparation(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){

			
			omega[i][j] = 0.50*((UwallsMR[i][j+1] - UwallsMR[i][j])/MESH.DeltasMP[i][j][1] - (VwallsMU[i+1][j] - VwallsMU[i][j])/MESH.DeltasMP[i][j][0]);
			S[i][j] = sqrt(2.0*omega[i][j]*omega[i][j]);
			fv2[i][j] = 1.0 - X[i][j]/(1.0 + X[i][j]*fv1[i][j]);
			Smodel[i][j] = S[i][j] + (vpresent[i][j]/(pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10))*fv2[i][j];

			
			//if(vpresent[i][j]/(Smodel[i][j]*pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10) < 10.0){
				r[i][j] = vpresent[i][j]/(Smodel[i][j]*pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10);
			//}
			//else{
			//	r[i][j] = 10.0;
			//}
			
			
			
		
			
			g[i][j] = r[i][j] + Cw2*(pow(r[i][j],6.0) - r[i][j]);
			
			fw[i][j] = g[i][j]*pow((1.0 + pow(Cw3,6.0))/(pow(g[i][j],6.0) + pow(Cw3,6.0) + 1e-10),1.0/6.0);
			
			
			ft2[i][j] = Ct3*exp(-Ct4*pow(X[i][j],2.0));
			
		}
	}

}

//Cálculo de la viscosidad con el modelo Spalart-Allmaras
void Solver::SpalartAllmaras(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){

			vfut[i][j] = (2.0*TimeBetta*vpresent[i][j] - (TimeBetta - 0.50)*vprevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*Fmupresent[i][j] - TimeBetta*Fmuprevious[i][j])/RhoPres[i][j];

			X[i][j] = vpresent[i][j]/(muBase[i][j]/RhoPres[i][j] + 1e-10);

			fv1[i][j] = pow(X[i][j],3.0)/(pow(X[i][j],3.0) + pow(Cv1,3.0) + 1e-10);
	
			muTurb[i][j] = RhoPres[i][j]*vfut[i][j]*fv1[i][j];
			
		}
	}

}