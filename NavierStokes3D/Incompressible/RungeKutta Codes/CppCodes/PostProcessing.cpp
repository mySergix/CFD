#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Mesher.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Solver.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/PostProcessing.h"

#define PI 3.141592653589793

#define DIRECTORIO "/home/sergiogus/Desktop/CFD/Incompressible3D/"

#define GP(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+2*HP)*(NY+2*HP)*(NZ+2*HP)*(dim)) //Global Index P Mesh
#define GU(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+1+2*HP)*(NY+2*HP)*(NZ+2*HP)*(dim)) //Global Index U Mesh
#define GV(i,j,k,dim) ((NY+1+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+1+2*HP)) + ((NX+2*HP)*(NY+1+2*HP)*(NZ+2*HP)*(dim)) //Global Index V Mesh
#define GW(i,j,k,dim) ((NY+2*HP)*(NZ+1+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) + ((NX+2*HP)*(NY+2*HP)*(NZ+1+2*HP)*(dim)) //Global Index W Mesh

//Local Index P Mesh
#define LPL(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Left Core
#define LPC(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Center Cores
#define LPR(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Right Core

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

PostProc::PostProc(Memory M1, ReadData R1, Mesher MESH, Solver S1){
	
    Problema = R1.ProblemNumericalData[0];

    //Datos necesarios para computación paralela
	Rank = MPI1.Rank;
	Procesos = MPI1.Procesos;

    NX = MESH.NX;
    NY = MESH.NY;
    NZ = R1.ProblemNumericalData[4];

    Reynolds = R1.ProblemPhysicalData[2];
    Rayleigh = R1.ProblemPhysicalData[3];

    Tleft = R1.ProblemPhysicalData[9]; 
	Tright = R1.ProblemPhysicalData[10];  

    //Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	Xcentroide = R1.GeometryData[3];
	Ycentroide = R1.GeometryData[4];

	Xcuadrado = R1.GeometryData[5];
	Ycuadrado = R1.GeometryData[6];

}

//Pasar los resultados de variables escalares a un archivo VTK en 3D
void PostProc::EscalarVTK3D(string Carpeta, string Variable, string NombreFile, double *PropertyMatrix, double *MC, int Nx, int Ny, int Nz){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<DIRECTORIO<<"ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<Nx<<"   "<<Ny<<"   "<<Nz<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<Nx*Ny*Nz<<"   "<<"double"<<endl;
	
	for(k = 0; k < Nz; k++){
		for(j = 0; j < Ny; j++){
			for(i = 0; i < Nx; i++){
				file<<MC[GP(i,j,k,0)]<<"   "<<MC[GP(i,j,k,1)]<<"   "<<MC[GP(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
	file<<"POINT_DATA"<<"   "<<Nx*Ny*Nz<<endl;
    file<<"SCALARS "<<Variable<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
    file<<endl;
	for(k = 0; k < Nz; k++){
		for(j = 0; j < Ny; j++){
			for(i = 0; i < Nx; i++){
				file<<PropertyMatrix[GP(i,j,k,0)]<<" ";
			}
		}
	}

    file.close();
}

//Pasar los resultados de variables vectoriales a un archivo VTK en 3D
void PostProc::VectorialVTK3D(Mesher MESH, string Carpeta, string Variable, string NombreFile, double *Field1, double *Field2, double *Field3, double *MC, int Nx, int Ny, int Nz){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<DIRECTORIO<<"ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<Nx<<"   "<<Ny<<"   "<<Nz<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<Nx*Ny*Nz<<"   "<<"double"<<endl;
	
	for(k = 0; k < Nz; k++){
		for(j = 0; j < Ny; j++){
			for(i = 0; i < Nx; i++){
				file<<MC[GP(i,j,k,0)]<<"   "<<MC[GP(i,j,k,1)]<<"   "<<MC[GP(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
    file<<"POINT_DATA"<<"   "<<(Nx)*(Ny)*(Nz)<<endl;
    file<<"VECTORS "<<Variable<<" double"<<endl;
    file<<endl;

	for(k = 0; k < Nz; k++){
		for(j = 0; j < Ny; j++){
			for(i = 0; i < Nx; i++){
				file<<0.50*(Field1[GU(i,j,k,0)] + Field1[GU(i+1,j,k,0)])<<" "<<0.50*(Field2[GV(i,j,k,0)] + Field2[GV(i,j+1,k,0)])<<" "<<0.50*(Field3[GW(i,j,k,0)] + Field3[GW(i,j,k+1,0)])<<endl;
			}
		}
	}
//
    file.close();

}

//Pasar los resultados del Driven Cavity a un TXT
void PostProc::Get_DrivenResults(Solver S1, Mesher MESH){
int i, j, k;
string Carpeta = "GnuPlotResults/Results/DrivenCavity/";
ofstream file;
string FileName;
stringstream InitialNameMP;
string FinalNameMP;
int RE = S1.Reynolds;
char Name[200];

	sprintf(Name, "DrivenCavityResults_RE_%d.txt", RE);

	FileName = Name;

	InitialNameMP<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMP = InitialNameMP.str();
    file.open(FinalNameMP.c_str());

	file<<"Simulación Driven Cavity"<<endl;		
	file<<endl;
	file<<"Reynolds: "<<RE<<endl;
	file<<"Tiempo de Simulación: "<<S1.Time<<endl;
	file<<"Tiempo de Cómputo total (s): "<<S1.RealTime<<endl;
	file<<"Número de Procesadores Utilizados: "<<S1.Procesos + 1<<endl;
	file<<"Número de Steps: "<<S1.Step<<endl;
	file<<"Número de Nodos en X: "<<NX<<endl;
	file<<"Número de Nodos en Y: "<<NY<<endl;
	file<<"Número de Nodos en Z: "<<NZ<<endl;

	file.close();


//Pasar a un TXT la línea de Velocidades U en la vertical de la cavidad
double MeanValue;
stringstream InitialNameVelocidadesU;
string FinalNameVelocidadesU;

	sprintf(Name, "VelocidadesUResults_RE_%d.txt", RE);

	FileName = Name;

	InitialNameVelocidadesU<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameVelocidadesU = InitialNameVelocidadesU.str();
    file.open(FinalNameVelocidadesU.c_str());

    for(j = 0; j < NY; j++){
		MeanValue = 0.0;

		for(k = 0; k < NZ; k++){
			MeanValue += S1.UGPRES[GU(NX/2,j,k,0)]*(MESH.DeltasMU[GU(NX/2,j,k,2)]/Zdominio);
		}		

		file<<MESH.MU[GU(NX/2,j,NZ/2,1)]<<"\t"<<MeanValue<<endl;
	}   

	file.close();

//Pasar a un TXT la línea de Velocidades V en la horizontal de la cavidad
stringstream InitialNameVelocidadesV;
string FinalNameVelocidadesV;

	sprintf(Name, "VelocidadesVResults_RE_%d.txt", RE);

	FileName = Name;

	InitialNameVelocidadesV<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameVelocidadesV = InitialNameVelocidadesV.str();
    file.open(FinalNameVelocidadesV.c_str());

    for(i = 0; i < NX; i++){
		MeanValue = 0.0;

		for(k = 0; k < NZ; k++){
			MeanValue += S1.VGPRES[GV(i,NY/2,k,0)]*(MESH.DeltasMV[GV(i,NY/2,k,2)]/Zdominio);
		}		

		file<<MESH.MV[GV(i,NY/2,NZ/2,0)]<<"\t"<<MeanValue<<endl;
	}   

	file.close();

}

//Pasar los resutlados del Differentially Heated a un TXT
void PostProc::Get_DifferentiallyResults(Solver S1, Mesher MESH){
int i, j, k;

int RE = S1.Reynolds;
int RA = S1.Rayleigh;


double MeanValue;

double *LocalNusselt;
LocalNusselt = new double [NY];

double NusseltMax, NusseltMin, NusseltMean;
double Y_NusseltMax, Y_NusseltMin;

	for(j = 0; j < NY; j++){
		MeanValue = 0.0;
		for(k = 0; k < NZ; k++){

			MeanValue += ((Xdominio/(Tleft - Tright))*((Tleft - S1.TGPRES[GP(0,j,k,0)])/MESH.DeltasMU[GU(0,j,k,0)]))*(MESH.DeltasMP[GP(0,j,k,2)]/Zdominio);

		}
		LocalNusselt[j] = MeanValue;

	}

	NusseltMax = 0.0;
	NusseltMin = 1000.0;
	NusseltMean = 0.0;

	for(j = 0; j < NY; j++){

		if(LocalNusselt[j] >= NusseltMax){
			NusseltMax = LocalNusselt[j];
			Y_NusseltMax = MESH.MP[GP(0,j,NZ/2,1)];
		}

		if(LocalNusselt[j] <= NusseltMin){
			NusseltMin = LocalNusselt[j];
			Y_NusseltMin = MESH.MP[GP(0,j,NZ/2,1)];
		}

		NusseltMean += LocalNusselt[j]*(MESH.MP[GP(0,j,NZ/2,1)]/Ydominio);

	}

double Umax, Vmax;
double Y_Umax, X_Vmax;

	double *LocalUmax;
	LocalUmax = new double [NY];

	double *LocalVmax;
	LocalVmax = new double [NX];

	for(j = 0; j < NY; j++){
		MeanValue = 0.0;

		for(k = 0; k < NZ; k++){
			MeanValue += S1.UGPRES[GU(NX/2,j,k,0)]*(MESH.MU[GU(NX/2,j,k,2)]/Zdominio);
		}

		LocalUmax[j] = MeanValue;
	}

	for(i = 0; i < NX; i++){
		MeanValue = 0.0;

		for(k = 0; k < NZ; k++){
			MeanValue += S1.VGPRES[GV(i,NY/2,k,0)]*(MESH.MV[GV(i,NY/2,k,2)]/Zdominio);
		}

		LocalVmax[j] = MeanValue;
	}


	Umax = 0.0;
	Vmax = 0.0;

	for(j = 0; j < NY; j++){

		if(LocalUmax[j] >= Umax){
			Umax = LocalUmax[j];
			Y_Umax = MESH.MU[GU(0,j,NZ/2,0)];
		}

	}

	for(i = 0; i < NX; i++){

		if(LocalVmax[i] >= Vmax){
			Vmax = LocalVmax[i];
			X_Vmax = MESH.MV[GV(i,NY/2,NZ/2,1)];
		}

	}
	
string Carpeta = "GnuPlotResults/Results/DifferentiallyHeated/";
ofstream file;
string FileName;
stringstream InitialNameMP;
string FinalNameMP;
char Name[200];

	sprintf(Name, "DifferentiallyHeatedResults_RA_%d.txt", RA);

	FileName = Name;

	InitialNameMP<<DIRECTORIO<<Carpeta<<FileName;

	FinalNameMP = InitialNameMP.str();
    file.open(FinalNameMP.c_str());

	file<<"Simulación Differentially Heated"<<endl;		
	file<<endl;
	file<<"Rayleigh: "<<RA<<endl;
	file<<"Prandtl: "<<S1.Prandtl<<endl;
	file<<"Tiempo de Simulación: "<<S1.Time<<endl;
	file<<"Tiempo de Cómputo Total: "<<S1.RealTime<<endl;
	file<<"Número de Procesadores Utilizados: "<<S1.Procesos + 1<<endl;
	file<<"Número de Steps: "<<S1.Step<<endl;
	file<<"Número de Nodos en X: "<<NX<<endl;
	file<<"Número de Nodos en Y: "<<NY<<endl;
	file<<"Número de Nodos en Z: "<<NZ<<endl;

	file<<endl;
	file<<"Resultados del problema:"<<endl;
	file<<endl;
	file<<"Número de Nusselt Máximo: "<<NusseltMax<<endl;
	file<<"Posición Vertical del Número de Nusselt Máximo: "<<Y_NusseltMax<<endl;
	file<<endl;
	file<<"Número de Nusselt Mínimo: "<<NusseltMin<<endl;
	file<<"Posición Vertical del Número de Nusselt Mínimo: "<<Y_NusseltMin<<endl;
	file<<endl;
	file<<"Número de Nusselt Medio: "<<NusseltMean<<endl;

	file<<endl;
	file<<"Velocidad U Máxima en la línea central vertical de la cavidad: "<<Umax<<endl;
	file<<"Posición de la Velocidad U Máxima en la línea central vertical de la cavidad: "<<Y_Umax<<endl;

	file<<endl;
	file<<"Velocidad V Máxima en la línea central horizontal de la cavidad: "<<Vmax<<endl;
	file<<"Posición de la Velocidad V Máxima en la línea central horizontal de la cavidad: "<<X_Vmax<<endl;

	file.close();

	delete LocalNusselt;
	delete LocalUmax;
	delete LocalVmax;

}

//Borrado de toda la memoria alocada para el mallador y el solver
void PostProc::DeleteEverything(Solver S1, Mesher MESH, ReadData R1){

    //Memoria de la lectura de datos de entrada
    delete R1.GeometryData; //Datos de la geometría del problema
	delete R1.ProblemNumericalData; //Datos numéricos del problema
	delete R1.ProblemData; //Datos del problema
	delete R1.ProblemPhysicalData; //Datos físicos sobre las condiciones del problema


    //Memoria del mallador

    //Matrices del Mallado del problema
	delete MESH.MP; //Coordenadas matriz de presión/temperatura
	delete MESH.MU; //Coordenadas matriz velocidad U
	delete MESH.MV; //Coordenadas matriz velocidad V
	delete MESH.MW; //Coordenadas matriz velocidad W

    //Matrices de distancias de volúmenes de control
	delete MESH.DeltasMP; //Deltas de la matriz de Presión/Temperatura
	delete MESH.DeltasMU; //Deltas de la matriz de velocidades (U)
	delete MESH.DeltasMV; //Deltas de la matriz de velocidades (V)
	delete MESH.DeltasMW; //Deltas de la matriz de velocidades (W)

    //Matrices de superficies de volúmenes de control
	delete MESH.SupMP; //Superficies de los volúmenes de la matriz de Presión/Temperatura
	delete MESH.SupMU; //Superficies de los volúmenes de la matriz de velocidades (U)
	delete MESH.SupMV; //Superficies de los volúmenes de la matriz de velocidades (V)
	delete MESH.SupMW; //Superficies de los volúmenes de la matriz de velocidades (V)

    //Matrices de volúmenes de los volúmenes de control
	delete MESH.VolMP; //Volúmenes de los volúmenes de la matriz de Presión/Temperatura
	delete MESH.VolMU; //Volúmenes de los volúmenes de la matriz de velocidades (U)
	delete MESH.VolMV; //Volúmenes de los volúmenes de la matriz de velocidades (V)
	delete MESH.VolMW; //Volúmenes de los volúmenes de la matriz de velocidades (W)

    //Matrices del Solver

		//Matrices locales de propiedades

		//Presion	
		delete S1.PLPRES;
		delete S1.PLFUT;
		delete S1.PLFUTsup;

		//Velocidad U
		delete S1.ULPAS;
		delete S1.ULPRES;
		delete S1.ULFUT;

		//Velocidad V
		delete S1.VLPAS;
		delete S1.VLPRES;
		delete S1.VLFUT;

		//Velocidad W
		delete S1.WLPAS;
		delete S1.WLPRES;
		delete S1.WLFUT;

		//Matrices de las contribuciones de cada una de las ecuaciones
		delete S1.ConvectiveU;
		delete S1.ConvectiveV;
		delete S1.ConvectiveW;

		delete S1.DiffusiveU;
		delete S1.DiffusiveV;
		delete S1.DiffusiveW;

		delete S1.UcontributionPast;
		delete S1.VcontributionPast;
		delete S1.WcontributionPast;

		delete S1.UcontributionPres;
		delete S1.VcontributionPres;
		delete S1.WcontributionPres;

		

		if(Problema == 2){ //Problema Differentially Heated

            delete S1.BoussinesqU;
	    	delete S1.BoussinesqV;
	    	delete S1.BoussinesqW; 

			//Temperatura
			delete S1.TLPAS;
			delete S1.TLPRES;
			delete S1.TLFUT;

			//Matrices de las contribuciones de cada una de las ecuaciones
			delete S1.ConvectiveT;

			delete S1.DiffusiveT;

			delete S1.TcontributionPast;
			delete S1.TcontributionPres;

		}
			
	
    if(Rank == 0){

		//Matrices globales de propiedades

		//Presion
		delete S1.PGPRES; //Presión P Global Presente
		delete S1.PGFUT; //Presión P Global Presente

		//Velocidad U
		delete S1.UGPRES; //Velocidad U Global Presente
		delete S1.UGFUT; //Velocidad U Global Futuro

		//Velocidad V
		delete S1.VGPRES; //Velocidad V Global Presente
		delete S1.VGFUT; //Velocidad V Global Futuro

		//Velocidad W
		delete S1.WGPRES; //Velocidad W Global Presente
		delete S1.WGFUT; //Velocidad W Global Futuro

		delete S1.PDT; //Posibles Delta T (Step Time)
		delete S1.PDiff; //Posibles Delta T (Step Time)

        if(Problema == 2){

            //Temperatura globales
			delete S1.TGPRES; //Temperatura T Global Presente
			delete S1.TGFUT; //Temperatura T Global Futuro

        }
    }

    //Matrices de los coeficientes de discretización de presión
	delete S1.aw;
	delete S1.ae;

	delete S1.as;
	delete S1.an;

	delete S1.ah;
	delete S1.at;

	delete S1.ap;
	delete S1.bp;

	//Matrices de condiciones de contorno

	if(Rank == 0){
		delete S1.Uleft;
		delete S1.Vleft;
		delete S1.Wleft;
		delete S1.TLEFT;
	}	

	if(Rank == Procesos - 1){
		delete S1.Uright;
		delete S1.Vright;
		delete S1.Wright;
		delete S1.TRIGHT;
	}

	//Parte Abajo
	delete S1.Ubot;
	delete S1.Vbot;
	delete S1.Wbot;
	delete S1.TBOT;

	//Parte Arriba
	delete S1.Utop;
	delete S1.Vtop;
	delete S1.Wtop;
	delete S1.TTOP;

	//Parte Here
	delete S1.Uhere;
	delete S1.Vhere;
	delete S1.Where;
	delete S1.There;

	//Parte There
	delete S1.Uthere;
	delete S1.Vthere;
	delete S1.Wthere;
	delete S1.Tthere;

}

//Ejecución del postprocesado de los resultados obtenidos
void PostProc::ExecutePostProcessing(Solver S1, Mesher MESH, Memory M1, ParPro MPI1){
char File1Name[300];
char File2Name[300];

    sprintf(File1Name, "MapaPresiones_Step_%d", S1.Step);
	sprintf(File2Name, "MapaVelocidades_Step_%d", S1.Step);

	if(Problema == 1){ //Problema Driven Cavity

		EscalarVTK3D("DrivenCavity/", "Presion", File1Name, S1.PGFUT, MESH.MP, NX, NY, NZ);
		VectorialVTK3D(MESH, "DrivenCavity/", "Velocidad", File2Name, S1.UGFUT, S1.VGFUT, S1.WGFUT, MESH.MP, NX, NY, NZ);
		Get_DrivenResults(MESH, S1.Time, S1.Step);

	}
	else if(Problema == 2){ //Problema Differentially Heated

		EscalarVTK3D("DifferentiallyHeated/", "Presion", File1Name, S1.PGFUT, MESH.MP, NX, NY, NZ);
		sprintf(File1Name, "MapaTemperaturas_Step_%d", S1.Step);
		EscalarVTK3D("DifferentiallyHeated/", "Temperatura", File1Name, S1.TGFUT, MESH.MP, NX, NY, NZ);

		VectorialVTK3D(MESH, "DifferentiallyHeated/", "Velocidad", File2Name, S1.UGFUT, S1.VGFUT, S1.WGFUT, MESH.MP, NX, NY, NZ);
		Get_DifferentiallyResults(MESH, S1.Time, S1.Step);
	}

}


