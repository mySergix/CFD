#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h>
#include <string>
#include <mpi.h>

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ReadData.h"
using namespace std;

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/InputData/"

//Constructor del lector de datos
ReadData::ReadData(Memory M1){

	GeometryData = M1.AllocateDouble(17,1,1); //Datos de la geometría del problema
	NumericalData = M1.AllocateInt(9,1,1); //Datos numéricos del problema
	ProblemData = M1.AllocateDouble(4,1,1); //Datos del problema
	HeatData = M1.AllocateDouble(4,1,1); //Datos físicos del problema sobre transferencia de calor
	ProblemPhysicalData = M1.AllocateDouble(7, 1, 1); //Datos físicos sobre las condiciones del problema

	SpalartAllmarasData = M1.AllocateDouble(11,1,1); //Datos del modelo Spalart-Allmaras
	CombustionData = M1.AllocateDouble(9,1,1); //Datos sobre la combustión de las especies químicas
	ViscosityData = M1.AllocateDouble(9,1,1); //Datos sobre la viscosidad de las especies químicas

	JannafO2 = M1.AllocateDouble(24,1,1); //Coeficientes función Cp(T) y H(T) del O2
	JannafH2 = M1.AllocateDouble(24,1,1); //Coeficientes función Cp(T) y H(T) del H2
	JannafH2O = M1.AllocateDouble(16,1,1); //Coeficientes función Cp(T) y H(T) del H2O
}

void ReadData::ReadArrays(string FileName, int TotalData, double *Array){
int i = 0;
stringstream InitialName;
string FinalName;

	InitialName<<DIRECTORIO<<FileName;
	FinalName=InitialName.str();

	ifstream Data(FinalName.c_str());

		if (Data){
        		string line;
        		while (getline(Data, line)){
        	 		istringstream iss(line);
					if(i < TotalData){
						if (iss >> Array[i]){ i++; }	
					}
        		 }
   	 	}
		
    	Data.close();

}

void ReadData::ReadInputs(){

	//Lectura datos en Arrays
	ReadArrays("GeometryData.txt", 17, GeometryData); //Input Datos Geometría del problema
	ReadArrays("HeatData.txt", 4, HeatData); //Input Datos Transferencia de calor del problema
	ReadArrays("ProblemPhysicalData.txt", 7, ProblemPhysicalData);	  //Input Datos físicos de las condiciones del problema
	ReadArrays("Spalart_Allmaras_Data.txt", 11, SpalartAllmarasData); //Input Datos modelo Spalart-Allmaras
	ReadArrays("CombustionData.txt", 9, CombustionData); //Input Datos modelo Spalart-Allmaras
	ReadArrays("SutherlandLaw_Data.txt", 9, ViscosityData); //Input Datos modelo Spalart-Allmaras

int i = 0;
string FileName;
stringstream InitialName1;
string FinalName1;

	FileName = "ProblemData.txt";

	InitialName1<<DIRECTORIO<<FileName;
	FinalName1=InitialName1.str();
	
	ifstream DatosProblema(FinalName1.c_str());

		if (DatosProblema){
        		string line;
        		while (getline(DatosProblema, line)){
        	 		istringstream iss(line);
					if(i < 9){
						if (iss >> NumericalData[i]){ i++; }	
					}
					else if(i >= 9 && i < 13){
						if (iss >> ProblemData[i-9]){ i++; }	
					}
					else if(i == 14){
						if (iss >> TipoProblema){ i++; }	
					}
					else if(i == 16){
						if (iss >> TipoTobera){ i++; }	
					}
					else{
						if (iss >> TypeMeshing){ i++; }	
					}
					
        		 }
   	 	}
			
    	DatosProblema.close();	


i = 0;
stringstream InitialName2;
string FinalName2;

	FileName = "JanafData.txt";

	InitialName2<<DIRECTORIO<<FileName;
	FinalName2=InitialName2.str();
	
	ifstream DatosJannaf(FinalName2.c_str());
	
		if (DatosJannaf){
        		string line;
        		while (getline(DatosJannaf, line)){
        	 		istringstream iss(line);
					if(i < 24){
						if (iss >> JannafO2[i]){ i++; }	
					}
					else if(i >= 24 && i < 48){
						if (iss >> JannafH2[i - 24]){ i++; }	
					}
					else{
						if (iss >> JannafH2O[i - 48]){ i++; }	
					}
        		 }
   	 	}
		
    	DatosJannaf.close();

}

