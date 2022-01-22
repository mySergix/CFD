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

#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ReadData.h"

using namespace std;

#define DIRECTORIO "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/InputData/"

//Constructor del lector de datos
ReadData::ReadData(Memory M1, int N){

	GeometryData = M1.AllocateDouble(3, 1, 1); //Datos de la geometría del problema
	NumericalData = M1.AllocateInt(4, 1, 1); //Datos numéricos del problema
	ProblemData = M1.AllocateDouble(2, 1, 1); //Datos del problema
	ProblemPhysicalData = M1.AllocateDouble(8, 1, 1); //Datos físicos sobre las condiciones del problema

	MeshStudyNX = M1.AllocateInt(N, 1, 1);
	MeshStudyNY = M1.AllocateInt(N, 1, 1);

	AttackAngle = M1.AllocateDouble(14, 1, 1);
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
	ReadArrays("GeometryData.txt", 3, GeometryData); //Input Datos Geometría del problema
	ReadArrays("PhysicalData.txt", 8, ProblemPhysicalData);	  //Input Datos físicos de las condiciones del problema
	ReadArrays("AttackAngleData.txt", 14, AttackAngle); 
	
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
					if(i < 4){
						if (iss >> NumericalData[i]){ i++; }	
					}
					else if(i >= 4 && i < 6){
						if (iss >> ProblemData[i-4]){ i++; }	
					}
					else{
						if (iss >> Problema){ i++; }	
					}
					
        		 }
   	 	}
			
    	DatosProblema.close();	

}

void ReadData::ReadMeshStudy(string FileName, int TotalData){
int a, b;

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
						if (iss >> MeshStudyNX[i] >> MeshStudyNY[i]){ i++; }	
					
        		 }
   	 	}
	
    	Data.close();

}
