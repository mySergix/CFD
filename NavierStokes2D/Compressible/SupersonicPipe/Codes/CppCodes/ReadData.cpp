#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <sstream>

#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Memory.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/ReadData.h"
using namespace std;

//Constructor del lector de datos
ReadData::ReadData(Memory M1){

		GeometryData = M1.AllocateDouble1D(14); //Datos de la geometría del problema
		NumericalData = M1.AllocateInt1D(9); //Datos numéricos del problema
		ProblemData = M1.AllocateDouble1D(4); //Datos del problema
		HeatData = M1.AllocateDouble1D(4); //Datos físicos del problema sobre transferencia de calor
		//double *BoundaryData; //Condiciones de contorno del problema
		SpalartAllmarasData = M1.AllocateDouble1D(11); //Datos sobre el modelo de turbulencia Spalart-Allmaras

		CombustionData = M1.AllocateDouble1D(9); //Datos químicos sobre la combustión y sus especies presentes

		JannafO2 = M1.AllocateDouble1D(24); //Coeficientes función Cp(T) y H(T) del O2
		JannafH2 = M1.AllocateDouble1D(24); //Coeficientes función Cp(T) y H(T) del H2
		JannafH2O = M1.AllocateDouble1D(16); //Coeficientes función Cp(T) y H(T) del H2O

		ViscosityData = M1.AllocateDouble1D(9); //Datos sobre la viscosidad de las especies termoquímicas
}

void ReadData::ReadInputs(){
int i = 0;
	
	ifstream DatosGeometria("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/GeometryData.txt");

		if (DatosGeometria){
        		string line;
        		while (getline(DatosGeometria, line)){
        	 		istringstream iss(line);
					if(i < 14){
						if (iss >> GeometryData[i]){ i++; }	
					}
        		 }
   	 	}
		
    	DatosGeometria.close();

i = 0;

	ifstream DatosProblema("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/ProblemData.txt");

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
	
	ifstream DatosCalor("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/HeatData.txt");

		if (DatosCalor){
        		string line;
        		while (getline(DatosCalor, line)){
        	 		istringstream iss(line);
					if(i < 4){
						if (iss >> HeatData[i]){ i++; }	
					}
        		 }
   	 	}
		
    	DatosCalor.close();


i = 0;
	
	ifstream DatosSpalartAllmaras("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/Spalart_Allmaras_Data.txt");

		if (DatosSpalartAllmaras){
        		string line;
        		while (getline(DatosSpalartAllmaras, line)){
        	 		istringstream iss(line);
					if(i < 11){
						if (iss >> SpalartAllmarasData[i]){ i++; }	
					}
        		 }
   	 	}
		
    	DatosSpalartAllmaras.close();


i = 0;

	ifstream DatosJannaf("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/JanafData.txt");

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


i = 0;

	ifstream DatosCombustion("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/CombustionData.txt");

		if (DatosCombustion){
        		string line;
        		while (getline(DatosCombustion, line)){
        	 		istringstream iss(line);
					if (iss >> CombustionData[i]){ i++; }	
        		 }
   	 	}
		
    	DatosCombustion.close();



i = 0;

	ifstream DatosViscosidad("/home/sergio/Desktop/TFG/RegimenTurbulento/InputData/SutherlandLaw_Data.txt");

		if (DatosViscosidad){
        		string line;
        		while (getline(DatosViscosidad, line)){
        	 		istringstream iss(line);
					if (iss >> ViscosityData[i]){ i++; }	
        		 }
   	 	}
		
    	DatosViscosidad.close();

}



