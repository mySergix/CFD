#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <sstream>

#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/ReadTXT.h"
using namespace std;

#define PI 3.141592653589793

ReadTXT::ReadTXT(){

		NumericalData = new int [8];
		PrecisionData = new double [3];
		GeometryData = new double [7];
		FluidData = new double [6];
		GravityData = new double [2];
		
		
		BoundaryData = new double [24];
		GradientData = new int [24];

}

void ReadTXT::ReadFromTXT(){
int i = 0;

	ifstream ArchivoDatos("/home/sergiogus/Desktop/CTTC/NavierStokes/Datos iniciales/ProblemData.txt");

		if (ArchivoDatos){
        		string line;
        		while (getline(ArchivoDatos, line)){

        	 		istringstream iss(line);
					if(i < 8){
						if (iss >> NumericalData[i]){ i++; }	
					}
					else if(i >= 8 && i < 11){
						if (iss >> PrecisionData[i - 8]){ i++; }
					}
					else if(i >= 11 && i < 18){
						if (iss >> GeometryData[i - 11]){ i++; }
					}
					else if(i >= 18 && i < 24){
						if (iss >> FluidData[i - 18]){ i++; }
					} 
					else if(i >= 24 && i < 26){
						if (iss >> GravityData[i - 24]){ i++; }
					}	
					else if(i == 26){
						if (iss >> Scheme){ i++; }
					}
					else{
						if (iss >> Scheme2){ i++; }
					}
        		 }
   	 	}
	
		
    	ArchivoDatos.close();

i = 0;
	
	ifstream ArchivoCondiciones("/home/sergiogus/Desktop/CTTC/NavierStokes/Datos iniciales/BoundaryConditions.txt");
	
		if (ArchivoCondiciones){
        		string line;
        		while (getline(ArchivoCondiciones, line)){

        	 		istringstream iss(line);
					if(i%2 == 0){
						if (iss >> BoundaryData[i/2]) { i++; }
					}
					else{
						if (iss >> GradientData[i/2]) { i++; }
					}	
        		 }
   	 	}

    	ArchivoCondiciones.close();
}

int ReadTXT::Get_NumericalData(int i){ return NumericalData[i]; }

double ReadTXT::Get_PrecisionData(int i){ return PrecisionData[i]; }

double ReadTXT::Get_GeometryData(int i){ return GeometryData[i]; }

double ReadTXT::Get_FluidData(int i){ return FluidData[i]; }

double ReadTXT::Get_GravityData(int i){ return GravityData[i]; }

double ReadTXT::Get_BoundaryData(int i){ return BoundaryData[i]; }
	
int ReadTXT::Get_GradientData(int i){ return GradientData[i]; }	

string ReadTXT::Get_ConvectiveScheme(){ return Scheme; }

string ReadTXT::Get_ConvectiveScheme2(){ return Scheme2; }


