//Prueba Bash Scripting

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

int main(){

//Problems data
	int *NumericalData;
	double *PrecisionData;
	double *GeometryData;
	double *FluidData;
	double *GravityData;
	string Scheme;
	string Scheme2;


	NumericalData = new int [8];
	PrecisionData = new double [3];
	GeometryData = new double [7];
	FluidData = new double [6];
	GravityData = new double [2];
		

int i = 0;

	ifstream ArchivoDatos("/home/sergiogus/Desktop/CTTC/BashScripting/ProblemData.txt");

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

double Reynolds = FluidData[2];
int RE = Reynolds;


	char Results[300];
	sprintf(Results, "/home/sergiogus/Desktop/CTTC/BashScripting/RE_%d/Archivo_%d.txt", RE, RE);

	FILE *fp1;
	fp1 = fopen(Results,"w");
		fprintf(fp1,"Número de Reynolds: %f \n",Reynolds);	
	fclose(fp1);

return 0;

}

