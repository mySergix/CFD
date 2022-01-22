//Prueba Bash Scripting

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main(){

int N = 50;

double Vector[N+1];

int i;
double I;
double n = N;

	for(i = 0; i < N+1; i++){
		I = i;
		Vector[i] = I*(10.0/n);
	}		

	FILE *fp1;
	fp1 = fopen("/home/sergiogus/Desktop/CTTC/BashScripting/Archivo.txt","w");
	for(i = 0; i < N+1; i++){
		fprintf(fp1,"%f \n",Vector[i]);	
    	}
	fclose(fp1);

return 0;

}

