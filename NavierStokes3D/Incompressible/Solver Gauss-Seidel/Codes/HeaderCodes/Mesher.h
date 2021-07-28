#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

using namespace std;

class Mesher{
	private:

	public:
		//Constructor de la clase
		Mesher(Memory, ReadData, ParPro, string);

		//Datos Numéricos del problema
		int Problema;

		int ProgresionLineal;

		int NX;
		int NY; 
		int NZ;

		int NX1;
		int NX2; 
		int NX3;

		int NY1;
		int NY2; 
		int NY3;

		int OptionX;
		int OptionY;
		int OptionZ;

		double SFX;
		double SFY;
		double SFZ;

		//Datos Geométricos del problma
		double Xdominio;
		double Ydominio;
		double Zdominio;

		double Xcentroide;
		double Ycentroide;

		double Xcuadrado;
		double Ycuadrado;

		//Parámetros de computación paralela
		int Rank;
		int Procesos;
		int Ix;
		int Fx;
		
		int HaloPressure;
		int HaloU;
		int HaloV; 
		int HaloW;
		
		int HP;
		
		string DIRECTORIO;
		
		//Matrices del Mallado del problema
		double *MP; //Coordenadas matriz de presión/temperatura
		double *MU; //Coordenadas matriz velocidad U
		double *MV; //Coordenadas matriz velocidad V
		double *MW; //Coordenadas matriz velocidad W

		//Matrices de distancias de volúmenes de control
		double *DeltasMP; //Deltas de la matriz de Presión/Temperatura
		double *DeltasMU; //Deltas de la matriz de velocidades (U)
		double *DeltasMV; //Deltas de la matriz de velocidades (V)
		double *DeltasMW; //Deltas de la matriz de velocidades (W)

		//Matrices de superficies de volúmenes de control
		double *SupMP; //Superficies de los volúmenes de la matriz de Presión/Temperatura
		double *SupMU; //Superficies de los volúmenes de la matriz de velocidades (U)
		double *SupMV; //Superficies de los volúmenes de la matriz de velocidades (V)
		double *SupMW; //Superficies de los volúmenes de la matriz de velocidades (V)

		//Matrices de volúmenes de los volúmenes de control
		double *VolMP; //Volúmenes de los volúmenes de la matriz de Presión/Temperatura
		double *VolMU; //Volúmenes de los volúmenes de la matriz de velocidades (U)
		double *VolMV; //Volúmenes de los volúmenes de la matriz de velocidades (V)
		double *VolMW; //Volúmenes de los volúmenes de la matriz de velocidades (W)

		void AllocateMemory(Memory); //Alojamiento de memoria para cada matriz
		void Get_Meshes(); //Creación de todas las mallas
		void Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
		void Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control
		void Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen
		
		void PrintTxt(); //Pasar los resultados de las mallas a un txt
		void MallaVTK3D(string, string, string, double*, int, int, int); //Pasar todos los resultados escalares en 2D a un archivo VTK
		

		void ExecuteMesher(Memory, ParPro); //Ejecutar todos los procesos del mallador
};
		
