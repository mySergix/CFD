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
		//Datos sobre el problema
			string Problema; //Problema (Canal/Canal con Cilindro/Perfil)
			
		//Datos sobre la geometría del problema
			double ChannelLength; //Longitud del canal
			double ChannelHeight; //Altura del canal
			double ChannelDepth; //Profundidad del canal

		//Datos del mallado
			int OpcionR; //Tipo de discretización en la dirección radial
			int OpcionA; //Tipo de discretización en la dirección axial
			double StretFactorX; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección axial
			double StretFactorY; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección radial

	public:
		//Constructor de la clase
		Mesher(Memory, ReadData, ParPro, int);

		//Datos para la computación en paralelo
		int Procesos;
		int Rank;
		int Ix;
		int Fx;

		//Densidad del mallado
		int NX; //Dirección axial total
		int NY; //Dirección radial
		
		//Matrices de coordenadas de discretización
		double *AxialCoord; //Matriz coordenada axial de la distribución

		//Caso mallado tipo Staggered
		double *MP; //Coordenadas matriz de presión/temperatura
		double *MU; //Coordenadas matriz velocidades axiales
		double *MR; //Coordenadas matriz velocidades radiales

		//Matrices de distancias de volúmenes de control
		double *DeltasMP; //Deltas X y R de la matriz de Presión/Temperatura
		double *DeltasMU; //Deltas X y R de la matriz de velocidades axiales (U)
		double *DeltasMR; //Deltas X y R de la matriz de velocidades radiales (V)

		//Matrices de superficies de los volúmenes de control
		double *SupMP;

		//Matriz de volúmenes de control
		double *VolMP;

		//Métodos de la clase del mallador
		void AllocateMemory(Memory); //Alojamiento de memoria para cada matriz

		void Get_AxialCoordinates(); //Cálculo del la posición axial de los nodos de velocidad axial	
		void Get_Mesh(); //Creación de todas las mallas
		void Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
		void Get_Surfaces(); //Cálculo de las superficies de los volúmenes de control	
		void Get_Volumes(); //Cálculo de los volúmenes de control

		void PrintTxt();
		void MallaVTK2D(string, string, string, double*, int, int);

		void ExecuteMesher(Memory, ParPro); //Ejecutar todos los procesos del mallador
};
		
