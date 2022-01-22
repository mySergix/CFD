#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
using namespace std;

class Mesher{
	private:
		//Datos sobre el problema
			string Problema; //Problema (Tobera/Tubería)
			string Tobera; //Tipo de tobera (NASA/RAO)

		//Datos sobre la geometría del problema

			//Selección problema tubería

			double PipeLength; //Longitud de la tubería
			double PipeDiameter; //Diametro de la tubería

			//Selección problema tobera

			double P1; //Separación secciones 1 y 2
			double P2; //Separación secciones 2 y 3
			double P3; //Separación secciones 3 y 4
			double P4; //Separación secciones 4 y 5
			double P5; //Separación secciones 5 y 6
			double P6; //Final geometría
			
			double Rinicio; //Radio inicial de la geometría
			double Rfinal; //Radio final de la geometría

			double Xinicio; //Coordenada axial inicial de la geometría
			double Xfinal; //Coordenada axial final de la geometría

			//Datos geometría sección convergente
			double L1; //Longitud de la cámara de combustión (Sección 1)
			double D1; //Diámetro de la cámara de combustión (Sección 1)

			double Rad1; //Radio Sección 2
			double Theta1; //Ángulo Sección 2

			double m; //Pendiente recta Sección 3

			double Rad2; //Radio Sección 4
			double Theta2; //Ángulo Sección 4

			//Datos geometría sección divergente
			double Dt; //Diámetro de la garganta de la tobera
			double Xt; //Coordenada axial de la garganta de la tobera
			double Rad3; //Diámetro Seccion 5 de la tobera
		
			double xW; //Coordenada axial del punto de tangencia sección divergente
			double rW; //Radio del punto de tangencia sección divergente
			double Theta3; //Ángulo del punto de tangencia sección divergente

			double xE; //Coordenada axial de la salida de la tobera
			double rE; //Radio de salida de la tobera
			double ThetaE; //Ángulo de salida de la tobera

		//Datos del mallado
			int OpcionR; //Tipo de discretización en la dirección radial
			int OpcionA; //Tipo de discretización en la dirección axial
			int OpcionSmooth; //Opción para la transición en la garganta
			int OpcionSide; //Opción para el número de nodos totales
			double StretFactorX; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección axial
			double StretFactorR; //Factor de estrechamiento de la discretización tangencial hiperbólica dirección radial

		//Opciones de pared del mallado
			int OpcionYplus;
			double Yplus;
			double DeltaYplus;

		//Características
			double Lcar; //Longitud característica
			double ThroatArea; //Área de la garganta de la tobera

	public:
		//Constructor de la clase
		Mesher(Memory, ReadData, Geometry);

		//Tipo del mallado utilizado
		string TypeMeshing; //Tipo de mallado utilizado (Collocated/Staggered)

		//Densidad del mallado
		int NAConv; //DIrección axial parte convergente
		int NADiv; //Dirección axial parte divergente
		int NA; //Dirección axial total
		int NRad; //Dirección radial
		int NAng; //Dirección angular (solo para postproceso)

		//Matrices de coordenadas de discretización
		double *AxialCoord; //Matriz coordenada axial de la distribución
		double *Radius; //Matriz radios en cada coordenada axial

		//Caso mallado tipo Staggered
		double ***MP; //Coordenadas matriz de presión/temperatura
		double ***MU; //Coordenadas matriz velocidades axiales
		double ***MR; //Coordenadas matriz velocidades radiales

		//Matrices de distancias de volúmenes de control
		double ***DeltasMP; //Deltas X y R de la matriz de Presión/Temperatura
		double ***DeltasMU; //Deltas X y R de la matriz de velocidades axiales (U)
		double ***DeltasMR; //Deltas X y R de la matriz de velocidades radiales (V)

		//Matrices de superficies de volúmenes de control
		double ***SupMP; //Superficies de los volúmenes de la matriz de Presión/Temperatura
		double ***SupMU; //Superficies de los volúmenes de la matriz de velocidades axiales (U)
		double **SupMR; //Superficies de los volúmenes de la matriz de velocidades radiales (V)

		//Matrices de volúmenes de los volúmenes de control
		double **VolMP; //Volúmenes de los volúmenes de la matriz de Presión/Temperatura
		double **VolMU; //Volúmenes de los volúmenes de la matriz de velocidades axiales (U)
		double *VolMR; //Volúmenes de los volúmenes de la matriz de velocidades radiales (V)
		
		//Matrices de los ángulos entre las velocidades y las superficies de los volúmenes de control
		double **AngleMR; //Ángulos entre las superficies de los v.c y las velocidades radiales (V)
		double **AngleMU; //Ángulos entre las velocidades axiales (U)....

		//Distancia mínima de cada nodo a la pared
		double **minDist; //Distancias mínimas de cada nodo a la pared mas cercana

		//Matrices de las superficies efectivas de los volúmenes de control
		double ***eSupMP; //Superficies efectiva de los volúmenes de control de Presión

		//Caso mallado tipo Collocated
		//Matrices de coordenadas de discretización
		double ***MRho; //Coordenadas de las matrices de densidad y velocidades

		//Matrices de distancias de volúmenes de control
		double ***DeltasMRho; //Deltas X y R de las matrices

		//Matrices de superficies de volúmenes de control
		double ***SupMRho; //Superficies de los volúmenes de las matrices

		//Matrices de volúmenes de los volúmenes de control
		double **VolMRho; //Volúmenes de los volúmenes de las matrices

		//Matrices de los ángulos entre las velocidades y las superficies de los volúmenes de control
		double ***Angle; //Ángulos de los volúmenes de control




		//Métodos de la clase del mallador
		void GetTotalAxialNodes(); //Cálculo del número de nodos totales en función de las opciones de discretización
		void GetTotalRadialNodes(); //Calcular el número total de nodos en la dirección radial

		void AllocateMemory(Memory); //Alojamiento de memoria para cada matriz
		void Get_AxialCoordinates(); //Cálculo del la posición axial de los nodos de velocidad axial
		void Get_Radius(Memory, Geometry); //Cálculo del radio de cada coordenada axial discretizada
		void Get_Mesh(); //Creación de todas las mallas
		void Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
		void Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control
		void Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen
		void Get_Angles();
		void Get_EffectiveSurfaces(); //Cálculo de las superficies efectivas de los volúmenes de control
		void Get_CharLength(); //Cálculo de la longitud característica de la geometría
		void ClosestDist(); //Cálculo de la distancia mínima de cad anodo a la pared más cercana

		void PrintWindow(); //Printar por pantalla

		int ReturnNA(){ return NA; }
		int ReturnNRad(){ return NRad; }

		void PrintTxt(); //Pasar los resultados de las mallas a un txt
		void EscalarVTK2D(string, string, double**); //Pasar todos los resultados escalares en 2D a un archivo VTK
		void EscalarVTK2D_Especial(string, string, double**);
		void EscalarVTK2Dv2(string, string);
		void EscalarVTK3D(Memory, string, string); //Pasar todos los resultados escalares en 3D a un archivo VTK

		void ExecuteMesher(Memory, Geometry); //Ejecutar todos los procesos del mallador
};
		
