#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
using namespace std;

class Geometry{
	private:
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
		
		//Datos método iterativo
		double Convergencia; //Convergencia para el método iterativo de los coeficientes (Tobera NASA)

		//Datos método límites
		double LimInfXE; //Límite inferior de la coordenada axial de salida
		double LimSupXE; //Límite superior de la coordenada axial de salida
		int AxialCoordIterations; //Número total de coordenadas axiales de salida a estudiar dentro del rango

		double LimInfRE; //Límite inferior de la coordenada radial de salida
		double LimSupRE; //Límite superior de la coordenada radial de salida
		int ExitRadiusIterations; //Número total de coordenadas radiales de salida a estudiar dentro del rango

		double LimInfThetaE; //Límite superior del ángulo de salida
		double LimSupThetaE; //Límite inferior del ángulo de salida
		int AngleIterations; //Número total de ángulos a estudiar dentro del rango

		double *XEpositions;
		double *UpperRE;
		double *LowerRE;

	public:

		//Constructor de la clase
		Geometry(Memory, ReadData);

		//Datos geometría sección divergente
		double *CoeficientesExpresion; //Coeficientes de la expresión matemática	
		double *CoeficientesExpresionSup; //Coeficientes de la expresión matemática (Valores supuestos)

		//Métodos de la clase 
		void GeometryNasa(Memory); //Obtener geometria Tobera NASA
		void GeometryRAO(Memory); //Obtener geometria Tobera RAO

		void NasaLimits(Memory M1); //Método de cálculo de los límites de la Tobera geometría NASA	
};
