#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
using namespace std;

//Malla 1 -> Presión/Temperatura
//Malla 2 -> V
//Malla 3 -> U

class Mallador{
	private:
	//Datos del problema
		double Xtotal;
		double Ytotal;
		double Wtotal;
		double Xcent;
		double Ycent;
		double Xsol;
		double Ysol;
		int NX;
		int NY;
		int OpcionX;
		int OpcionY;
		int Eleccion;
		int Problema;
		int NXTOT;
		int NYTOT;
		int Progresion;
		double StretFactor;

	//Matrices de las 3 mallas
		double ***M1;
		double ***M2;
		double ***M3;

	//Matrices de cálculos
		double **DeltasYP;
		double **DeltasXP;
		double **DeltasYU;
		double **DeltasXU;
		double **DeltasYV;
		double **DeltasXV;
		
		double **SurfaceXU;
		double **SurfaceYU;
		double **SurfaceXV;
		double **SurfaceYV;
		double **SurfaceXP;
		double **SurfaceYP;

	//Cálculos geométricos
		double _X1, _X2, _X3;
		double _Y1, _Y2, _Y3;
		int NX1, NX2, NX3;
		int NY1, NY2, NY3;	
		double Delta1, Delta2, Delta3, Delta4;

	public:
	//Métodos de la clase
		Mallador(ReadTXT);
		
		int Get_NodosX();
		int Get_NodosY();
		int Get_NX1();
		int Get_NX2();
		int Get_NX3();
		int Get_NY1();
		int Get_NY2();
		int Get_NY3();

		double Get_Xtotal();
		double Get_Ytotal();

		double M1_1(int, int, int);
		double M1_2(int, int, int);
		double M1_3(int, int, int);

		double DeltaXP(int, int);
		double DeltaYP(int, int);
		double DeltaXU(int, int);
		double DeltaYU(int, int);
		double DeltaXV(int, int);
		double DeltaYV(int, int);

		double SupXU(int, int);
		double SupYU(int, int);
		double SupXV(int, int);
		double SupYV(int, int);
		double SupXP(int, int);
		double SupYP(int, int);
		
		double **Allocate2D(int, int);
		double ***Allocate3D(int, int);

		void Get_Mesh();
		void Get_Walls();
		void Get_Surfaces();
		void PrintTxt();

		void EjecucionMesh();
				
};




