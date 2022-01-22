#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>

#include "/home_nobck/sergiogus/libSLU/SLU.h"
#include "/home_nobck/sergiogus/SuiteSparse/UMFPACK/Include/umfpack.h"

using namespace std;

#define PI 3.141592653589793



class CalculoSolver{
	private:
	//Valores iniciales de los mapas de propiedades
		
		int NXTOT;
		int NYTOT;
		int NX1, NX2, NX3;
		int NY1, NY2, NY3;	
		int Problema;
		double Xtotal;
		double Ytotal;
		double Wtotal;
		double Xsolido;
		double Ysolido;
		double mu;
		double Rho;
		double Re;
		double Velocidad;
		double ConvergenciaStep;
		double ConvergenciaSteady;
		double gx;
		double gy;
		string Esquema;
		string Esquema2;
		double Producto;
		double Beta;
		double Alpha;
		double To;
		double Rayleigh;
		double Cp;
		double k;
		double Pr;
		
		double Ubot, Utop, Uleft, Uright;
		double Vbot, Vtop, Vleft, Vright;
		double Tbot, Ttop, Tleft, Tright;

		int GUbot, GUtop, GUleft, GUright;
		int GVbot, GVtop, GVleft, GVright;
		int GTbot, GTtop, GTleft, GTright;

		double UbotSol, UtopSol, UleftSol, UrightSol;
		double VbotSol, VtopSol, VleftSol, VrightSol;
		double TbotSol, TtopSol, TleftSol, TrightSol;

		int GUbotSol, GUtopSol, GUleftSol, GUrightSol;
		int GVbotSol, GVtopSol, GVleftSol, GVrightSol;
		int GTbotSol, GTtopSol, GTleftSol, GTrightSol;

		int StepsPantalla;
		int StepsGnuplot;

		//Coeficientes de lift y drag
		double DragCoefficient;
		double LiftCoefficient;

		//Calculo del Strouhal
		double ValorInstante;
		double Periodo;
		double Tiempo;

		//Posicion de comprobación del perfil con RE 100
		int PosicionI;
		int PosicionJ;

		//Matrices de propiedades
		double **P;
		double **Pant;
		double **Ps;

		double **Upres;
		double **Ufut;

		double **Vpres;
		double **Vfut;

		double **FuAnt;
		double **FuPres; 

		double **FvAnt;
		double **FvPres;

		double **Tpres;
		double **Tfut;
		
		double **FTant;
		double **FTpres;

		//Matrices de coeficientes de discretización (calcular presión)
		double **ae;
		double **as;
		double **an;
		double **aw;
		double **ap;
		double **bp;

		//Matrices para aplicar el TDMA
		double **P1;
		double **Q1;
		double **P2;
		double **Q2;
		double **bPtdma;

		//Vectores para aplicar las condiciones de contorno
//----------------------------------------------------------------------------------------------------------------------------------
		//Problemas 1 y 2
		double *utop;
		double *ubot;
		double *uleft;
		double *uright;

		double *vtop;
		double *vbot;
		double *vleft;
		double *vright;

		double *ttop;
		double *tbot;
		double *tleft;
		double *tright;
//----------------------------------------------------------------------------------------------------------------------------------
		//Problema 3
		double *utopSol;
		double *ubotSol;
		double *uleftSol;
		double *urightSol;

		double *vtopSol;
		double *vbotSol;
		double *vleftSol;
		double *vrightSol;

		double *ttopSol;
		double *tbotSol;
		double *tleftSol;
		double *trightSol;

	public:
		//Métodos de la clase

		//Constructor
		CalculoSolver(ReadTXT, SLU, Mallador);

		//Calculo del mapa de presiones
		void InitialFields();
		void BoundaryConditions(Mallador);
		double InterpolacionCDS(double, double, double, double, double);
		double Interpolacion(double, double, double, double, double, double, double, double, double, double, string);
		double Fu(Mallador, int, int);
		double Fv(Mallador, int, int);
		void CalculoBouss(Mallador);
		void CalculoF(Mallador, int);
		void CalculoFinicial(Mallador, int);
		double Predictora(Mallador, int, int, double);
		void CalculoCoeficientes(Mallador);
		void CalculoBP(Mallador, double);
		void CalculoPresionTDMA(double);
		void CalculoVelocidades(Mallador, double);
		void Actualizar();
		double StepTime(Mallador);
		double Parada();
		void EjecucionSolver(SLU, Mallador);
		void TXT(Mallador);
		void ShowWindow();
		void ShowWindow2();
		void CalculoCoefficients(Mallador);
		void PosicionCoordenadas(Mallador);
		void ComprobacionPerfil(Mallador);
		void Inflow(Mallador);
		void MeanSpeed(Mallador);
		void CoeficientesFluctuantes(Mallador, double);

		void CalculoStrouhal(double, double);
		void PrintVTK(Mallador);
		void PrintGnuPlot(int , double, Mallador);
		void VelocidadesVTK(double**, double**, Mallador, string, string);
		void PresionesVTK(double**, Mallador, string, string);

		//Calculo del mapa de temperaturas
		double BoussinesqU(Mallador, int, int, int);
		double BoussinesqV(Mallador, int, int, int);
		void CalculoFTinicial(Mallador, int);
		double FTpresente(Mallador, int, int, int);
		void CalculoFT(Mallador, int);
		void CalculoTemperaturas(double);
		void CalculoNusselt(Mallador);
		
};
#endif

