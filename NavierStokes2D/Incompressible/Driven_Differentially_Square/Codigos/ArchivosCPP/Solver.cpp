#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <cassert>

#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/ReadTXT.h"
#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/Mesh.h"
#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/Gnuplot.h"
#include "/home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosHEADER/Solver.h"

using namespace std;


#define PI 3.141592653589793


CalculoSolver::CalculoSolver(ReadTXT R1, SLU C1, Mallador M1){

	Problema = R1.Get_NumericalData(0);

	NXTOT = M1.Get_NodosX();
	NYTOT = M1.Get_NodosY();

	if(Problema == 3){
		NX1 = M1.Get_NX1();
		NX2 = M1.Get_NX2();
		NX3 = M1.Get_NX3();

		NY1 = M1.Get_NY1();
		NY2 = M1.Get_NY2();
		NY3 = M1.Get_NY3();
	}

	Utop = R1.Get_BoundaryData(0);
	Uleft = R1.Get_BoundaryData(1);
	Uright = R1.Get_BoundaryData(2);
	Ubot = R1.Get_BoundaryData(3);
	
	Vtop = R1.Get_BoundaryData(4);
	Vleft = R1.Get_BoundaryData(5);
	Vright = R1.Get_BoundaryData(6);
	Vbot = R1.Get_BoundaryData(7);

	Ttop = R1.Get_BoundaryData(8);
	Tleft = R1.Get_BoundaryData(9);
	Tright = R1.Get_BoundaryData(10);
	Tbot = R1.Get_BoundaryData(11);

	Xtotal = R1.Get_GeometryData(0);
	Ytotal = R1.Get_GeometryData(1);
	Wtotal = R1.Get_GeometryData(2);

	Xsolido = R1.Get_GeometryData(5);
	Ysolido = R1.Get_GeometryData(6);

	Rho = R1.Get_FluidData(0);
	
	gx = R1.Get_GravityData(0);
	gy = R1.Get_GravityData(1);

	
	Velocidad = R1.Get_FluidData(1);
	Re = R1.Get_FluidData(2);
		
	Rayleigh = R1.Get_FluidData(3);
	Cp = R1.Get_FluidData(4);
	Pr = R1.Get_FluidData(5);		
		
	if(abs(Tleft + Tright)/2.0 != 0.0){ To = abs(Tleft + Tright)/2; }	
	else{ To = abs(Tbot + Ttop)/2; }
		
	Beta = 1.0/To;	
	Producto = (pow(Rho,2)*abs(gy)*pow(Xtotal,3)*Beta*abs(Tleft-Tright)*Pr)/Rayleigh;

	if(Problema == 1){
		mu = (Rho*Velocidad*Xtotal)/Re;
	}
	else if(Problema == 2){
		mu = sqrt(Producto);
	}
	else if(Problema == 3){
		mu = (Rho*Velocidad*Xsolido)/(Re);
	}
	k = (Cp*mu)/Pr;

	ConvergenciaStep = R1.Get_PrecisionData(1);
	ConvergenciaSteady = R1.Get_PrecisionData(2);
	
	StepsPantalla = R1.Get_NumericalData(6);
	StepsGnuplot = R1.Get_NumericalData(7);

	Esquema = R1.Get_ConvectiveScheme();
	Esquema2 = R1.Get_ConvectiveScheme2();

	GUtop = R1.Get_GradientData(0);
	GUleft = R1.Get_GradientData(1);
	GUright = R1.Get_GradientData(2);
	GUbot = R1.Get_GradientData(3);
	
	GVtop = R1.Get_GradientData(4);
	GVleft = R1.Get_GradientData(5);
	GVright = R1.Get_GradientData(6);
	GVbot = R1.Get_GradientData(7);

	GTtop = R1.Get_GradientData(8);
	GTleft = R1.Get_GradientData(9);
	GTright = R1.Get_GradientData(10);
	GTbot = R1.Get_GradientData(11);
	
	P = M1.Allocate2D(NXTOT, NYTOT);
	Pant = M1.Allocate2D(NXTOT, NYTOT);
	Ps = M1.Allocate2D(NXTOT, NYTOT);

	Vpres = M1.Allocate2D(NXTOT, NYTOT+1);
	Vfut = M1.Allocate2D(NXTOT, NYTOT+1);

	Upres = M1.Allocate2D(NXTOT+1, NYTOT);
	Ufut = M1.Allocate2D(NXTOT+1, NYTOT);

	FuAnt = M1.Allocate2D(NXTOT, NYTOT);
	FuPres = M1.Allocate2D(NXTOT, NYTOT); 

	FvAnt = M1.Allocate2D(NXTOT, NYTOT);
	FvPres = M1.Allocate2D(NXTOT, NYTOT);

	Tpres = M1.Allocate2D(NXTOT, NYTOT);
	Tfut = M1.Allocate2D(NXTOT, NYTOT);

	FTant = M1.Allocate2D(NXTOT, NYTOT);
	FTpres = M1.Allocate2D(NXTOT, NYTOT);

	aw = M1.Allocate2D(NXTOT, NYTOT);
	ae = M1.Allocate2D(NXTOT, NYTOT);
	an = M1.Allocate2D(NXTOT, NYTOT);
	as = M1.Allocate2D(NXTOT, NYTOT);
	ap = M1.Allocate2D(NXTOT, NYTOT);
	bp = M1.Allocate2D(NXTOT, NYTOT);

	P1 = M1.Allocate2D(NXTOT, NYTOT);
	Q1 = M1.Allocate2D(NXTOT, NYTOT);

	P2 = M1.Allocate2D(NXTOT, NYTOT);
	Q2 = M1.Allocate2D(NXTOT, NYTOT);

	bPtdma = M1.Allocate2D(NXTOT, NYTOT);

	utop = new double [NXTOT+1];
	ubot = new double [NXTOT+1];
	uleft = new double [NYTOT];
	uright = new double [NYTOT];

	vtop = new double [NXTOT];
	vbot = new double [NXTOT];
	vleft = new double [NYTOT+1];
	vright = new double [NYTOT+1];

	ttop = new double [NXTOT];
	tbot = new double [NXTOT];
	tleft = new double [NYTOT];
	tright = new double [NYTOT];

//---------------------------------------------------------------------------------------------------------------------------------
	//Condiciones de contorno para el problema 3

	if(Problema == 3){
		//Valores fijados de las propiedades
		UtopSol = R1.Get_BoundaryData(12);
		UleftSol = R1.Get_BoundaryData(13);
		UrightSol = R1.Get_BoundaryData(14);
		UbotSol = R1.Get_BoundaryData(15);
	
		VtopSol = R1.Get_BoundaryData(16);
		VleftSol = R1.Get_BoundaryData(17);
		VrightSol = R1.Get_BoundaryData(18);
		VbotSol = R1.Get_BoundaryData(19);

		TtopSol = R1.Get_BoundaryData(20);
		TleftSol = R1.Get_BoundaryData(21);
		TrightSol = R1.Get_BoundaryData(22);
		TbotSol = R1.Get_BoundaryData(23);

		//Gradientes de las propiedades
		GUtopSol = R1.Get_GradientData(12);
		GUleftSol = R1.Get_GradientData(13);
		GUrightSol = R1.Get_GradientData(14);
		GUbotSol = R1.Get_GradientData(15);
	
		GVtopSol = R1.Get_GradientData(16);
		GVleftSol = R1.Get_GradientData(17);
		GVrightSol = R1.Get_GradientData(18);
		GVbotSol = R1.Get_GradientData(19);

		GTtopSol = R1.Get_GradientData(20);
		GTleftSol = R1.Get_GradientData(21);
		GTrightSol = R1.Get_GradientData(22);
		GTbotSol = R1.Get_GradientData(23);
	
		utopSol = new double [NX2 + 1];
		ubotSol = new double [NX2 + 1];
		uleftSol = new double [NY2];
		urightSol = new double [NY2];

		vtopSol = new double [NX2];
		vbotSol = new double [NX2];
		vleftSol = new double [NY2 + 1];
		vrightSol = new double [NY2 + 1];

		ttopSol = new double [NX2];
		tbotSol = new double [NX2];
		tleftSol = new double [NY2];
		trightSol = new double [NY2];
	}

	//Tiempo de inicio para calcular el número de Strouhal
	Tiempo = 0.0;

}

void CalculoSolver::InitialFields(){
int i, j;
double I, nytot;
nytot = NYTOT;		
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			P[i][j] = 0.0;
			Ps[i][j] = 0.0;
			Tpres[i][j] = Tleft + (i/nytot)*(Tright-Tleft);;
		}
	}

	for(i = 0; i < NXTOT+1; i++){
		for(j = 0; j < NYTOT; j++){
			Upres[i][j] = 0.001;
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT+1; j++){
			Vpres[i][j] = 0.001;
		}
	}
}

void CalculoSolver::BoundaryConditions(Mallador M1){
int i, j;
double J;
double nytot = NYTOT;

	for(i = 0; i < NXTOT; i++){
		if(GVtop == 0){ vtop[i] = Vtop; }
		else{ vtop[i] = Vpres[i][NYTOT-1]; }

		if(GVbot == 0){ vbot[i] = Vbot; }
		else{ vbot[i] = Vpres[i][1]; }

		if(GTtop == 0){ ttop[i] = Ttop; }
		else{ ttop[i] = Tpres[i][NYTOT-1]; }

		if(GTbot == 0){ tbot[i] = Tbot; }
		else{ tbot[i] = Tpres[i][0]; }
	}

	for(j = 0; j < NYTOT; j++){
		J = j;
		if(GUleft == 0){ 
			if(Problema == 3){ uleft[j] = Velocidad*Uleft*(1 - (pow(M1.M1_3(0,j,1) - Ytotal/2,2))/(pow(Ytotal/2,2))); }
			else{ uleft[j] = Uleft;}
		}
		else{ uleft[j] = Upres[0][j]; }

		if(GUright == 0){ uright[j] = Uright; }
		else{ uright[j] = Upres[NXTOT-1][j]; }

		if(GTleft == 0){ tleft[j] = Tleft; }
		else{ tleft[j] = Tpres[0][j]; }

		if(GTright == 0){ tright[j] = Tright;}
		else{ tright[j] = Tpres[NXTOT-1][j]; }
	}

	for(i = 0; i < NXTOT+1; i++){
		if(GUtop == 0){ utop[i] = Utop; }
		else{ utop[i] = Upres[i][NYTOT-1]; }

		if(GUbot == 0){ ubot[i] = Ubot; }
		else{ ubot[i] = Upres[i][0]; }
	}

	for(j = 0; j < NYTOT+1; j++){
		if(GVleft == 0){ vleft[j] = Vleft; }
		else{ vleft[j] = Vpres[0][j]; }

		if(GVright == 0){ vright[i] = Vright; }
		else{ vright[i] = Vpres[NXTOT-1][j]; }
	}

	if(Problema == 3){

		for(i = 0; i < NX2; i++){
			if(GVtopSol == 0){ vtopSol[i] = VtopSol; }
			else{ vtopSol[i] = Vpres[NX1 + i][NY1 + NY2 + 2]; }

			if(GVbotSol == 0){ vbotSol[i] = VbotSol; }
			else{ vbotSol[i] = Vpres[NX1 + i][NY1]; }

			if(GTtopSol == 0){ ttopSol[i] = TtopSol; }
			else{ ttopSol[i] = Tpres[NX1 + i][NY1 + NY2 + 1]; }

			if(GTbotSol == 0){ tbotSol[i] = TbotSol; }
			else{ tbotSol[i] = Tpres[NX1 + i][NY1]; }
		}

		for(j = 0; j < NY2; j++){
			if(GUleftSol == 0){ uleftSol[j] = UleftSol;}
			else{ uleftSol[j] = Upres[NX1][NY1 + j]; }

			if(GUrightSol == 0){ urightSol[j] = UrightSol; }
			else{ urightSol[j] = Upres[NX1 + NX2 + 2][NY1 + j]; }

			if(GTleftSol == 0){ tleftSol[j] = TleftSol; }
			else{ tleftSol[j] = Tpres[NX1-1][NY1+j]; }

			if(GTrightSol == 0){ trightSol[j] = TrightSol;}
			else{ trightSol[j] = Tpres[NX1+NX2][NY1+j]; }
		}

		for(i = 0; i < NX2+1; i++){
			if(GUtopSol == 0){ utopSol[i] = UtopSol; }
			else{ utopSol[i] = Upres[NX1+i][NY1+NY2+1]; }

			if(GUbotSol == 0){ ubotSol[i] = UbotSol; }
			else{ ubotSol[i] = Upres[NX1+i][NY1]; }
		}

		for(j = 0; j < NY2+1; j++){
			if(GVleftSol == 0){ vleftSol[j] = VleftSol; }
			else{ vleftSol[j] = Vpres[NX1-1][NY1+j]; }

			if(GVrightSol == 0){ vrightSol[i] = VrightSol; }
			else{ vrightSol[i] = Vpres[NX1+NX2][NY1+j]; }
		}

	}

}


double CalculoSolver::InterpolacionCDS(double Coordenada1, double Valor1, double Coordenada2, double Valor2, double CoordenadaPunto){

double Valor;
	
	Valor = Valor1 + ((CoordenadaPunto - Coordenada1)/(Coordenada2 - Coordenada1))*(Valor2 - Valor1);

	return Valor;
}

//Funcion para elegir el esquema de interpolacion en las caras
double CalculoSolver::Interpolacion(double Coordenada1, double Phi1, double Coordenada2, double Phi2, double Coordenada3, double Phi3, double Coordenada4, double Phi4, double CoordenadaCara, double Velocidad, string Esquema){
	
	double PHI;

	double XD;
	double PhiD;
	double XC;
	double PhiC;
	double XU;
	double PhiU;

	if (Velocidad < 0.0 || (Phi1 == 0.0 && Coordenada1 == 0.0)){
		XD = Coordenada2;
		PhiD = Phi2;
		XC = Coordenada3;
		PhiC = Phi3;
		XU = Coordenada4;
		PhiU = Phi4;
	}
	else if(Velocidad > 0.0 || (Phi4 == 0.0 && Coordenada4 == 0.0)){
		XD = Coordenada3;
		PhiD = Phi3;
		XC = Coordenada2;
		PhiC = Phi2;
		XU = Coordenada1;
		PhiU = Phi1;
	}

	//Adimensionalizacion
	double PhC;
	double XadC;
	double Xade;

	PhC = (PhiC - PhiU)/(PhiD - PhiU);

	XadC = (XC - XU)/(XD - XU);

	Xade = (CoordenadaCara - XU)/(XD - XU);

	//Evaluacion
	double Phf;
	if (PhiD == PhiU){
		PHI = PhiD;
	}
	else{
	if(Esquema == "CDS"){
		Phf = ((Xade - XadC)/(1 - XadC)) + ((Xade - 1)/(XadC - 1))*PhC;	
	}
	else if(Esquema == "UDS"){
		Phf = PhC;	
	}
	else if(Esquema == "SUDS"){
		Phf = (Xade/XadC)*PhC;
	}
	else if(Esquema == "QUICK"){
		Phf = Xade + (((Xade*(Xade - 1))/(XadC*(XadC-1))))*(PhC - XadC);
	}

	else if(Esquema == "SMART"){
		if(PhC > 0 && PhC < XadC/3){
			Phf = -((Xade*(1 - 3*XadC + 2*Xade))/(XadC*(XadC - 1)))*PhC;
		}
		else if(PhC > XadC/6 && PhC < (XadC/Xade)*(1 + Xade - XadC)){
			Phf = ((Xade*(Xade - XadC))/(1 - XadC)) + ((Xade*(Xade - 1))/(XadC*(XadC - 1)))*PhC;
		}
		else if(PhC > (XadC/Xade)*(1 + Xade - XadC) && PhC < 1){
			Phf = 1;
		}
		else{
			Phf = PhC;
		}
	}

	//Dimensionalizacion
	PHI = PhiU + (PhiD - PhiU)*Phf;

	}

	return PHI;
}

double CalculoSolver::Fu(Mallador M1, int i, int j){

double FU;

double Parte1;
double ue, uw, un, us;
double ue_pred, uw_pred, un_pred, us_pred, vs_pred, vn_pred;

double vs;
double vn;

double Parte2;
	
if(i == 1){
	uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

	uw = Interpolacion(0, 0, M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema2);

	ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

	if(j == 0){
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));
	us = ubot[i];
	un = Interpolacion(M1.M1_2(i,j,1), ubot[i], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);
	vs = vbot[i];
	vn = Interpolacion(M1.M1_3(i-1,j+1,0), vleft[j], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);
	
	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - ubot[i])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}
	else if(j == 1){
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(0, vleft[j], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}
	else if(j == NYTOT-2){
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);

	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], Ytotal, utop[i], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);

	vn = Interpolacion(0, vleft[j], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == NYTOT-1){
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), us_pred, Esquema);
	un = utop[i];
	vs = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = vtop[i];

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(utop[i] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else{
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(0, vleft[j], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}
}

else if(i == NXTOT-1){
	uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

	uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), uright[j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

	ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], 0, 0, M1.M1_1(i,j,0), ue_pred, Esquema2);

	if(j == 0){
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = ubot[i];
	un = Interpolacion(0, ubot[i], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);
	vs = vbot[i];
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], Xtotal, vright[j+1], M1.M1_3(i,j,0), vn_pred, Esquema); 

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - ubot[i])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == 1){
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], Xtotal, vright[j], M1.M1_3(i,j,0), vn_pred, Esquema); 
	
	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == NYTOT-2){
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], Ytotal, utop[i], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], Xtotal, vright[j], M1.M1_3(i,j,0), vn_pred, Esquema);

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == NYTOT-1){
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), us_pred, Esquema);
	un = utop[i];
	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = vtop[i];

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(utop[i] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else{
	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], Xtotal, vright[j], M1.M1_3(i,j,0), vn_pred, Esquema);

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}
}

else{

	if(j == 0){
	uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

	uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

	ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = ubot[i];
	un = Interpolacion(0, ubot[i], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);
	vs = vbot[i];
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);
	
	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - ubot[i])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == 1){
	uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

	uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

	ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == NYTOT-2){
	uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

	uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

	ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
	vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
	un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], Ytotal, utop[i], M1.M1_2(i,j+1,1), un_pred, Esquema);
	
	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);
	
	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else if(j == NYTOT-1){
	uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

	uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

	ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

	us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));

	us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), us_pred, Esquema);
	un = utop[i];
	vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
	vn = vtop[i];

	Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
		(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
		(utop[i] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
		(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	}

	else{
		
			if(Problema == 1 || Problema == 2){
			uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
			ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

			uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

			ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

			us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
			un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
			vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
			vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

			us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
			un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

			vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
			vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

			Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
				(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
				(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
				(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
		}
		else if(Problema == 3){
			if(i >= NX1 && i <= NX1+NX2){
				uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

				uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

				ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);
				if(j == NY1-1){
					us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));

					us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j+1,1), ubotSol[i-NX1], M1.M1_2(i,j,1), us_pred, Esquema);
					un = ubotSol[i-NX1];
					vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
					vn = vbotSol[i-NX1];

					Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
						(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
						(ubotSol[i-NX1] - Upres[i][j])*(2*M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
						(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
	
				}
				else if(j == NY1-2){
					us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
					vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
					vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

					us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
					un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+2,1), ubotSol[i-NX1], M1.M1_2(i,j+1,1), un_pred, Esquema);
	
					vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
					vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);
	
					Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
						(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
						(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
						(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
				}
				else if(j == NY1+NY2){
					un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
					vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

					us = utopSol[i-NX1];
					un = Interpolacion(M1.M1_2(i,j,1), utopSol[i-NX1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);
					vs = vtopSol[i-NX1];
					vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);
	
					Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
						(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
						(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
						(Upres[i][j] - utopSol[i-NX1])*(2*M1.SupXU(i,j)/M1.DeltaYV(i,j));
				}
				else if(j == NY1+NY2+1){
					us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
					vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
					vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

					us = Interpolacion(M1.M1_2(i,j-1,1), utopSol[i-NX1], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
					un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

					vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
					vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

					Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
						(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
						(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
						(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
				}
				else{
					us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
					vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
					vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

					us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
					un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

					vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
					vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

					Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
						(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
						(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
						(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
				}
			}
			else if(i == NX1-1 && j >= NY1 && j < NY1+NY2){
				uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

				uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), urightSol[j-NY1], M1.M1_1(i-1,j,0), uw_pred, Esquema);

				ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], 0, 0, M1.M1_1(i,j,0), ue_pred, Esquema2);

				us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
				un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
				vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
				vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

				us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
				un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

				vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0), vs_pred, Esquema);
				vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i+1,j,0), vleftSol[j-NY1], M1.M1_3(i,j,0), vn_pred, Esquema);

				Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
					(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
					(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
					(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
			}
			else if(i == NX1+NX2+1 && j >= NY1 && j < NY1+NY2){
				uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

				uw = Interpolacion(0, 0, M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema2);

				ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);

				us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
				un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
				vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
				vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

				us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
				un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

				vs = Interpolacion(M1.M1_3(i-1,j,0), vleftSol[j-NY1], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
				vn = Interpolacion(M1.M1_3(i-1,j,0), vleftSol[j-NY1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

				Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
					(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
					(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
					(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
				
			}
			else{
				uw_pred = InterpolacionCDS(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_1(i-1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i,j,0));

				uw = Interpolacion(M1.M1_3(i-2,j,0), Upres[i-2][j], M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_1(i-1,j,0), uw_pred, Esquema);

				ue = Interpolacion(M1.M1_3(i-1,j,0), Upres[i-1][j], M1.M1_3(i,j,0), Upres[i][j], M1.M1_3(i+1,j,0), Upres[i+1][j], M1.M1_3(i+2,j,0), Upres[i+2][j], M1.M1_1(i,j,0), ue_pred, Esquema);			
				
				us_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
				un_pred = InterpolacionCDS(M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_3(i,j,0));
				vs_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j,1));
				vn_pred = InterpolacionCDS(M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j+1,1));

				us = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), us_pred, Esquema);
				un = Interpolacion(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_3(i,j+2,1), Upres[i][j+2], M1.M1_2(i,j+1,1), un_pred, Esquema);

				vs = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vs_pred, Esquema);
				vn = Interpolacion(M1.M1_2(i-2,j+1,0), Vpres[i-2][j+1], M1.M1_2(i-1,j+1,0), Vpres[i-1][j+1], M1.M1_2(i,j+1,0), Vpres[i][j+1], M1.M1_2(i+1,j+1,0), Vpres[i+1][j+1], M1.M1_3(i,j,0), vn_pred, Esquema);

				Parte2 = (Upres[i+1][j] - Upres[i][j])*(M1.SupYU(i,j)/M1.DeltaXP(i,j)) - 
					(Upres[i][j] - Upres[i-1][j])*(M1.SupYU(i,j)/M1.DeltaXP(i-1,j)) + 
					(Upres[i][j+1] - Upres[i][j])*(M1.SupXU(i,j)/M1.DeltaYV(i,j+1)) - 
					(Upres[i][j] - Upres[i][j-1])*(M1.SupXU(i,j)/M1.DeltaYV(i,j));
			}
		}
	}
}

	Parte1 = ue*ue*M1.SupYU(i,j) - uw*uw*M1.SupYU(i,j) + vn*un*M1.SupXU(i,j) - vs*us*M1.SupXU(i,j);

	FU = (-Parte1 + (mu/Rho)*Parte2)/(M1.SupXU(i,j)*M1.DeltaYU(i,j));


	return FU;

}

double CalculoSolver::Fv(Mallador M1, int i, int j){

double FV;

double Parte1;
double ve, vw, vs, vn;
double uw, ue;

double ve_pred, vw_pred, vs_pred, vn_pred, uw_pred, ue_pred;

double Parte2;

if(j == 1){
	vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
	vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

	vs = Interpolacion(0, 0, 0, vbot[i], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema2);
	vn = Interpolacion(0, vbot[i], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

	if(i == 0){
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));

	vw = vleft[j];
	ve = Interpolacion(0, vleft[j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);
	uw = uleft[j];
	ue = Interpolacion(0, ubot[i+1], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - vleft[j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == 1){
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(0, ubot[i], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == NXTOT-2){
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], Xtotal, vright[j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(0, ubot[i], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == NXTOT-1){
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = vright[j];
	uw = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = uright[j];

	Parte2 = (vright[j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else{
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(0, ubot[i], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(0, ubot[i], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);


	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}
}

else if(j == NYTOT-1){
	vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
	vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

	vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], Ytotal, vtop[i], 0, 0, M1.M1_1(i,j,1), vn_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], Ytotal, vtop[i], M1.M1_1(i,j-1,1), vs_pred, Esquema);

	if(i == 0){
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));

	vw = vleft[j];
	ve = Interpolacion(0, vleft[j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);
	uw = uleft[j];
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], Ytotal, utop[i+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - vleft[j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == 1){
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], Ytotal, utop[i], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == NXTOT-2){
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], Xtotal, vright[j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], Ytotal, utop[i], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == NXTOT-1){
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = vright[j];
	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = uright[j];

	Parte2 = (vright[j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else{
	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], Ytotal, utop[i], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], Ytotal, utop[i+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}
}

else{
	
	if(i == 0){
	vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
	vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

	vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);


	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));

	vw = vleft[j];
	ve = Interpolacion(0, vleft[j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);
	uw = uleft[j];
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - vleft[j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == 1){
	vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
	vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

	vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);


	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(0, vleft[j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else if(i == NXTOT-2){
	vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
	vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

	vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);


	ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
	ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], Xtotal, vright[j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

	Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}
	else if(i == NXTOT-1){
	vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
	vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

	vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

	vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);


	vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
	uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

	vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], Xtotal, vright[j], M1.M1_3(i,j,0), vw_pred, Esquema);
	ve = vright[j];
	uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
	ue = uright[j];

	Parte2 = (vright[j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
		(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
		(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
		(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
	}

	else{
		
			if(Problema == 1 || Problema == 2){
			vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
			vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

			vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

			vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);


			ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
			ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
			vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
			uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

			vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
			ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

			uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
			ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

			Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
				(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
				(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
				(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
		}
		else if(Problema == 3){
			if(j >= NY1 && j <= NY1+NY2){
				vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
				vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

				vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

				vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);
				if(i == NX1-1){
					vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

					vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i+1,j,0), vleftSol[j-NY1], M1.M1_3(i,j,0), vw_pred, Esquema);
					ve = vleftSol[j-NY1];
					uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
					ue = uleftSol[j-NY1];

					Parte2 = (vleftSol[j] - Vpres[i][j])*(2*M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
						(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
						(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
						(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
				}
				else if(i == NX1-2){
					ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
					ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
					vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

					vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
					ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+2,j,0), vleftSol[j-NY1], M1.M1_3(i+1,j,0), ve_pred, Esquema);

					uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
					ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

					Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
						(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
						(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
						(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
				}
				else if(i == NX1+NX2){
					ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
					ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));

					vw = vrightSol[j-NY1];
					ve = Interpolacion(M1.M1_3(i,j,0), vrightSol[j-NY1], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);
					uw = urightSol[j-NY1];
					ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

					Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
						(Vpres[i][j] - vrightSol[j-NY1])*(2*M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
						(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
						(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
				}
				else if(i == NX1+NX2+1){
					ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
					ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
					vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

					vw = Interpolacion(M1.M1_3(i-1,j,0), vrightSol[j-NY1], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
					ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

					uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
					ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

					Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
						(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
						(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
						(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
				}
				else{
					ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
					ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
					vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
					uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

					vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
					ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

					uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
					ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

					Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
						(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
						(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
						(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
				}
			}
			else if(j == NY1-1 && i >= NX1 && i < NX1+NX2){
				vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
				vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

				vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_3(i,j+1,1), vbotSol[i-NX1], 0, 0, M1.M1_1(i,j,1), vn_pred, Esquema);

				vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_3(i,j+1,1), vbotSol[i-NX1], M1.M1_1(i,j-1,1), vs_pred, Esquema);

				ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
				vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
				uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

				vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
				ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

				uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_2(i,j+1,1), ubotSol[i-NX1], M1.M1_2(i,j,1), uw_pred, Esquema);
				ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_2(i,j+1,1), ubotSol[i+1-NX1], M1.M1_2(i,j,1), ue_pred, Esquema);

				Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
					(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
					(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
					(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));

			}
			else if(j == NY1+NY2+1 && i >= NX1 && i < NX1+NX2){
				vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
				vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

				vs = Interpolacion(0, 0, M1.M1_2(i,j-1,1), vtopSol[i-NX1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema2);
				vn = Interpolacion(M1.M1_2(i,j-1,1), vtopSol[i-NX1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

				ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
				vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
				uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

				vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
				ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

				uw = Interpolacion(M1.M1_2(i,j-1,1), utopSol[i-NX1], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
				ue = Interpolacion(M1.M1_2(i,j-1,1), utopSol[i+1-NX1], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);


				Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
					(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
					(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
					(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
			}
			else{
				vs_pred = InterpolacionCDS(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_1(i,j-1,1));
				vn_pred = InterpolacionCDS(M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j,1));

				vn = Interpolacion(M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_2(i,j+2,1), Vpres[i][j+2], M1.M1_1(i,j,1), vn_pred, Esquema);

				vs = Interpolacion(M1.M1_2(i,j-2,1), Vpres[i][j-2], M1.M1_2(i,j-1,1), Vpres[i][j-1], M1.M1_2(i,j,1), Vpres[i][j], M1.M1_2(i,j+1,1), Vpres[i][j+1], M1.M1_1(i,j-1,1), vs_pred, Esquema);


				ve_pred = InterpolacionCDS(M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i+1,j,0));
				ue_pred = InterpolacionCDS(M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j,0));
				vw_pred = InterpolacionCDS(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_3(i,j,0));
				uw_pred = InterpolacionCDS(M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Vpres[i][j], M1.M1_2(i,j,1));

				vw = Interpolacion(M1.M1_2(i-2,j,0), Vpres[i-2][j], M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_3(i,j,0), vw_pred, Esquema);
				ve = Interpolacion(M1.M1_2(i-1,j,0), Vpres[i-1][j], M1.M1_2(i,j,0), Vpres[i][j], M1.M1_2(i+1,j,0), Vpres[i+1][j], M1.M1_2(i+2,j,0), Vpres[i+2][j], M1.M1_3(i+1,j,0), ve_pred, Esquema);

				uw = Interpolacion(M1.M1_3(i,j-2,1), Upres[i][j-2], M1.M1_3(i,j-1,1), Upres[i][j-1], M1.M1_3(i,j,1), Upres[i][j], M1.M1_3(i,j+1,1), Upres[i][j+1], M1.M1_2(i,j,1), uw_pred, Esquema);
				ue = Interpolacion(M1.M1_3(i+1,j-2,1), Upres[i+1][j-2], M1.M1_3(i+1,j-1,1), Upres[i+1][j-1], M1.M1_3(i+1,j,1), Upres[i+1][j], M1.M1_3(i+1,j+1,1), Upres[i+1][j+1], M1.M1_2(i,j,1), ue_pred, Esquema);

				Parte2 = (Vpres[i+1][j] - Vpres[i][j])*(M1.SupYV(i,j)/M1.DeltaXU(i+1,j)) - 
					(Vpres[i][j] - Vpres[i-1][j])*(M1.SupYV(i,j)/M1.DeltaXU(i,j)) + 
					(Vpres[i][j+1] - Vpres[i][j])*(M1.SupXV(i,j)/M1.DeltaYP(i,j)) - 
					(Vpres[i][j] - Vpres[i][j-1])*(M1.SupXV(i,j)/M1.DeltaYP(i,j-1));
			}

		}
	}
	
}


Parte1 = ue*ve*M1.SupYV(i,j) - uw*vw*M1.SupYV(i,j) + vn*vn*M1.SupXV(i,j) - vs*vs*M1.SupXV(i,j);

FV = (-Parte1 + (mu/Rho)*Parte2)/(M1.SupXV(i,j)*M1.DeltaYV(i,j));

return FV;

}

double CalculoSolver::BoussinesqU(Mallador M1, int i, int j, int Step){
double Temperatura;
double Resultado;
double GX = gx;
double B = Beta;
double T = To;

	if(i > 1 && i < NXTOT-2){
		Temperatura = Interpolacion(M1.M1_1(i-2,j,0), Tpres[i-2][j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
	}
	else if(i == 1){
		Temperatura = Interpolacion(0, tleft[j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
	}
	else if(i == NXTOT-1){
		Temperatura = Interpolacion(M1.M1_1(i-2,j,0), Tpres[i-2][j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], Xtotal, tright[j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
	}
	
	Resultado = GX*(1 - B*(Temperatura - T));
	
	return Resultado;
}

double CalculoSolver::BoussinesqV(Mallador M1, int i, int j, int Step){
double Temperatura;
double Resultado;
double GY = gy;
double B = Beta;
double T = To;

	if(j > 1 && j < NYTOT-1){
		Temperatura = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
	}
	else if(j == 1){
		Temperatura = Interpolacion(0, tbot[i], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
	}
	else if(j == NYTOT-1){
		Temperatura = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], Ytotal, ttop[i], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
	}
	
	Resultado = GY*(1 - B*(Temperatura - T));
	
	return Resultado;
}

void CalculoSolver::CalculoF(Mallador M1, int Step){
double Temperatura;

		for(int j = 0; j < NYTOT; j++){
			for(int i = 1; i < NXTOT; i++){
				FuAnt[i][j] = FuPres[i][j];	
				FuPres[i][j] = Fu(M1, i, j) + BoussinesqU(M1, i, j, Step);	
			}
		}
		for(int i = 0; i < NXTOT; i++){
			for(int j = 1; j < NYTOT; j++){
				FvAnt[i][j] = FvPres[i][j];
				FvPres[i][j] = Fv(M1, i, j) + BoussinesqV(M1, i, j, Step);	
			}
		}
	
}
void CalculoSolver::CalculoFinicial(Mallador M1, int Step){
double T;
	
	for(int i = 1; i < NXTOT; i++){
		for(int j = 0; j < NYTOT; j++){
			FuPres[i][j] = Fu(M1, i, j) + BoussinesqU(M1, i, j, Step);
				
		}
	}
	for(int i = 0; i < NXTOT; i++){
		for(int j = 1; j < NYTOT; j++){
			FvPres[i][j] = Fv(M1, i, j) + BoussinesqV(M1, i, j, Step);	
		}
	}

}
double CalculoSolver::Predictora(Mallador M1, int i, int j, double DeltaT){

double Predict;
double Ue, Uw, Vs, Vn;

		if (i>0 && i < NXTOT-1 && j > 0 && j < NYTOT-1){
			if(Problema == 1 || Problema == 2){
				Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
				Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
				Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
				Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
			}
			else if(Problema == 3){
				if(i >= NX1 && i < NX1+NX2 && j >= NY1 && j < NY1+NY2){
					Ue = 0.0;
					Uw = 0.0;
					Vs = 0.0;
					Vn = 0.0;
				}
				else{
					if(i == NX1-1 && j >= NY1 && j < NY1+NY2){
						Ue = uleftSol[j-NY1];
						Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
						Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
						Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
					}
					else if(i == NX1+NX2 && j >= NY1 && j < NY1+NY2){
						Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
						Uw = urightSol[j-NY1];
						Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
						Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
					}
					else if(j == NY1+NY2 && i >= NX1 && i < NX1+NX2){
						Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
						Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
						Vs = vtopSol[i-NX1];
						Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
					}
					else if(j == NY1-1 && i >= NX1 && i < NX1+NX2){
						Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
						Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
						Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
						Vn = vbotSol[i-NX1];
					}
					else{
						Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
						Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
						Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
						Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
					}
				}
			}
		}
		if(i == 0){
			if(j == 0){
				Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
				Uw = uleft[j];
				Vs = vbot[i];
				Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
			}
			else if(j == NYTOT-1){
				Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
				Uw = uleft[j];
				Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
				Vn = vtop[i];
			}
			else{
				Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
				Uw = uleft[j];
				Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
				Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
			}
		}
		if(i == NXTOT-1){
			if(j == 0){
				Ue = uright[j];
				Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
				Vs = vbot[i];
				Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
			}
			else if(j == NYTOT-1){
				Ue = uright[j];
				Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
				Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
				Vn = vtop[i];
			}
			else{
				Ue = uright[j];
				Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
				Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
				Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
			}
		}
		if(j == 0){
			if(i > 0 && i < NXTOT-1){
				Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
				Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
				Vs = vbot[i];
				Vn = Vpres[i][j+1] + DeltaT*(1.5*FvPres[i][j+1] - 0.5*FvAnt[i][j+1]);
			}
		}
		if(j == NYTOT-1){
			if(i > 0 && i <NXTOT-1){
				Ue = Upres[i+1][j] + DeltaT*(1.5*FuPres[i+1][j] - 0.5*FuAnt[i+1][j]);
				Uw = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]);
				Vs = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]);
				Vn = vtop[i];
			}
		}
	
	
	Predict = Ue*M1.SupYP(i,j) - Uw*M1.SupYP(i,j) - Vs*M1.SupXP(i,j) + Vn*M1.SupXP(i,j);
	return Predict;
}

void CalculoSolver::CalculoCoeficientes(Mallador M1){

int i, j;

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
				if(i > 0 && i < NXTOT-1 && j > 0 && j < NYTOT-1){
					if(Problema == 1 || Problema == 2){
						aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
						ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
						as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
						an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
					}
					else if(Problema == 3){
						if(i == NX1-1 && j >= NY1 && j < NY1+NY2){
							aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
							ae[i][j] = 0.0;
							as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
							an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
						}
						else if(i == NX1+NX2 && j >= NY1 && j < NY1+NY2){
							aw[i][j] = 0.0;
							ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
							as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
							an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
						}
						else if(j == NY1+NY2 && i >= NX1 && i < NX1+NX2){
							aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
							ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
							as[i][j] = 0.0;
							an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
						}
						else if(j == NY1-1 && i >= NX1 && i < NX1+NX2){
							aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
							ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
							as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
							an[i][j] = 0.0;
						}
						else{
							aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
							ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
							as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
							an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
						}
					}
				}
				if(i == 0){
					aw[i][j] = 0.0;
					if(j == 0){
						ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
						an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
						as[i][j] = 0.0;
					}
					else if(j == NYTOT-1){
						ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
						as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
						an[i][j] = 0.0;
					}
					else{
						ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
						as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
						an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
					}
				}
				else if(i == NXTOT-1){
					ae[i][j] = 0.0;
					if(j == 0){
						aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
						an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
						as[i][j] = 0.0;
					}
					else if(j == NYTOT-1){
						aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
						as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
						an[i][j] = 0.0;
					}
					else{
						aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
						as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
						an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
					}
				}
				else{
					if(j == 0){
						as[i][j] = 0;
						ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
						aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
						an[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j+1);
					}
					else if(j == NYTOT-1){
						as[i][j] = M1.SupXP(i,j)/M1.DeltaYV(i,j);
						ae[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i+1,j);
						aw[i][j] = M1.SupYP(i,j)/M1.DeltaXU(i,j);
						an[i][j] = 0;
					}
				}

			ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j];

		}
	}
}

void CalculoSolver::CalculoBP(Mallador M1, double DeltaT){
int i, j;

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			bp[i][j] = -(Rho/DeltaT)*Predictora(M1, i, j, DeltaT);
		}
	}
}

void CalculoSolver::CalculoVelocidades(Mallador M1, double DeltaT){

int i,j;
double T;

//Velocidades en dirección horizontal (U)
	for(i = 0; i < NXTOT+1; i++){
		for(j = 0; j < NYTOT; j++){
			if(Problema == 1 || Problema == 2){
				if(i > 0 && i < NXTOT){
					Ufut[i][j] = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]) - (DeltaT/(M1.DeltaXU(i,j)*Rho))*(P[i][j] - P[i-1][j]);
				}
				else if(i == 0){
					Ufut[i][j] = uleft[j];
				}
				else if(i == NXTOT){
					Ufut[i][j] = uright[j];
				}
			}
			else if(Problema == 3){
				if (i >= NX1 && i <= NX1 + NX2 && j >= NY1 && j < NY1 + NY2){ Ufut[i][j] = 0.0;}
				else{
					if(i == 0){ Ufut[i][j] = uleft[j];}
					else if(i == NXTOT){ Ufut[i][j] = uright[j];}
					else if(i == NX1 && j >= NY1 && j < NY1+NY2){ Ufut[i][j] = uleftSol[j-NY1];}
					else if(i == NX1+NX2 && j >= NY1 && j < NY1+NY2){ Ufut[i][j] = urightSol[j-NY1];}
					else{
						Ufut[i][j] = Upres[i][j] + DeltaT*(1.5*FuPres[i][j] - 0.5*FuAnt[i][j]) - (DeltaT/(M1.DeltaXU(i,j)*Rho))*(P[i][j] - P[i-1][j]);
					}
				}


			}
		}
	}

//Velocidades en dirección vertical (V)
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT+1; j++){
			if(Problema == 1 || Problema == 2){
				if(j > 0 && j < NYTOT){
					Vfut[i][j] = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]) - (DeltaT/(M1.DeltaYV(i,j)*Rho))*(P[i][j] - P[i][j-1]);
				}
				else if(j == 0){
					Vfut[i][j] = vbot[i];
				}
				else if(j == NYTOT){
					Vfut[i][j] = vtop[i];
				}
			}
			else if(Problema == 3){
				if (i >= NX1 && i < NX1+NX2 && j > NY1 && j < NY1 + NY2){ Vfut[i][j] = 0.0; }
				else{
					if(j == 0){ Vfut[i][j] = vbot[i];}
					else if(j == NYTOT){ Vfut[i][j] = vtop[i];}
					else if(j == NY1 && i >= NX1 && i < NX1+NX2){ Vfut[i][j] = vbotSol[i-NX1];}
					else if(j == NY1+NY2 && i >= NX1 && i < NX1+NX2){ Vfut[i][j] = vtopSol[i-NX1];}
					else{
						Vfut[i][j] = Vpres[i][j] + DeltaT*(1.5*FvPres[i][j] - 0.5*FvAnt[i][j]) - (DeltaT/(M1.DeltaYV(i,j)*Rho))*(P[i][j] - P[i][j-1]);
					}
				}
			}
		}
	}
}

void CalculoSolver::CalculoPresionTDMA(double Fr){

int i, j;
double MaxDif = 2*ConvergenciaStep;


	while(MaxDif >= ConvergenciaStep){
		for(j = 0; j < NYTOT; j++){
			bPtdma[0][j] = bp[0][j] + ae[0][j]*P[1][j];
			bPtdma[NXTOT-1][j] = bp[NXTOT-1][j] + aw[NXTOT-1][j]*P[NXTOT-2][j];
			for(i = 1; i < NXTOT-1; i++){
				bPtdma[i][j] = bp[i][j] + aw[i][j]*P[i-1][j] + ae[i][j]*P[i+1][j];
			}
		}
	
		for(i = 0; i < NXTOT; i++){
			P1[i][0] = an[i][0]/ap[i][0];
			Q1[i][0] = bPtdma[i][0]/ap[i][0];

			P1[i][NYTOT-1] = 0;
			Q1[i][NYTOT-1] = bPtdma[i][NYTOT-1]/ap[i][NYTOT-1];
			for(j = 1; j < NYTOT-1; j++){
				P1[i][j] = (an[i][j])/(ap[i][j] - as[i][j]*P1[i][j-1]);
				Q1[i][j] = (bPtdma[i][j] + as[i][j]*Q1[i][j-1])/(ap[i][j] - as[i][j]*P1[i][j-1]);
			}
		}

		for(i = 0; i < NXTOT; i++){
			P[i][NYTOT-1] = Q1[i][NYTOT-1];
			for (j = NYTOT-2; j >= 0; j--){
				P[i][j] = Ps[i][j] + Fr*(P1[i][j]*P[i][j+1] + Q1[i][j] - Ps[i][j]);
			}
		}


		for(i = 0; i < NXTOT; i++){
			bPtdma[i][0] = bp[i][0] + an[i][0]*P[i][1];
			bPtdma[i][NYTOT-1] = bp[i][NYTOT-1] + as[i][NYTOT-1]*P[i][NYTOT-2];
			for(j = 1; j < NYTOT-1; j++){
				bPtdma[i][j] = bp[i][j] + as[i][j]*P[i][j-1] + an[i][j]*P[i][j+1];
			}
		}
		
		for(j = 0; j < NYTOT; j++){
			P2[0][j] = ae[0][j]/ap[0][j];
			Q2[0][j] = bPtdma[0][j]/ap[0][j];

			P2[NXTOT-1][j] = 0;
			Q2[NXTOT-1][j] = bPtdma[NXTOT-1][j]/ap[NXTOT-1][j];
			for(i = 1; i < NXTOT-1; i++){
				P2[i][j] = (ae[i][j])/(ap[i][j] - aw[i][j]*P2[i-1][j]);
				Q2[i][j] = (bPtdma[i][j] + aw[i][j]*Q2[i-1][j])/(ap[i][j] - aw[i][j]*P2[i-1][j]);
			}
		}

		for(j = 0; j < NYTOT; j++){
			P[NXTOT-1][j] = Q2[NXTOT-1][j];
			for(i = NXTOT-2; i>=0; i--){
				P[i][j] = Ps[i][j] + Fr*(P2[i][j]*P[i+1][j] + Q2[i][j] - Ps[i][j]);
			}
		}

		MaxDif = 0.0;
		for(i = 0; i < NXTOT; i++){
			for(j = 0; j < NYTOT; j++){
				if(abs(P[i][j] - Ps[i][j]) >= MaxDif){
					MaxDif = abs(P[i][j] - Ps[i][j]);
					
				}
				Ps[i][j] = P[i][j];
			}
		}
	
	}



}


double CalculoSolver::StepTime(Mallador M1){
int i, j;

double DeltaTD = 1000;
double DeltaTC = 1000;
double DeltaTT = 1000;

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			if(abs((0.20*M1.DeltaXP(i,j))/(0.5*(Upres[i][j] + Upres[i+1][j]))) < DeltaTD){
				DeltaTD = abs((0.20*M1.DeltaXP(i,j))/(0.5*(Upres[i][j] + Upres[i+1][j])));
			}
			if(abs((0.20*M1.DeltaYP(i,j))/(0.5*(Vpres[i][j] + Vpres[i][j+1]))) < DeltaTD){
				DeltaTD = abs((0.20*M1.DeltaYP(i,j))/(0.5*(Vpres[i][j] + Vpres[i][j+1])));
			}	
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			if(abs((0.10*Rho*pow(M1.DeltaXP(i,j),2))/mu) < DeltaTC){
				DeltaTC = abs((0.10*Rho*pow(M1.DeltaXP(i,j),2))/mu);
			}
			if(abs((0.10*Rho*pow(M1.DeltaYP(i,j),2))/mu) < DeltaTC){
				DeltaTC = abs((0.10*Rho*pow(M1.DeltaYP(i,j),2))/mu);
			}	
		}
	}
	
	if(Problema == 2){
		for(i = 0; i < NXTOT; i++){
			for(j = 0; j < NYTOT; j++){
				if(abs((0.05*Rho*Cp*pow(M1.DeltaXP(i,j),2))/k) < DeltaTT){
					DeltaTT = abs((0.05*Rho*Cp*pow(M1.DeltaXP(i,j),2))/k);
				}
				if(abs((0.05*Rho*Cp*pow(M1.DeltaYP(i,j),2))/k) < DeltaTT){
					DeltaTT = abs((0.05*Rho*Cp*pow(M1.DeltaYP(i,j),2))/k);
				}	
			}
		}
	}


	if(DeltaTD < DeltaTC && DeltaTD < DeltaTT){
		return DeltaTD;
	}
	else if(DeltaTC < DeltaTD && DeltaTC < DeltaTT){
		return DeltaTC;
	}
	else{
		return DeltaTT;
	}
}
void CalculoSolver::Actualizar(){

int i, j;

	for(i = 0; i < NXTOT+1; i++){
		for(j = 0; j < NYTOT; j++){
			Upres[i][j] = Ufut[i][j];
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT+1; j++){
			Vpres[i][j] = Vfut[i][j];
		}
	}
	
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			Pant[i][j] = P[i][j];
			Tpres[i][j] = Tfut[i][j];
		}
	}	
}

double CalculoSolver::Parada(){

double MaxDifParar = 0.0;
double Diferencia;
int i, j;

	if(Problema == 2){
		for(i = 0; i < NXTOT; i++){
			for(j = 0; j < NYTOT; j++){
				if(abs(Tfut[i][j] - Tpres[i][j]) >= MaxDifParar){
					MaxDifParar = abs(Tfut[i][j] - Tpres[i][j]);
					Diferencia = Tfut[i][j] - Tpres[i][j];
				}
			}
		}
	}

	for(i = 0; i < NXTOT+1; i++){
		for(j = 0; j < NYTOT; j++){
			if(abs(Ufut[i][j] - Upres[i][j]) >= MaxDifParar){
				MaxDifParar = abs(Ufut[i][j] - Upres[i][j]);
				Diferencia = Ufut[i][j] - Upres[i][j];
			}
		}
	}

	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT+1; j++){
			if(abs(Vfut[i][j] - Vpres[i][j]) >= MaxDifParar){
				MaxDifParar = abs(Vfut[i][j] - Vpres[i][j]);
				Diferencia = Vfut[i][j] - Vpres[i][j];
			}
		}
	}

	return Diferencia;
}



void CalculoSolver::TXT(Mallador M1){
int i, j;

	FILE *fp1;
	fp1 = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Presion/Presion.txt","w");
	
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			fprintf(fp1,"%f %f %f \n", M1.M1_1(i,j,0), M1.M1_1(i,j,1), P[i][j]);
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);

	FILE *fpVel;
	fpVel = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Velocidades/Velocidades.txt","w");
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			if(Problema == 1 || Problema == 2){
				fprintf(fpVel,"%f %f %f %f \n", M1.M1_1(i,j,0), M1.M1_1(i,j,1), 0.5*(Ufut[i][j] + Ufut[i+1][j]), 0.5*(Vfut[i][j] + Vfut[i][j+1]));
			}
			else if(Problema == 3){
				if(i >= NX1 && i < NX1+NX2 && j > NY1 && j < NY1+NY2){ }
				else{ fprintf(fpVel,"%f %f %f %f \n", M1.M1_1(i,j,0), M1.M1_1(i,j,1), 0.5*(Ufut[i][j] + Ufut[i+1][j]), 0.5*(Vfut[i][j] + Vfut[i][j+1])); }	
			}
		}
		fprintf(fpVel, "\n");
	}
	fclose(fpVel);

	FILE *fpT;
	fpT = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Temperatura/Temperaturas.txt","w");
	
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j <NYTOT; j++){
			if(Problema == 1 || Problema == 2){
				fprintf(fpT,"%f %f %f \n", M1.M1_1(i,j,0), M1.M1_1(i,j,1), Tfut[i][j]);
			}
			else if(Problema == 3){
				if(i >= NX1+1 && i < NX1+NX2 && j >= NY1 && j < NY1+NY2){ }
				else{ fprintf(fpT,"%f %f %f \n", M1.M1_1(i,j,0), M1.M1_1(i,j,1), Tfut[i][j]); }	
			}
		}
		fprintf(fpT, "\n");
	}
	fclose(fpT);

}

double CalculoSolver::FTpresente(Mallador M1, int i, int j, int Step){

double Parte1;
double Parte2;

double Te, Tw, Ts, Tn;
double FT;

	if(i == 0){
		Te = Interpolacion(0, tleft[j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_1(i+2,j,0), Tpres[i+2][j], M1.M1_3(i+1,j,0), Upres[i+1][j], Esquema);
		Tw = tleft[j];
		if(j == 0){
			Ts = tbot[i];
			Tn = Interpolacion(0, tbot[i], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - tleft[j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) -
				(Tpres[i][j] - tbot[i])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j));	

		}
		else if(j == 1){
			Ts = Interpolacion(0, tbot[i], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - tleft[j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));	

		}
		else if(j == NYTOT-1){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], Ytotal, ttop[i], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = ttop[i];

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - tleft[j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j)) +
				(ttop[i] - Tpres[i][j])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));
		}	
		else if(j == NYTOT-2){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], Ytotal, ttop[i], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);
			
			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - tleft[j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));
		}
		else{
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);
			
			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - tleft[j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));
		}
	}		
	else if(i == NXTOT-1){
		Te = tright[j];
		Tw = Interpolacion(M1.M1_1(i-2,j,0), Tpres[i-2][j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], Xtotal, tright[j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
		if(j == 0){
			Ts = tbot[i];
			Tn = Interpolacion(0, tbot[i], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);
				
			Parte2 = (tright[j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) -
				(Tpres[i][j] - tbot[i])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j));	
				
		}
		else if(j == 1){
			Ts = Interpolacion(0, tbot[i], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);

			Parte2 = (tright[j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));	

		}
		else if(j == NYTOT-1){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], Ytotal, ttop[i], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = ttop[i];

			Parte2 = (tright[j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) -	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j)) + 
				(ttop[i] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j));
		}	
		else if(j == NYTOT-2){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], Ytotal, ttop[i], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);
			
			Parte2 = (tright[j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));

		}
		else{
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);
		
			Parte2 = (tright[j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) -
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));
		}
	}
	else{
		if(i == 1){
			Te = Interpolacion(M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_1(i+2,j,0), Tpres[i+2][j], M1.M1_3(i+1,j,0), Upres[i+1][j], Esquema);
			Tw = Interpolacion(0, tleft[j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
		}
		else if(i == NXTOT-2){
			Te = Interpolacion(M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], Xtotal, tright[j], M1.M1_3(i+1,j,0), Upres[i+1][j], Esquema);
			Tw = Interpolacion(M1.M1_1(i-2,j,0), Tpres[i-2][j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
		}
		else if(i > 1 && i < NXTOT-2){
			Te = Interpolacion(M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_1(i+2,j,0), Tpres[i+2][j], M1.M1_3(i+1,j,0), Upres[i+1][j], Esquema);
			Tw = Interpolacion(M1.M1_1(i-2,j,0), Tpres[i-2][j], M1.M1_1(i-1,j,0), Tpres[i-1][j], M1.M1_1(i,j,0), Tpres[i][j], M1.M1_1(i+1,j,0), Tpres[i+1][j], M1.M1_3(i,j,0), Upres[i][j], Esquema);
		}
		if(j == 0){
			Ts = tbot[i];
			Tn = Interpolacion(0, tbot[i], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) -
				(Tpres[i][j] - tbot[i])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));	
				
		}
		else if(j == 1){
			Ts = Interpolacion(0, tbot[i], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));	

		}
		else if(j == NYTOT-1){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], Ytotal, ttop[i], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = ttop[i];

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) + 
				(ttop[i] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) -	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));
		}	
		else if(j == NYTOT-2){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], Ytotal, ttop[i], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);
			
			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));

		}
		else if(j > 1 && j < NYTOT-2){
			Ts = Interpolacion(M1.M1_1(i,j-2,1), Tpres[i][j-2], M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_2(i,j,1), Vpres[i][j], Esquema);
			Tn = Interpolacion(M1.M1_1(i,j-1,1), Tpres[i][j-1], M1.M1_1(i,j,1), Tpres[i][j], M1.M1_1(i,j+1,1), Tpres[i][j+1], M1.M1_1(i,j+2,1), Tpres[i][j+2], M1.M1_2(i,j+1,1), Vpres[i][j+1], Esquema);

			Parte2 = (Tpres[i+1][j] - Tpres[i][j])/(M1.DeltaXU(i+1,j)*M1.DeltaXP(i,j)) -
				(Tpres[i][j] - Tpres[i-1][j])/(M1.DeltaXU(i,j)*M1.DeltaXP(i,j)) +
				(Tpres[i][j+1] - Tpres[i][j])/(M1.DeltaYV(i,j+1)*M1.DeltaYP(i,j)) - 	
				(Tpres[i][j] - Tpres[i][j-1])/(M1.DeltaYP(i,j)*M1.DeltaYV(i,j));
		}
	}	

	Parte1 = Upres[i+1][j]*Te*M1.SupYP(i,j) + Vpres[i][j+1]*Tn*M1.SupXP(i,j) - Upres[i][j]*Tw*M1.SupYP(i,j) - Vpres[i][j]*Ts*M1.SupXP(i,j);
	
	FT = (k/(Rho*Cp))*Parte2 - Parte1/(M1.SupXP(i,j)*M1.DeltaYP(i,j));
	
	return FT;
}

void CalculoSolver::CalculoFT(Mallador M1, int Step){
int i, j;
		for(i = 0; i < NXTOT; i++){
			for(j = 0; j < NYTOT; j++){
				FTant[i][j] = FTpres[i][j];		
				FTpres[i][j] = FTpresente(M1, i, j, Step);
				
			}
		}
}

void CalculoSolver::CalculoFTinicial(Mallador M1, int Step){
int i, j;
	for(i = 0; i < NXTOT; i++){
		for(j = 0; j < NYTOT; j++){
			FTpres[i][j] = FTpresente(M1, i, j, Step);	
		}
	}
}

void CalculoSolver::CalculoTemperaturas(double DeltaT){

int i, j;
		for(i = 0; i < NXTOT; i++){
			for(j = 0; j < NYTOT; j++){
				Tfut[i][j] = Tpres[i][j] + DeltaT*(1.5*FTpres[i][j] - 0.5*FTant[i][j]);
			}
		}	
}

void CalculoSolver::CalculoNusselt(Mallador M1){
int i, j;
double NusseltMax, NusseltMin, NusseltMedio;
double Umax = 0.0;
int ImaxU = 0;
double Vmax = 0.0;
int ImaxV = 0;
double Sumatorio = 0;
int ImaxNu = 0;
int IminNu = 0;
double *NusseltLocal;
NusseltLocal = new double [NYTOT];

	if(Tleft > Tright && Tleft != 0.0){
		for(j = 0; j < NYTOT; j++){
			NusseltLocal[j] = ((Xtotal/(Tleft - Tright)))*((Tleft - Tpres[0][j])/M1.DeltaXU(0,j));
		}	
	}
	

	NusseltMax = 0.0;
	NusseltMin = 1e3;
	
		for(i = 0; i < NXTOT; i++){
			if(NusseltLocal[i] > NusseltMax){
				NusseltMax = NusseltLocal[i];
				ImaxNu = i;
			}
			if(NusseltLocal[i] < NusseltMin){
				NusseltMin = NusseltLocal[i];
				IminNu = i;
			}
			Sumatorio = Sumatorio + NusseltLocal[i];
		}
	double nytot = NYTOT;
		NusseltMedio = Sumatorio/nytot;

		for(j = 0; j < NYTOT; j++){
			if(Upres[NXTOT/2][j] > Umax){ 
				Umax = Upres[NXTOT/2][j];
				ImaxU = j;
			}
			
		}

		for(i = 0; i < NXTOT; i++){
			if(Vpres[i][NYTOT/2] > Vmax){ 
				Vmax = Vpres[i][NYTOT/2];
				ImaxV = i;
			}
			
		}

				FILE *fp1;
				fp1 = fopen("/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosFinales/DifferentialHeated/RA_1E6/ResultadosNusselt.txt","w");
					fprintf(fp1,"Número de Rayleigh: %f \n", Rayleigh);
					fprintf(fp1, "\n");
					fprintf(fp1,"Velocidad U máxima (m/s): %f \n", (Xtotal/k)*Umax);
					fprintf(fp1,"Posicion Y de la velocidad U máxima (m): %f \n", M1.M1_3(NXTOT/2,ImaxU,1));
					fprintf(fp1, "\n");
					fprintf(fp1,"Velocidad V máxima (m/s): %f \n", (Xtotal/k)*Vmax);
					fprintf(fp1,"Posicion X de la velocidad V máxima (m): %f \n", M1.M1_2(ImaxV,NYTOT/2,0));
					fprintf(fp1, "\n");
					fprintf(fp1,"Número Nusselt medio: %f \n", NusseltMedio);
					fprintf(fp1, "\n");
					fprintf(fp1,"Número de Nusselt máximo: %f \n", NusseltMax);
					fprintf(fp1,"Posición Y del número de Nusselt máximo (m): %f \n", M1.M1_1(1,ImaxNu,1));
					fprintf(fp1, "\n");
					fprintf(fp1,"Número de Nusselt mínimo: %f \n", NusseltMin);
					fprintf(fp1,"Posición Y del número de Nusselt mínimo (m): %f \n", M1.M1_1(1,IminNu,1));
					
				fclose(fp1);

	
delete NusseltLocal;

}

void CalculoSolver::CalculoCoefficients(Mallador M1){
int i,j;
int REYNOLDS = Re;
double Strouhal = Xsolido/(Periodo*Uleft*Velocidad);

	char Results[300];
	sprintf(Results, "/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosFinales/SquareCylinder/RE_%d/Resultados.txt", REYNOLDS);

	FILE *fp1;
		fp1 = fopen(Results,"w");
			fprintf(fp1,"Número de Reynolds de la simulación: %f \n", Re);
			fprintf(fp1,"Coeficiente de drag: %f \n", DragCoefficient);
			fprintf(fp1,"Coeficiente de lift: %f \n", LiftCoefficient);
			fprintf(fp1,"Periodo de los vórtices: %f s \n", Periodo);
			fprintf(fp1,"Frecuencia de los vórtices: %f hz \n", 1/Periodo);
			fprintf(fp1,"Número de Strouhal de la simulación: %f \n", Strouhal);
			fprintf(fp1,"\n");
			fprintf(fp1,"Datos de la simulación: \n");
			fprintf(fp1,"Número de nodos en X: %d \n", NXTOT);
			fprintf(fp1,"Número de nodos en Y: %d \n", NYTOT);
			fprintf(fp1,"\n");			
	fclose(fp1);
}

void CalculoSolver::CoeficientesFluctuantes(Mallador M1, double TotalTime){

int i, j;
double FuerzaFrontal1 = 0.0; 
double FuerzaTrasera1 = 0.0;
double CortanteUp = 0.0;
double CortanteDown = 0.0;
double CortanteLeft = 0.0;
double CortanteRight = 0.0;
double FuerzaNetaX;

	for(j = 0; j < NY2; j++){
		FuerzaFrontal1 = FuerzaFrontal1 + P[NX1-1][NY1+j]*M1.SupYP(NX1-1,NY1+j);
		FuerzaTrasera1 = FuerzaTrasera1 + P[NX1+NX2][NY1+j]*M1.SupYP(NX1+NX2,NY1+j);
	}
	
	for(i = 0; i < NX2; i++){
		CortanteUp = CortanteUp + 2*(M1.SupXU(NX1+i,NY1+NY2)*mu*Ufut[NX1+i][NY1+NY2])/M1.DeltaYV(NX1+i,NY1+NY2);
		CortanteDown = CortanteDown + 2*(M1.SupXU(NX1+i,NY1-1)*mu*Ufut[NX1+i][NY1-1])/M1.DeltaYV(NX1+i,NY1);
		
	}

	FuerzaNetaX = FuerzaFrontal1 - FuerzaTrasera1 + CortanteUp + CortanteDown;

	DragCoefficient = FuerzaNetaX/(0.5*Rho*pow(Velocidad,2)*Wtotal*Ysolido);

double FuerzaDown = 0.0; 
double FuerzaUp = 0.0;
double ShearLeft = 0.0;
double ShearRight = 0.0;
double FuerzaNetaY;

	for(i = 0; i < NX2; i++){
		FuerzaDown = FuerzaDown + P[NX1+i][NY1-1]*M1.SupXP(NX1+i,NY1-1);
		FuerzaUp = FuerzaUp + P[NX1+i][NY1+NY2]*M1.SupXP(NX1+i,NY1+NY2);
	}

	for(j = 0; j < NY2; j++){
		ShearRight = ShearRight + (2*mu*Vfut[NX1+NX2][NY1+j]*M1.SupYP(NX1+NX2,NY1+j))/M1.DeltaXU(NX1+NX2,NY1+j);
		ShearLeft = ShearLeft + (2*mu*Vfut[NX1-1][NY1+j]*M1.SupYP(NX1-1,NY1+j))/M1.DeltaXU(NX1-1,NY1+j);
	}

	FuerzaNetaY = FuerzaDown - FuerzaUp + ShearLeft + ShearRight;

	LiftCoefficient = FuerzaNetaY/(0.5*Rho*pow(Velocidad,2)*Wtotal*Xsolido);
	
	int REYNOLDS = Re;
	char Results[300];
	sprintf(Results, "/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosFinales/SquareCylinder/RE_%d/DragCoefficient.txt", REYNOLDS);

	FILE *fp1;
		fp1 = fopen(Results,"a");
			fprintf(fp1,"%f \t %f \n", TotalTime, DragCoefficient);
					
	fclose(fp1);
	
	char Results2[300];
	sprintf(Results2, "/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosFinales/SquareCylinder/RE_%d/LiftCoefficient.txt", REYNOLDS);

	FILE *fp2;
		fp2 = fopen(Results2,"a");
			fprintf(fp2,"%f \t %f \n", TotalTime, LiftCoefficient);
					
	fclose(fp2);
	
}
void CalculoSolver::CalculoStrouhal(double DeltaT, double TotalTime){
double Tolerancia = 1e-6;


	if(TotalTime > 30.0){
		if(abs(Ufut[NX1 + NX2 + NX3/3][NY1 + NY2/2] - ValorInstante) <= Tolerancia){
			Periodo = Tiempo;
			Tiempo = 0.0;
			ValorInstante = Ufut[NX1 + NX2 + NX3/3][NY1 + NY2/2];
			
		}
		else{
			Tiempo = Tiempo + DeltaT;
			if(Tiempo > 25.0){
				ValorInstante = Ufut[NX1 + NX2 + NX3/3][NY1 + NY2/2];
				Tiempo = 0.0;
			}
		}
	}

}


void CalculoSolver::PosicionCoordenadas(Mallador M1){
int i, j;
double DiferenciaX = Xtotal;
double DiferenciaY = Ytotal;

		for(i = NX1 + NX2-1; i < NXTOT; i++){
			if(abs(M1.M1_2(i,1,0) - 23) <= DiferenciaX){
				DiferenciaX = abs(M1.M1_2(i,1,0) - 23.0);
				PosicionI = i;
			}
		}

		for(j = NY1; j < NY1+NY2; j++){
			if(abs(M1.M1_2(1,j,1) - 4) <= DiferenciaY){
				DiferenciaY = abs(M1.M1_2(1,j,1) - 4.0);
				PosicionJ = j;
			}
		}
}
void CalculoSolver::ComprobacionPerfil(Mallador M1){
int RE = Re;

	if(Vpres[PosicionI][PosicionJ] < 0.0 && Vfut[PosicionI][PosicionJ] > 0.0 && Re == 100){

		char Velocity[100];
		sprintf(Velocity, "PerfilVelocidades_RE_%d", RE);
	
		char Pressure[100];
		sprintf(Pressure, "PerfilPresiones_RE_%d", RE);

		VelocidadesVTK(Upres, Vpres, M1, "Velocidades", Velocity); 
		PresionesVTK(P, M1, "Presiones", Pressure);	
	
	}
}

void CalculoSolver::PrintVTK(Mallador M1){
int RE = Re;

	char Velocity[100];
	sprintf(Velocity, "ResultadosVelocidades_RE_%d", RE);
	
	char Pressure[100];
	sprintf(Pressure, "ResultadosPresiones_RE_%d", RE);

	VelocidadesVTK(Upres, Vpres, M1, "Velocidades", Velocity); 
	PresionesVTK(P, M1, "Presiones", Pressure);

}

void CalculoSolver::PrintGnuPlot(int Step, double TotalTime, Mallador M1){
int RE = Re;

//Clase de Gnuplot
Gnuplot G1;

	if(Step%StepsGnuplot == 0){
		TXT(M1);
		char Directorio[200];
		char RangoX[200];
		char RangoY[200];
		char Title[200];
			
		G1("set terminal pngcairo size 1920,1080");
		if(Problema == 1){
			sprintf(Directorio, "set output \'/home/sergiogus/Desktop/CTTC/NavierStokes/Imagenes/DrivenCavity/RE_%d/VelocidadesStep_%d.png", RE, Step);	
		}
		else if(Problema == 2){
			sprintf(Directorio, "set output \'/home/sergiogus/Desktop/CTTC/NavierStokes/Imagenes/DifferentiallyHeated/RE_%d/VelocidadesStep_%d.png", RE, Step);
		}
		else if(Problema == 3){
			sprintf(Directorio, "set output \'/home/sergiogus/Desktop/CTTC/NavierStokes/Imagenes/SquareCylinder/RE_%d/VelocidadesStep_%d.png", RE, Step);
		}
 			
 		G1(Directorio);

		sprintf(RangoX, "set xrange [0:%f]", Xtotal);
		G1(RangoX);

		sprintf(RangoY, "set yrange [0:%f]", Ytotal);
		G1(RangoY);

		G1("set xlabel 'Coordenada X'");
		G1("set ylabel 'Coordenada Y'");

		if(Problema == 1 || Problema == 3){
		sprintf(Title, "set title 'Mapa de velocidades RE %d, tiempo %f'", RE, TotalTime);
		}
		else if(Problema == 2){
		sprintf(Title, "set title 'Mapa de velocidades Ra %f, tiempo %f'",Rayleigh, TotalTime);
		}
		G1(Title);

 		G1("plot \'/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Velocidades/Velocidades.txt\' using 1:2:3:4 with vectors filled head lw 0.2");
			
		if(Problema == 2){
			char Contour[200];
			G1("set terminal pngcairo size 1920,1080");
 			sprintf(Contour, "set output \'/home/sergiogus/Desktop/CTTC/NavierStokes/Imagenes/DifferentiallyHeated/MapaTemperaturasStep_%d.png", Step);	
 			G1(Contour);
 			G1(RangoX);
 			G1(RangoY);

			G1("set title 'Mapa de temperaturas Differential Heated problem' ");
			G1("set sample 11; set isosamples 11");
			G1("set pm3d map");
			G1("set palette");
			G1("set colorbox");
			G1("set lmargin 0");
			G1("set pm3d flush begin");
 			G1("splot '/home/sergiogus/Desktop/CTTC/NavierStokes/Resultados/Temperatura/Temperaturas.txt' u 1:2:3");
		}
	}

}

void CalculoSolver::VelocidadesVTK(double **phi1, double **phi2, Mallador M1, string name, string Doc){
        ofstream file;
        stringstream name_t;
        string name_f;

	if(Problema == 1) name_t<<"/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosParaview/DrivenCavity/"<<Doc<<".vtk";

        else if(Problema == 2)  name_t<<"/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosParaview/DifferentiallyHeated/"<<Doc<<".vtk";

        else name_t<<"/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosParaview/SquareCylinder/"<<Doc<<".vtk";

	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<name<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET RECTILINEAR_GRID"<<endl;
        file<<"DIMENSIONS"<<"   "<<NXTOT<<"   "<<NYTOT<<"   "<<1<<endl;
        file<<endl;
        file<<"X_COORDINATES"<<"   "<<NXTOT<<"   "<<"double"<<endl;
        for(int i = 0; i < NXTOT; i++){
                file<<M1.M1_1(i,0,0)<<"   ";
        }
        file<<endl;
        file<<"Y_COORDINATES"<<"   "<<NYTOT<<"   "<<"double"<<endl;
        for(int j = 0;j < NYTOT;j++){
                file<<M1.M1_1(0,j,1)<<"   ";
        }
        file<<endl;
        file<<"Z_COORDINATES"<<"   "<<1<<"   "<<"double"<<endl;
        file<<0.0<<endl;
        file<<endl;
        file<<"POINT_DATA"<<"   "<<NXTOT*NYTOT<<endl;
        file<<"VECTORS "<<name<<" double"<<endl;
        file<<endl;

        for(int j = 0; j < NYTOT; j++){
                for(int i = 0; i < NXTOT; i++){
		        file<<0.5*(phi1[i][j] + phi1[i+1][j])<<"   "<<0.5*(phi2[i][j] + phi2[i][j+1])<<"   "<<"0.0"<<endl;
                }
        }
        file.close();
}

void CalculoSolver::PresionesVTK(double **phi, Mallador M1, string name, string Doc){
        ofstream file;
        stringstream name_t;
        string name_f;

	if(Problema == 1) name_t<<"/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosParaview/DrivenCavity/"<<Doc<<".vtk";

        else if(Problema == 2)  name_t<<"/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosParaview/DifferentiallyHeated/"<<Doc<<".vtk";

        else name_t<<"/home/sergiogus/Desktop/CTTC/NavierStokes/ResultadosParaview/SquareCylinder/"<<Doc<<".vtk";

	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<name<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET RECTILINEAR_GRID"<<endl;
         file<<"DIMENSIONS"<<"   "<<NXTOT<<"   "<<NYTOT<<"   "<<1<<endl;
        file<<endl;
        file<<"X_COORDINATES"<<"   "<<NXTOT<<"   "<<"double"<<endl;
        for(int i = 0; i < NXTOT; i++){
                file<<M1.M1_1(i,0,0)<<"   ";
        }
        file<<endl;
        file<<"Y_COORDINATES"<<"   "<<NYTOT<<"   "<<"double"<<endl;
        for(int j = 0;j < NYTOT;j++){
                file<<M1.M1_1(0,j,1)<<"   ";
        }
        file<<endl;
        file<<"Z_COORDINATES"<<"   "<<1<<"   "<<"double"<<endl;
        file<<0<<endl;
        file<<endl;
        file<<"POINT_DATA"<<"   "<<NXTOT*NYTOT<<endl;
        file<<"SCALARS "<<name<<" double"<<endl;
        file<<"LOOKUP_TABLE"<<"   "<<name<<endl;
        file<<endl;
        for(int j = 0; j < NYTOT; j++){
                for(int i = 0; i < NXTOT; i++){
            		file<<phi[i][j]<<endl;
                }
        }
        file.close();
}
