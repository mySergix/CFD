#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Memory.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/ReadData.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Geometry.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Mesher.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Solver.h"

using namespace std;

#define PI 3.141592653589793

//Constructor del mallador
Solver::Solver(Memory M1, ReadData R1, Mesher MESH){
		
	//Datos del problema
	NA = MESH.ReturnNA();
	NRad = MESH.ReturnNRad();

	TypeMesh = R1.Get_MeshingType(); //Tipo de mallado seleccionado (Collocated/Staggered)
	EsquemaAmplio = "CDS";
	Problema = R1.Get_ProblemType(); //Problema (Tobera/Tubería)

	Pr = 0.71; //Número de Prandtl

	n = 0.7;
	B = 1.0;

	PipeDiameter = 2.0;
	ReynoldsM = 3181.0;
	MachM = 1.40; //Número de Mach
	Twall = 298.15; //K

	VelocidadSonido = sqrt(1.40*287.0*Twall);
	uM = MachM*VelocidadSonido;

	B = 1.0;

	muW = B*pow(Twall,n);
	
	RhoM = (ReynoldsM*muW)/(uM*0.50*PipeDiameter);
	
	TimeBetta = 0.5; //Adams-Bashforth
	Convergencia = 1e-6;
	
	Ratio = 1.0;
	
	Cmu = 0.09;
	I = 0.05;
	k_modelo = 1.50*pow(uM,2.0)*pow(I,2.0);
	Epsilon = (Cmu*pow(k_modelo,1.50))/(0.1*(PipeDiameter/2.0));
	vt = (Cmu*pow(k_modelo,2.0))/Epsilon;
	

	Tref = 273.15;
 	muRef = 1.76e-5;
 	Sref = 110.4;

 //	muW = muRef*pow(Twall/Tref,1.50)*((Tref + Sref)/(Twall + Sref));
   


	//Propiedades termoquímicas de las especies presentes
	FracO2 = 0.82099; //Fracción de O2 en la mezcla de gases
	FracH2 = 0.135592; //Fracción de H2 en la mezcla de gases
	FracH2O = 0.0434175; //Fracción de H2O en la mezcla de gases

	MmO2 = R1.CombustionData[0]; //Masa molar del O2
	MmH2 = R1.CombustionData[3]; //Masa molar del H2	
	MmH2O = R1.CombustionData[6]; //Masa molar del H2O

	RO2 = R1.CombustionData[1]; //Constante de los gases ideales (R) del O2
	RH2 = R1.CombustionData[4]; //Constante de los gases ideales (R) del H2
	RH2O = R1.CombustionData[7]; //Constante de los gases ideales (R) del H2O

	HfO2 = R1.CombustionData[2]; //Entalpía de formación del O2 a 298.15 K 
	HfH2 = R1.CombustionData[5]; //Entalpía de formación del H2 a 298.15 K 
	HfH2O = R1.CombustionData[8]; //Entalpía de formación del H2O a 298.15 K 

	//Variables relacionadas con las viscosidad de las especies termoquímicas y la ley de Sutherland
	mu0O2 = R1.ViscosityData[0]; //Viscosidad dinámica base del O2
	TsO2 = R1.ViscosityData[1]; //Ts (Constante) del O2 (K)
	ToO2 = R1.ViscosityData[2]; //Temperatura To del O2 (K)

	mu0H2 = R1.ViscosityData[3]; //Viscosidad dinámica base del H2
	TsH2 = R1.ViscosityData[4]; //Ts (Constante) del H2 (K)
	ToH2 = R1.ViscosityData[5]; //Temperatura To del H2 (K)
 
	mu0H2O = R1.ViscosityData[6]; //Viscosidad dinámica base del H2O
	TsH2O = R1.ViscosityData[7]; //Ts (Constante) del H2O (K)
	ToH2O = R1.ViscosityData[8]; //Temperatura To del H2O (K)

	//Constantes del modelo de turbulencia Spalart-Allmaras
	Sigma = (2.0/3.0); //Sigma

	Cb1 = 0.1355; //Cb1
	Cb2 = 0.622; //Cb2
	
	K_VK = 0.41; //Constante de Von-Kárman
		
	Cw1 = Cb1/pow(K_VK,2.0) + (1.0 + Cb2)/Sigma; //Cw1
	Cw2 = 0.30; //Cw2
	Cw3 = 2.0; //Cw3

	Cv1 = 7.1; //Cv1

	Ct1 = 1.0; //Ct1
	Ct2 = 2.0; //Ct2
	Ct3 = 1.1; //Ct3
	Ct4 = 2.0; //Ct4
/*
	Sigma = R1.SpalartAllmarasData[0]; //Sigma

	Cb1 = R1.SpalartAllmarasData[1]; //Cb1
	Cb2 = R1.SpalartAllmarasData[2]; //Cb2
	
	K_VK = R1.SpalartAllmarasData[3]; //Constante de Von-Kárman
		
	Cw1 = Cb1/pow(K_VK,2.0) + (1.0 + Cb2)/Sigma; //Cw1
	Cw2 = R1.SpalartAllmarasData[4]; //Cw2
	Cw3 = R1.SpalartAllmarasData[5]; //Cw3

	Cv1 = R1.SpalartAllmarasData[6]; //Cv1

	Ct1 = R1.SpalartAllmarasData[7]; //Ct1
	Ct2 = R1.SpalartAllmarasData[8]; //Ct2
	Ct3 = R1.SpalartAllmarasData[9]; //Ct3
	Ct4 = R1.SpalartAllmarasData[10]; //Ct4*/

}

//Alojamiento de memoria para las matrices necesarias
void Solver::AllocateMatrix(Memory M1){
	
	JanafH2 = M1.AllocateDouble2D(3, 8); //Coeficientes función Cp(T) y H(T) del H2
	JanafO2 = M1.AllocateDouble2D(3, 8); //Coeficientes función Cp(T) y H(T) del O2
	JanafH2O = M1.AllocateDouble2D(2, 8); //Coeficientes función Cp(T) y H(T) del H2O
	
	//Matrices de los mapas de propiedades del problema
	
	//Densidad
	RhoPrevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de densidades previo
	RhoPres = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de densidades presente
	RhoFut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de densidades futuro

	//Presión
	Pprevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de presiones previo
	P = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de presiones presente

	//Velocidad axial (U)
	Uprevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades axiales (U) previo
	Upres = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades axiales (U) presente
	Ufut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades axiales (U) futuro

	//Velocidad radial (V)
	Vprevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades radiales (V) previo
	Vpres = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades radiales (V) presente
	Vfut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades radiales (V) futuro

	//Velocidad tangencial (W)
	Wprevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades radiales (V) previo
	Wpres = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades radiales (V) presente
	Wfut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de velocidades radiales (V) futuro

	//Temperatura
	Tpres = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de temperaturas presente
	Tfut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de temperaturas futuro

	//Energía Interna
	eprevious = M1.AllocateDouble2D(NA, NRad); 
	epres = M1.AllocateDouble2D(NA, NRad); 
	efut = M1.AllocateDouble2D(NA, NRad); 

	//Entalpía
	Hprevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de entalpías previo
	Hpres = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de entalpías presente
	Hfut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de entalpías futuro

	//Viscosidad turbulenta (cinemática)
	vprevious = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de viscosidades turbulentas previo
	vpresent = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de viscosidades turbulentas presente
	vfut = M1.AllocateDouble2D(NA, NRad); //Matriz del mapa de viscosidades turbulentas futuro

	//Matrices del modelo Spalart-Allmaras
	muTotal = M1.AllocateDouble2D(NA, NRad); //Viscosidad dinámica total
	muWallsMU = M1.AllocateDouble2D(NA+1, NRad); //Viscosidad dinámica en las paredes de los volúmenes de control (Nodos U)
	muWallsMR = M1.AllocateDouble2D(NA, NRad+1); //Viscosidad dinámica en las paredes de los volúmenes de control (Nodos R)

	nuWallsMU = M1.AllocateDouble2D(NA+1, NRad); //Viscosidad cinemática en las paredes de los volúmenes de control (Nodos U)
	nuWallsMR = M1.AllocateDouble2D(NA, NRad+1); //Viscosidad cinemática en las paredes de los volúmenes de control (Nodos R)

	UwallsMU = M1.AllocateDouble2D(NA+1, NRad); //Velocidad axial en las paredes de los volúmenes de control (Nodos U)
	UwallsMR = M1.AllocateDouble2D(NA, NRad+1); //Velocidad axial en las paredes de los volúmenes de control (Nodos R)

	VwallsMU = M1.AllocateDouble2D(NA+1, NRad); //Velocidad radial en las paredes de los volúmenes de control (Nodos U)
	VwallsMR = M1.AllocateDouble2D(NA, NRad+1); //Velocidad radial en las paredes de los volúmenes de control (Nodos R)

	WwallsMU = M1.AllocateDouble2D(NA+1, NRad); //Velocidad tangencial en las paredes de los volúmenes de control (Nodos U)
	WwallsMR = M1.AllocateDouble2D(NA, NRad+1); //Velocidad tangencial en las paredes de los volúmenes de control (Nodos R)

	eWallsMU = M1.AllocateDouble2D(NA+1, NRad); //Energia en las paredes de los volúmenes de control (Nodos U)
	eWallsMR = M1.AllocateDouble2D(NA, NRad+1); //Energia en las paredes de los volúmenes de control (Nodos R)

	//Matrices de cálculo de contribuciones a las propiedades

	//Densidad
	FRhoPrevious = M1.AllocateDouble2D(NA, NRad); 
	FRhoPresent = M1.AllocateDouble2D(NA, NRad);

	//Velocidades
	FuPrevious = M1.AllocateDouble2D(NA, NRad);
	FuPresent = M1.AllocateDouble2D(NA, NRad);

	FvPrevious = M1.AllocateDouble2D(NA, NRad);
	FvPresent = M1.AllocateDouble2D(NA, NRad);

	FwPrevious = M1.AllocateDouble2D(NA, NRad);
	FwPresent = M1.AllocateDouble2D(NA, NRad);

	//Viscosidad
	Fmuprevious = M1.AllocateDouble2D(NA, NRad);
	Fmupresent = M1.AllocateDouble2D(NA, NRad);

	//Energia Interna
	Fepresent = M1.AllocateDouble2D(NA, NRad);
	Feprevious = M1.AllocateDouble2D(NA, NRad);

	//Entalpía
	FHpresent = M1.AllocateDouble2D(NA, NRad);
	FHprevious = M1.AllocateDouble2D(NA, NRad);

	EnergyDifusive = M1.AllocateDouble2D(NA, NRad); //Término difusivo de la ecuación de conservación de la energía
	EnergyPressureTerm = M1.AllocateDouble2D(NA, NRad); //Término de presión de la ecuación de la energía
	EnergyConvective = M1.AllocateDouble2D(NA, NRad); //Término convectivo de la ecuación de energía
	EnergyViscous = M1.AllocateDouble2D(NA, NRad); //Término de disipación viscosa en la ecuación de energía

	Divergence = M1.AllocateDouble2D(NA, NRad); //Matriz de divergencias de la velocidad en cada nodo (Ecuación de energía, término viscoso)
	DivergenceMU = M1.AllocateDouble2D(NA+1, NRad);
	DivergenceMR = M1.AllocateDouble2D(NA, NRad+1);

	GradU_DxMU = M1.AllocateDouble2D(NA + 1, NRad); //Gradiente de velocidad U respecto de X en las caras de los volúmenes de control
	GradU_DyMU = M1.AllocateDouble2D(NA + 1, NRad); //Gradiente de velocidad U respecto de Y en las caras de los volúmenes de control

	GradV_DxMU = M1.AllocateDouble2D(NA + 1, NRad); //Gradiente de velocidad V respecto de X en las caras de los volúmenes de control
	GradV_DyMU = M1.AllocateDouble2D(NA + 1, NRad); //Gradiente de velocidad V respecto de Y en las caras de los volúmenes de control

	GradU_DxMR = M1.AllocateDouble2D(NA, NRad + 1); //Gradiente de velocidad U respecto de X en las caras de los volúmenes de control
	GradU_DyMR = M1.AllocateDouble2D(NA, NRad + 1); //Gradiente de velocidad U respecto de Y en las caras de los volúmenes de control

	GradV_DxMR = M1.AllocateDouble2D(NA, NRad + 1); //Gradiente de velocidad V respecto de X en las caras de los volúmenes de control
	GradV_DyMR = M1.AllocateDouble2D(NA, NRad + 1); //Gradiente de velocidad V respecto de Y en las caras de los volúmenes de control

	GradW_DxMU = M1.AllocateDouble2D(NA + 1, NRad); //Gradiente de velocidad U respecto de X en las caras de los volúmenes de control
	GradW_rMR = M1.AllocateDouble2D(NA, NRad + 1);

	K = M1.AllocateDouble2D(NA, NRad); //Conductividad térmica en los nodos
	KwallsMU = M1.AllocateDouble2D(NA+1, NRad); //Conductividad térmica en las paredes de los volúmenes de control (Nodos U)
	KwallsMV = M1.AllocateDouble2D(NA, NRad+1); //Conductividad térmica en las paredes de los volúmenes de control (Nodos R)

	Cp = M1.AllocateDouble2D(NA, NRad);
 
	//Presión
	PxGradient = M1.AllocateDouble2D(NA, NRad); //Gradiente de presión dirección Axial
	PyGradient = M1.AllocateDouble2D(NA, NRad); //Gradiente de presión dirección Radial

	PwallsU = M1.AllocateDouble2D(NA+1, NRad); //Presión en las paredes de lo volúmenes de control (Nodos U)
	PwallsR = M1.AllocateDouble2D(NA, NRad+1); //Presión en las paredes de lo volúmenes de control (Nodos V)

	MomentumDifusiveU = M1.AllocateDouble2D(NA, NRad); //Término difusivo de la ecuación de cantidad de movimiento (Velocidad axial U)
	MomentumDifusiveV = M1.AllocateDouble2D(NA, NRad); //Término difusivo de la ecuación de cantidad de movimiento (Velocidad raidal V)
	MomentumDifusiveW = M1.AllocateDouble2D(NA, NRad); //Término difusivo de la ecuación de cantidad de movimiento (Velocidad tangencial V)

	MomentumConvectiveU = M1.AllocateDouble2D(NA, NRad); //Término convective de la ecuación de cantidad de movimiento (Velocidad axial U)
	MomentumConvectiveV = M1.AllocateDouble2D(NA, NRad); //Término convective de la ecuación de cantidad de movimiento (Velocidad raidal V)
	MomentumConvectiveW = M1.AllocateDouble2D(NA, NRad); //Término convective de la ecuación de cantidad de movimiento (Velocidad tangencial V)

	//Viscosidad
	muBase = M1.AllocateDouble2D(NA, NRad);
	muTotal = M1.AllocateDouble2D(NA, NRad);
	muTurb = M1.AllocateDouble2D(NA, NRad);

	//Matrices de variables del modelo Spalart-Allmaras
	ft1 = M1.AllocateDouble2D(NA, NRad);
	r = M1.AllocateDouble2D(NA, NRad);
	g = M1.AllocateDouble2D(NA, NRad);
	ft2 = M1.AllocateDouble2D(NA, NRad);
	fw = M1.AllocateDouble2D(NA, NRad);
	X = M1.AllocateDouble2D(NA, NRad);
	S = M1.AllocateDouble2D(NA, NRad);
	Smodel = M1.AllocateDouble2D(NA, NRad);
	omega = M1.AllocateDouble2D(NA, NRad);
	fv1 = M1.AllocateDouble2D(NA, NRad);
	fv2 = M1.AllocateDouble2D(NA, NRad);

	vprevious = M1.AllocateDouble2D(NA, NRad);
	vpresent = M1.AllocateDouble2D(NA, NRad);
	vfut = M1.AllocateDouble2D(NA, NRad);

	Fmuprevious = M1.AllocateDouble2D(NA, NRad);
	Fmupresent = M1.AllocateDouble2D(NA, NRad);

	SA_Termino1 = M1.AllocateDouble2D(NA, NRad); //Término 1 de la ecuación de Spalart-Allmaras
	SA_Termino2 = M1.AllocateDouble2D(NA, NRad); //Término 2 de la ecuación de Spalart-Allmaras
	SA_Termino3 = M1.AllocateDouble2D(NA, NRad); //Término 3 de la ecuación de Spalart-Allmaras
	SA_Termino4 = M1.AllocateDouble2D(NA, NRad); //Término 4 de la ecuación de Spalart-Allmaras
	SA_Termino5 = M1.AllocateDouble2D(NA, NRad); //Término 5 de la ecuación de Spalart-Allmaras
	SA_Termino6 = M1.AllocateDouble2D(NA, NRad); //Término 6 de la ecuación de Spalart-Allmaras

	TauRR_mu = M1.AllocateDouble2D(NA+1, NRad);
	TauRR_mr = M1.AllocateDouble2D(NA, NRad+1);

	TauRZ_mu = M1.AllocateDouble2D(NA+1, NRad); 
	TauRZ_mr = M1.AllocateDouble2D(NA, NRad+1); 
		
	TauZZ_mu = M1.AllocateDouble2D(NA+1, NRad);
	TauZZ_mr = M1.AllocateDouble2D(NA, NRad+1);

	RhoWallsMU = M1.AllocateDouble2D(NA+1, NRad); //Densidad en las paredes de los volúmenes de control (Nodos U)
	RhoWallsMR = M1.AllocateDouble2D(NA, NRad+1);; //Densidad en las paredes de los volúmenes de control (Nodos R)

	//Matrices y arrays de condiciones de contorno
	
	//Condiciones de contorno del mapa de densidades (Rho)
	RhoLeft = M1.AllocateDouble1D(NRad); //Densidad pared izquierda
	RhoRight = M1.AllocateDouble1D(NRad); //Densidad pared derecha
	
	RhoUp = M1.AllocateDouble1D(NA); //Densidad pared superior
	RhoDown = M1.AllocateDouble1D(NA); //Densidad pared inferior
	
	//Condiciones de contorno del mapa de velocidades axiales (U)
	Uleft = M1.AllocateDouble1D(NRad); //Velocidad axial pared izquierda
	Uright = M1.AllocateDouble1D(NRad); //Velocidad axial pared derecha
	
	Uup = M1.AllocateDouble1D(NA); //Velocidad axial pared superior
	Udown = M1.AllocateDouble1D(NA); //Velocidad axial pared inferior
 
	//Condiciones de contorno del mapa de velocidades radiales (V)
	Vleft = M1.AllocateDouble1D(NRad); //Velocidad radial pared izquierda
	Vright = M1.AllocateDouble1D(NRad); //Velocidad radial pared derecha
	
	Vup = M1.AllocateDouble1D(NA); //Velocidad radial pared superior
	Vdown = M1.AllocateDouble1D(NA); //Velocidad radial pared inferior

	//Condiciones de contorno del mapa de velocidades tangenciales (W)
	Wleft = M1.AllocateDouble1D(NRad); //Velocidad tangencial pared izquierda
	Wright = M1.AllocateDouble1D(NRad); //Velocidad tangencial pared derecha
	
	Wup = M1.AllocateDouble1D(NA); //Velocidad tangencial pared superior
	Wdown = M1.AllocateDouble1D(NA); //Velocidad tangencial pared inferior

	//Condiciones de contorno del mapa de temperaturas (T)
	Tleft = M1.AllocateDouble1D(NRad); //Temperaturas pared izquierda
	Tright = M1.AllocateDouble1D(NRad); //Temperaturas pared derecha

	Tup = M1.AllocateDouble1D(NA); //Temperaturas pared superior
	Tdown = M1.AllocateDouble1D(NA); //Temperaturas pared inferior

	//Condiciones de contorno del mapa de presiones (P)
	Pleft = M1.AllocateDouble1D(NRad); //Presión pared izquierda
	Pright = M1.AllocateDouble1D(NRad); //Presión axial pared derecha
	
	Pup = M1.AllocateDouble1D(NA); //Presión axial pared superior
	Pdown = M1.AllocateDouble1D(NA); //Presión axial pared inferior

	//Condiciones de contorno del mapa de viscosidades dinámicas
	muLeft = M1.AllocateDouble1D(NRad); //Viscosidad dinámica pared izquierda
	muRight = M1.AllocateDouble1D(NRad); //Viscosidad dinámica pared derecha

	muUp = M1.AllocateDouble1D(NA); //Viscosidad dinámica pared superior
	muDown = M1.AllocateDouble1D(NA); //Viscosidad dinámica pared inferior

	//Condiciones de contorno del mapa de viscosidades dinámicas
	KLeft = M1.AllocateDouble1D(NRad); //Conductividad térmica pared izquierda
	KRight = M1.AllocateDouble1D(NRad); //Conductividad térmica pared derecha

	KUp = M1.AllocateDouble1D(NA); //Conductividad térmica pared superior
	KDown = M1.AllocateDouble1D(NA); //Conductividad térmica pared inferior

	//Condiciones de contorno del mapa de energías internas
	eleft = M1.AllocateDouble1D(NRad); //Energía Interna pared izquierda
	eright = M1.AllocateDouble1D(NRad); //Energía Interna pared derecha

	eup = M1.AllocateDouble1D(NA); //Energía Interna pared superior
	edown = M1.AllocateDouble1D(NA); //Energía Interna pared inferior

	//Condiciones de contorno del mapa de entalpías
	Hleft = M1.AllocateDouble1D(NRad); //Entalpía pared izquierda
	Hright = M1.AllocateDouble1D(NRad); //Entalpía pared derecha

	Hup = M1.AllocateDouble1D(NA); //Entalpía pared superior
	Hdown = M1.AllocateDouble1D(NA); //Entalpía pared inferior

	//Condiciones de contorno del mapa de viscosidades cinemáticas turbulentas
	vleft = M1.AllocateDouble1D(NRad); //Viscosidad cinemática turbulenta pared izquierda
	vright = M1.AllocateDouble1D(NRad); //Viscosidad cinemática turbulenta pared derecha

	vup = M1.AllocateDouble1D(NA); //Viscosidad cinemática turbulenta pared arriba
	vdown = M1.AllocateDouble1D(NA); //Viscosidad cinemática turbulenta pared abajo

	//Condiciones de contorno del mapa de viscosidades dinámicas base
	muBaseLeft = M1.AllocateDouble1D(NRad); //Viscosidad dinámica base pared izquierda
	muBaseRight = M1.AllocateDouble1D(NRad); //Viscosidad dinámica base pared derecha

	muBaseDown = M1.AllocateDouble1D(NA); //Viscosidad dinámica base pared abajo
	muBaseUp = M1.AllocateDouble1D(NA); //Viscosidad dinámica base pared arriba


	//Valores medios de las propiedades del flujo
	DensidadMedia = M1.AllocateDouble1D(NRad);
	TemperaturaMedia = M1.AllocateDouble1D(NRad);
	PresionMedia = M1.AllocateDouble1D(NRad);

	DensidadMediaTemporal = M1.AllocateDouble1D(NRad);
	TemperaturaMediaTemporal = M1.AllocateDouble1D(NRad);
	PresionMediaTemporal = M1.AllocateDouble1D(NRad);

	DensidadMediaTemporalDato = M1.AllocateDouble1D(NRad);
	TemperaturaMediaTemporalDato = M1.AllocateDouble1D(NRad);
	PresionMediaTemporalDato = M1.AllocateDouble1D(NRad);

	//Código;
	//1 -> West (aW)
	//2 -> East (aE)
	//3 -> South (aS)
	//4 -> North (aN)
	//5 -> Central (aP)
	
}

//Inicialización de los campos de temperaturas
void Solver::InitializeFields(Mesher MESH){
int i, j;
double I, J;
double nrad = NRad;
double na = NA;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			
			J = j;
			//Campos de Temperaturas
			Tpres[i][j] = Twall;
	
			//Campos de densidades
			RhoPrevious[i][j] = RhoM;
			RhoPres[i][j] = RhoM;
			
			//Campos de velocidades axiales (U)
			Uprevious[i][j] = uM;
			//*((pow(PipeDiameter/2.0,2.0) - pow(MESH.MP[i][j][1],2.0))/pow(PipeDiameter/2.0,2.0));
			Upres[i][j] = uM;
			//*((pow(PipeDiameter/2.0,2.0) - pow(MESH.MP[i][j][1],2.0))/pow(PipeDiameter/2.0,2.0));
			
			//Campos de velocidades radiales (V)
			Vprevious[i][j] = 0.0;
			Vpres[i][j] = 0.0;

			//Campos de energías internas
			eprevious[i][j] = (Cp[i][j]/1.40)*Tpres[i][j];
			epres[i][j] = (Cp[i][j]/1.40)*Tpres[i][j];

			//Campos de entalpías
			Hprevious[i][j] = Cp[i][j]*Tpres[i][j];
			Hpres[i][j] = Cp[i][j]*Tpres[i][j];

			muBase[i][j] = B*pow(Tpres[i][j],n);
			//(RhoPres[i][j]*Upres[i][j]*PipeDiameter)/(ReynoldsM);

			//Campos de viscosidades cinemáticas
			vprevious[i][j] = 0.0;
			vpresent[i][j] = 0.0;
			
			//Matrices de variables del modelo Spalart-Allmaras
			X[i][j] = 0.0;
			fv1[i][j] = 0.0;
			fv2[i][j] = 0.0;
			r[i][j] = 0.0;
			g[i][j] = 0.0;
			fw[i][j] = 0.0;
			ft2[i][j] = 0.0;
			omega[i][j] = 0.0;
			S[i][j] = 0.0;
			Smodel[i][j] = 0.0;
		
		}
	}		
}


//Asignación de temperaturas a las condiciones de contorno
void Solver::UpdateBoundaryConditions(Mesher MESH){
int i, j;
double I, J;
double nrad = NRad;

	for(i = 0; i < NA; i++){
		//Densidad parte arriba
		RhoUp[i] = RhoPres[i][NRad-1];

		//Densidad parte abajo
		RhoDown[i] = RhoPres[i][0];

		//Presión parte abajo
		Pdown[i] = P[i][0];

		//Presión parte arriba
		Pup[i] = P[i][NRad-1];

		//Velocidad radial (V) parte ariba
		Vup[i] = 0.0;

		//Velocidad axial (U) parte arriba
		Uup[i] = 0.0;

		//Velocidad radial (V) parte abajo
		Vdown[i] = Vpres[i][0];

		//Velocidad axial (U) parte abajo
		Udown[i] = Upres[i][0];
	
		//Velocidad tangencial (W) parte abajo
		Wdown[i] = Wpres[i][0];

		//Temperatura parte arriba
		Tup[i] = Twall;

		//Temperatura parte abajo
		Tdown[i] = Tpres[i][0];

		//Viscosidad parte abajo
		muDown[i] = muTotal[i][0];

		//Viscosidad parte arriba
		muUp[i] = muTotal[i][NRad-1];

		//Conductividad térmica parte abajo
		KDown[i] = K[i][0];

		//Conductividad térmica parte arriba
		KUp[i] = K[i][NRad-1];

		//Energía Interna parte abajo
		edown[i] = epres[i][0];

		//Energía Interna parte arriba
		eup[i] = Cp[i][NRad-1]*Tup[i]/1.40;

		//Entalpía parte abajo
		Hdown[i] = Hpres[i][0];

		//Entalpía parte arriba
		Hup[i] = Hpres[i][NRad-1];

		//Viscosidad turbulenta parte abajo
		vdown[i] = vpresent[i][0];
		
		//Viscosidad turbulenta parte arriba
		vup[i] = 0.0;

		//Viscosidad dinámica base pared abajo
		muBaseDown[i] = muBase[i][0]; 

		//Viscosidad dinámica base pared arriba
		muBaseUp[i] = muBase[i][NRad-1]; 

	}

	for(j = 0; j < NRad; j++){
		//Densidad parte izquierda
		RhoLeft[j] = RhoM;
		
		//Densidad parte derecha
		RhoRight[j] = RhoPres[NA-1][j];

		//Presión parte izquierda
		Pleft[j] = P[0][j];

		//Presión parte derecha
		Pright[j] = P[NA-1][j];

		//Velocidad radial (V) parte izquierda
		Vleft[j] = 0.0;

		//Velocidad axial (U) parte izquierda
		Uleft[j] = uM;
		
		//Velocidad radial (V) parte derecha
		Vright[j] = Vpres[NA-1][j];

		//Velocidad axial (U) parte derecha
		Uright[j] = Upres[NA-1][j];
	
		//Velocidad tangencial (W) parte derecha
		Wright[j] = Wpres[NA-1][j];

		//Temperatura parte izquierda
		Tleft[j] = Twall;

		//Temperatura parte derecha
		Tright[j] = Tpres[NA-1][j];

		//Viscosidad parte izquierda
		muLeft[j] = muTotal[0][j];

		//Viscosidad parte derecha
		muRight[j] = muTotal[NA-1][j];

		//Conductividad térmica parte izquierda
		KLeft[j] = K[0][j];

		//Conductividad térmica parte derecha
		KRight[j] = K[NA-1][j];

		//Energía Interna parte izquierda
		eleft[j] = Tleft[j]*Cp[0][j]/1.40;

		//Energía Interna parte derecha
		eright[j] = epres[NA-1][j];
		
		//Entalpía parte izquierda
		Hleft[j] = Hpres[0][j];

		//Entalpía parte derecha
		Hright[j] = Hpres[NA-1][j];

		//Viscosidad turbulenta parte izquierda
		vleft[j] = 4.0*B*pow(Twall,n)/RhoM;
		//
		//1.0*B*pow(Twall,n)/RhoM;

		//4.0*B*pow(Twall,n)/RhoM;
		//
		//

		//Viscosidad turbulenta parte derecha
		vright[j] = vpresent[NA-1][j];
		//1.0*B*pow(Twall,n)/RhoM;
		
		//Viscosidad dinámica base pared izquierda
		muBaseLeft[j] = muBase[0][j]; 

		//Viscosidad dinámica base pared derecha
		muBaseRight[j] = muBase[NA-1][j]; 

	}
	
}

void Solver::Get_CpAir(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Cp[i][j] = 1000.0;
		}
	}

}


//Cálculo de la constante ideal de la mezcla de gases
void Solver::Get_Rideal(){ Rideal = 287.0; }
//Rideal = FracO2*RO2 + FracH2*RH2 + FracH2O*RH2O; } // (J/KgK) 

void Solver::Get_AirViscosity(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			muBase[i][j] = B*pow(Tpres[i][j],n);
			//(RhoPres[i][j]*sqrt(pow(Upres[i][j],2.0) + pow(Vpres[i][j],2.0))*PipeDiameter)/ReynoldsM;
		}
	}

}

//Cálculo del DeltaT para el siguiente Step
double Solver::StepTime(Mesher MESH){
int i, j;
double DeltaT = 1000.0;
double Tpar = 0.05;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){

			//CFL Temperatura
			if(abs((Tpar*RhoPres[i][j]*Cp[i][j]*pow(MESH.DeltasMP[i][j][0],2.0))/(K[i][j]+1e-10)) < DeltaT){
				DeltaT = abs((Tpar*RhoPres[i][j]*Cp[i][j]*pow(MESH.DeltasMP[i][j][0],2.0))/(K[i][j]+1e-10));
			}
			if(abs((Tpar*RhoPres[i][j]*Cp[i][j]*pow(MESH.DeltasMP[i][j][1],2.0))/(K[i][j]+1e-10)) < DeltaT){
				DeltaT = abs((Tpar*RhoPres[i][j]*Cp[i][j]*pow(MESH.DeltasMP[i][j][1],2.0))/(K[i][j]+1e-10));
			}

			//CFL Velocidades (U + V)
			if(abs((Tpar*MESH.DeltasMP[i][j][0])/(VelocidadSonido + sqrt(pow(Upres[i][j],2.0) + 1e-10))) < DeltaT){
				DeltaT = abs((Tpar*MESH.DeltasMP[i][j][0])/(VelocidadSonido + sqrt(pow(Upres[i][j],2.0) + 1e-10 )));
			}
			if(abs((Tpar*MESH.DeltasMP[i][j][1])/(VelocidadSonido + abs(Vpres[i][j]) + 1e-10)) < DeltaT){
				DeltaT = abs((Tpar*MESH.DeltasMP[i][j][1])/(VelocidadSonido + abs(Vpres[i][j]) + 1e-10));
			}

			//CFL Difusivo
			if((Tpar*RhoPres[i][j]*pow(MESH.DeltasMP[i][j][0],2.0))/(muTotal[i][j]+1e-10) < DeltaT){
				DeltaT = (Tpar*RhoPres[i][j]*pow(MESH.DeltasMP[i][j][0],2))/(muTotal[i][j]+1e-10);
			}
			if((Tpar*RhoPres[i][j]*pow(MESH.DeltasMP[i][j][1],2))/(muTotal[i][j]+1e-10) < DeltaT){
				DeltaT = (Tpar*RhoPres[i][j]*pow(MESH.DeltasMP[i][j][1],2))/(muTotal[i][j]+1e-10);
			}	
		}
	}
	
	return DeltaT;
	
}

//Seteo inicial de las F
void Solver::InitialF(){
int i, j;
	
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			FRhoPrevious[i][j] = 0.0;
			FuPrevious[i][j] = 0.0;
			FvPrevious[i][j] = 0.0;
			Feprevious[i][j] = 0.0;
			Fmuprevious[i][j] = 0.0;
		}
	}
	
}

//Cálculo de las velocidades en las paredes de los volúmenes de control
void Solver::Get_VelocityWalls(Mesher MESH){
int i, j;

	//Nodos U
	for(j = 0; j < NRad; j++){
		UwallsMU[0][j] = Uleft[j];
		UwallsMU[NA][j] = Uright[j];

		VwallsMU[0][j] = Vleft[j];
		VwallsMU[NA][j] = Vright[j];

 		for(i = 1; i < NA; i++){
			UwallsMU[i][j] = 0.50*(Upres[i][j] + Upres[i-1][j]);
			VwallsMU[i][j] = 0.50*(Vpres[i][j] + Vpres[i-1][j]);
		}
	}

	//Nodos R
	for(i = 0; i < NA; i++){
		UwallsMR[i][0] = Udown[i];
		UwallsMR[i][NRad] = Uup[i];

		VwallsMR[i][0] =  Vdown[i];
		VwallsMR[i][NRad] = Vup[i];

		for(j = 1; j < NRad; j++){
			UwallsMR[i][j] = 0.50*(Upres[i][j] + Upres[i][j-1]);
			
			VwallsMR[i][j] = 0.50*(Vpres[i][j] + Vpres[i][j-1]);
			
		}
	}
}

void Solver::Get_RhoWalls(Mesher MESH){
int i, j;

	//Nodos U
	for(j = 0; j < NRad; j++){
		RhoWallsMU[0][j] = RhoLeft[j];	
		RhoWallsMU[NA][j] = RhoRight[j];
		
		for(i = 1; i < NA; i++){
			RhoWallsMU[i][j] = sqrt(RhoPres[i][j])*sqrt(RhoPres[i-1][j]);
		}
	}

	//Nodos R
	for(i = 0; i < NA; i++){
		RhoWallsMR[i][0] = RhoDown[i];
		RhoWallsMR[i][NRad] = RhoUp[i];
	
		for(j = 1; j < NRad; j++){
			RhoWallsMR[i][j] = sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j-1]);
		}
	}
}

//Calculo de las FRho (Densidad, Conservación de masa)
void Solver::Get_FRho(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			FRhoPresent[i][j] = -(1.0/MESH.VolMP[i][j])*(
								MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*RhoWallsMU[i][j]
							  + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*RhoWallsMU[i+1][j]
							  + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*RhoWallsMR[i][j]
							  + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*RhoWallsMR[i][j+1]
							  );
		}
	}	
}

//Cálculo de las presiones en las paredes de los volúmenes de control
void Solver::Get_Pwalls(Mesher MESH){
int i, j;

	//Nodos U 
	for(j = 0; j < NRad; j++){
		PwallsU[0][j] = RhoLeft[j]*Rideal*Tleft[j];
		PwallsU[NA][j] = RhoRight[j]*Rideal*Tright[j];
		for(i = 1; i < NA; i++){
			PwallsU[i][j] = 0.50*(P[i][j] + P[i-1][j]);
		}
	}

	//Nodos V
	for(i = 0; i < NA; i++){
		PwallsR[i][0] = RhoDown[i]*Rideal*Tdown[i];
		PwallsR[i][NRad] = Pup[i];
		for(j = 1; j < NRad; j++){
			PwallsR[i][j] =  0.50*(P[i][j] + P[i][j-1]);
		}
	}

}
//Cálculo de los gradientes de presión en cada nodo en ambas direciones (U y V)
void Solver::Get_Pgradients(Mesher MESH){
int i, j;

	//Centro Gradiente Axial
	for(j = 0; j < NRad; j++){
		PxGradient[0][j] = 0.0;
		//PxGradient[NA-1][j] = 0.0;
		for(i = 1; i < NA; i++){
			PxGradient[i][j] = (1.0/MESH.DeltasMP[i][j][0])*(PwallsU[i+1][j] - PwallsU[i][j]);
		}
	}

	//Centro Gradiente Radial
	for(i = 0; i < NA; i++){
		PyGradient[i][0] = 0.0;
		PyGradient[i][NRad-1] = 0.0;
		for(j = 1; j < NRad-1; j++){
			PyGradient[i][j] = (1.0/MESH.VolMP[i][j])*(MESH.SupMP[i][j][2]*PwallsR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]) + MESH.SupMP[i][j][3]*PwallsR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1])) - P[i][j]/MESH.MP[i][j][1];			
		}
	}
}

//Cálculo de las viscosidades dinámicas en las paredes de los volúmenes de control
void Solver::Get_muWalls(Mesher MESH){
int i, j;

	//Nodos U
	for(j = 0; j < NRad; j++){
		muWallsMU[0][j] = muTotal[0][j];
		muWallsMU[NA][j] = muTotal[NA-1][j];
		for(i = 1; i < NA; i++){
			muWallsMU[i][j] = 0.50*(muTotal[i-1][j] + muTotal[i][j]);
		}

	}

	//Nodos R
	for(i = 0; i < NA; i++){
		muWallsMR[i][0] = muTotal[i][0];	
		muWallsMR[i][NRad] = muTotal[i][NRad-1];
		for(j = 1; j < NRad; j++){	
			muWallsMR[i][j] =  0.50*(muTotal[i][j] + muTotal[i][j-1]);
		}

	}

}

//Cálculo de las viscosidades cinemáticas en las paredes de los volúmenes de control
void Solver::Get_nuWalls(Mesher MESH){
int i,j;

	//Nodos U
	for(j = 0; j < NRad; j++){
		nuWallsMU[0][j] = muBase[0][j]/RhoLeft[j] + vleft[j];	
		nuWallsMU[NA][j] = muBase[NA-1][j]/RhoRight[j] + vright[j];
		for(i = 1; i < NA; i++){
			nuWallsMU[i][j] = 0.50*(muBase[i][j]/RhoPres[i][j] + muBase[i-1][j]/RhoPres[i-1][j]) + 0.50*(vpresent[i][j] + vpresent[i-1][j]);
		}

	}

	//Nodos R
	for(i = 0; i < NA; i++){
		nuWallsMR[i][0] = muBase[i][0]/RhoDown[i] + vdown[i];	
		nuWallsMR[i][NRad] = muBase[i][NRad-1]/RhoUp[i] + vup[i];
		for(j = 1; j < NRad; j++){	
			nuWallsMR[i][j] = 0.50*(muBase[i][j]/RhoPres[i][j] + muBase[i][j-1]/RhoPres[i][j-1]) + 0.50*(vpresent[i][j] + vpresent[i][j-1]);
		}

	}
}

//Cálculo de los gradientes de velocidades en las caras de los volúmenes de control
void Solver::Get_VelocityGradients(Mesher MESH){
int i, j;
	
	//Nodos U, Gradiente X
	for(j = 0; j < NRad; j++){
		GradU_DxMU[0][j] = (1.0/MESH.DeltasMU[0][j][0])*(Upres[0][j] - Uleft[j]);
		GradU_DxMU[NA][j] = (1.0/MESH.DeltasMU[NA][j][0])*(Uright[j] - Upres[NA-1][j]);

		GradV_DxMU[0][j] = (1.0/MESH.DeltasMU[0][j][0])*(Vpres[0][j] - Vleft[j]);
		GradV_DxMU[NA][j] = (1.0/MESH.DeltasMU[NA][j][0])*(Vright[j] - Vpres[NA-1][j]);

		for(i = 1; i < NA; i++){
			GradU_DxMU[i][j] = (1.0/MESH.DeltasMU[i][j][0])*(Upres[i][j] - Upres[i-1][j]);

			GradV_DxMU[i][j] = (1.0/MESH.DeltasMU[i][j][0])*(Vpres[i][j] - Vpres[i-1][j]);
		}
	}

	//Nodos U, Gradiente Y
	for(i = 1; i < NA; i++){
		GradU_DyMU[i][0] = (1.0/MESH.DeltasMU[i][0][1])*(0.50*(UwallsMU[i][0] + UwallsMU[i][1]) - 0.50*(Udown[i] + Udown[i-1]));
		GradU_DyMU[i][NRad-1] = (1.0/MESH.DeltasMU[i][NRad-1][1])*(0.50*(Uup[i] + Uup[i-1]) - 0.50*(UwallsMU[i][NRad-1] + UwallsMU[i][NRad-2]));

		GradV_DyMU[i][0] = (1.0/MESH.DeltasMU[i][0][1])*(0.50*(VwallsMU[i][0] + VwallsMU[i][1]) - 0.50*(Vdown[i] + Vdown[i-1]));
		GradV_DyMU[i][NRad-1] = (1.0/MESH.DeltasMU[i][NRad-1][1])*(0.50*(Vup[i] + Vup[i-1]) - 0.50*(VwallsMU[i][NRad-1] + VwallsMU[i][NRad-2]));

		for(j = 1; j < NRad-1; j++){
			GradU_DyMU[i][j] = (1.0/MESH.DeltasMU[i][j][1])*(0.50*(UwallsMU[i][j] + UwallsMU[i][j+1]) - 0.50*(UwallsMU[i][j] + UwallsMU[i][j-1]));

			GradV_DyMU[i][j] = (1.0/MESH.DeltasMU[i][j][1])*(0.50*(VwallsMU[i][j] + VwallsMU[i][j+1]) - 0.50*(VwallsMU[i][j] + VwallsMU[i][j-1]));
		}	
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		GradU_DyMU[0][j] = (1.0/MESH.DeltasMU[0][j][1])*(0.50*(UwallsMU[0][j] + UwallsMU[0][j+1]) - 0.50*(UwallsMU[0][j] + UwallsMU[0][j-1]));

		GradV_DyMU[0][j] = (1.0/MESH.DeltasMU[0][j][1])*(0.50*(VwallsMU[0][j] + VwallsMU[0][j+1]) - 0.50*(VwallsMU[0][j] + VwallsMU[0][j-1]));

		//Parte derecha
		GradU_DyMU[NA][j] = (1.0/MESH.DeltasMU[NA][j][1])*(0.50*(UwallsMU[NA][j] + UwallsMU[NA][j+1]) - 0.50*(UwallsMU[NA][j] + UwallsMU[NA][j-1]));

		GradV_DyMU[NA][j] = (1.0/MESH.DeltasMU[NA][j][1])*(0.50*(VwallsMU[NA][j] + VwallsMU[NA][j+1]) - 0.50*(VwallsMU[NA][j] + VwallsMU[NA][j-1]));

	}

	//Esquina abajo izquierda
	GradU_DyMU[0][0] = (1.0/MESH.DeltasMU[0][0][1])*(0.50*(UwallsMU[0][0] + UwallsMU[0][1]) - Udown[0]);
	GradV_DyMU[0][0] = (1.0/MESH.DeltasMU[0][0][1])*(0.50*(VwallsMU[0][0] + VwallsMU[0][1]) - Vdown[0]);

	//Esquina arriba izquierda
	GradU_DyMU[0][NRad-1] = (1.0/MESH.DeltasMU[0][NRad-1][1])*(Uup[0] - 0.50*(UwallsMU[0][NRad-1] + UwallsMU[0][NRad-2]));
	GradV_DyMU[0][NRad-1] = (1.0/MESH.DeltasMU[0][NRad-1][1])*(Vup[0] - 0.50*(VwallsMU[0][NRad-1] + VwallsMU[0][NRad-2]));

	//Esquina abajo derecha
	GradU_DyMU[NA][0] = (1.0/MESH.DeltasMU[NA][0][1])*(0.50*(UwallsMU[NA][0] + UwallsMU[NA][1]) - Udown[NA-1]);
	GradV_DyMU[NA][0] = (1.0/MESH.DeltasMU[NA][0][1])*(0.50*(VwallsMU[NA][0] + VwallsMU[NA][1]) - Vdown[NA-1]);

	//Esquina arriba derecha
	GradU_DyMU[NA][NRad-1] = (1.0/MESH.DeltasMU[NA][NRad-1][1])*(Uup[NA-1] - 0.50*(UwallsMU[NA][NRad-1] + UwallsMU[NA][NRad-2]));
	GradV_DyMU[NA][NRad-1] = (1.0/MESH.DeltasMU[NA][NRad-1][1])*(Vup[NA-1] - 0.50*(VwallsMU[NA][NRad-1] + VwallsMU[NA][NRad-2]));



	//Nodos R, Gradiente Y
	for(i = 0; i < NA; i++){
		GradU_DyMR[i][0] = (1.0/MESH.DeltasMR[i][0][1])*(Upres[i][0] - Udown[i]);
		GradU_DyMR[i][NRad] = (1.0/MESH.DeltasMR[i][NRad][1])*(Uup[i] - Upres[i][NRad-1]);
		
		GradV_DyMR[i][NRad] = (1.0/MESH.DeltasMR[i][NRad][1])*(Vup[i] - Vpres[i][NRad-1]);
		GradV_DyMR[i][0] = (1.0/MESH.DeltasMR[i][0][1])*(Vpres[i][0] - Vdown[i]);
		for(j = 1; j < NRad; j++){
			GradU_DyMR[i][j] = (1.0/MESH.DeltasMR[i][j][1])*(Upres[i][j] - Upres[i][j-1]);

			GradV_DyMR[i][j] = (1.0/MESH.DeltasMR[i][j][1])*(Vpres[i][j] - Vpres[i][j-1]);
		}	
	}
  
	//Nodos R, Gradiente X
	for(i = 1; i < NA-1; i++){	
		GradU_DxMR[i][0] = (1.0/MESH.DeltasMR[i][0][0])*(0.50*(UwallsMR[i+1][0] + UwallsMR[i][0]) - 0.50*(UwallsMR[i-1][0] + UwallsMR[i][0]));
		GradU_DxMR[i][NRad] = (1.0/MESH.DeltasMR[i][NRad][0])*(0.50*(UwallsMR[i+1][NRad] + UwallsMR[i][NRad]) - 0.50*(UwallsMR[i-1][NRad] + UwallsMR[i][NRad]));

		GradV_DxMR[i][0] = (1.0/MESH.DeltasMR[i][0][0])*(0.50*(VwallsMR[i+1][0] + VwallsMR[i][0]) - 0.50*(VwallsMR[i-1][0] + VwallsMR[i][0]));
		GradV_DxMR[i][NRad] = (1.0/MESH.DeltasMR[i][NRad][0])*(0.50*(VwallsMR[i+1][NRad] + VwallsMR[i][NRad]) - 0.50*(VwallsMR[i-1][NRad] + VwallsMR[i][NRad]));

		for(j = 1; j < NRad-1; j++){
			GradU_DxMR[i][j] = (1.0/MESH.DeltasMR[i][j][0])*(0.50*(UwallsMR[i+1][j] + UwallsMR[i][j]) - 0.50*(UwallsMR[i-1][j] + UwallsMR[i][j]));

			GradV_DxMR[i][j] = (1.0/MESH.DeltasMR[i][j][0])*(0.50*(VwallsMR[i+1][j] + VwallsMR[i][j]) - 0.50*(VwallsMR[i-1][j] + VwallsMR[i][j]));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		GradU_DxMR[0][j] = (1.0/MESH.DeltasMR[0][j][0])*(0.50*(UwallsMR[1][j] + UwallsMR[0][j]) - 0.50*(Uleft[j] + Uleft[j-1]));
		GradV_DxMR[0][j] = (1.0/MESH.DeltasMR[0][j][0])*(0.50*(VwallsMR[1][j] + VwallsMR[0][j]) - 0.50*(Vleft[j] + Vleft[j-1]));

		//Parte derecha
		GradU_DxMR[NA-1][j] = (1.0/MESH.DeltasMR[NA-1][j][0])*(0.50*(Uright[j] + Uright[j-1]) - 0.50*(UwallsMR[NA-2][j] + UwallsMR[NA-1][j]));
		GradV_DxMR[NA-1][j] = (1.0/MESH.DeltasMR[NA-1][j][0])*(0.50*(Vright[j] + Vright[j-1]) - 0.50*(VwallsMR[NA-2][j] + VwallsMR[NA-1][j]));

	}

	//Esquina abajo izquierda
	GradU_DxMR[0][0] = (1.0/MESH.DeltasMR[0][0][0])*(0.50*(UwallsMR[0][0] + UwallsMR[1][0]) - Uleft[0]);
	GradV_DxMR[0][0] = (1.0/MESH.DeltasMR[0][0][0])*(0.50*(VwallsMR[0][0] + VwallsMR[1][0]) - Vleft[0]);

	//Esquina arriba izquierda
	GradU_DxMR[0][NRad] = (1.0/MESH.DeltasMR[0][NRad-1][0])*(0.50*(UwallsMR[0][NRad-1] + UwallsMR[1][NRad-1]) - UwallsMU[0][NRad-1]);
	GradV_DxMR[0][NRad] = (1.0/MESH.DeltasMR[0][NRad-1][0])*(0.50*(VwallsMR[0][NRad-1] + VwallsMR[1][NRad-1]) - VwallsMU[0][NRad-1]);

	//Esquina abajo derecha
	GradU_DxMR[NA-1][0] = (1.0/MESH.DeltasMR[NA-1][0][0])*(0.50*(UwallsMU[NA][0] + Udown[NA-1]) - 0.50*(UwallsMU[NA-1][0] + Udown[NA-1]));
	GradV_DxMR[NA-1][0] = (1.0/MESH.DeltasMR[NA-1][0][0])*(0.50*(VwallsMU[NA][0] + Vdown[NA-1]) - 0.50*(VwallsMU[NA-1][0] + Vdown[NA-1]));

	//Esquina arriba derecha
	GradU_DxMR[NA-1][NRad] = (1.0/MESH.DeltasMR[NA-1][NRad][0])*(0.50*(UwallsMU[NA][NRad-1] + Uup[NA-1]) - 0.50*(UwallsMU[NA-1][NRad-1] + Uup[NA-1]));

	GradV_DxMR[NA-1][NRad] = (1.0/MESH.DeltasMR[NA-1][NRad][0])*(0.50*(VwallsMU[NA][NRad-1] + Vup[NA-1]) - 0.50*(VwallsMU[NA-1][NRad-1] + Vup[NA-1]));
	
}

//Cálculo de la divergencia de la velocidad en cada nodo (Término viscoso ecuación de energía)
void Solver::Get_Divergence(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Divergence[i][j] = (1.0/MESH.VolMP[i][j])*(
							   MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))
							 + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))
							 + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
							 + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
							 );
		}
	}

	//Nodos U
	//Centro
	for(j = 0; j < NRad; j++){
		DivergenceMU[0][j] = GradU_DxMU[0][j] + GradV_DyMU[0][j] + VwallsMU[0][j]/MESH.MU[0][j][1];
		DivergenceMU[NA][j] = GradU_DxMU[NA][j] + GradV_DyMU[NA][j] + VwallsMU[NA][j]/MESH.MU[NA][j][1];	
		for(i = 1; i < NA; i++){
			DivergenceMU[i][j] = 0.50*(Divergence[i][j] + Divergence[i-1][j]);
		}
	}

	//Nodos R
	//Centro
	for(i = 0; i < NA; i++){	
		DivergenceMR[i][NRad] = Divergence[i][NRad-1];	
		for(j = 1; j < NRad-1; j++){
			DivergenceMR[i][j] = 0.50*(Divergence[i][j] + Divergence[i][j-1]);
		}
	}

}

void Solver::Get_Stresses(Mesher MESH){
int i, j;

	for(i = 0; i < NA+1; i++){
		for(j = 0; j < NRad; j++){
			TauRR_mu[i][j] = -muWallsMU[i][j]*(2.0*GradV_DyMU[i][j] - (2.0/3.0)*DivergenceMU[i][j]);

			TauRZ_mu[i][j] = -muWallsMU[i][j]*(GradU_DyMU[i][j] + GradV_DxMU[i][j]);

			TauZZ_mu[i][j] = -muWallsMU[i][j]*(2.0*GradU_DxMU[i][j] - (2.0/3.0)*DivergenceMU[i][j]);
		}
	}

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad+1; j++){
			TauRR_mr[i][j] = -muWallsMR[i][j]*(2.0*GradV_DyMR[i][j] - (2.0/3.0)*DivergenceMR[i][j]);

			TauRZ_mr[i][j] = -muWallsMR[i][j]*(GradU_DyMR[i][j] + GradV_DxMR[i][j]);

			TauZZ_mr[i][j] = -muWallsMR[i][j]*(2.0*GradU_DxMR[i][j] - (2.0/3.0)*DivergenceMR[i][j]);
		}
	}
}

//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveU(Mesher MESH){
int i, j;

for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumDifusiveU[i][j] = -(1.0/MESH.VolMP[i][j])*(
									  MESH.SupMP[i][j][0]*(TauZZ_mu[i][j]*cos(PI + MESH.AngleMU[i][j]) + TauRZ_mu[i][j]*sin(PI + MESH.AngleMU[i][j]))
									+ MESH.SupMP[i][j][1]*(TauZZ_mu[i+1][j]*cos(MESH.AngleMU[i+1][j]) + TauRZ_mu[i+1][j]*sin(MESH.AngleMU[i+1][j]))  
									+ MESH.SupMP[i][j][2]*(TauZZ_mr[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + TauRZ_mr[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
									+ MESH.SupMP[i][j][3]*(TauZZ_mr[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + TauRZ_mr[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
									);
		}
	}

}


//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveV(Mesher MESH){
int i, j;
	
	for(i = 0; i < NA; i++){
		for(j = 1; j < NRad; j++){
			MomentumDifusiveV[i][j] = -(1.0/MESH.VolMP[i][j])*( 
									+ MESH.SupMP[i][j][0]*(TauRZ_mu[i][j]*cos(PI + MESH.AngleMU[i][j]) + TauRR_mu[i][j]*sin(PI + MESH.AngleMU[i][j]))
									+ MESH.SupMP[i][j][1]*(TauRZ_mu[i+1][j]*cos(MESH.AngleMU[i+1][j]) + TauRR_mu[i+1][j]*sin(MESH.AngleMU[i+1][j]))  
									+ MESH.SupMP[i][j][2]*(TauRZ_mr[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + TauRR_mr[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
									+ MESH.SupMP[i][j][3]*(TauRZ_mr[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + TauRR_mr[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
									)
									- (2.0*muTotal[i][j]*Vpres[i][j])/pow(MESH.MP[i][j][1],2.0)
									+ (2.0/3.0)*(muTotal[i][j]/MESH.MP[i][j][1])*Divergence[i][j]
									;
		}
	}
}

//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
void Solver::Get_MomentumConvectiveU(Mesher MESH){
int i, j;

	//Centro
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumConvectiveU[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*RhoWallsMU[i][j]*UwallsMU[i][j]
									  + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*RhoWallsMU[i+1][j]*UwallsMU[i+1][j]
									  + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*RhoWallsMR[i][j]*UwallsMR[i][j]
									  + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*RhoWallsMR[i][j+1]*UwallsMR[i][j+1]
									  );
		}
	}
									  
}

//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomentumConvectiveV(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumConvectiveV[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*RhoWallsMU[i][j]*VwallsMU[i][j]
									  + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*RhoWallsMU[i+1][j]*VwallsMU[i+1][j]
									  + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*RhoWallsMR[i][j]*VwallsMR[i][j]
									  + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*RhoWallsMR[i][j+1]*VwallsMR[i][j+1]
									  );
		}
	}
	
}

//Cálculo de la conductividad térmica de los gases en cada nodo y paredes
void Solver::Get_K(Mesher MESH){
int i, j;
double Predict;

	//Cálculo de K en los nodos
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			K[i][j] = (Cp[i][j]*muTotal[i][j])/Pr;
		}
	}

	//Cálculo de K en las paredes de los volúmenes de control	
	//Nodos U
	for(j = 0; j < NRad; j++){
		KwallsMU[0][j] = KLeft[j]; //Parte izquierda
		KwallsMU[NA][j] = KRight[j]; //Parte derecha
		for(i = 1; i < NA; i++){
			KwallsMU[i][j] = 0.50*(K[i][j] + K[i-1][j]); //Centro
		}
	}
	
	//Nodos V
	for(i = 0; i < NA; i++){
		KwallsMV[i][0] = K[i][0]; //Parte abajo
		KwallsMV[i][NRad] = KUp[i]; //Parte arriba
		for(j = 1; j < NRad; j++){
			KwallsMV[i][j] = 0.50*(K[i][j] + K[i][j-1]); //Centro
		}
	}
}

//Cálculo del término difusivo de la ecuación de conservación de la energía
void Solver::Get_EnergyDifusive(Mesher MESH){
int i, j;
	
	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			EnergyDifusive[i][j] = (1.0/MESH.VolMP[i][j])*(
								   (KwallsMU[i][j]*(Tpres[i][j] - Tpres[i-1][j])*(MESH.SupMP[i][j][0])*cos(PI + MESH.AngleMU[i][j]))/MESH.DeltasMU[i][j][0]
								 + (KwallsMU[i+1][j]*(Tpres[i+1][j] - Tpres[i][j])*MESH.SupMP[i][j][1]*cos(MESH.AngleMU[i+1][j]))/MESH.DeltasMU[i+1][j][0] 
								 + (KwallsMV[i][j]*(Tpres[i][j] - Tpres[i][j-1])*MESH.SupMP[i][j][2]*sin(1.50*PI + MESH.AngleMR[i][j]))/MESH.DeltasMR[i][j][1] 
								 + (KwallsMV[i][j+1]*(Tpres[i][j+1] - Tpres[i][j])*MESH.SupMP[i][j][3]*sin(0.50*PI + MESH.AngleMR[i][j+1]))/MESH.DeltasMR[i][j+1][1]
								 );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		EnergyDifusive[0][j] = (1.0/MESH.VolMP[0][j])*(
							   KwallsMU[0][j]*(Tpres[0][j] - Tleft[j])*(MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*cos(PI + MESH.AngleMU[0][j]) 
							 + KwallsMU[1][j]*(Tpres[1][j] - Tpres[0][j])*(MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*cos(MESH.AngleMU[1][j]) 
							 + KwallsMV[0][j]*(Tpres[0][j] - Tpres[0][j-1])*(MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*sin(1.50*PI + MESH.AngleMR[0][j]) 
							 + KwallsMV[0][j+1]*(Tpres[0][j+1] - Tpres[0][j])*(MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*sin(0.50*PI + MESH.AngleMR[0][j+1])
							 );

		//Parte derecha
		EnergyDifusive[NA-1][j] = (1.0/MESH.VolMP[NA-1][j])*(
								  KwallsMU[NA-1][j]*(Tpres[NA-1][j] - Tpres[NA-2][j])*(MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*cos(PI + MESH.AngleMU[NA-1][j]) 
								+ KwallsMU[NA][j]*(Tright[j] - Tpres[NA-1][j])*(MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*cos(MESH.AngleMU[NA][j]) 
								+ KwallsMV[NA-1][j]*(Tpres[NA-1][j] - Tpres[NA-1][j-1])*(MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*sin(1.50*PI + MESH.AngleMR[NA-1][j]) 
								+ KwallsMV[NA-1][j+1]*(Tpres[NA-1][j+1] - Tpres[NA-1][j])*(MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])
								);
	}
	
	for(i = 1; i < NA-1; i++){
		//Parte abajo
		EnergyDifusive[i][0] = (1.0/MESH.VolMP[i][0])*(
							   KwallsMU[i][0]*(Tpres[i][0] - Tpres[i-1][0])*(MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*cos(PI + MESH.AngleMU[i][0]) 
							 + KwallsMU[i+1][0]*(Tpres[i+1][0] - Tpres[i][0])*(MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*cos(MESH.AngleMU[i+1][0]) 
							 + KwallsMV[i][0]*(Tpres[i][0] - Tdown[i])*(MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*sin(1.50*PI + MESH.AngleMR[i][0]) 
							 + KwallsMV[i][1]*(Tpres[i][1] - Tpres[i][0])*(MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*sin(0.50*PI + MESH.AngleMR[i][1])
							 );

		//Parte arriba
		EnergyDifusive[i][NRad-1] = (1.0/MESH.VolMP[i][NRad-1])*(
								    KwallsMU[i][NRad-1]*(Tpres[i][NRad-1] - Tpres[i-1][NRad-1])*(MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*cos(PI + MESH.AngleMU[i][NRad-1]) 
								  + KwallsMU[i+1][NRad-1]*(Tpres[i+1][NRad-1] - Tpres[i][NRad-1])*(MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*cos(MESH.AngleMU[i+1][NRad-1]) 
								  + KwallsMV[i][NRad-1]*(Tpres[i][NRad-1] - Tpres[i][NRad-2])*(MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[i][NRad-1]) 
								  + KwallsMV[i][NRad]*(Tup[i] - Tpres[i][NRad-1])*(MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*sin(0.50*PI + MESH.AngleMR[i][NRad])
								  );
	}
	
	//Esquina abajo izquierda
	EnergyDifusive[0][0] = (1.0/MESH.VolMP[0][0])*(
						   KwallsMU[0][0]*(Tpres[0][0] - Tleft[0])*(MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*cos(PI + MESH.AngleMU[0][0]) 
						 + KwallsMU[1][0]*(Tpres[1][0] - Tpres[0][0])*(MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*cos(MESH.AngleMU[1][0]) 
						 + KwallsMV[0][0]*(Tpres[0][0] - Tdown[0])*(MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*sin(1.50*PI + MESH.AngleMR[0][0]) 
						 + KwallsMV[0][1]*(Tpres[0][1] - Tpres[0][0])*(MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*sin(0.50*PI + MESH.AngleMR[0][1])
						 );

	//Esquina arriba izquierda
	EnergyDifusive[0][NRad-1] = (1.0/MESH.VolMP[0][NRad-1])*(
						  	    KwallsMU[0][NRad-1]*(Tpres[0][NRad-1] - Tleft[NRad-1])*(MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*cos(PI + MESH.AngleMU[0][NRad-1]) 
						 	  + KwallsMU[1][NRad-1]*(Tpres[1][NRad-1] - Tpres[0][NRad-1])*(MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*cos(MESH.AngleMU[1][NRad-1]) 
						 	  + KwallsMV[0][NRad-1]*(Tpres[0][NRad-1] - Tpres[0][NRad-2])*(MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]) 
						 	  + KwallsMV[0][NRad]*(Tup[0] - Tpres[0][NRad-1])*(MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*sin(0.50*PI + MESH.AngleMR[0][NRad])
						 	  );

	//Esquina abajo derecha
	EnergyDifusive[NA-1][0] = (1.0/MESH.VolMP[NA-1][0])*(
							  KwallsMU[NA-1][0]*(Tpres[NA-1][0] - Tpres[NA-2][0])*(MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*cos(PI + MESH.AngleMU[NA-1][0]) 
							+ KwallsMU[NA][0]*(Tright[0] - Tpres[NA-1][0])*(MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*cos(MESH.AngleMU[NA][0]) 
							+ KwallsMV[NA-1][0]*(Tpres[NA-1][0] - Tdown[NA-1])*(MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*sin(1.50*PI + MESH.AngleMR[NA-1][0]) 
							+ KwallsMV[NA-1][1]*(Tpres[NA-1][1] - Tpres[NA-1][0])*(MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1])
							);

	//Esquina arriba derecha
	EnergyDifusive[NA-1][NRad-1] = (1.0/MESH.VolMP[NA-1][NRad-1])*(
								   KwallsMU[NA-1][NRad-1]*(Tpres[NA-1][NRad-1] - Tpres[NA-2][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) 
								 + KwallsMU[NA][NRad-1]*(Tright[NRad-1] - Tpres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*cos(MESH.AngleMU[NA][NRad-1]) 
								 + KwallsMV[NA-1][NRad-1]*(Tpres[NA-1][NRad-1] - Tpres[NA-1][NRad-2])*(MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) 
								 + KwallsMV[NA-1][NRad]*(Tup[NA-1] - Tpres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])
								 );
}

void Solver::Get_eWalls(Mesher MESH){
int i, j;

		//Nodos U
		for(j = 0; j < NRad; j++){
			eWallsMU[0][j] = eleft[j];
			eWallsMU[NA][j] = eright[j];
			for(i = 1; i < NA; i++){
					eWallsMU[i][j] = sqrt(epres[i][j])*sqrt(epres[i-1][j]);
			}
		}

		//Nodos R
		for(i = 0; i < NA; i++){
			eWallsMR[i][0] = edown[i];
			eWallsMR[i][NRad] = eup[i];
			for(j = 1; j < NRad; j++){
				eWallsMR[i][j] = sqrt(epres[i][j])*sqrt(epres[i][j-1]);
			}
		}

}

//Cálculo del término convectivo de la energía
void Solver::Get_EnergyConvective(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			EnergyConvective[i][j] = -(1.0/MESH.VolMP[i][j])*(
									 MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*RhoWallsMU[i][j]*eWallsMU[i][j]
								   + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*RhoWallsMU[i+1][j]*eWallsMU[i+1][j]
								   + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*RhoWallsMR[i][j]*eWallsMR[i][j]
								   + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*RhoWallsMR[i][j+1]*eWallsMR[i][j+1]
								   );
		}
	}

}

//Cálculo del término de presión en la ecuación de conservación de la energía
void Solver::Get_EnergyPressureTerm(Mesher MESH){
int i, j;
		
		//Centro
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad; j++){
				EnergyPressureTerm[i][j] = -P[i][j]*(
										 + (1.0/MESH.DeltasMP[i][j][0])*(UwallsMU[i+1][j] - UwallsMU[i][j])
										 + (1.0/MESH.DeltasMP[i][j][1])*(VwallsMR[i][j+1] - VwallsMR[i][j])
										 );
			}
		}

	
}

//Cálculo del término viscoso de la ecuaciń de conservación de la energía
void Solver::Get_EnergyViscousTerm(Mesher MESH){
int i, j;


	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			EnergyViscous[i][j] = 
							    + 0.50*(TauRR_mr[i][j] + TauRR_mr[i][j+1])*(1.0/MESH.DeltasMP[i][j][1])*(VwallsMR[i][j+1] - VwallsMR[i][j])
							    + 0.50*(TauZZ_mu[i][j] + TauZZ_mu[i+1][j])*(1.0/MESH.DeltasMP[i][j][0])*(UwallsMU[i+1][j] - UwallsMU[i][j])
							    + 0.50*(TauRZ_mr[i][j] + TauRZ_mr[i][j+1])*((1.0/MESH.DeltasMP[i][j][0])*(VwallsMU[i+1][j] - VwallsMU[i][j]) + (1.0/MESH.DeltasMP[i][j][1])*(UwallsMR[i][j+1] - UwallsMR[i][j]))
							    ;
		}
	}
	
}

//Cálculo del término 1 de la ecuación de Spalart-Allmaras
void Solver::SA_Term1(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino1[i][j] = -(1.0/MESH.VolMP[i][j])*(
								  MESH.SupMP[i][j][0]*0.50*(vpresent[i][j] + vpresent[i-1][j])*0.50*(RhoPres[i][j] + RhoPres[i-1][j])*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))
								+ MESH.SupMP[i][j][1]*0.50*(vpresent[i][j] + vpresent[i+1][j])*0.50*(RhoPres[i][j] + RhoPres[i+1][j])*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))
								+ MESH.SupMP[i][j][2]*0.50*(vpresent[i][j] + vpresent[i][j-1])*0.50*(RhoPres[i][j] + RhoPres[i][j-1])*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
								+ MESH.SupMP[i][j][3]*0.50*(vpresent[i][j] + vpresent[i][j+1])*0.50*(RhoPres[i][j] + RhoPres[i][j+1])*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino1[0][j] = -(1.0/MESH.VolMP[0][j])*(
								MESH.SupMP[0][j][0]*vleft[j]*RhoLeft[j]*(UwallsMU[0][j]*cos(PI + MESH.AngleMU[0][j]) + VwallsMU[0][j]*sin(PI + MESH.AngleMU[0][j]))
								+ MESH.SupMP[0][j][1]*0.50*(vpresent[0][j] + vpresent[1][j])*0.50*(RhoPres[0][j] + RhoPres[1][j])*(UwallsMU[1][j]*cos(MESH.AngleMU[1][j]) + VwallsMU[1][j]*sin(MESH.AngleMU[1][j]))
								+ MESH.SupMP[0][j][2]*0.50*(vpresent[0][j] + vpresent[0][j-1])*0.50*(RhoPres[0][j] + RhoPres[0][j-1])*(UwallsMR[0][j]*cos(1.50*PI + MESH.AngleMR[0][j]) + VwallsMR[0][j]*sin(1.50*PI + MESH.AngleMR[0][j]))
								+ MESH.SupMP[0][j][3]*0.50*(vpresent[0][j] + vpresent[0][j+1])*0.50*(RhoPres[0][j] + RhoPres[0][j+1])*(UwallsMR[0][j+1]*cos(0.50*PI + MESH.AngleMR[0][j+1]) + VwallsMR[0][j+1]*sin(0.50*PI + MESH.AngleMR[0][j+1]))
							  );

		//Parte derecha
		SA_Termino1[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(
								MESH.SupMP[NA-1][j][0]*0.50*(vpresent[NA-1][j] + vpresent[NA-2][j])*0.50*(RhoPres[NA-1][j] + RhoPres[NA-2][j])*(UwallsMU[NA-1][j]*cos(PI + MESH.AngleMU[NA-1][j]) + VwallsMU[NA-1][j]*sin(PI + MESH.AngleMU[NA-1][j]))
								+ MESH.SupMP[NA-1][j][1]*vright[j]*RhoRight[j]*(UwallsMU[NA][j]*cos(MESH.AngleMU[NA][j]) + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][j]))
								+ MESH.SupMP[NA-1][j][2]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1])*0.50*(RhoPres[NA-1][j] + RhoPres[NA-1][j-1])*(UwallsMR[NA-1][j]*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + VwallsMR[NA-1][j]*sin(1.50*PI + MESH.AngleMR[NA-1][j]))
								+ MESH.SupMP[NA-1][j][3]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1])*0.50*(RhoPres[NA-1][j] + RhoPres[NA-1][j+1])*(UwallsMR[NA-1][j+1]*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + VwallsMR[NA-1][j+1]*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))
							  );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino1[i][0] = -(1.0/MESH.VolMP[i][0])*(
								MESH.SupMP[i][0][0]*0.50*(vpresent[i][0] + vpresent[i-1][0])*0.50*(RhoPres[i][0] + RhoPres[i-1][0])*(UwallsMU[i][0]*cos(PI + MESH.AngleMU[i][0]) + VwallsMU[i][0]*sin(PI + MESH.AngleMU[i][0]))
								+ MESH.SupMP[i][0][1]*0.50*(vpresent[i][0] + vpresent[i+1][0])*0.50*(RhoPres[i][0] + RhoPres[i+1][0])*(UwallsMU[i+1][0]*cos(MESH.AngleMU[i+1][0]) + VwallsMU[i+1][0]*sin(MESH.AngleMU[i+1][0]))
								+ MESH.SupMP[i][0][2]*vdown[i]*RhoDown[i]*(UwallsMR[i][0]*cos(1.50*PI + MESH.AngleMR[i][0]) + VwallsMR[i][0]*sin(1.50*PI + MESH.AngleMR[i][0]))
								+ MESH.SupMP[i][0][3]*0.50*(vpresent[i][0] + vpresent[i][1])*0.50*(RhoPres[i][0] + RhoPres[i][1])*(UwallsMR[i][1]*cos(0.50*PI + MESH.AngleMR[i][1]) + VwallsMR[i][1]*sin(0.50*PI + MESH.AngleMR[i][1]))
							  );

		//Parte arriba
		SA_Termino1[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(
								MESH.SupMP[i][NRad-1][0]*0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1])*0.50*(RhoPres[i][NRad-1] + RhoPres[i-1][NRad-1])*(UwallsMU[i][NRad-1]*cos(PI + MESH.AngleMU[i][NRad-1]) + VwallsMU[i][NRad-1]*sin(PI + MESH.AngleMU[i][NRad-1]))
								+ MESH.SupMP[i][NRad-1][1]*0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1])*0.50*(RhoPres[i][NRad-1] + RhoPres[i+1][NRad-1])*(UwallsMU[i+1][NRad-1]*cos(MESH.AngleMU[i+1][NRad-1]) + VwallsMU[i+1][NRad-1]*sin(MESH.AngleMU[i+1][NRad-1]))
								+ MESH.SupMP[i][NRad-1][2]*0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2])*0.50*(RhoPres[i][NRad-1] + RhoPres[i][NRad-2])*(UwallsMR[i][NRad-1]*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + VwallsMR[i][NRad-1]*sin(1.50*PI + MESH.AngleMR[i][NRad-1]))
								+ MESH.SupMP[i][NRad-1][3]*vup[i]*RhoUp[i]*(UwallsMR[i][NRad]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + VwallsMR[i][NRad]*sin(0.50*PI + MESH.AngleMR[i][NRad]))
							  );
	}

	//Esquina abajo izquierda
	SA_Termino1[0][0] = -(1.0/MESH.VolMP[0][0])*(
								MESH.SupMP[0][0][0]*vleft[0]*RhoLeft[0]*(UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0]) + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0]))
								+ MESH.SupMP[0][0][1]*0.50*(vpresent[0][0] + vpresent[1][0])*0.50*(RhoPres[0][0] + RhoPres[1][0])*(UwallsMU[1][0]*cos(MESH.AngleMU[1][0]) + VwallsMU[1][0]*sin(MESH.AngleMU[1][0]))
								+ MESH.SupMP[0][0][2]*vdown[0]*RhoDown[0]*(UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))
								+ MESH.SupMP[0][0][3]*0.50*(vpresent[0][0] + vpresent[0][1])*0.50*(RhoPres[0][0] + RhoPres[0][1])*(UwallsMR[0][1]*cos(0.50*PI + MESH.AngleMR[0][1]) + VwallsMR[0][1]*sin(0.50*PI + MESH.AngleMR[0][1]))
							  );

	//Esquina arriba izquierda
	SA_Termino1[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
								MESH.SupMP[0][NRad-1][0]*vleft[NRad-1]*RhoLeft[NRad-1]*(UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))
								+ MESH.SupMP[0][NRad-1][1]*0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1])*0.50*(RhoPres[0][NRad-1] + RhoPres[1][NRad-1])*(UwallsMU[1][NRad-1]*cos(MESH.AngleMU[1][NRad-1]) + VwallsMU[1][NRad-1]*sin(MESH.AngleMU[1][NRad-1]))
								+ MESH.SupMP[0][NRad-1][2]*0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2])*0.50*(RhoPres[0][NRad-1] + RhoPres[0][NRad-2])*(UwallsMR[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + VwallsMR[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))
								+ MESH.SupMP[0][NRad-1][3]*vup[0]*RhoUp[0]*(UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad]))
							  );

	//Esquina abajo derecha
	SA_Termino1[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								MESH.SupMP[NA-1][0][0]*0.50*(vpresent[NA-1][0] + vpresent[NA-2][0])*0.50*(RhoPres[NA-1][0] + RhoPres[NA-2][0])*(UwallsMU[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0]) + VwallsMU[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0]))
								+ MESH.SupMP[NA-1][0][1]*vright[0]*RhoRight[0]*(UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0]) + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0]))
								+ MESH.SupMP[NA-1][0][2]*vdown[NA-1]*RhoDown[NA-1]*(UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))
								+ MESH.SupMP[NA-1][0][3]*0.50*(vpresent[NA-1][0] + vpresent[NA-1][1])*0.50*(RhoPres[NA-1][0] + RhoPres[NA-1][1])*(UwallsMR[NA-1][1]*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + VwallsMR[NA-1][1]*sin(0.50*PI + MESH.AngleMR[NA-1][1]))
							  );

	//Esquina arriba derecha
	SA_Termino1[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
								MESH.SupMP[NA-1][NRad-1][0]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1])*0.50*(RhoPres[NA-1][NRad-1] + RhoPres[NA-2][NRad-1])*(UwallsMU[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + VwallsMU[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1]))
								+ MESH.SupMP[NA-1][NRad-1][1]*vright[NRad-1]*RhoRight[NRad-1]*(UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + VwallsMU[NA][NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))
								+ MESH.SupMP[NA-1][NRad-1][2]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2])*0.50*(RhoPres[NA-1][NRad-1] + RhoPres[NA-1][NRad-2])*(UwallsMR[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + VwallsMR[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))
								+ MESH.SupMP[NA-1][NRad-1][3]*vup[NA-1]*RhoUp[NA-1]*(UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))
							  );
}

//Cálculo del término 3 de la ecuación de Spalart-Allmaras
void Solver::SA_Term2(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			SA_Termino2[i][j] = Cb1*Smodel[i][j]*RhoPres[i][j]*vpresent[i][j];
		}
	}
}

//Cálculo del término 2 de la ecuación de Spalart-Allmaras
void Solver::SA_Term3(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino3[i][j] = (1.0/(Sigma*MESH.VolMP[i][j]))*(
								(MESH.SupMP[i][j][0]/MESH.DeltasMU[i][j][0])*muWallsMU[i][j]*(vpresent[i][j] - vpresent[i-1][j])*cos(PI + MESH.AngleMU[i][j])
							  + (MESH.SupMP[i][j][1]/MESH.DeltasMU[i+1][j][0])*muWallsMU[i+1][j]*(vpresent[i+1][j] - vpresent[i][j])*cos(MESH.AngleMU[i+1][j])
							  + (MESH.SupMP[i][j][2]/MESH.DeltasMR[i][j][1])*muWallsMR[i][j]*(vpresent[i][j] - vpresent[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j])
							  + (MESH.SupMP[i][j][3]/MESH.DeltasMR[i][j+1][1])*muWallsMR[i][j+1]*(vpresent[i][j+1] - vpresent[i][j])*sin(0.50*PI + MESH.AngleMR[i][j+1])	
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino3[0][j] = (1.0/(Sigma*MESH.VolMP[0][j]))*(
								(MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*muWallsMU[0][j]*(vpresent[0][j] - vleft[j])*cos(PI + MESH.AngleMU[0][j])
							  + (MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*muWallsMU[1][j]*(vpresent[1][j] - vpresent[0][j])*cos(MESH.AngleMU[1][j])
							  + (MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*muWallsMR[0][j]*(vpresent[0][j] - vpresent[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j])
							  + (MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*muWallsMR[0][j+1]*(vpresent[0][j+1] - vpresent[0][j])*sin(0.50*PI + MESH.AngleMR[0][j+1])	
							  );

		//Parte derecha
		SA_Termino3[NA-1][j] = (1.0/(Sigma*MESH.VolMP[NA-1][j]))*(
								(MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*muWallsMU[NA-1][j]*(vpresent[NA-1][j] - vpresent[NA-2][j])*cos(PI + MESH.AngleMU[i][j])
							  + (MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*muWallsMU[NA][j]*(vright[j] - vpresent[NA-1][j])*cos(MESH.AngleMU[NA][j])
							  + (MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*muWallsMR[NA-1][j]*(vpresent[NA-1][j] - vpresent[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])
							  + (MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*muWallsMR[NA-1][j+1]*(vpresent[NA-1][j+1] - vpresent[NA-1][j])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])	
							  );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino3[i][0] = (1.0/(Sigma*MESH.VolMP[i][0]))*(
								(MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*muWallsMU[i][0]*(vpresent[i][0] - vpresent[i-1][0])*cos(PI + MESH.AngleMU[i][0])
							  + (MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*muWallsMU[i+1][0]*(vpresent[i+1][0] - vpresent[i][0])*cos(MESH.AngleMU[i+1][0])
							  + (MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*muWallsMR[i][0]*(vpresent[i][0] - vdown[i])*sin(1.50*PI + MESH.AngleMR[i][0])
							  + (MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*muWallsMR[i][1]*(vpresent[i][1] - vpresent[i][0])*sin(0.50*PI + MESH.AngleMR[i][1])	
							  );

		//Parte arriba
		SA_Termino3[i][NRad-1] = (1.0/(Sigma*MESH.VolMP[i][NRad-1]))*(
								(MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*muWallsMU[i][NRad-1]*(vpresent[i][NRad-1] - vpresent[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1])
							  + (MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*muWallsMU[i+1][NRad-1]*(vpresent[i+1][NRad-1] - vpresent[i][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1])
							  + (MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*muWallsMR[i][NRad-1]*(vpresent[i][NRad-1] - vpresent[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])
							  + (MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*muWallsMR[i][NRad]*(vup[i] - vpresent[i][NRad-1])*sin(0.50*PI + MESH.AngleMR[i][NRad])	
							  );
	}

	//Esquina abajo izquierda
	SA_Termino3[0][0] = (1.0/(Sigma*MESH.VolMP[0][0]))*(
								(MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*muWallsMU[0][0]*(vpresent[0][0] - vleft[0])*cos(PI + MESH.AngleMU[0][0])
							  + (MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*muWallsMU[1][0]*(vpresent[1][0] - vpresent[0][0])*cos(MESH.AngleMU[1][0])
							  + (MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*muWallsMR[0][0]*(vpresent[0][0] - vdown[0])*sin(1.50*PI + MESH.AngleMR[0][0])
							  + (MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*muWallsMR[0][1]*(vpresent[0][1] - vpresent[0][0])*sin(0.50*PI + MESH.AngleMR[0][1])	
							  );

	//Esquina arriba izquierda
	SA_Termino3[0][NRad-1] = (1.0/(Sigma*MESH.VolMP[0][NRad-1]))*(
								(MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*muWallsMU[0][NRad-1]*(vpresent[0][NRad-1] - vleft[NRad-1])*cos(PI + MESH.AngleMU[0][NRad-1])
							  + (MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*muWallsMU[1][NRad-1]*(vpresent[1][NRad-1] - vpresent[0][NRad-1])*cos(MESH.AngleMU[1][NRad-1])
							  + (MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*muWallsMR[0][NRad-1]*(vpresent[0][NRad-1] - vpresent[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])
							  + (MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*muWallsMR[0][NRad]*(vup[0] - vpresent[0][NRad-1])*sin(0.50*PI + MESH.AngleMR[0][NRad])	
							  );

	//Esquina abajo derecha
	SA_Termino3[NA-1][0] = (1.0/(Sigma*MESH.VolMP[NA-1][0]))*(
								(MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*muWallsMU[NA-1][0]*(vpresent[NA-1][0] - vpresent[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0])
							  + (MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*muWallsMU[NA][0]*(vright[0] - vpresent[NA-1][0])*cos(MESH.AngleMU[NA][0])
							  + (MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*muWallsMR[NA-1][0]*(vpresent[NA-1][0] - vdown[NA-1])*sin(1.50*PI + MESH.AngleMR[NA-1][0])
							  + (MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*muWallsMR[NA-1][1]*(vpresent[NA-1][1] - vpresent[NA-1][0])*sin(0.50*PI + MESH.AngleMR[NA-1][1])	
							  );

	//Esquina arriba derecha
	SA_Termino3[NA-1][NRad-1] = (1.0/(Sigma*MESH.VolMP[NA-1][NRad-1]))*(
								(MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*muWallsMU[NA-1][NRad-1]*(vpresent[NA-1][NRad-1] - vpresent[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1])
							  + (MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*muWallsMU[NA][NRad-1]*(vright[NRad-1] - vpresent[NA-1][NRad-1])*cos(MESH.AngleMU[NA][NRad-1])
							  + (MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*muWallsMR[NA-1][NRad-1]*(vpresent[NA-1][NRad-1] - vpresent[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])
							  + (MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*muWallsMR[NA-1][NRad]*(vup[NA-1] - vpresent[NA-1][NRad-1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])	
							  );

}

//Cálculo del término 4 de la ecuación de Spalart-Allmaras
void Solver::SA_Term4(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino4[i][j] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[i][j] + vpresent[i+1][j]) - 0.50*(vpresent[i][j] + vpresent[i-1][j]))/MESH.DeltasMP[i][j][0],2.0)
							  + pow((0.50*(vpresent[i][j] + vpresent[i][j+1]) - 0.50*(vpresent[i][j] + vpresent[i][j-1]))/MESH.DeltasMP[i][j][1],2.0)		
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino4[0][j] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[0][j] + vpresent[1][j]) - vleft[j])/MESH.DeltasMP[0][j][0],2.0)
							  + pow((0.50*(vpresent[0][j] + vpresent[0][j+1]) - 0.50*(vpresent[0][j] + vpresent[0][j-1]))/MESH.DeltasMP[0][j][1],2.0)		
							  );

		//Parte derecha
		SA_Termino4[NA-1][j] = ((Cb2)/Sigma)*(
								pow((vright[j] - 0.50*(vpresent[NA-1][j] + vpresent[NA-2][j]))/MESH.DeltasMP[NA-1][j][0],2.0)
							  + pow((0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1]) - 0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1]))/MESH.DeltasMP[i][j][1],2.0)		
							  );

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino4[i][0] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[i][0] + vpresent[i+1][0]) - 0.50*(vpresent[i][0] + vpresent[i-1][0]))/MESH.DeltasMP[i][0][0],2.0)
							  + pow((0.50*(vpresent[i][0] + vpresent[i][1]) - vdown[i])/MESH.DeltasMP[i][0][1],2.0)		
							  );

		//Parte arriba
		SA_Termino4[i][NRad-1] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1]) - 0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1]))/MESH.DeltasMP[i][NRad-1][0],2.0)
							  + pow((vup[i] - 0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2]))/MESH.DeltasMP[i][NRad-1][1],2.0)		
							  );
	
	}

	//Esquina abajo izquierda
	SA_Termino4[0][0] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[0][0] + vpresent[1][0]) - vleft[0])/MESH.DeltasMP[0][0][0],2.0)
							  + pow((0.50*(vpresent[0][0] + vpresent[0][1]) - vdown[0])/MESH.DeltasMP[0][0][1],2.0)		
							  );

	//Esquina arriba izquierda
	SA_Termino4[0][NRad-1] = ((Cb2)/Sigma)*(
								pow((0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1]) - vleft[NRad-1])/MESH.DeltasMP[0][NRad-1][0],2.0)
							  + pow((vup[0] - 0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2]))/MESH.DeltasMP[0][NRad-1][1],2.0)		
							  );

	//Esquina abajo derecha
	SA_Termino4[NA-1][0] = ((Cb2)/Sigma)*(
								pow((vright[0] - 0.50*(vpresent[NA-1][0] + vpresent[NA-2][0]))/MESH.DeltasMP[NA-1][0][0],2.0)
							  + pow((0.50*(vpresent[NA-1][0] + vpresent[NA-1][1]) - vdown[NA-1])/MESH.DeltasMP[NA-1][0][1],2.0)		
							  );

	//Esquina arriba derecha
	SA_Termino4[NA-1][NRad-1] = ((Cb2)/Sigma)*(
								pow((vright[NRad-1] - 0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1]))/MESH.DeltasMP[NA-1][NRad-1][0],2.0)
							  + pow((vup[NA-1] - 0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2]))/MESH.DeltasMP[NA-1][NRad-1][1],2.0)		
							  );

}

//Cálculo del término 5 de la ecuación de Spalart-Allmaras
void Solver::SA_Term5(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			SA_Termino5[i][j] = -Cw1*RhoPres[i][j]*fw[i][j]*pow((vpresent[i][j])/MESH.minDist[i][j],2.0);
		}
	}

}

//Cálculo de todas las contribuciones del modelo Spalart-Allmaras
void Solver::FmuSA(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Fmupresent[i][j] = 
							 + SA_Termino1[i][j] 
							 + SA_Termino2[i][j] 
							 + SA_Termino3[i][j] 
							 + SA_Termino4[i][j] 
							 + SA_Termino5[i][j]					
							 ;
		}
	}

}

//Inicialización matrices del modelo Spalart Allmaras
void Solver::InitialSpalartAllmaras(){
int i,j;	
	
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad; j++){
				X[i][j] = vpresent[i][j]/(muBase[i][j]/RhoPres[i][j] + 1e-10);
				fv1[i][j] = pow(X[i][j],3.0)/(pow(X[i][j],3.0) + pow(Cv1,3.0) + 1e-10);

			}
		}
}

//Cálculo de variables y matrices necesarias para aplicar el modelo Spalart-Allmaras
void Solver::SpalartAllmarasPreparation(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){

			
			omega[i][j] = 0.50*((UwallsMR[i][j+1] - UwallsMR[i][j])/MESH.DeltasMP[i][j][1] - (VwallsMU[i+1][j] - VwallsMU[i][j])/MESH.DeltasMP[i][j][0]);
			S[i][j] = sqrt(2.0*omega[i][j]*omega[i][j]);
			fv2[i][j] = 1.0 - X[i][j]/(1.0 + X[i][j]*fv1[i][j]);
			Smodel[i][j] = S[i][j] + (vpresent[i][j]/(pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10))*fv2[i][j];

			
			//if(vpresent[i][j]/(Smodel[i][j]*pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10) < 10.0){
				r[i][j] = vpresent[i][j]/(Smodel[i][j]*pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10);
			//}
			//else{
			//	r[i][j] = 10.0;
			//}
			
			
			
		
			
			g[i][j] = r[i][j] + Cw2*(pow(r[i][j],6.0) - r[i][j]);
			
			fw[i][j] = g[i][j]*pow((1.0 + pow(Cw3,6.0))/(pow(g[i][j],6.0) + pow(Cw3,6.0) + 1e-10),1.0/6.0);
			
			
			ft2[i][j] = Ct3*exp(-Ct4*pow(X[i][j],2.0));
			
		}
	}

}

//Cálculo de la viscosidad con el modelo Spalart-Allmaras
void Solver::SpalartAllmaras(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){

			vfut[i][j] = (2.0*TimeBetta*vpresent[i][j] - (TimeBetta - 0.50)*vprevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*Fmupresent[i][j] - TimeBetta*Fmuprevious[i][j])/RhoPres[i][j];

			X[i][j] = vpresent[i][j]/(muBase[i][j]/RhoPres[i][j] + 1e-10);

			fv1[i][j] = pow(X[i][j],3.0)/(pow(X[i][j],3.0) + pow(Cv1,3.0) + 1e-10);
	
			muTurb[i][j] = RhoPres[i][j]*vfut[i][j]*fv1[i][j];
			
		}
	}

}

//Cálculo de la contribución total de momentum (Fu/Fv)
void Solver::Calculo_Fu_Fv(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			FuPresent[i][j] = MomentumDifusiveU[i][j]
							- PxGradient[i][j]
							+ MomentumConvectiveU[i][j]
							;
					
			FvPresent[i][j] = MomentumDifusiveV[i][j]
							- PyGradient[i][j]
							+ MomentumConvectiveV[i][j]
							;	
		}
	}
}

//Cálculo del mapa de densidades futuro
void Solver::Get_Densities(){
int i, j;
	
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			RhoFut[i][j] = (2.0*TimeBetta*RhoPres[i][j] - (TimeBetta - 0.50)*RhoPrevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*FRhoPresent[i][j] - TimeBetta*FRhoPrevious[i][j]);
		}	
	}		
}

//Cálculo del mapa de presiones con la Ley de los Gases Ideales
void Solver::Get_Pressure(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			P[i][j] = RhoPres[i][j]*Rideal*Tpres[i][j];
		}
	}

}

//Cálculo del mapa de velocidades futuro
void Solver::Get_Velocities(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 1; j < NRad; j++){

			//Cálculo de velocidades axiales (U)
			Ufut[i][j] = (2.0*TimeBetta*Upres[i][j] - (TimeBetta - 0.50)*Uprevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FuPresent[i][j]/RhoPres[i][j]) - TimeBetta*(FuPrevious[i][j]/RhoPrevious[i][j]));
			
			//Cálculo de velocidades radiales (V)
			Vfut[i][j] = (2.0*TimeBetta*Vpres[i][j] - (TimeBetta - 0.50)*Vprevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FvPresent[i][j]/RhoPres[i][j]) - TimeBetta*(FvPrevious[i][j]/RhoPrevious[i][j]));
			
		}
		Ufut[i][0] = (2.0*TimeBetta*Upres[i][0] - (TimeBetta - 0.50)*Uprevious[i][0])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*(FuPresent[i][0]/RhoPres[i][0]) - TimeBetta*(FuPrevious[i][0]/RhoPrevious[i][0]));
		Vfut[i][0] = 0.0;
	}

}

//Cálculo de la contribución total de energía (Fe)
void Solver::Get_Fe(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Fepresent[i][j] =
							+ EnergyConvective[i][j]
						    + EnergyDifusive[i][j]
						    + EnergyPressureTerm[i][j]
					        - EnergyViscous[i][j]
							;
		}
	}

}

//Cálculo del mapa de energías internas específicas
void Solver::Get_InternalEnergy(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			efut[i][j] = (2.0*TimeBetta*epres[i][j] - (TimeBetta - 0.50)*eprevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*Fepresent[i][j] - TimeBetta*Feprevious[i][j])/RhoPres[i][j];
		}

	}	
}

//Cálculo del mapa de entalpías específicas
void Solver::Get_SpecEnthalpy(){
int i,j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Hfut[i][j] = efut[i][j] + P[i][j]/RhoFut[i][j];
		}
	}
}

//Cálculo del mapa de temperaturas futuro
void Solver::Get_Temperatures(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Tfut[i][j] = Tpres[i][j] + (efut[i][j] - epres[i][j])/(Cp[i][j]/1.40 + 1e-10);
		}
	}
}

//Cálculo de la diferencia de resultados entre Steps
double Solver::Get_Stop(){
int i, j;
double MaximaDiferencia = 0.0;

	
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad; j++){
				//Mapa de densidades
				if(abs(RhoFut[i][j] - RhoPres[i][j])/(RhoPres[i][j]+1e-10) >= MaximaDiferencia){ MaximaDiferencia = abs(RhoFut[i][j] - RhoPres[i][j])/(RhoPres[i][j]+1e-10); }
				
				//Mapa de temperaturas
				if(abs(Tfut[i][j] - Tpres[i][j])/(Tpres[i][j] + 1e-10) >= MaximaDiferencia){ MaximaDiferencia = abs(Tfut[i][j] - Tpres[i][j])/(Tpres[i][j] + 1e-10); }	
	
				//Mapa de velocidades axiales (U)
				if(abs(Ufut[i][j] - Upres[i][j])/(abs(Upres[i][j])+1e-10) >= MaximaDiferencia){ MaximaDiferencia = abs(Ufut[i][j] - Upres[i][j])/(abs(Upres[i][j])+1e-10); }

				//Mapa de velocidades radiales (V)
				if(abs(Vfut[i][j] - Vpres[i][j])/(abs(Vpres[i][j])+1e-10) >= MaximaDiferencia){ MaximaDiferencia = abs(Vfut[i][j] - Vpres[i][j])/(abs(Vpres[i][j])+1e-10); }
			}
		}

	
	return MaximaDiferencia;

}

//Update de las F
void Solver::UpdateF(){
int i, j;

	for(i = NA-1; i >= 0; i--){
		for(j = NRad-1; j >= 0; j--){
			//Densidad
			FRhoPrevious[i][j] = FRhoPresent[i][j];

			//Velocidad axial (U)
			FuPrevious[i][j] = FuPresent[i][j];

			//Velocidad radial (V)
			FvPrevious[i][j] = FvPresent[i][j];

			//Energía Interna
			Feprevious[i][j] = Fepresent[i][j];

			//Viscosidades turbulentas
			Fmuprevious[i][j] = Fmupresent[i][j];
		}
	}

}

void Solver::UpdatePropertiesFields(){
int i, j;
	
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			//Densidades
			RhoPrevious[i][j] = RhoPres[i][j];
			RhoPres[i][j] = RhoFut[i][j];

			//Velocidades axiales (U)
			Uprevious[i][j] = Upres[i][j];
			Upres[i][j] = Ufut[i][j];

			//Velocidades radiales (V)
			Vprevious[i][j] = Vpres[i][j];
			Vpres[i][j] = Vfut[i][j];

			//Presiones
			Pprevious[i][j] = P[i][j];

			//Temperaturas
			Tpres[i][j] = Tfut[i][j];

			//Energías Internas
			eprevious[i][j] = epres[i][j];
			epres[i][j] = efut[i][j];

			//Entalpías
			//Hprevious[i][j] = Hpres[i][j];
			//Hpres[i][j] = Hfut[i][j];

			//Viscosidades turbulentas
			vprevious[i][j] = vpresent[i][j];
			vpresent[i][j] = vfut[i][j];

		}
	}
	
	
	

}

//Sumar las contribuciones de todas las viscosidades
void Solver::Get_TotalMU(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			muTotal[i][j] = muBase[i][j] + muTurb[i][j];
		}
	}



}
//Pasar a un .vtk los resultados de campos vectoriales
void Solver::VectorialVTK(Mesher MESH, double **MatrizX, double **MatrizY, string Documento, string NombrePropiedad){
int i, j;

        ofstream file;
        stringstream name_t;
        string name_f;
	
	name_t<<"/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/Paraview2/"<<Problema<<"_"<<Documento<<".vtk";
	
	name_f=name_t.str();
        file.open(name_f.c_str());

        file<<"# vtk DataFile Version 2.0"<<endl;
        file<<NombrePropiedad<<endl;
        file<<"ASCII"<<endl;
        file<<endl;
        file<<"DATASET STRUCTURED_GRID"<<endl;
        file<<"DIMENSIONS"<<"   "<<NA<<"   "<<NRad<<"   "<<1<<endl;
        file<<endl;
        file<<"POINTS"<<"   "<<NA*NRad<<"   "<<"double"<<endl;
	
	for(j = 0; j < NRad; j++){
		for(i = 0; i < NA; i++){
			file<<MESH.MP[i][j][0]<<"   "<<MESH.MP[i][j][1]<<"   "<<0.0<<endl;
			
		}
		
	}
       
        file<<endl;
        file<<"POINT_DATA"<<"   "<<NA*NRad<<endl;
        file<<"VECTORS "<<NombrePropiedad<<" double"<<endl;
        file<<endl;

       for(j = 0; j < NRad; j++){
                for(i = 0; i < NA; i++){
		        file<<MatrizX[i][j]/VelocidadSonido<<"   "<<MatrizY[i][j]/VelocidadSonido<<"   "<<"0.0"<<endl;
                }
        }
        file.close();
}

void Solver::ComputeParameters(Mesher MESH, int Step){
int i, j;
double ReynoldsTau, MachTau;
int LimInf, LimSup;
double SumVel = 0.0;
double MeanVel, MeanRho;
double ReynoldsMean;

double RhoM = 0.0;

RhoW = 0.0;
TauW = 0.0;

LimInf = 0.80*NA;
LimSup = 0.95*NA;
	
	TiempoEstadistico += DeltaT;

	for(j = 0; j < NRad; j++){
		DensidadMedia[j] = 0.0;
		TemperaturaMedia[j] = 0.0;
		PresionMedia[j] = 0.0;
	}


	

	for(j = 0; j < NRad; j++){
		for(i = LimInf; i <= LimSup; i++){
			DensidadMedia[j] = DensidadMedia[j] + (RhoPres[i][j]*MESH.DeltasMP[i][j][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
			TemperaturaMedia[j] = TemperaturaMedia[j] + (Tpres[i][j]*MESH.DeltasMP[i][j][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
			PresionMedia[j] = PresionMedia[j] + (P[i][j]*MESH.DeltasMP[i][j][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
		}
	}

	for(j = 0; j < NRad; j++){
			DensidadMediaTemporal[j] = DensidadMediaTemporal[j] + (DensidadMedia[j]*DeltaT);
			TemperaturaMediaTemporal[j] = TemperaturaMediaTemporal[j] + (TemperaturaMedia[j]*DeltaT);
			PresionMediaTemporal[j] = PresionMediaTemporal[j] + (PresionMedia[j]*DeltaT);
	}

	
	if(Step%1000 == 0){

	for(i = LimInf; i <= LimSup; i++){
		RhoW = RhoW + (RhoPres[i][NRad-1]*MESH.DeltasMP[i][NRad-1][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
		TauW = TauW + (muTotal[i][NRad-1]*((Upres[i][NRad-1]/MESH.DeltasMR[i][NRad][1])*MESH.DeltasMP[i][NRad-1][0]))/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);

	}

	for(j = 0; j < NRad; j++){
			SumVel = SumVel + (Upres[LimSup+1][j]*MESH.DeltasMP[LimSup+1][j][1])/VelocidadSonido;
	}

	for(j = 0; j < NRad; j++){
		RhoM = RhoM + (RhoPres[LimSup+1][j]*MESH.DeltasMP[LimSup+1][j][1]);
	}
	MeanRho = RhoM/(0.50*PipeDiameter);
	MeanVel = SumVel/(PipeDiameter/2.0);

	ReynoldsMean = (MeanRho*MeanVel*VelocidadSonido*0.50*PipeDiameter)/muW;
	uT = sqrt(abs(TauW)/RhoW);
	
	ReynoldsTau = (RhoW*uT*PipeDiameter)/(2.0*muW);
	MachTau = uT/sqrt(1.40*Rideal*Twall);


		for(j = 0; j < NRad; j++){
			DensidadMediaTemporalDato[j] = DensidadMediaTemporal[j]/TiempoEstadistico;
			TemperaturaMediaTemporalDato[j] = TemperaturaMediaTemporal[j]/TiempoEstadistico;
			PresionMediaTemporalDato[j] = PresionMediaTemporal[j]/TiempoEstadistico;
	}

	FILE *fp1;
		fp1 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/ResultadosParametros.txt","w");
			fprintf(fp1,"ReynoldsMean: %f \n", ReynoldsMean);
			fprintf(fp1,"ReynoldsTau: %f \n", ReynoldsTau);
			fprintf(fp1,"MachTau: %f \n", MachTau);
			fprintf(fp1,"Limite Inferior: %d \n", LimInf);
			fprintf(fp1,"Limite Superior: %d \n", LimSup);
			fprintf(fp1,"RhoW: %f \n", RhoW);
			fprintf(fp1,"Mean cross velocity: %f \n", MeanVel);
			fprintf(fp1,"TauW: %f \n", TauW);
			fprintf(fp1,"uT: %f \n", uT);
			fprintf(fp1,"TiempoEstadistico: %f \n", TiempoEstadistico);
					
	fclose(fp1);

	FILE *fp2;
		fp2 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/MeanDensity.txt","w");
			for(j = 0; j < NRad; j++){
				fprintf(fp2,"%f \t %f \n", MESH.MP[LimInf][NRad-1-j][1]/(PipeDiameter/2.0), DensidadMediaTemporalDato[j]/RhoW);
			}
	fclose(fp2);

	FILE *fp3;
		fp3 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/MeanTemperature.txt","w");
			for(j = 0; j < NRad; j++){
				fprintf(fp3,"%f \t %f \n", MESH.MP[LimInf][NRad-1-j][1]/(PipeDiameter/2.0), TemperaturaMediaTemporalDato[j]/Twall);
			}
	fclose(fp3);

	FILE *fp4;
		fp4 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/MeanPressure.txt","w");
			for(j = 0; j < NRad; j++){
				fprintf(fp4,"%f \t %f \n", MESH.MP[LimInf][NRad-1-j][1]/(PipeDiameter/2.0), PresionMediaTemporalDato[j]/PresionMediaTemporalDato[NRad-1]);
			}
	fclose(fp4);
	}
}


//Ejecución de todos los procesos del solver
void Solver::ExecuteSolver(Memory M1, ReadData R1, Mesher MESH){
int i, j;
double Time = 0.0;
double MaxDif = 2.0*Convergencia;
int Step = 0;
char Nombre[300];
int Alerta = 0;


	//Alojamiento de la memoria de todas las matrices utilizadas
	AllocateMatrix(M1);
	
	//Cálculos termoquímicos previos a la simulación
	Get_Rideal(); //Cálculo de la constante ideal de la mezcla de gases
	Get_CpAir();
	InitializeFields(MESH); //Pasar todos los coeficientes termoqímicos de las especies a sus matrices
	InitialF(); //Seteo inicial de las F
	InitialSpalartAllmaras();
	Get_Pressure(); //Cálculo del mapa de presiones con la Ley de los Gases Ideales

	sprintf(Nombre, "MapaPresiones_Step_%d", Step);
	MESH.EscalarVTK2D("Presiones", Nombre, P);

	sprintf(Nombre, "MapaTemperaturas_Step_%d", Step);
	MESH.EscalarVTK2D("Temperaturas", Nombre, Tpres);

	sprintf(Nombre, "MapaVelocidades_Step_%d", Step);
	VectorialVTK(MESH, Upres, Vpres, Nombre, "Velocidades");

	sprintf(Nombre, "MapaDensidades_Step_%d", Step);
	MESH.EscalarVTK2D("Densidades", Nombre, RhoPres);

	while(MaxDif >= Convergencia){
		Step += 1;
		UpdateBoundaryConditions(MESH);
		
		Get_AirViscosity();
	//	Get_nuWalls(MESH);
		
		Get_VelocityWalls(MESH); //Cálculo de las velocidades en las paredes de los volúmenes de control
	
		//Cálculo del modelo Spalart-Allmaras
	/*	SpalartAllmarasPreparation(MESH);
		SA_Term1(MESH);
		SA_Term2();
		SA_Term3(MESH);
		SA_Term4(MESH);
		SA_Term5(MESH);
		FmuSA();
		SpalartAllmaras();*/

		Get_TotalMU();
		Get_muWalls(MESH); //Cálculo de las conductividades térmicas en las paredes de los volúmenes de control 
		Get_VelocityGradients(MESH);
		Get_Divergence(MESH);
		Get_Stresses(MESH);
		
		DeltaT = StepTime(MESH);

		Time += DeltaT;

		
		//Cálculo de la densidad	
		Get_RhoWalls(MESH);
		Get_FRho(MESH); //Cálculo de las contribuciones de la ecuación de continuidad
		Get_Densities(); //Cálculo del mapa de densidades futuro
	
		//Cálculo de la presión
		Get_Pressure(); //Cálculo del mapa de presiones con la Ley de los Gases Ideales
		Get_Pwalls(MESH);
		Get_Pgradients(MESH); //Cálculo de los gradientes de presión en cada nodo en ambas direciones (U y V)
		
		
		
		//Cálculo de las velocidades
		Get_MomemtumDifusiveU(MESH); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Axial U)
		//Get_MomemtumDifusiveV(MESH); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
		Get_MomentumConvectiveU(MESH); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
		Get_MomentumConvectiveV(MESH); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)
		Calculo_Fu_Fv(); //Cálculo de la contribución total de momentum (Fu/Fv)
		Get_Velocities(); //Cálculo del mapa de velocidades futuro

		//Cálculo de la tempratura
		Get_K(MESH);
		Get_EnergyDifusive(MESH);
		Get_eWalls(MESH);
		Get_EnergyConvective(MESH);
		Get_EnergyViscousTerm(MESH);
		Get_EnergyPressureTerm(MESH);
		Get_Fe();
		Get_InternalEnergy();
		Get_Temperatures();

		if(Step > 5000){ComputeParameters(MESH, Step); }
		
		//Estudio de la convergencia
		if(Step%1000 == 0){
			MaxDif = Get_Stop();
			cout<<"Step: "<<Step<<", DeltaT: "<<DeltaT<<", Tiempo: "<<Time<<", MaxDif: "<<MaxDif<<endl;

			sprintf(Nombre, "MapaPresiones_Step_%d", Step);
			MESH.EscalarVTK2D("Presiones", Nombre, P);

			sprintf(Nombre, "MapaTurbulencia_Step_%d", Step);
			MESH.EscalarVTK2D("ViscosidadTurbulenta", Nombre, muTurb);

			sprintf(Nombre, "MapaDensidades_Step_%d", Step);
			MESH.EscalarVTK2D("Densidades", Nombre, RhoFut);

			sprintf(Nombre, "MapaTemperaturas_Step_%d", Step);
			MESH.EscalarVTK2D("Temperaturas", Nombre, Tfut);

			sprintf(Nombre, "MapaVelocidades_Step_%d", Step);
			VectorialVTK(MESH, Ufut, Vfut, Nombre, "Velocidades");

			
		}

		UpdateF();
		UpdatePropertiesFields();

		
		
	}
	


}

