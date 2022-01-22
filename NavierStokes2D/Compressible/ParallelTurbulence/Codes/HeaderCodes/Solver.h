#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <stdio.h>
#include <mpi.h>

using namespace std;

class Solver{
	private:
		//Densidad del mallado
		int NR; //Dirección radial
		int NA; //Dirección axial

		string TypeMesh; //Tipo de mallado seleccionado (Collocated/Staggered)
		string EsquemaAmplio;
		string Problema;

		//Datos y variables de computación paralela
		int Rank;
		int Procesos;
		int Ix;
		int Fx;
		int Halo;

		//Datos Geométricos del problema
		double PipeLength;   //Longitud de la tubería
		double PipeDiameter; //Diametro de la tubería

		double H;

		//Condiciones físicas del problema
		double BulkRe;
		double BulkM;
		double Pr;
		double Twall;
		double Gamma;
		double Rideal;
		double CpBase;
		
		double SoundSpeed;
		double uM;

		double B;
		double n;
		double muW;
		double RhoM;
		double S;
		
		double TimeBetta; //Esquema de integración temporal
		double DeltaT;
		double Convergencia;
		double MaxDifference;

		//Variable de interpolación
		double Valor;

		//Variables de la función del esquema convectivo utilizado
		double PHI, XD, PhiD, XC, PhiC, XU, PhiU;
		
		//Variables de la función del esquema convectivo utilizado de adimensionalizacion
		double PhC, Xadc, Xade;
		
		//Variables de la función del esquema convectivo utilizado en evaluacion
		double Phf;
	
		double uT;
		double RhoW;
		double TauW;


	public:
		Solver(Memory, ReadData, ParPro, Mesher);

		//Matrices globales de propiedades
		double *UglobalFut;
		double *VglobalFut;
		double *TglobalFut;
		double *RhoGlobalFut;
		double *PresionGlobal;

		double *UglobalPres;
		double *VglobalPres;
		double *TglobalPres;
		double *RhoGlobalPres;

		//Matrices de propiedades locales
		//Densidad
		double *RhoLocalPrev;
		double *RhoLocalPres;
		double *RhoLocalFut;

		//Presión
		double *Pressure;

		//Velocidad Axial
		double *UlocalPrev;
		double *UlocalPres;
		double *UlocalFut;

		//Velocidad Radial
		double *VlocalPrev;
		double *VlocalPres;
		double *VlocalFut;

		//Energía Interna
		double *ElocalPrev;
		double *ElocalPres;
		double *ElocalFut;

		//Temperatura
		double *TlocalPres;
		double *TlocalFut;

		//Matrices de propiedades en las paredes de los volúmenes de control
		//Densidad
		double *RhoWallsMU;
		double *RhoWallsMR;

		//Presión
		double *PwallsMU;
		double *PwallsMR;

		//Velocidad Axial
		double *UwallsMU;
		double *UwallsMR;
		
		double *UwallsMUglobal;

		//Velocidad Radial
		double *VwallsMU;
		double *VwallsMR;

		//Conductividad Térmica
		double *KwallsMU;
		double *KwallsMR;

		//Energía Interna
		double *EwallsMU;
		double *EwallsMR;

		//Viscosidad Dinámica
		double *muTotalMU;
		double *muTotalMR;

		//Matrices de cálculo de contribuciones a las propiedades
		//Densidad
		double *FRhoPrev;
		double *FRhoPres;

		//Velocidad Axial
		double *FUPrev;
		double *FUPres;

		//Velocidad Radial
		double *FVPrev;
		double *FVPres;

		//Energía Interna
		double *FEpres;
		double *FEprev;

		//Términos de las ecuaciones de Navier-Stokes
		//Momentum Axial
		double *MomentumConvectiveAxial;
		double *MomentumDiffusiveAxial;
		double *PressureGradientAxial;

		//Momentum Radial
		double *MomentumConvectiveRadial;
		double *MomentumDiffusiveRadial;
		double *PressureGradientRadial;

		//Energía
		double *EnergyConvective;	//Término convectivo de la ecuación de energía
		double *EnergyDiffusive;	//Término difusivo de la ecuación de conservación de la energía
		double *EnergyPressureTerm; //Término de presión de la ecuación de la energía	
		double *EnergyViscous;		//Término de disipación viscosa en la ecuación de energía

		//Matrices de Esfuerzos Viscosos
		double *TauRR;
		double *TauRZ;
		double *TauZZ;
		double *TauThetaTheta;

		double *TauRZglobal;

		double *TauRR_mu;
		double *TauRR_mr;

		double *TauRZ_mu;
		double *TauRZ_mr;
		
		double *TauZZ_mu;
		double *TauZZ_mr;
		
		//Matrices locales de propiedades físicas del fluido
		//Propiedades térmicas
		double *Cp;		//Calor específico de los gases en cada nodo
		double *K;		//Conductividad térmica en los nodos

		//Viscosidad
		double *muBase;
		double *muTurb;
		double *muTotal; //Viscosidad dinámica total en cada uno de los nodos

		//Matrices adicionales necesarias
		double *Divergence; //Matriz de divergencias de la velocidad en cada nodo 
		double *DivergenceMU;
		double *DivergenceMR;
		
		double *PDT; //Possible Delta T

		//Matrices y arrays de condiciones de contorno
		//Condiciones de contorno del mapa de densidades
		double *RhoLeft; //Densidad pared izquierda
		double *RhoRight; //Densidad pared derecha
	
		double *RhoDown; //Densidad pared abajo
		double *RhoUp; //Densidad pared superior

		//Condiciones de contorno del mapa de velocidades axiales (U)
		double *Uleft; //Velocidad axial pared izquierda
		double *Uright; //Velocidad axial pared derecha
	
		double *Udown;
		double *Uup; //Velocidad axial pared superior
 
		//Condiciones de contorno del mapa de velocidades radiales (V)
		double *Vleft; //Velocidad radial pared izquierda
		double *Vright; //Velocidad radial pared derecha
	
		double *Vdown;
		double *Vup; //Velocidad radial pared superior
	
		//Condiciones de contorno del mapa de energías internas
		double *Eleft;	//Energía Interna pared izquierda
		double *Eright; //Energía Interna pared derecha

		double *Edown;
		double *Eup;   //Energía Interna pared superior

		//Condiciones de contorno del mapa de temperaturas (T)
		double *Tleft; //Temperaturas pared izquierda
		double *Tright; //Temperaturas pared derecha

		double *Tdown;
		double *Tup; //Temperaturas pared superior
		
		//Condiciones de contorno del mapa de viscosidades dinámicas
		double *Pleft; //Presión pared izquierda
		double *Pright; //Presión pared derecha

		double *Pdown;
		double *Pup; //Presión pared superior
		
		//Condiciones de contorno del mapa de viscosidades dinámicas
		double *muLeft; //Viscosidad dinámica pared izquierda
		double *muRight; //Viscosidad dinámica pared derecha

		double *muDown;
		double *muUp; //Viscosidad dinámica pared superior
		
		//Condiciones de contorno del mapa de condutividades térmicas
		double *KLeft; //Conductividad térmica pared izquierda
		double *KRight; //Conductividad térmica pared derecha

		double *Kdown;
		double *KUp; //Conductividad térmica pared superior
	
		//Esquema convectivo

		//Función de interpolación lineal
		inline double Interpolacion(double C1, double V1, double C2, double V2, double CPunto){
			return V1 + ((CPunto - C1)/(C2 - C1))*(V2 - V1); }
	
		//Función de aplicación del esquema convectivo
		double EsquemaConvectivo(double, double, double, double, double, double, double, double, double, double, string);

		void AllocateMatrix(Memory);
		
		void InitializeFields(Mesher); //Pasar todos los coeficientes termoqímicos de las especies a sus matrices
		void InitialF(); //Seteo inicial de las Contribuciones de las ecuaciones

		void UpdateBoundaryConditions(Mesher);

		void Get_FluidViscosity();
		void Get_TotalViscosity();
		void Get_StepTime(Mesher, ParPro);

		void Get_VelocityWalls(Mesher); //Cálculo de las velocidades en las paredes de los volúmenes de control
		void Get_DensityWalls(Mesher);

		void Get_DensityConvective(Mesher);
		void Get_Densities(); //Cálculo del mapa de densidades futuro

		void Get_Pressure(); 	//Cálculo del mapa de presiones con la Ley de los Gases Ideales	
		void Get_PressureWalls(Mesher); //Cálculo de las presiones en las paredes de los volúmenes de control

		void Get_MomentumConvectiveAxial(Mesher); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
		void Get_MomentumConvectiveRadial(Mesher);	  //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)

		void Get_Stresses(Mesher, ParPro);
		void Get_MomemtumDifusiveAxial(Mesher); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Axial U)
		void Get_MomemtumDifusiveRadial(Mesher); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)

		void Get_PressureGradients(Mesher);	//Cálculo de los gradientes de presión en cada nodo en ambas direciones (U y V)
				
		void Get_Velocities(); //Cálculo del mapa de velocidades futuro

		void Get_K(Mesher, ParPro); //Cálculo de la conductividad térmica de los gases en cada nodo y paredes
		void Get_EnergyDiffusive(Mesher); //Cálculo del término difusivo de la ecuación de conservación de la energía

		void Get_EnergyWalls(Mesher);
		void Get_EnergyConvective(Mesher); //Cálculo del término convectivo de la ecuación de energía

		void Get_EnergyPressureTerm(Mesher); //Cálculo del término de presión en la ecuación de la energía
		void Get_EnergyViscousTerm(Mesher); //Cálculo del término de disipación viscosa en la ecuación de energía

		void Get_Temperatures(); //Cálculo del mapa de temperaturas futuro

		void Get_Stop(ParPro);

		void EscalarVTK2D(Mesher, string, string, string, double*, double*, int, int); //Pasar los resultados a un archivo VTK en 2D
		void VectorialVTK2D(Mesher, string, string, string, double*, double*, double*, int, int); //Pasar a un .vtk los resultados de campos vectoriales
	    
		void UpdateParameters();
		
		void ExecuteSolver(Memory, ParPro, Mesher);
		
};
