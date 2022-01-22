#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
using namespace std;

class Solver{
	private:
		//Densidad del mallado
		int NRad; //Dirección radial
		int NA; //Dirección axial

		string TypeMesh; //Tipo de mallado seleccionado (Collocated/Staggered)
		string EsquemaAmplio;
		string Problema;

		double TimeBetta; //Esquema de integración temporal
		double DeltaT;
		double Convergencia; 
 		double PipeDiameter;

 		double Tref;
 		double muRef;
 		double Sref;

		double Pr; //Número de Prandtl
		double ReynoldsM;
		double MachM;
		double VelocidadSonido;
		double uM;

		double B;
		double n;
		double muW;
		double RhoM;
		double mu; //Viscosidad dinámica

		double Cmu;
		double k_modelo;
		double I;
		double Epsilon;
		double vt;

		//Constantes del modelo Spalart-Allmaras
		double Cb1; //Cb1
		double Cb2; //Cb2

		double Sigma; //Sigma
		double K_VK; //Constante de Von-Kárman
		
		double Cw1; //Cw1
		double Cw2; //Cw2
		double Cw3; //Cw3

		double Cv1; //Cv1

		double Ct1; //Ct1
		double Ct2; //Ct2
		double Ct3; //Ct3
		double Ct4; //Ct4
		
		//Variable de interpolación
		double Valor;

		//Variables de la función del esquema convectivo utilizado
		double PHI, XD, PhiD, XC, PhiC, XU, PhiU;
		
		//Variables de la función del esquema convectivo utilizado de adimensionalizacion
		double PhC, Xadc, Xade;
		
		//Variables de la función del esquema convectivo utilizado en evaluacion
		double Phf;
		
		//Variables relacionadas con las propiedades termoquímicas de los gases
		double MmH2; //Masa molar del H2
		double MmO2; //Masa molar del O2
		double MmH2O; //Masa molar del H2O
		
		double FracH2; //Fracción de H2 en la mezcla de gases
		double FracO2; //Fracción de O2 en la mezcla de gases
		double FracH2O; //Fracción de H2O en la mezcla de gases

		double RO2; //Constante de los gases ideales (R) del O2
		double RH2; //Constante de los gases ideales (R) del H2
		double RH2O; //Constante de los gases ideales (R) del H2O
		double Rideal; //Constante de los gases ideales (R) de la mezcla de gases

		double HfO2; //Entalpía de formación del O2 a 298.15 K 
		double HfH2; //Entalpía de formación del H2 a 298.15 K 
		double HfH2O; //Entalpía de formación del H2O a 298.15 K 

		//Variables relacionadas con las viscosidad de las especies termoquímicas y la ley de Sutherland
		double mu0O2; //Viscosidad dinámica base del O2
		double TsO2; //Ts (Constante) del O2 (K)
		double ToO2; //Temperatura To del O2 (K)

		double mu0H2; //Viscosidad dinámica base del H2
		double TsH2; //Ts (Constante) del H2 (K)
		double ToH2; //Temperatura To del H2 (K)
 
		double mu0H2O; //Viscosidad dinámica base del H2O
		double TsH2O; //Ts (Constante) del H2O (K)
		double ToH2O; //Temperatura To del H2O (K)

		double uT;
		double RhoW;
		double TauW;

		double Twall;
		double Ratio;

	public:
		Solver(Memory, ReadData, Mesher);

		//Matrices de propiedades termoquímicas de las especies
		double **JanafH2; //Coeficientes función Cp(T) y H(T) del H2
		double **JanafO2; //Coeficientes función Cp(T) y H(T) del O2
		double **JanafH2O; //Coeficientes función Cp(T) y H(T) del H2O

		double **Cp; //Calor específico de los gases en cada nodo

		//Matrices del modelo Spalart-Allmaras
		double **muTotal; //Viscosidad dinámica total en cada uno de los nodos
		double **muWallsMU; //Viscosidad dinámica en las paredes de los volúmenes de control (Nodos U)
		double **muWallsMR; //Viscosidad dinámica en las paredes de los volúmenes de control (Nodos R)

		double **nuWallsMU; //Viscosidad cinemática en las paredes de los volúmenes de control (Nodos U)
		double **nuWallsMR; //Viscosidad cinemática en las paredes de los volúmenes de control (Nodos R)

		double **UwallsMU; //Velocidad axial en las paredes de los volúmenes de control (Nodos U)
		double **UwallsMR; //Velocidad axial en las paredes de los volúmenes de control (Nodos R)

		double **VwallsMU; //Velocidad radial en las paredes de los volúmenes de control (Nodos U)
		double **VwallsMR; //Velocidad radial en las paredes de los volúmenes de control (Nodos R)

		double **WwallsMU; //Velocidad tangencial en las paredes de los volúmenes de control (Nodos U)
		double **WwallsMR; //Velocidad tangencial en las paredes de los volúmenes de control (Nodos R)

		double **RhoWallsMU; //Densidad en las paredes de los volúmenes de control (Nodos U)
		double **RhoWallsMR; //Densidad en las paredes de los volúmenes de control (Nodos R)

		double **eWallsMU; //Energia en las paredes de los volúmenes de control (Nodos U)
		double **eWallsMR; //Energia en las paredes de los volúmenes de control (Nodos R)

		//Matrices de los mapas de propiedades del problema
	
		//Densidad
		double **RhoPrevious; //Matriz del mapa de densidades previo
		double **RhoPres; //Matriz del mapa de densidades presente
		double **RhoFut; //Matriz del mapa de densidades futuro

		//Presión
		double **Pprevious; //Matriz del mapa de presiones previo
		double **P; //Matriz del mapa de presiones presente

		//Velocidad axial (U)
		double **Uprevious; //Matriz del mapa de velocidades axiales (U) previo
		double **Upres; //Matriz del mapa de velocidades axiales (U) presente
		double **Ufut; //Matriz del mapa de velocidades axiales (U) futuro

		//Velocidad radial (V)
		double **Vprevious; //Matriz del mapa de velocidades radiales (V) previo
		double **Vpres; //Matriz del mapa de velocidades radiales (V) presente
		double **Vfut; //Matriz del mapa de velocidades radiales (V) futuro

		//Velocidad tangencial (W)
		double **Wprevious; //Matriz del mapa de velocidades tangneciales (W) previo
		double **Wpres; //Matriz del mapa de velocidades tangneciales (W) presente
		double **Wfut; //Matriz del mapa de velocidades tangneciales (W) futuro

		//Temperatura
		double **Tpres; //Matriz del mapa de temperaturas presente
		double **Tfut; //Matriz del mapa de temperaturas futuro

		//Energía Interna
		double **eprevious;
		double **epres;
		double **efut;

		//Entalpía
		double **Hprevious; //Matriz del mapa de entalpías previo
		double **Hpres; //Matriz del mapa de entalpías presente
		double **Hfut; //Matriz del mapa de entalpías futuro

		//Viscosidad turbulenta (cinemática)
		double **vprevious; //Matriz del mapa de viscosidades turbulentas previo
		double **vpresent; //Matriz del mapa de viscosidades turbulentas presente
		double **vfut; //Matriz del mapa de viscosidades turbulentas futuro

		//Matrices de cálculo de contribuciones a las propiedades

		double **MatrixNueva;
		
		//Densidad
		double **FRhoPrevious;
		double **FRhoPresent;

		//Velocidades
		double **FuPrevious;
		double **FuPresent;

		double **FvPrevious;
		double **FvPresent;

		double **FwPrevious;
		double **FwPresent;

		//Viscosidad
		double **Fmuprevious;
		double **Fmupresent;

		//Energía Interna
		double **Fepresent;
		double **Feprevious;

		//Entalpía
		double **FHpresent;
		double **FHprevious;

		//Temperaturas
		double **EnergyDifusive; //Término difusivo de la ecuación de conservación de la energía
		double **EnergyPressureTerm; //Término de presión de la ecuación de la energía
		double **EnergyConvective; //Término convectivo de la ecuación de energía
		double **EnergyViscous; //Término de disipación viscosa en la ecuación de energía

		double **Divergence; //Matriz de divergencias de la velocidad en cada nodo (Eq energía, término viscoso)
		double **DivergenceMU; //Matriz de divergencias de la velocidad en cada nodo (Eq energía, término viscoso)
		double **DivergenceMR; //Matriz de divergencias de la velocidad en cada nodo (Eq energía, término viscoso)

		double **GradU_DxMU; //Gradiente de velocidad U respecto de X en las caras de los volúmenes de control
		double **GradU_DyMU; //Gradiente de velocidad U respecto de Y en las caras de los volúmenes de control

		double **GradV_DxMU; //Gradiente de velocidad V respecto de X en las caras de los volúmenes de control
		double **GradV_DyMU; //Gradiente de velocidad V respecto de Y en las caras de los volúmenes de control

		double **GradU_DxMR; //Gradiente de velocidad U respecto de X en las caras de los volúmenes de control
		double **GradU_DyMR; //Gradiente de velocidad U respecto de Y en las caras de los volúmenes de control

		double **GradV_DxMR; //Gradiente de velocidad V respecto de X en las caras de los volúmenes de control
		double **GradV_DyMR; //Gradiente de velocidad V respecto de Y en las caras de los volúmenes de control

		double **GradW_DxMU;
		double **GradW_rMR;

		double **K; //Conductividad térmica en los nodos
		double **KwallsMU; //Conductividad térmica en las paredes de los volúmenes de control (Nodos U)
		double **KwallsMV; //Conductividad térmica en las paredes de los volúmenes de control (Nodos R)

		//Presión
		double **PxGradient; //Gradiente de presión dirección Axial
		double **PyGradient; //Gradiente de presión dirección Radial

		double **PwallsU; //Presión en las paredes de lo volúmenes de control (Nodos U)
		double **PwallsR; //Presión en las paredes de lo volúmenes de control (Nodos V)

		double **MomentumDifusiveU; //Término difusivo de la ecuación de cantidad de movimiento (Velocidad axial U)
		double **MomentumDifusiveV; //Término difusivo de la ecuación de cantidad de movimiento (Velocidad raidal V)
		double **MomentumDifusiveW; //Término difusivo de la ecuación de cantidad de movimiento (Velocidad tangencial W)

		double **MomentumConvectiveU; //Término convectivo de la ecuación de cantidad de movimiento (Velocidad axial U)
		double **MomentumConvectiveV; //Término convectivo de la ecuación de cantidad de movimiento (Velocidad raidal V)
		double **MomentumConvectiveW; //Término convectivo de la ecuación de cantidad de movimiento (Velocidad tangencial W)

		//Viscosidad
		double **muBase;
		double **muTurb;

		double **vwallsU;
		double **vwallsR;

		double **muBasewallsU;
		double **muBasewallsR;

		//Matrices de variables del modelo Spalart-Allmaras
		double **ft1;
		double **r;
		double **g;
		double **ft2;
		double **X;
		double **S;
		double **Smodel;
		double **omega;
		double **fw;
		double **fv1;
		double **fv2;

		double **SA_Termino1; //Término 1 de la ecuación de Spalart-Allmaras
		double **SA_Termino2; //Término 2 de la ecuación de Spalart-Allmaras
		double **SA_Termino3; //Término 3 de la ecuación de Spalart-Allmaras
		double **SA_Termino4; //Término 4 de la ecuación de Spalart-Allmaras
		double **SA_Termino5; //Término 5 de la ecuación de Spalart-Allmaras
		double **SA_Termino6; //Término 6 de la ecuación de Spalart-Allmaras

		double **TauRR_mu;
		double **TauRR_mr;

		double **TauRZ_mu; 
		double **TauRZ_mr; 

		double **TauZZ_mu;
		double **TauZZ_mr;

		double **TauRR;
		double **TauRZ;
		double **TauZZ;

		//Matrices y arrays de condiciones de contorno

		//Condiciones de contorno del mapa de densidades
		double *RhoLeft; //Presión pared izquierda
		double *RhoRight; //Presión pared derecha
	
		double *RhoUp; //Presión pared superior
		double *RhoDown; //Presión pared inferior

		//Condiciones de contorno del mapa de velocidades axiales (U)
		double *Uleft; //Velocidad axial pared izquierda
		double *Uright; //Velocidad axial pared derecha
	
		double *Uup; //Velocidad axial pared superior
		double *Udown; //Velocidad axial pared inferior
 
		//Condiciones de contorno del mapa de velocidades radiales (V)
		double *Vleft; //Velocidad radial pared izquierda
		double *Vright; //Velocidad radial pared derecha
	
		double *Vup; //Velocidad radial pared superior
		double *Vdown; //Velocidad radial pared inferior

		//Condiciones de contorno del mapa de velocidades tangenciales (W)
		double *Wleft; //Velocidad radial pared izquierda
		double *Wright; //Velocidad radial pared derecha
	
		double *Wup; //Velocidad radial pared superior
		double *Wdown; //Velocidad radial pared inferior

		//Condiciones de contorno del mapa de temperaturas (T)
		double *Tleft; //Temperaturas pared izquierda
		double *Tright; //Temperaturas pared derecha

		double *Tup; //Temperaturas pared superior
		double *Tdown; //Temperaturas pared inferior

		//Condiciones de contorno del mapa de viscosidades dinámicas
		double *Pleft; //Presión pared izquierda
		double *Pright; //Presión pared derecha

		double *Pup; //Presión pared superior
		double *Pdown; //Presión pared inferior

		//Condiciones de contorno del mapa de viscosidades dinámicas
		double *muLeft; //Viscosidad dinámica pared izquierda
		double *muRight; //Viscosidad dinámica pared derecha

		double *muUp; //Viscosidad dinámica pared superior
		double *muDown; //Viscosidad dinámica pared inferior

		//Condiciones de contorno del mapa de condutividades térmicas
		double *KLeft; //Conductividad térmica pared izquierda
		double *KRight; //Conductividad térmica pared derecha

		double *KUp; //Conductividad térmica pared superior
		double *KDown; //Conductividad térmica pared inferior

		//Condiciones de contorno del mapa de energías internas
		double *eleft; //Energía Interna pared izquierda
		double *eright; //Energía Interna pared derecha

		double *eup; //Energía Interna pared superior
		double *edown; //Energía Interna pared inferior

		//Condiciones de contorno del mapa de entalpías
		double *Hleft; //Entalpía pared izquierda
		double *Hright; //Entalpía pared derecha

		double *Hup; //Entalpía pared superior
		double *Hdown; //Entalpía pared inferior

		//Condiciones de contorno del mapa de viscosidades cinemáticas turbulentas
		double *vleft; //Viscosidad cinemática turbulenta pared izquierda
		double *vright; //Viscosidad cinemática turbulenta pared derecha

		double *vup; //Viscosidad cinemática turbulenta pared arriba
		double *vdown; //Viscosidad cinemática turbulenta pared abajo

		//Condiciones de contorno del mapa de viscosidades dinámicas base
		double *muBaseLeft; //Viscosidad dinámica base pared izquierda
		double *muBaseRight; //Viscosidad dinámica base pared derecha

		double *muBaseDown; //Viscosidad dinámica base pared abajo
		double *muBaseUp; //Viscosidad dinámica base pared arriba

		//Matrices de Valores medios del flujo
		double *TotalMeanTemperature;
		double *TotalMeanDensity;
		double *TotalMeanPressure;
		double *TotalMeanStress;

		double *MeanTemperature;
		double *MeanDensity;
		double *MeanPressure;
		double *MeanStress;

		double *MeanTemperatureInstant;
		double *MeanDensityInstant;
		double *MeanPressureInstant;
		double *MeanStressInstant;

		double TiempoEstadistico;
		
		//Esquema convectivo

		//Función de interpolación lineal
		inline double Interpolacion(double C1, double V1, double C2, double V2, double CPunto){
			return V1 + ((CPunto - C1)/(C2 - C1))*(V2 - V1); }
	
		//Función de aplicación del esquema convectivo
		double EsquemaConvectivo(double, double, double, double, double, double, double, double, double, double, string);

		//Cálculos termoquímicos
		void GetJanaf(ReadData); //Pasar todos los coeficientes termoqímicos de las especies a sus matrices
		void Get_CpCombustion(); //Calculo del Cp de cada uno de los nodos del dominio
		void Get_Enthalpy(); //Cálculo de la entalpía de los gases para una determinada temperatura
		void Get_Rideal(); //Cálculo de la constante ideal de la mezcla de gases

		void Get_CpAir();

		void InitialBoundary();
		void InitializeFields(Mesher); //Pasar todos los coeficientes termoqímicos de las especies a sus matrices
		void UpdateBoundaryConditions(Mesher);

		void InitialF(); //Seteo inicial de las F
		void UpdateF();
		void Get_RhoWalls(Mesher);
		void Get_muWalls(Mesher); //Cálculo de las viscosidades dinámicas en las paredes de los volúmenes de control 
		void Get_nuWalls(Mesher); //Cálculo de las viscosidades cinemáticas en las paredes de los volúmenes de control 
		void Get_VelocityWalls(Mesher); //Cálculo de las velocidades en las paredes de los volúmenes de control
		void Get_K(Mesher); //Cálculo de la conductividad térmica de los gases en cada nodo y paredes

		void Get_AirViscosity();
		void Get_BaseViscosity(); //Cálculo de la viscosidad base de cada nodo con la ley de Sutherland
		void Get_Densities(); //Cálculo del mapa de densidades futuro
		void Get_FRho(Mesher);
		
		void Get_Pgradients(Mesher); //Cálculo de los gradientes de presión en cada nodo en ambas direciones (U y V)
		void Get_Stresses(Mesher);
		void Get_MomemtumDifusiveU(Mesher); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Axial U)
		void Get_MomemtumDifusiveV(Mesher); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
		void Get_MomemtumDifusiveW(Mesher); //Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Tangencial W)

		void Get_MomentumConvectiveU(Mesher); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
		void Get_MomentumConvectiveV(Mesher); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)
		void Get_MomentumConvectiveW(Mesher); //Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Tangencial W)

		void Get_Pwalls(Mesher); //Cálculo de las presiones en las paredes de los volúmenes de control
		void Get_EnergyPressureTerm(Mesher); //Cálculo del término de presión en la ecuación de la energía
		void Get_EnergyDifusive(Mesher); //Cálculo del término difusivo de la ecuación de conservación de la energía
		void Get_EnergyConvective(Mesher); //Cálculo del término convectivo de la ecuación de energía
		void Get_EnergyViscousTerm(Mesher); //Cálculo del término de disipación viscosa en la ecuación de energía

		void Get_eWalls(Mesher);

		void Get_Divergence(Mesher); //Cálculo de la divergencia de la velocidad en cada nodo (Término viscoso eq energía)
		void Get_VelocityGradients(Mesher); //Cálculo de los gradientes de velocidades en cada una de las direcciones
		void Get_InternalEnergy(); //Cálculo del mapa de energías internas específicas
		void Get_SpecEnthalpy(); //Cálculo del mapa de entalpías específicas
		void Get_Fe(); //Cálculo de la contribución total de energía (Fe)
		void Get_Temperatures(); //Cálculo del mapa de temperaturas futuro

		void Calculo_Fu_Fv(); //Cálculo de la contribución total de momentum (Fu/Fv)
		void Get_Velocities(); //Cálculo del mapa de velocidades futuro

		void Get_Pressure(); //Cálculo del mapa de presiones con la Ley de los Gases Ideales

		//Métodos modelo de turbulencia Spalart-Allmaras
		void Get_vwalls(Mesher);
		void Get_muBasewalls(Mesher);
		void SA_Term1(Mesher); //Cálculo del término 1 de la ecuación de Spalart-Allmaras
		void SA_Term2(); //Cálculo del término 2 de la ecuación de Spalart-Allmaras
		void SA_Term3(Mesher); //Cálculo del término 3 de la ecuación de Spalart-Allmaras
		void SA_Term4(Mesher); //Cálculo del término 4 de la ecuación de Spalart-Allmaras
		void SA_Term5(Mesher); //Cálculo del término 5 de la ecuación de Spalart-Allmaras
		void SA_Term6(Mesher); //Cálculo del término 6 de la ecuación de Spalart-Allmaras

		void FmuSA(); //Cálculo de todas las contribuciones del modelo Spalart-Allmaras
		
		void SpalartAllmarasPreparation(Mesher); //Cálculo de variables y matrices necesarias para aplicar el modelo Spalart-Allmaras
		void InitialSpalartAllmaras();
		void SpalartAllmaras(); //Cálculo de la viscosidad con el modelo Spalart-Allmaras


		double Get_Stop();
		void UpdatePropertiesFields();
		void UpdatePeriodicConditions(Mesher);

		double StepTime(Mesher);
		void AllocateMatrix(Memory);
		
		void Get_TotalMU();
		
		void VectorialVTK(Mesher, double**, double**, string, string); //Pasar a un .vtk los resultados de campos vectoriales

		void ComputeParameters(Mesher, int);
		void ComputeMeanVariables(Mesher);

		void ExecuteSolver(Memory, ReadData, Mesher);

		

		

		
};
