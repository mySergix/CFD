#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

class Solver{
	private:

		//Problema a simular 
		string Problema; //(Canal/Canal con Cilindro/Perfil)

		//Densidad del mallado
		int NX; //Dirección horizontal
		int NY; //Dirección vertical

		//Datos para la computación en paralelo
		int Ix;
		int Fx;
		int Procesos;
		int Halo;
		int Rank;
		
		//Datos de la geometría
		double ChannelLength; //Longitud del canal
		double ChannelHeight; //Altura del canal
		double ChannelDepth; //Profundidad del canal

		double CylinderRadius;

		//Datos Físicos del problema
		double Uref;
		double Vref;
		double Tref;
		double RhoRef;
		double Pref;
		double Cp;
		double Rideal;
		double Gamma;

		double Circulation;
		double PhiCilindro;

		//Datos sobre el solver iterativo
		double ConvergenciaGlobal; 
		double ConvergenciaGS;

		double MaxDiffGS;
		double MaxDiffGlobal;

		double Cd;
		double Cl;

		double AnguloAtaque;
		
		double NumericalCirculation;

		double *PDiff;
	
	public:
		Solver(Memory, ReadData, ParPro, Mesher, int);

		//Matrices de los mapas de propiedades del problema
	
		//Matrices globales de propiedades
		double *PCirc;

		//Step Presente
		double *PhiGlobalPres;
		double *UglobalPres;
		double *VglobalPres;
		double *TglobalPres;
		double *RhoGlobalPres;
		double *PglobalPres;

		//Step Futuro
		double *PhiGlobalFut;
		double *UglobalFut;
		double *VglobalFut;
		double *TglobalFut;
		double *RhoGlobalFut;
		double *PglobalFut;

		//Matrices locales de propiedades

		int *DensityIndicator;
		int *GlobalDensityIndicator;

		//Streamline
		double *PhiLocalPres;
		double *PhiLocalFut;

		double *PhiLocalSup;

		//Velocidad Horizontal
		double *UlocalPres;
		double *UlocalFut;

		//Velocidad Vertical
		double *VlocalPres;
		double *VlocalFut;

		//Temperatura
		double *TlocalPres;
		double *TlocalFut;

		//Densidad
		double *RhoLocalPres;
		double *RhoLocalFut;

		//Presión
		double *PlocalPres;
		double *PlocalFut;
		
		//Matrices de valores de las propiedades en las paredes de los volúmenes de control
		double *UwallsMR;
		double *VwallsMU;

		double *UwallsMR_Global;
		double *VwallsMU_Global;

		double *RhoWallsMU;
		double *RhoWallsMR;

		//Matrices de cálculo de contribuciones a las propiedades

		double *aw;
		double *ae;
		double *as;
		double *an;
		double *ap;

		double *awGlobal;
		double *aeGlobal;
		double *asGlobal;
		double *anGlobal;
		double *apGlobal;
		
		//Matrices y arrays de condiciones de contorno

		//Condiciones de contorno del mapa de valor de streamlines
		double *PhiLeft; //Presión pared izquierda
		double *PhiRight; //Presión pared derecha
	
		double *PhiUp; //Presión pared superior
		double *PhiDown; //Presión pared inferior

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

		//Métodos de la clase Solver
		void AllocateMatrix(Memory);

		void InitializeFields(Mesher); //Pasar todos los coeficientes termoqímicos de las especies a sus matrices
		void UpdateBoundaryConditions(Mesher);

		void Get_DensityIndicatorCylinder(Mesher);
		void Get_DensityIndicatorAirfoil(Mesher);
		
		void Get_DensityWalls(Mesher, ParPro);
		void Get_Coeficients(Mesher); //Cálculo de los coeficientes de discretización

		void Get_Streamlines(ParPro, Mesher);
		void Get_MaxDifGS(ParPro);
		void Get_Velocities(Mesher, ParPro); //Cálculo del mapa de velocidades futuro
		void Get_Temperatures(); //Cálculo del mapa de temperaturas futuro
		void Get_Density();
		void Get_Pressure();

		void Get_Stop(ParPro);

		void UpdatePropertiesFields();
		void Get_Circulation(Mesher, ParPro);
		void Get_AeroCoefficients(Mesher);

		void Send_Coefficients(Mesher, ParPro);

		void EscalarVTK2D(Mesher, string, string, string, double*, double*, int, int);
		void EscalarVTK2DInt(Mesher, string, string, string, double*, int*, int, int);
		void VectorialVTK2D(Mesher, string, string, string, double*, double*, double*, int, int);

		void ExecuteSolver(Memory, ReadData, ParPro, Mesher);
		
};
