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
		int Problema;

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

		//Datos Físicos del problema
		double Peclet;
		double Uref;
		double Vref;
		double RhoRef;
		double Alpha;
		double AlphaRad;

		double PhiIz;
		double PhiDer;
		double PhiTop;
		double PhiBotton;

		string EsquemaLargo;
		string EsquemaCorto;

		double FR;
		
		double DeltaT;

		//Datos sobre el solver iterativo
		double ConvergenciaGlobal; 
		double ConvergenciaGS;

		double MaxDiffGS;
		double MaxDiffGlobal;

		double *PDiff;
	
	public:
		Solver(Memory, ReadData, ParPro, Mesher, int);

		//Matrices de los mapas de propiedades del problema
	
		//Matrices globales de propiedades

		//Step Presente
		double *PhiGlobalPres;
		double *UglobalPres;
		double *VglobalPres;

		//Step Futuro
		double *PhiGlobalFut;
		double *UglobalFut;
		double *VglobalFut;

		double *PDT;

		//Matrices locales de propiedades

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

		//Matrices de valores de las propiedades en las paredes de los volúmenes de control
		double *UwallsMU;
		double *VwallsMR;

		double *PhiWallsMU;
		double *PhiWallsMR;

		//Matrices de cálculo de contribuciones a las propiedades

		double *Dw;
		double *De;
		double *Ds;
		double *Dn;

		double *Cw;
		double *Ce;
		double *Cs;
		double *Cn;

		double *aw;
		double *ae;
		double *as;
		double *an;
		double *ap;
		
		double *bp;

		//Matrices y arrays de condiciones de contorno

		//Condiciones de contorno del mapa de valor de streamlines
		double *PhiLeft; //Presión pared izquierda
		double *PhiRight; //Presión pared derecha
	
		double *PhiUp; //Presión pared superior
		double *PhiDown; //Presión pared inferior

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

		//Matrices para verificar el caso analítico
		double *CoordenadasAnalitico;
		double *PhiAnalitico;

		double *PhiReal;

		double *MinAbajo;
		double *MinArriba;

		double *PhiAbajo;
		double *PhiArriba;

		//Métodos de la clase Solver
		void AllocateMatrix(Memory);
		void Get_AnalyticalResults();
		double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);
		void Get_StepTime(Mesher, ParPro);

		void InitializeFields(Mesher); //Pasar todos los coeficientes termoqímicos de las especies a sus matrices
		void UpdateBoundaryConditions(Mesher);

		void Get_PhiWalls(Mesher, ParPro);
		void Get_PhiWalls2(Mesher, ParPro);
		void Get_Diffusive(Mesher);
		void Get_Convective(Mesher);

		void Get_Coeficients(Mesher); //Cálculo de los coeficientes de discretización
		void Get_BetaCoefficient(Mesher);

		void Get_StreamlinesFR(ParPro);
		void Get_Streamlines(ParPro);
		void Get_MaxDifGS(ParPro);
		void Get_Velocities(Mesher); //Cálculo del mapa de velocidades futuro

		void Get_Stop(ParPro);

		void UpdatePropertiesFields();

		void Get_NumericalResults(Mesher, string, string);
		void RelativeError(Mesher, string, string);
		void EscalarVTK2D(Mesher, string, string, string, double*, double*, int, int);
		void VectorialVTK2D(Mesher, string, string, string, double*, double*, double*, int, int);
		void Get_SymetricOutlet(Mesher, string, string);

		void ExecuteSolver(Memory, ReadData, ParPro, Mesher, int);
		
};
