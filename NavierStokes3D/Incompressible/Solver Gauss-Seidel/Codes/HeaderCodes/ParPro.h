#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

class ParPro{
	private:
		
		
	public:

		//Datos generales del problema
		int NX;
		int NY;
		int NZ;
		
		//Datos y variables de computación paralela
		int Rank;
		int Procesos;
		int Ix;
		int Fx;
		int Halo;
		int pix, pfx;

		int HP;

		//Constructor de la clase
		ParPro(ReadData);
		
		//Métodos de la clase
		void Rango();
		void Processes();
		void Initial_WorkSplit(int, int&, int&);
		void Get_Worksplit(int, int, int, int&, int&);

		void CommunicateDataLP(double*, double*, int, int);
		void CommunicateDataLU(double*, double*, int, int);
		void CommunicateDataLV(double*, double*, int, int);
		void CommunicateDataLW(double*, double*, int, int);
		
		void SendMatrixToZeroMP(double*, double*, int, int, int, int, int, int);
		void SendMatrixToZeroMU(double*, double*, int, int, int, int, int, int);
		void SendMatrixToZeroMV(double*, double*, int, int, int, int, int, int);
		void SendMatrixToZeroMW(double*, double*, int, int, int, int, int, int);
		
		void SendDataToZero(double, double*);
		void SendDataToAll(double, double&);

		void SendMatrixToZeroBP(double*, double*, int , int , int , int , int , int );
		void SendMatrixToAll_Pressure(double*, double*, int, int, int, int, int, int);

		void Execute();

		int Get_Rank(){ return Rank; }
		int Get_Processes(){ return Procesos; }
		int Get_Halo(){ return Halo; }
		int Get_Ix(){ return Ix; }
		int Get_Fx(){ return Fx; }

};
