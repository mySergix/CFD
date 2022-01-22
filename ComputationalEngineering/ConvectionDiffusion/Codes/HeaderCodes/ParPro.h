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

		//Datos y variables de computación paralela
		int Rank;
		int Procesos;
		int Ix;
		int Fx;
		int Halo;
		int pix, pfx;

		//Constructor de la clase
		ParPro(ReadData, int);
		
		//Métodos de la clase
		void Rango();
		void Processes();
		void Initial_WorkSplit(int, int&, int&);
		void Get_Worksplit(int, int, int, int&, int&);

		void SendData(double*, int, int);
		void ReceiveData(double*, int, int);
		void SendMatrixToZero(double*, double*, int, int, int, int, int);
		void SendDataToZero(double, double*);
		void SendDataToAll(double, double&);
		
		void Execute();

		int Get_Rank(){ return Rank; }
		int Get_Processes(){ return Procesos; }
		int Get_Halo(){ return Halo; }
		int Get_Ix(){ return Ix; }
		int Get_Fx(){ return Fx; }

};
