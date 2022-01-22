#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cstdio> 
#include <time.h>
#include <chrono>

using namespace std;

//Declaración de funciones a utilizar
double *AllocateDouble(int NX, int NY, int Dim);
int *AllocateInt(int NX, int NY, int Dim);
double Get_TimeStep(int N, double Re, double Courant);
void InitialConditions(int N, double* UpresRe, double* UpresIm, double* DiffusivePast, double* ConvectivePast);
double Get_TotalViscosity(string CFD, int N, int k, double Ck, double* UpresRe, double* UpresIm, double Reynolds);
void Get_Diffusive(double* DiffusivePres, double Viscosity, int k, double* UpresRe, double* UpresIm);
void Get_Convective(double* ConvectivePres, int N, int k, double* UpresRe, double* UpresIm);
void Get_Velocity(int k, string Esquema, double* DiffusivePres, double* ConvectivePres, double* DiffusivePast, double* ConvectivePast, double DeltaT, double* UpresRe, double* UpresIm, double* UfutRe, double* UfutIm);
double Get_StopError(int N, double* UpresRe, double* UpresIm, double* UfutRe, double* UfutIm);
void Get_Update(int N, double* UpresRe, double* UpresIm, double* UfutRe, double* UfutIm, double* DiffusivePres, double* ConvectivePres, double* DiffusivePast, double* ConvectivePast);
void Get_FinalResults(double Re, int N, double* UfutRe, double* UfutIm, string CFD, string Esquema);
void Get_SlopeData(int N);

int main(){

     auto start = std::chrono::high_resolution_clock::now();

//Declaración de variables del problema

    //Variables físicas
    double Re = 500.0; //Número de Reynolds

    double Viscosity;

    //Variables numéricas
    int N = 200; //Número de modos a calcular
    string Esquema = "Euler"; //Esquema de integración a usar (Euler/AdamBashforth)
    string CFD = "DNS"; //Tipo de simulaición (DNS/LES) (Para la viscosidad)
    double Convergencia = 1e-6; //Convergencia simulación
    double Divergencia = 1e10; //Divergencia simulación
    double Ck = 0.4223; //Constante de Kolmogorov (LES)
    double Courant = 0.005; //Número de Courant para el CFL

    //Declaración de Arrays
    double *UpresRe; //Step presente parte real
    double *UpresIm; //Step presente parte imaginaria

    double *UfutRe; //Step futuro parte real
    double *UfutIm; //Step futuro parte imaginaria

    double *DiffusivePres; //Término difusivo 
    double *ConvectivePres; //Término convectivo 

    double *DiffusivePast; //Término difusivo 
    double *ConvectivePast; //Término convectivo 

    //Alojamiento de memoria de los Arrays
    UpresRe = AllocateDouble(N, 1, 1); //Step presente parte real
    UpresIm = AllocateDouble(N, 1, 1); //Step presente parte imaginaria

    UfutRe = AllocateDouble(N, 1, 1); //Step futuro parte real
    UfutIm = AllocateDouble(N, 1, 1); //Step futuro parte imaginaria

    DiffusivePres = AllocateDouble(2, 1, 1); //Término difusivo
    ConvectivePres = AllocateDouble(2, 1, 1); //Término convectivo 

    DiffusivePast = AllocateDouble(2, 1, 1); //Término difusivo
    ConvectivePast = AllocateDouble(2, 1, 1); //Término convectivo 

    //Declaración de variables necesarias para el problema
    double DeltaT = Get_TimeStep(N, Re, Courant);

    //Set of the Initial Condition of the Velocity Arrays
    InitialConditions(N, UpresRe, UpresIm, DiffusivePast, ConvectivePast);

    //Variables para la simulación
    int NITER = 0;
    double Difference = 2.0*Convergencia;
    int k;

    //Simulation Process
    while(Difference >= Convergencia && Difference <= Divergencia){

        //Velocidades del modo inicial (1)
        UpresRe[0] = 1.0;
        UpresIm[0] = 0.0;

        UfutRe[0] = 1.0;
        UfutIm[0] = 0.0;

        for(k = 1; k < N; k++){

            //Cálculo de la viscosidad
            Viscosity = Get_TotalViscosity(CFD, N, k, Ck, UpresRe, UpresIm, Re);
            
            //Cálculo del término difusivo
            Get_Diffusive(DiffusivePres, Viscosity, k, UpresRe, UpresIm);
            
            //Cálculo del término convectivo
            Get_Convective(ConvectivePres, N, k, UpresRe, UpresIm);

            //Cálculo de las velocidades
            Get_Velocity(k, Esquema, DiffusivePres, ConvectivePres, DiffusivePast, ConvectivePast, DeltaT, UpresRe, UpresIm, UfutRe, UfutIm);

        }
    
        //Cálculo del error relativo entre Steps
        Difference = Get_StopError(N, UpresRe, UpresIm, UfutRe, UfutIm);

        //Update de los Arrays
        Get_Update(N, UpresRe, UpresIm, UfutRe, UfutIm, DiffusivePres, ConvectivePres, DiffusivePast, ConvectivePast);

        NITER += 1;

        if(NITER%1000 == 0){
            cout<<"Iteración: "<<NITER<<", Error Relativo: "<<Difference<<endl;
        }

    }

    cout<<"Final de la simulación."<<endl;
    cout<<"Iteraciones Totales: "<<NITER<<endl;
    cout<<"Error Relativo Final: "<<Difference<<endl;

    Get_FinalResults(Re, N, UfutRe, UfutIm, CFD, Esquema);
    Get_SlopeData(N);

    // Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}

//Memoria dinámica matriz (double)
double *AllocateDouble(int NX, int NY, int Dim){
double *M1;

	M1 = new double [NX*NY*Dim];				
	return M1;
}

//Memoria dinámica matriz (int)
int *AllocateInt(int NX, int NY, int Dim){ 
int *M1;

	M1 = new int [NX*NY*Dim];				
	return M1;
}

//Cálculo del time Step
double Get_TimeStep(int N, double Re, double Courant){
double DeltaT;

    DeltaT = Courant*(Re/pow(N,2.0));

    return DeltaT;

}

//Set de las condiciones iniciales de las variables
void InitialConditions(int N, double* UpresRe, double* UpresIm, double* DiffusivePast, double* ConvectivePast){
int i;
double I;

    for(i = 0; i < N; i++){
        I = i;
        UpresRe[i] = 1.0/(I + 1.0);
        UpresIm[i] = 0.0;
    }

    DiffusivePast[0] = 0.0;
    DiffusivePast[1] = 0.0;

    ConvectivePast[0] = 0.0;
    ConvectivePast[1] = 0.0;

}

//Cálculo de la viscosidad 
double Get_TotalViscosity(string CFD, int N, int k, double Ck, double* UpresRe, double* UpresIm, double Reynolds){
double K = k;
double m = 2.0;
double Ekn;
double Viscosity = 1.0/Reynolds;

    if(CFD == "LES"){
        Ekn = pow(UpresRe[N-1],2.0) + pow(UpresIm[N-1],2.0);
        Viscosity += (0.31*((5.0 - m)/(m + 1.0))*sqrt(3.0 - m)*pow(Ck, -1.5))*(pow(Ekn/N, 0.50))*(1 + 34.5*exp(-3.03*(N/K)));
    }
   
    return Viscosity;
}

//Cálculo del término difusivo de la ecuación de Momentum
void Get_Diffusive(double* DiffusivePres, double Viscosity, int k, double* UpresRe, double* UpresIm){

	DiffusivePres[0] = Viscosity*pow(k + 1,2.0)*UpresRe[k]; //Parte Real Difusivo (+1 porque empieza en 0 el array)
	DiffusivePres[1] = Viscosity*pow(k + 1,2.0)*UpresIm[k]; //Parte Imaginaria Difusivo

}

//Cálculo del término convectivo de la ecuación de Momentum
void Get_Convective(double* ConvectivePres, int N, int k, double* UpresRe, double* UpresIm){
int p, q;
int ip, iq;

    ConvectivePres[0] = 0.0;
    ConvectivePres[1] = 0.0;
	
	for (int p = - N; p <= N; p++) {
		q = k - p + 1;

		ip = abs(p) - 1; 
		iq = abs(q) - 1; 

		if (q >= - N && q <= N) { 
			if (q == 0 || p == 0) {

                ConvectivePres[0] += 0.0;
                ConvectivePres[1] += 0.0;

            }
			else if (q < 0) {

                ConvectivePres[0] += + q*(+ UpresRe[ip]*UpresIm[iq] - UpresRe[iq]*UpresIm[ip]);
				ConvectivePres[1] += + q*(UpresRe[ip]*UpresRe[iq] + UpresIm[ip]*UpresIm[iq]);     

			}
			else if (p < 0) {

                ConvectivePres[0] += + q*(- UpresRe[ip]*UpresIm[iq] + UpresRe[iq]*UpresIm[ip]);
				ConvectivePres[1] += + q*(UpresRe[ip]*UpresRe[iq] + UpresIm[ip]*UpresIm[iq]);     
              
			}
			else {

				ConvectivePres[0] += + q*(- UpresRe[ip]*UpresIm[iq] - UpresRe[iq]*UpresIm[ip]);
				ConvectivePres[1] += + q*(UpresRe[ip]*UpresRe[iq] - UpresIm[ip]*UpresIm[iq]);     

			}
		}
	}

}

void Get_Velocity(int k, string Esquema, double* DiffusivePres, double* ConvectivePres, double* DiffusivePast, double* ConvectivePast, double DeltaT, double* UpresRe, double* UpresIm, double* UfutRe, double* UfutIm){

    if(Esquema == "Euler"){
        UfutRe[k] = UpresRe[k] + DeltaT*(- ConvectivePres[0] - DiffusivePres[0]);
        UfutIm[k] = UpresIm[k] + DeltaT*(- ConvectivePres[1] - DiffusivePres[1]);
    }
    else if(Esquema == "AdamBashforth"){
        UfutRe[k] = UpresRe[k] + DeltaT*(1.50*(- ConvectivePres[0] - DiffusivePres[0]) - 0.50*(- ConvectivePast[0] - DiffusivePast[0]));
        UfutIm[k] = UpresIm[k] + DeltaT*(1.50*(- ConvectivePres[1] - DiffusivePres[1]) - 0.50*(- ConvectivePast[1] - DiffusivePast[1]));
    }
    
}

double Get_StopError(int N, double* UpresRe, double* UpresIm, double* UfutRe, double* UfutIm){
int k;
double Upres, Ufut;
double MaxDiff = 0.0;

    for(k = 0; k < N; k++){

        Upres = sqrt(pow(UpresRe[k],2.0) + pow(UpresIm[k],2.0));
        Ufut = sqrt(pow(UfutRe[k],2.0) + pow(UfutIm[k],2.0));

        if(abs((Ufut - Upres)/(Upres + 1e-10)) >= MaxDiff){
            MaxDiff = abs((Ufut - Upres)/(Upres + 1e-10));
        }
    }

    return MaxDiff;

}

void Get_Update(int N, double* UpresRe, double* UpresIm, double* UfutRe, double* UfutIm, double* DiffusivePres, double* ConvectivePres, double* DiffusivePast, double* ConvectivePast){
int k;

    for(k = 0; k < N; k++){
        UpresRe[k] = UfutRe[k];
        UpresIm[k] = UfutIm[k];

        DiffusivePast[0] = DiffusivePres[0];
        DiffusivePast[1] = DiffusivePres[1];

        ConvectivePast[0] = ConvectivePres[0];
        ConvectivePast[1] = ConvectivePres[1];
    }

}

void Get_FinalResults(double Re, int N, double* UfutRe, double* UfutIm, string CFD, string Esquema){
int k;
double Ek;
char Directorio[500];
int RE = Re;

    sprintf(Directorio,"/home/sergiogus/Desktop/ComputationalEngineering/BurgersEquation/Results/Espectro_RE_%d_N_%d_.txt", RE, N);

    FILE *fp1;
	fp1 = fopen(Directorio,"w");

        for(k = 0; k < N; k++){
            Ek = pow(UfutRe[k],2.0) + pow(UfutIm[k],2.0);

            fprintf(fp1,"%d \t %.12f \n", k + 1, Ek);	
        }
							
	fclose(fp1);

}

void Get_SlopeData(int N){
int k;
double Exponente = -2.0;
int Exp = Exponente;
double K;
char Directorio2[500];

     sprintf(Directorio2,"/home/sergiogus/Desktop/ComputationalEngineering/BurgersEquation/Results/Slope%d.txt", abs(Exp));

    FILE *fp2;
	fp2 = fopen(Directorio2,"w");

        for(k = 0; k < N; k++){
            K = k;
            fprintf(fp2,"%d \t %f \n", k + 1, pow(K, Exponente));	
        }
							
	fclose(fp2);
}