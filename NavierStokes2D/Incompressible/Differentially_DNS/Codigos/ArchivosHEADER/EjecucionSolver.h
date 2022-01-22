#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>

using namespace std;

void CalculoSolver::EjecucionSolver(SLU C1, Mallador M1){
double DeltaT;
double TotalTime = 0.0;
double MaxDifParada = 1.0;
int Step = 0;
int RE = Re;
int logRA = log10(Rayleigh);

//Inicialización de las variables
InitialFields();  

//Creación matriz de coeficientes "a"
C1.constSLU(NXTOT, NYTOT);
CalculoCoeficientes(M1);
C1.changeA(ap, ae, aw, an, as, -1);

//Condiciones y cálculos iniciales
BoundaryConditions(M1);
CalculoVelocidades(M1, DeltaT);
CalculoFinicial(M1);
if(Problema == 2){ CalculoFTinicial(M1, Step); }

	while(abs(MaxDifParada) >= ConvergenciaSteady){
		//Calculo del Time Step
		DeltaT = StepTime(M1);
		TotalTime = TotalTime + DeltaT;
		Step = Step + 1;
		
		//Calculo de las presiones y velocidades
		BoundaryConditions(M1);
		CalculoF(M1);
		CalculoBP(M1, DeltaT);
		C1.changeB(bp);
		C1.calcP();
		C1.ReturnP(P);
		CalculoVelocidades(M1, DeltaT);

		//Comprobacion perfil de velocidades Reynolds 100
		if(Problema == 3 && Re == 100.0){
			ComprobacionPerfil(M1);
		}
		
		//Calculo de las temperaturas
		if(Problema == 2){
			CalculoFT(M1, Step);
			CalculoTemperaturas(DeltaT);
		}

		//Calculo de los coeficientes adimensionales
		if(Problema == 3){ CoeficientesFluctuantes(M1, TotalTime); }
		if(Problema == 3 && Re >= 60){ CalculoStrouhal(DeltaT, TotalTime); }
		
		

		//Mostrar por pantalla los resultados cada X steps
		if(Step%StepsPantalla == 0){
			MaxDifParada = Parada();
			
			
			
			//Pasar los resultados a formato VTK para usarlos en Paraview
			PrintVTK(M1, Step);

			//Cálculo coeficientes adimensionales no fluctuantes
			if(Problema == 3){
				CalculoCoefficients(M1);
			}
		
			//Mostrar los resultados por pantalla en ese instante de simulación
			if(Problema == 1 || Problema == 3){
				cout<<"Step: "<<Step<<" Salto: "<<MaxDifParada<<" Tiempo: "<<TotalTime<<" Periodo: "<<Periodo<<" RE: "<<Re<<endl;	
			}
			else if(Problema == 2){
				if(TotalTime > 50.0){
					ResultadosNusselt(M1, TotalTime);
				}
				cout<<"Step: "<<Step<<" Salto: "<<MaxDifParada<<" Tiempo: "<<TotalTime<<" Periodo: "<<Periodo<<" RA: "<<Rayleigh<<endl;	
			}
		}
		
		//PrintGnuPlot(Step, TotalTime, M1);

		//Calculo número de Nusselt
			if(Problema == 2){
				AverageFlowTurbulence(M1, TotalTime, DeltaT);
			}

		Actualizar();

		if(Problema == 3){
			if(abs(TotalTime - 30.0) <= 1e-1 && Re >= 60){ ValorInstante = Ufut[NX1 + NX2 + NX3/3][NY1 + NY2/2]; }
		}
	}
	
	//Pasar los resultados de velocidad, presión y temperatura a TXT
	TXT(M1);
}



