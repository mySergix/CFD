#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;

#define PI 3.141592653589793
#define max(a,b) ((a)>(b)?(a):(b));

//VELOCIDADES EJERCICIO 4

//Opcion -> Ejercicio (1-4)
double U(double X, double Y, double Angulo, double Uo, double Opcion){
	
		double VelocidadX;
		if (Opcion == 4){
			VelocidadX = 2*Y*(1-pow(X,2));
		}
		else{
			VelocidadX = Uo*cos(Angulo);
		}
		
		return VelocidadX;
}
//Opcion -> Ejercicio (1-4)
double V(double X, double Y, double Angulo, double Uo, double Opcion){

		double VelocidadY;
		if (Opcion == 4){
			VelocidadY = -2*X*(1-pow(Y,2));
		}
		else{
			VelocidadY = Uo*sin(Angulo);
		}
		
		return VelocidadY;
}


//Funcion para elegir el esquema de interpolacion en las caras
double Interpolacion(double Phi1, double Coordenada1, double Phi2, double Coordenada2, double Phi3, double Coordenada3, double Phi4, double Coordenada4, double CoordenadaCara, double Flujo, string Esquema){
	
	double PHI;

	double XD;
	double PhiD;
	double XC;
	double PhiC;
	double XU;
	double PhiU;

	if (Flujo < 0){
		XD = Coordenada2;
		PhiD = Phi2;
		XC = Coordenada3;
		PhiC = Phi3;
		XU = Coordenada4;
		PhiU = Phi4;
	}
	else{
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

	if (PhiD == PhiU){
		PHI = PhiD;
	}
	else{
	PhC = (PhiC - PhiU)/(PhiD - PhiU);

	XadC = (XC - XU)/(XD - XU);

	Xade = (CoordenadaCara - XU)/(XD - XU);

	//Evaluacion
	double Phf;

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
		Phf = Xade + ((Xade*(Xade - 1))/(XadC*(XadC-1)))*(PhC - XadC);
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

double MediaVelocidades(double Vel1, double Vel2){
	double Velocidad = (Vel1 + Vel2)/2;
	return Velocidad;
}
//Condiciones de contorno para el problema 4

double Cond1(double CoordenadaX, double Alpha){
	double PHI = 1 + tanh((2*CoordenadaX + 1)*Alpha);
	return PHI;
}

int main(){

int i, j;
int N = 20;
int Problema;
string Esquema;
string Esquema2;
double Convergencia = 1e-10;
double MaxDif;
cout<<"Introduce el número del problema"<<endl;
cin>>Problema;
cout<<"Elige un esquema de interpolacion (nodos centrales) (UDS, CDS, SUDS, QUICK, SMART)"<<endl;
cin>>Esquema;
cout<<"Elige un esquema de interpolacion nodos laterales (UDS, CDS)"<<endl;
cin>>Esquema2;

double Angulo;
int NX;
int NY;
if (Problema == 1){

	Angulo = 0;
	NX = N;
	NY = N;
}
else if (Problema == 2){
	Angulo = -PI/2;
	NX = N;
	NY = N;
}
else if (Problema == 3 || Problema == 4){
	Angulo = (PI/4);
	NX = N;
	NY = N;
}

double X = 2.0;
double Y = 1.0;


double Uo = 1;
double Rho = 10;
double Coordenadas[NX+2][NY+2][2];
double u[NX+2][NY+2][2];
double Phi[NX+2][NY+2];
double PhiS[NX+2][NY+2];

double DeltaX;
double DeltaY;

//Valores para las condiciones de contorno
double PHI0 = 5;
double PHIL = 10;
//Problema 3
//PHI1 = PHI0
//PHI2 = PHIL

double Alpha = 10; //Problema 4

//double **Ce;
//double **Cw;
//double **Cn;
//double **Cs;


double Ce[NX+1][NY+1];
double Cw[NX+1][NY+1];
double Cs[NX+1][NY+1];
double Cn[NX+1][NY+1];

//double **aE;
//double **aW;
//double **aS;
//double **aN;
//double **aP;
//double **bP;

double aE[NX+1][NY+1];
double aW[NX+1][NY+1];
double aS[NX+1][NY+1];
double aN[NX+1][NY+1];
double aP[NX+1][NY+1];
double bP[NX+1][NY+1];
//Coordenadas de los nodos
	// 1 -> Direccion X
	// 2 -> Direccion Y
if(Problema != 4){
DeltaX = X/NX;
DeltaY = Y/NY;
	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){
			if (i == 0){
				Coordenadas[i][j][1] = 0;
				if(j == 0){
					Coordenadas[i][j][2] = 0;
				}
				else if(j == NY+1){
					Coordenadas[i][j][2] = Y;
				}
				else{
					Coordenadas[i][j][2] = DeltaY/2 + DeltaY*(j-1);
				}
			}
			else if(i == NX+1){
				Coordenadas[i][j][1] = X;
				if(j == 0){
					Coordenadas[i][j][2] = 0;
				}
				else if(j == NY+1){
					Coordenadas[i][j][2] = Y;
				}
				else{
					Coordenadas[i][j][2] = DeltaY/2 + DeltaY*(j-1);
				}
			}
			if (j==0 && i > 0 && i < NX+1){
				Coordenadas[i][j][2] = 0;
				Coordenadas[i][j][1] = DeltaX/2 + DeltaX*(i-1);
			}
			else if(j==NY+1 && i > 0 && i < NX+1){
				Coordenadas[i][j][2] = Y;
				Coordenadas[i][j][1] = DeltaX/2 + DeltaX*(i-1);
			}
			if (i > 0 && i < NX+1 && j > 0 && j < NY+1){
				Coordenadas[i][j][1] = DeltaX/2 + DeltaX*(i-1);
				Coordenadas[i][j][2] = DeltaY/2 + DeltaY*(j-1);
			}
			
		}
	}
}else{
DeltaX = 2.0/NX;
DeltaY = 1.0/NY;
	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){
			if (i == 0){
				Coordenadas[i][j][1] = -1;
				if(j == 0){
					Coordenadas[i][j][2] = 0;
				}
				else if(j == NY+1){
					Coordenadas[i][j][2] = 1;
				}
				else{
					Coordenadas[i][j][2] = DeltaY/2 + DeltaY*(j-1);
				}
			}
			else if(i == NX+1){
				Coordenadas[i][j][1] = 1;
				if(j == 0){
					Coordenadas[i][j][2] = 0;
				}
				else if(j == NY+1){
					Coordenadas[i][j][2] = 1;
				}
				else{
					Coordenadas[i][j][2] = DeltaY/2 + DeltaY*(j-1);
				}
			}
			if (j==0 && i > 0 && i < NX+1){
				Coordenadas[i][j][2] = 0;
				Coordenadas[i][j][1] = -1 + DeltaX/2 + DeltaX*(i-1);
			}
			else if(j==NY+1 && i > 0 && i < NX+1){
				Coordenadas[i][j][2] = 1;
				Coordenadas[i][j][1] = -1 + DeltaX/2 + DeltaX*(i-1);
			}
			if (i > 0 && i < NX+1 && j > 0 && j < NY+1){
				Coordenadas[i][j][1] = -1 + DeltaX/2 + DeltaX*(i-1);
				Coordenadas[i][j][2] = DeltaY/2 + DeltaY*(j-1);
			}
			
		}
	}
}	
	


//Velocidades en los nodos

		for (i=0;i<=NX+1;i++){
			for (j=0;j<=NY+1;j++){
				u[i][j][1] = U(Coordenadas[i][j][1],Coordenadas[i][j][2],Angulo,Uo,Problema);
				u[i][j][2] = V(Coordenadas[i][j][1],Coordenadas[i][j][2],Angulo,Uo,Problema);
			}
		}

//Asignacion del mapa inicial de propiedades
		for (i=0;i<=NX+1;i++){
			for(j=0;j<=NY+1;j++){
				
				Phi[i][j] = i + j;
				PhiS[i][j] = Phi[i][j];
			}
			
		}

		//Condiciones de contorno
		if (Problema == 1 || Problema == 2){
			for (j=0;j<=NY+1;j++){
				Phi[0][j] = PHI0;
				Phi[NX+1][j] = PHIL;
			}
			for (i=0;i<=NX+1;i++){
				Phi[i][0] = Phi[i][1];
				Phi[i][NY+1] = Phi[i][NY];
			}
			
		}
		if (Problema == 3){
			for (j=0;j<=NY+1;j++){
				Phi[0][j] = PHI0;
				Phi[NX+1][j] = PHIL;
			}
			for (i=0;i<=NX+1;i++){
				Phi[i][0] = PHIL;
				Phi[i][NY+1] = PHI0;
			}
		} 
		if (Problema == 4){
			for (i=0;i<=NX+1;i++){
				if(Coordenadas[i][0][1] < 0){
					Phi[i][0] = Cond1(Coordenadas[i][0][1],Alpha);
				}else{
					Phi[i][0] = Phi[i][1];
				}
				Phi[i][NY+1] = 1 - tanh(Alpha);
			}
			for (j=0;j<=NY+1;j++){
				Phi[0][j] = 1 - tanh(Alpha);
				Phi[NX+1][j] = 1 - tanh(Alpha);
			}
		}

		
		
		//Calculo de los terminos difusivo y convectivo
		//Terminos difusivos
		double De[NX+1][NY+1];
		double Dw[NX+1][NY+1];
		double Ds[NX+1][NY+1];
		double Dn[NX+1][NY+1];

		for (i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){
				De[i][j] = DeltaY/DeltaX;
				Dw[i][j] = DeltaY/DeltaX;
				Ds[i][j] = DeltaX/DeltaY;
				Dn[i][j] = DeltaX/DeltaY;
			}
		}

		//Cabales masicos
		double Me[NX+1][NY+1];
		double Mw[NX+1][NY+1];
		double Ms[NX+1][NY+1];
		double Mn[NX+1][NY+1];
		for(i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){
				Me[i][j] = -Rho*MediaVelocidades(u[i+1][j][1],u[i][j][1])*DeltaY;
				Mw[i][j] = Rho*MediaVelocidades(u[i-1][j][1],u[i][j][1])*DeltaY;
				Ms[i][j] = Rho*MediaVelocidades(u[i][j-1][2],u[i][j][2])*DeltaX;
				Mn[i][j] = -Rho*MediaVelocidades(u[i][j+1][2],u[i][j][2])*DeltaX;
			}
		}

		//Terminos convectivos
		for(i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){
				Ce[i][j] = Me[i][j];
				Cw[i][j] = Mw[i][j];
				Cs[i][j] = Ms[i][j];
				Cn[i][j] = Mn[i][j];
			}
		}

		

		//Calculo de los valores de la propiedad Phi
		MaxDif = 2*Convergencia;
		int NITER = 0;
		while (MaxDif >= Convergencia){
		
		//Condiciones de contorno
		if (Problema == 1 || Problema == 2){
			for (j=0;j<=NY+1;j++){
				Phi[0][j] = PHI0;
				Phi[NX+1][j] = PHIL;
			}
			for (i=0;i<=NX+1;i++){
				Phi[i][0] = Phi[i][1];
				Phi[i][NY+1] = Phi[i][NY];
			}
			
		}
		if (Problema == 3){
			for (j=0;j<=NY+1;j++){
				Phi[0][j] = PHI0;
				Phi[NX+1][j] = PHIL;
			}
			for (i=0;i<=NX+1;i++){
				Phi[i][0] = PHIL;
				Phi[i][NY+1] = PHI0;
			}
		} 
		if (Problema == 4){
			for (i=0;i<=NX+1;i++){
				if(Coordenadas[i][0][1] < 0){
					Phi[i][0] = Cond1(Coordenadas[i][0][1],Alpha);
				}else{
					Phi[i][0] = Phi[i][1];
				}
				Phi[i][NY+1] = 1 - tanh(Alpha);
			}
			for (j=0;j<=NY+1;j++){
				Phi[0][j] = 1 - tanh(Alpha);
				Phi[NX+1][j] = 1 - tanh(Alpha);
			}
		}



		//Calculo de los coeficientes de discretizacion
	if(Problema != 4){
		for (i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){
				//Nodos centrales
				if (i != 1 && i != NX && j != 1 && j != NY){
					aE[i][j] = De[i][j];
					aW[i][j] = Dw[i][j];	
					aS[i][j] = Ds[i][j];
					aN[i][j] = Dn[i][j];
					aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]);
					bP[i][j] = + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema) 
					+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
					+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema) 
					+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
				}
				//Columna izquierda
				if (i==1){
					aW[i][j] = 0;	
					aE[i][j] = De[i][j];
					//Esquina inferior
					if(j==1){
						aS[i][j] = 0;
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i-1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema) 
						+ Mw[i][j]*Interpolacion(0, 0, Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema2) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema) 
						+ Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
					}
				
					//Esquina superior
					else if(j==NY){
						aS[i][j] = Ds[i][j];
						aN[i][j] = 0;
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i-1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema) 
						+ Mw[i][j]*Interpolacion(0, 0, Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema2) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], 0, 0, (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema2) 
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}

					//Resto
					else{
						aS[i][j] = Ds[i][j];
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i-1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)  
						+ Mw[i][j]*Interpolacion(0, 0, Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema2) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema) 
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}
				}

				//Columna derecha
				if (i==NX){
					aE[i][j] = 0;
					aW[i][j] = Dw[i][j];		
					//Esquina inferior
					if(j==1){
						aS[i][j] = 0;
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i+1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], 0, 0, (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema2) 
						 + Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						 + Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)  
						 + Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
					}
					//Esquina superior
					else if(j==NY){
						aS[i][j] = Ds[i][j];
						aN[i][j] = 0;
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i+1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], 0, 0, (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema2)  
						+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], 0, 0, (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema2) 
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}
					//Resto
					else{
						aS[i][j] = Ds[i][j];
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i+1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], 0, 0, (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema2)  
						+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)  
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}
				}
				//Fila inferior
				if (j==1 && i > 1 && i < NX){
					aE[i][j] = De[i][j];
					aW[i][j] = Dw[i][j];		
					aS[i][j] = 0;
					aN[i][j] = Dn[i][j];
					aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]);
					bP[i][j] = + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)
					+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
					+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)   
					+ Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
				}
				//Fila superior
				if (j==NY && i > 1 && i < NX){
					aE[i][j] = De[i][j];
					aW[i][j] = Dw[i][j];		
					aS[i][j] = Ds[i][j];
					aN[i][j] = 0;
					aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]);
					bP[i][j] = + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)
					+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
					+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], 0, 0, (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema2) 
					+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
				}
			}
		}
	}
	else{
		for (i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){
				//Nodos centrales
				if (i != 1 && i != NX && j != 1 && j != NY){
					aE[i][j] = De[i][j];
					aW[i][j] = Dw[i][j];	
					aS[i][j] = Ds[i][j];
					aN[i][j] = Dn[i][j];
					aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]);
					bP[i][j] = + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema) 
					+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
					+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema) 
					+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
				}
				//Columna izquierda
				if (i==1){
					aW[i][j] = 0;	
					aE[i][j] = De[i][j];
					//Esquina inferior
					if(j==1){
						aS[i][j] = 0;
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX) + 2*(DeltaX/DeltaY);
						bP[i][j] = 2*(DeltaY/DeltaX)*(1 - tanh(Alpha)) + 2*(DeltaX/DeltaY)*Phi[i][j-1] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema) 
						+ Mw[i][j]*Interpolacion(0, 0, Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema2) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema) 
						+ Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
					}
				
					//Esquina superior
					else if(j==NY){
						aS[i][j] = Ds[i][j];
						aN[i][j] = 0;
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX) + 2*(DeltaX/DeltaY);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i-1][j] + 2*(DeltaX/DeltaY)*Phi[i][j+1] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema) 
						+ Mw[i][j]*Interpolacion(0, 0, Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema2) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], 0, 0, (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema2) 
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}

					//Resto
					else{
						aS[i][j] = Ds[i][j];
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i-1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)  
						+ Mw[i][j]*Interpolacion(0, 0, Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema2) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema) 
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}
				}

				//Columna derecha
				if (i==NX){
					aE[i][j] = 0;
					aW[i][j] = Dw[i][j];		
					//Esquina inferior
					if(j==1){
						aS[i][j] = 0;
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i+1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], 0, 0, (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema2) 
						 + Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						 + Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)  
						 + Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
					}
					//Esquina superior
					else if(j==NY){
						aS[i][j] = Ds[i][j];
						aN[i][j] = 0;
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX) + 2*(DeltaX/DeltaY);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i+1][j] + 2*(DeltaX/DeltaY)*Phi[i][j+1] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], 0, 0, (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema2)  
						+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], 0, 0, (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema2) 
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}
					//Resto
					else{
						aS[i][j] = Ds[i][j];
						aN[i][j] = Dn[i][j];
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaY/DeltaX);
						bP[i][j] = 2*(DeltaY/DeltaX)*Phi[i+1][j] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], 0, 0, (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema2)  
						+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)  
						+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
					}
				}
				//Fila inferior
				if (j==1 && i > 1 && i < NX){
					aE[i][j] = De[i][j];
					aW[i][j] = Dw[i][j];		
					aS[i][j] = 0;
					aN[i][j] = Dn[i][j];
					if (Coordenadas[i][j][1] <= 0){
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaX/DeltaY);
						bP[i][j] = 2*(DeltaX/DeltaY)*Phi[i][j-1] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)
						+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)   
						+ Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
				
					}
					else{
						aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]);
						bP[i][j] = Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)
						+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
						+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], Phi[i][j+2], Coordenadas[i][j+2][2], (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema)   
						+ Ms[i][j]*Interpolacion(0, 0, Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema2);
				
					}
					
				}
				//Fila superior
				if (j==NY && i > 1 && i < NX){
					aE[i][j] = De[i][j];
					aW[i][j] = Dw[i][j];		
					aS[i][j] = Ds[i][j];
					aN[i][j] = 0;
					aP[i][j] = aE[i][j] + aS[i][j] + aN[i][j] + aW[i][j] + (Mw[i][j] + Me[i][j] + Ms[i][j] + Mn[i][j]) + 2*(DeltaX/DeltaY);
					bP[i][j] = 2*(DeltaX/DeltaY)*Phi[i][j+1] + Me[i][j]*Interpolacion(Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], Phi[i+2][j], Coordenadas[i+2][j][1], (Coordenadas[i][j][1] + Coordenadas[i+1][j][1])/2, Me[i][j], Esquema)
					+ Mw[i][j]*Interpolacion(Phi[i-2][j], Coordenadas[i-2][j][1], Phi[i-1][j], Coordenadas[i-1][j][1], Phi[i][j], Coordenadas[i][j][1], Phi[i+1][j], Coordenadas[i+1][j][1], (Coordenadas[i][j][1] + Coordenadas[i-1][j][1])/2, Mw[i][j], Esquema) 
					+ Mn[i][j]*Interpolacion(Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], 0, 0, (Coordenadas[i][j+1][2] + Coordenadas[i][j][2])/2, Mn[i][j], Esquema2) 
					+ Ms[i][j]*Interpolacion(Phi[i][j-2], Coordenadas[i][j-2][2], Phi[i][j-1], Coordenadas[i][j-1][2], Phi[i][j], Coordenadas[i][j][2], Phi[i][j+1], Coordenadas[i][j+1][2], (Coordenadas[i][j-1][2] + Coordenadas[i][j][2])/2, Ms[i][j], Esquema);
				}
			}
		}
	}
			//Calculo de los valores
		for (i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){

				Phi[i][j] = (aW[i][j]*Phi[i-1][j] + aE[i][j]*Phi[i+1][j] + aS[i][j]*Phi[i][j-1] + aN[i][j]*Phi[i][j+1] + bP[i][j])/aP[i][j];

			}
		}

		//Comprobacion de la convergencia
		MaxDif = 0;
		for (i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){

				if(abs(Phi[i][j] - PhiS[i][j]) >= MaxDif){
					MaxDif = abs(Phi[i][j] - PhiS[i][j]);
				}
			}
		}

		for (i=1;i<=NX;i++){
			for(j=1;j<=NY;j++){

				PhiS[i][j] = Phi[i][j];
			}
		}
		//cout<<MaxDif<<endl;
		NITER = NITER + 1;
		}
	//Fin while

	cout<<"Número de iteraciones: "<<NITER<<endl;
	cout<<endl;
	double P;
	P = Rho*Uo*X;
	if(Problema != 4){

		cout<<"Solución numérica"<<endl;
		for(j=NY+1;j>=0;j--){
			for(i=0;i<=NX+1;i++){
				
				cout<<Phi[i][j]<<", ";
			}
			cout<<"; "<<endl;
		}

		double PhiAnalitica[NX+1][NY+1];
		
		if(Problema == 1){
			for (i=0;i<=NX+1;i++){
				for(j=0;j<=NY+1;j++){
					
					PhiAnalitica[i][j] = PHI0 + (PHIL - PHI0)*((exp((P*Coordenadas[i][j][1])/X) - 1)/(exp(P) - 1));
				}
			}
		}
		else if(Problema == 2){
			for (i=0;i<=NX+1;i++){
				for(j=0;j<=NY+1;j++){
					
					PhiAnalitica[i][j] = PHI0 + ((PHIL - PHI0)/X)*Coordenadas[i][j][1];
				}
			}
		}

		cout<<endl;
		cout<<"Solución analitica"<<endl;
		for(j=NY+1;j>=0;j--){
			for(i=0;i<=NX+1;i++){
				cout<<PhiAnalitica[i][j]<<", ";
			}
			cout<<endl;
		}
	}
	else{
		cout<<"Valor Phi:"<<endl;
		for(i=0;i<=NX+1;i++){
			if(Coordenadas[i][0][1] >= 0){
				cout<<Phi[i][1]<<", ";
			}
			
		}
		cout<<endl;
		cout<<"Posición X:"<<endl;
		for(i=0;i<=NX+1;i++){
			if(Coordenadas[i][0][1] >= 0){
			cout<<Coordenadas[i][0][1]<<", ";
			}
		}
		cout<<endl;
	
	}
	
	cout<<endl;
	cout<<"Numero de Peclet: "<<P<<endl;
	cout<<"Delta X: "<<DeltaX<<endl;
	cout<<"Delta Y: "<<DeltaY<<endl;
	cout<<"Esquemas: "<<Esquema<<" + "<<Esquema2<<endl;
	cout<<"Valor convergencia: "<<MaxDif<<endl;

	
	return 0;
}