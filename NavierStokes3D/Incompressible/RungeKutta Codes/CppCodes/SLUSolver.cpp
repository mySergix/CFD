#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>
#include "umfpack.h"

using namespace std;

#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Mesher.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Solver.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/SLUSolver.h"


#define PI 3.141592653589793
#define ZERO 1e-12

#define DIRECTORIO "/home/sergiogus/Desktop/CFD/Incompressible3D/"

//Local Index P Mesh
#define LPL(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i)+HP) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Left Core
#define LPC(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Center Cores
#define LPR(i,j,k,dim) ((NY+2*HP)*(NZ+2*HP))*((i) - Ix + Halo) + (((j)+HP) + ((k)+HP)*(NY+2*HP)) //Right Core

//Coeficientes A
#define LAG(i,j,k,dim) ((NY*NZ)*((i))) + ((j) + (k)*NY)
#define LAL(i,j,k,dim) ((NY*NZ)*((i) - Ix)) + ((j) + (k)*NY)

#define I(row) row/(NY*NZ)
#define K(row,i) (( row - (((i) + 1)*(NY*NZ)) )/NY) - 1 
#define J(row,i) (( row - (((i) + 1)*(NY*NZ)) )%NY) - 1 

SLUSolver::SLUSolver(Solver S1){

    NX = S1.NX;
    NY = S1.NY;
    NZ = S1.NZ;

    b = new double [NX*NY*NZ];
    m_Ap = new int [N+1];

};

void SLUSolver::Build_Ap_vector(int N, Solver S1){
int i, j, k;
int n;

    for(i = 0; i < N + 1; i++){
        Ap[i] = 0;
    }

    for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){
            for(k = 0; k < NZ; k++){
                
                //Coeficiente ap
                if(abs(S1.ap[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j,k,0);
                    Ap[n]++;
                }

                //Coeficiente aw
                if(abs(S1.aw[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i - 1,j,k,0);
                    Ap[n]++;
                }

                //Coeficiente ae
                if(abs(S1.ae[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i + 1,j,k,0);
                    Ap[n]++;
                }   

                //Coeficiente as
                if(abs(S1.as[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j - 1,k,0);
                    Ap[n]++;
                }

                //Coeficiente an
                if(abs(S1.an[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j + 1,k,0);
                    Ap[n]++;
                }

            }
        }
    }

    for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){

            //Centro
            for(k = 1; k < NZ - 1; k++){

                //Coeficiente ah
                if(abs(S1.ah[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j,k - 1,0);
                    Ap[n]++;
                }

                //Coeficiente at
                if(abs(S1.at[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j,k + 1,0);
                    Ap[n]++;
                }

            }

            //Parte Here
            if(abs(S1.ah[LAG(i,j,0,0)]) > ZERO){
                n = LAG(i,j,NZ - 1,0);
                Ap[n]++;
            }
            if(abs(S1.at[LAG(i,j,0,0)]) > ZERO){
                n = LAG(i,j,1,0);
                Ap[n]++;
            }


            //Parte There
            if(abs(S1.ah[LAG(i,j,NZ - 1,0)]) > ZERO){
                n = LAG(i,j,NZ - 2,0);
                Ap[n]++;
            }

            //Coeficiente at
            if(abs(S1.at[LAG(i,j,NZ - 1,0)]) > ZERO){
                n = LAG(i,j,0,0);
                Ap[n]++;
            }

        }
    }

 
    //Suma de los Ap ascendentemente
    for(i = 2; i < N + 1; i++){
        Ap[i] += Ap[i-1];
    } 

}

void SLUSolver::Build_Ai_Ax_vectors(int N, Solver S1){
int i, j, k;
int n;
int contador = 0;

    for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){
            for(k = 0; k < NZ; k++){

                //Coeficiente ap
                if(abs(S1.ap[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j,k,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.ap[LAG(i,j,k,0)];
                    contador++;
                }

                //Coeficiente aw
                if(abs(S1.aw[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i - 1,j,k,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.aw[LAG(i,j,k,0)];
                    contador++;
                }

                //Coeficiente ae
                if(abs(S1.ae[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i + 1,j,k,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.ae[LAG(i,j,k,0)];
                    contador++;
                }   

                //Coeficiente as
                if(abs(S1.as[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j - 1,k,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.as[LAG(i,j,k,0)];
                    contador++;
                }

                //Coeficiente an
                if(abs(S1.an[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j + 1,k,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.an[LAG(i,j,k,0)];
                    contador++;
                }

            }
        }
    }

    for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){

            //Centro
            for(k = 1; k < NZ - 1; k++){

                //Coeficiente ah
                if(abs(S1.ah[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j,k - 1,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.ah[LAG(i,j,k,0)];
                    contador++;
                }

                //Coeficiente at
                if(abs(S1.at[LAG(i,j,k,0)]) > ZERO){
                    n = LAG(i,j,k + 1,0);
                    Ai[contador] = n;
                    Ax[contador] = S1.at[LAG(i,j,k,0)];
                    contador++;
                }

            }

            //Parte Here

            //Coeficiente ah
            if(abs(S1.ah[LAG(i,j,0,0)]) > ZERO){
                n = LAG(i,j,NZ - 1,0);
                Ai[contador] = n;
                Ax[contador] = S1.ah[LAG(i,j,0,0)];
                contador++;
            }

            //Coeficiente at
            if(abs(S1.at[LAG(i,j,0,0)]) > ZERO){
                n = LAG(i,j,1,0);
                Ai[contador] = n;
                Ax[contador] = S1.at[LAG(i,j,0,0)];
                contador++;
            }

            //Parte There

            //Coeficiente ah
            if(abs(S1.ah[LAG(i,j,NZ - 1,0)]) > ZERO){
                n = LAG(i,j,NZ - 2,0);
                Ai[contador] = n;
                Ax[contador] = S1.ah[LAG(i,j,NZ - 1,0)];
                contador++;
            }

            //Coeficiente at
            if(abs(S1.at[LAG(i,j,NZ - 1,0)]) > ZERO){
                n = LAG(i,j,0,0);
                Ai[contador] = n;
                Ax[contador] = S1.at[LAG(i,j,NZ - 1,0)];
                contador++;
            }

        }
    }

    //Comprobación
    if(contador!=Ap[N]){
	    cout<<"Error: Solver: size(Ap)!=size(Ai/Ax). Program failed."<<'\n';
	}


}

void SLUSolver::Build_BP_vector(double *BPmatrix){
int i, j, k;

    for(i = 0; i < NX; i++){
        for(k = 0; k < NZ; k++){
            for(j = 0; j < NY; j++){

                b[LAG(i,j,k,0)] = BPmatrix[LAG(i,j,k,0)];

            }
        }
    }
    
}

void SLUSolver::Get_LaplacianMatrixA(Solver S1, double *& null, void *& Symbolic, void *& Numeric){
int N;

    N = (NX)*(NY)*(NZ);
    Ap = new int [N + 1];

    //Construccion vector Ap
    Build_Ap_vector(N, S1);

    Ai = new int [Ap[N]];
    Ai = new double [Ap[N]];

    //Contruccion vectores Ai y Ax
    Build_Ai_Ax_vectors(N, S1);

    (void) umfpack_di_symbolic (N, N, m_Ap, m_Ai, m_Ax, &Symbolic, null, null);
	(void) umfpack_di_numeric (m_Ap, m_Ai, m_Ax, Symbolic, &Numeric, null, null);
	umfpack_di_free_symbolic (&Symbolic);

}

void SLUSolver::Get_Pressure(double *x, void *& Numeric, double *& null){

    umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);

}