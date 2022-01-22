//Problemas de Navier-Stokes

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

#include "/home/sergiogus/Desktop/CTTC/Turbulencia/Codigos/ArchivosHEADER/ReadTXT.h"
#include "/home/sergiogus/Desktop/CTTC/Turbulencia/Codigos/ArchivosHEADER/Mesh.h"
#include "/home_nobck/sergiogus/libSLU/SLU.h"
#include "/home_nobck/sergiogus/SuiteSparse/UMFPACK/Include/umfpack.h"
#include "/home/sergiogus/Desktop/CTTC/Turbulencia/Codigos/ArchivosHEADER/Solver.h"
#include "/home/sergiogus/Desktop/CTTC/Turbulencia/Codigos/ArchivosHEADER/EjecucionSolver.h"

using namespace std;

//#define PI 3.141592653589793

int main(){

//Data reading class
ReadTXT R1;

//Read all the data from TXT files
R1.ReadFromTXT();


//Mesh creator class
Mallador M1(R1);

//Creation of the meshes
M1.EjecucionMesh();

//Matrix solver class
SLU C1;

//Solver class
CalculoSolver S1(R1, C1, M1);

//Problem resolution
S1.EjecucionSolver(C1, M1);

return 0;

}



