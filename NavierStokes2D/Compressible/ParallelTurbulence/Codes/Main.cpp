#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Memory.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ReadData.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/ParPro.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Geometry.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Mesher.h"
#include "/home_nobck/sergiogus/ParallelTurbulence/Codes/HeaderCodes/Solver.h"
using namespace std;

#define DIRECTORIO "/home_nobck/sergiogus/ParallelTurbulence/"

#define G(i,j,dim) (((j) + (i)*NR) + NX*NR*(dim)) //Global Index
#define LNH(i,j,dim) (((j) + ((i) - Ix + 1)*NR)) //Local Index No Halo + NY*(Fx - Ix + 2*Halo)*dim
#define LSH(i,j,dim) (((j) + ((i) - Ix)*NY) + NY*(Fx-Ix + 2*Halo)*dim) //Local Index Si Halo

int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs();

ParPro MPI1(R1);
MPI1.Execute();

Geometry G1(M1, R1);

G1.GeometryNasa(M1);

Mesher MESH(M1, R1, G1, MPI1);
MESH.ExecuteMesher(M1, G1, MPI1);

printf("Proceso %d, Initial X: %d, Final X: %d \n", MESH.Rank, MESH.Ix, MESH.Fx);

Solver S1(M1, R1, MPI1, MESH);
S1.ExecuteSolver(M1, MPI1, MESH);

MPI_Finalize();
return 0;

}
