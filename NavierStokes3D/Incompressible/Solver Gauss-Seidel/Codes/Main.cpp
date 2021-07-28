#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include "HeaderCodes/Memory.h"
#include "HeaderCodes/ReadData.h"
#include "HeaderCodes/ParPro.h"
#include "HeaderCodes/Mesher.h"
#include "HeaderCodes/PostProcessing.h"
#include "HeaderCodes/Solver.h"

using namespace std;

#define DIRECTORIO "/home_nobck/sergiogus/Incompressible3D/"

int main(int argc, char* argv[]){

int i;

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1, DIRECTORIO);
R1.ReadInputs();

ParPro MPI1(R1);
MPI1.Execute();

Mesher MESH(M1, R1, MPI1, DIRECTORIO);
MESH.ExecuteMesher(M1, MPI1);

MPI_Barrier(MPI_COMM_WORLD);

PostProcessing POST1(M1, R1, MESH, DIRECTORIO);

Solver S1(M1, R1, MPI1, MESH, POST1, DIRECTORIO);

S1.ExecuteSolver(M1, MPI1, MESH, POST1);

MPI_Finalize();

return 0;

}
