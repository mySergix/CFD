#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Mesher.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/Solver.h"
#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/SLUSolver.h"

//#include "/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/HeaderCodes/PostProcessing.h"

using namespace std;

#define DIRECTORIO "/home/sergiogus/Desktop/CFD/Incompressible3D/"

int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs();

ParPro MPI1(R1);
MPI1.Execute();

Mesher MESH(M1, R1, MPI1);
MESH.ExecuteMesher(M1, MPI1);

MPI_Barrier(MPI_COMM_WORLD);

Solver S1(M1, R1, MPI1, MESH);
S1.ExecuteSolver(M1, MPI1, MESH);

//PostProc P1(M1, R1, MESH, S1);
//P1.ExecutePostProcessing(S1, MESH, M1, MPI1);

//P1.DeleteEverything(S1, MESH, R1);

MPI_Finalize();

return 0;

}
