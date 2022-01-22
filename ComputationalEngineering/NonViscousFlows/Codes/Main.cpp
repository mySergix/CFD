#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Memory.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ReadData.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/ParPro.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Mesher.h"
#include "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/HeaderCodes/Solver.h"

using namespace std;

#define DIRECTORIO "/home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/"

int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
cout << endl;
int N = 11;
int i = 0;

Memory M1;

ReadData R1(M1, N);
R1.ReadInputs();

    ParPro MPI1(R1, i);
    MPI1.Execute();

    Mesher MESH(M1, R1, MPI1, i);
    MESH.ExecuteMesher(M1, MPI1);

    printf("Proceso %d, Initial X: %d, Final X: %d \n", MESH.Rank, MESH.Ix, MESH.Fx);

    Solver S1(M1, R1, MPI1, MESH, i);
    S1.ExecuteSolver(M1, R1, MPI1, MESH);

MPI_Finalize();
return 0;

}
