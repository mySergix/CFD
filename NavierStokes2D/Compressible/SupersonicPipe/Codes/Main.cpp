#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Memory.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/ReadData.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Geometry.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Mesher.h"
#include "/home/sergio/Desktop/TFG/RegimenTurbulento/Codes/HeaderCodes/Solver.h"

using namespace std;

int main(){

Memory M1;

ReadData R1(M1);

R1.ReadInputs();

Geometry G1(M1, R1);

//G1.NasaLimits(M1);

G1.GeometryNasa(M1);

Mesher MESH(M1, R1, G1);

MESH.ExecuteMesher(M1, G1);

Solver S1(M1, R1, MESH);

S1.ExecuteSolver(M1, R1, MESH);

return 0;

}
