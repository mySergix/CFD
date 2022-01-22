//Cálculo de las conductividades térmicas en las paredes de los volúmenes de control
void Solver::Get_muWalls(Mesher MESH){
int i, j;
double Predict;

	//Nodos U
	for(j = 0; j < NRad; j++){
		muWallsMU[0][j] = muLeft[j];
		
		Predict = Interpolacion(MESH.MP[0][j][0], Upres[0][j], MESH.MP[1][j][0], Upres[1][j], MESH.MU[1][j][0]);
		muWallsMU[1][j] = EsquemaConvectivo(MESH.MU[0][j][0], muLeft[j], MESH.MP[0][j][0], muTotal[0][j], MESH.MP[1][j][0], muTotal[1][j], MESH.MP[2][j][0], muTotal[2][j], MESH.MU[1][j][0], Predict, EsquemaAmplio);

		
		muWallsMU[NA][j] = muRight[j];
	
		Predict = Interpolacion(MESH.MP[NA-2][j][0], Upres[NA-2][j], MESH.MP[NA-1][j][0], Upres[NA-1][j], MESH.MU[NA-1][j][0]);
		muWallsMU[NA-1][j] = EsquemaConvectivo(MESH.MP[NA-3][j][0], muTotal[NA-3][j], MESH.MP[NA-2][j][0], muTotal[NA-2][j], MESH.MP[NA-1][j][0], muTotal[NA-1][j], MESH.MU[NA][j][0], muRight[j], MESH.MU[NA-1][j][0], Predict, EsquemaAmplio);

		for(i = 2; i < NA-1; i++){
			Predict = Interpolacion(MESH.MP[i-1][j][0], Upres[i-1][j], MESH.MP[i][j][0], Upres[i][j], MESH.MU[i][j][0]);
			muWallsMU[i][j] = EsquemaConvectivo(MESH.MP[i-2][j][0], muTotal[i-2][j], MESH.MP[i-1][j][0], muTotal[i-1][j], MESH.MP[i][j][0], muTotal[i][j], MESH.MP[i+1][j][0], muTotal[i+1][j], MESH.MU[i][j][0], Predict, EsquemaAmplio);
		}

	}

	//Nodos R
	for(i = 0; i < NA; i++){
		muWallsMR[i][0] = muDown[i];

		Predict = Interpolacion(MESH.MP[i][0][1], Vpres[i][0], MESH.MP[i][1][1], Vpres[i][1], MESH.MR[i][1][1]);
		muWallsMR[i][1] = EsquemaConvectivo(MESH.MR[i][0][1], muDown[i], MESH.MP[i][0][1], muTotal[i][0], MESH.MP[i][1][1], muTotal[i][1], MESH.MP[i][2][1], muTotal[i][2], MESH.MR[i][1][1], Predict, EsquemaAmplio);
		
		Predict = Interpolacion(MESH.MP[i][NRad-2][1], Vpres[i][NRad-2], MESH.MP[i][NRad-1][1], Vpres[i][NRad-1], MESH.MR[i][NRad-1][1]);
		muWallsMR[i][NRad-1] = EsquemaConvectivo(MESH.MP[i][NRad-3][1], muTotal[i][NRad-3], MESH.MP[i][NRad-2][1], muTotal[i][NRad-2], MESH.MP[i][NRad-1][1], muTotal[i][NRad-1], MESH.MR[i][NRad][1], muUp[i], MESH.MR[i][NRad-1][1], Predict, EsquemaAmplio);
		
		muWallsMR[i][NRad] = muUp[i];

		for(j = 2; j < NRad-1; j++){
			Predict = Interpolacion(MESH.MP[i][j-1][1], Vpres[i][j-1], MESH.MP[i][j][1], Vpres[i][j], MESH.MR[i][j][1]);
			muWallsMR[i][j] = EsquemaConvectivo(MESH.MP[i][j-2][1], muTotal[i][j-2], MESH.MP[i][j-1][1], muTotal[i][j-1], MESH.MP[i][j][1], muTotal[i][j], MESH.MP[i][j+1][1], muTotal[i][j+1], MESH.MR[i][j][1], Predict, EsquemaAmplio);
		}

	}

}