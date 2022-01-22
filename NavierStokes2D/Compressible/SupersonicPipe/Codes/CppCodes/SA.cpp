//Cálculo del término 1 de la ecuación de Spalart-Allmaras
void Solver::SA_Term1(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino1[i][j] = -(1.0/MESH.VolMP[i][j])*(MESH.SupMP[i][j][0]*0.50*(vpresent[i][j] + vpresent[i-1][j])*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j])) + MESH.SupMP[i][j][1]*0.50*(vpresent[i][j] + vpresent[i+1][j])*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(MESH.AngleMU[i+1][j])) + MESH.SupMP[i][j][2]*0.50*(vpresent[i][j] + vpresent[i][j-1])*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j])) + MESH.SupMP[i][j][3]*0.50*(vpresent[i][j] + vpresent[i][j+1])*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1])));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino1[0][j] = -(1.0/MESH.VolMP[0][j])*(MESH.SupMP[0][j][0]*vleft[j]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j]) + Vleft[j]*sin(PI + MESH.AngleMU[0][j])) + MESH.SupMP[0][j][1]*0.50*(vpresent[0][j] + vpresent[1][j])*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(MESH.AngleMU[1][j])) + MESH.SupMP[0][j][2]*0.50*(vpresent[0][j] + vpresent[0][j-1])*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j])) + MESH.SupMP[0][j][3]*0.50*(vpresent[0][j] + vpresent[0][j+1])*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1])));

		//Parte derecha
		SA_Termino1[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(MESH.SupMP[NA-1][j][0]*0.50*(vpresent[NA-1][j] + vpresent[NA-2][j])*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j])) + MESH.SupMP[NA-1][j][1]*vright[j]*(Uright[j]*cos(MESH.AngleMU[NA][j]) + Vright[j]*sin(MESH.AngleMU[NA][j])) + MESH.SupMP[NA-1][j][2]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1])*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])) + MESH.SupMP[NA-1][j][3]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1])*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])));

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino1[i][0] = -(1.0/MESH.VolMP[i][0])*(MESH.SupMP[i][0][0]*0.50*(vpresent[i][0] + vpresent[i-1][0])*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0])) + MESH.SupMP[i][0][1]*0.50*(vpresent[i][0] + vpresent[i+1][0])*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(MESH.AngleMU[i+1][0])) + MESH.SupMP[i][0][2]*vdown[i]*(Udown[i]*cos(1.50*PI + MESH.AngleMR[i][j]) + Vdown[i]*sin(1.50*PI + MESH.AngleMR[i][0])) + MESH.SupMP[i][0][3]*0.50*(vpresent[i][0] + vpresent[i][1])*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1])));

		//Parte arriba
		SA_Termino1[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(MESH.SupMP[i][NRad-1][0]*0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1])*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1])) + MESH.SupMP[i][NRad-1][1]*0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1])*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1])) + MESH.SupMP[i][NRad-1][2]*0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2])*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])) + MESH.SupMP[i][NRad-1][3]*vup[i]*(Uup[i]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + Vup[i]*sin(0.50*PI + MESH.AngleMR[i][NRad])));

	}

	//Esquina abajo izquierda
	SA_Termino1[0][0] = -(1.0/MESH.VolMP[0][0])*(MESH.SupMP[0][0][0]*vleft[0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0]) + Vleft[0]*sin(PI + MESH.AngleMU[0][0])) + MESH.SupMP[0][0][1]*0.50*(vpresent[0][0] + vpresent[1][0])*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(MESH.AngleMU[1][0])) + MESH.SupMP[0][0][2]*vdown[0]*(Udown[0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vdown[0]*sin(1.50*PI + MESH.AngleMR[0][0])) + MESH.SupMP[0][0][3]*0.50*(vpresent[0][0] + vpresent[0][1])*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1])));

	//Esquina arriba izquierda
	SA_Termino1[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(MESH.SupMP[0][NRad-1][0]*vleft[NRad-1]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1])) + MESH.SupMP[0][NRad-1][1]*0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1])*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1])) + MESH.SupMP[0][NRad-1][2]*0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2])*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])) + MESH.SupMP[0][NRad-1][3]*vup[0]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad])));

	//Esquina abajo derecha
	SA_Termino1[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(MESH.SupMP[NA-1][0][0]*0.50*(vpresent[NA-1][0] + vpresent[NA-2][0])*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0])) + MESH.SupMP[NA-1][0][1]*vright[0]*(Uright[0]*cos(MESH.AngleMU[NA][0]) + Vright[0]*sin(MESH.AngleMU[NA][0])) + MESH.SupMP[NA-1][0][2]*vdown[NA-1]*(Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vdown[NA-1]*sin(1.50*PI + MESH.AngleMR[NA-1][0])) + MESH.SupMP[NA-1][0][3]*0.50*(vpresent[NA-1][0] + vpresent[NA-1][1])*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1])));

	//Esquina arriba derecha
	SA_Termino1[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1])*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1])) + MESH.SupMP[NA-1][NRad-1][1]*vright[NRad-1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1])) + MESH.SupMP[NA-1][NRad-1][2]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2])*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])) + MESH.SupMP[NA-1][NRad-1][3]*vup[NA-1]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])));

}


//Cálculo del término 2 de la ecuación de Spalart-Allmaras
void Solver::SA_Term2(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino2[i][j] = (1.0/(MESH.VolMP[i][j]*Sigma))*(MESH.SupMP[i][j][0]*(0.50*(muBase[i][j]/(RhoPres[i][j] + 1e-10) + muBase[i-1][j]/(RhoPres[i-1][j] + 1e-10)) + 0.50*(vpresent[i][j] + vpresent[i-1][j]))*((vpresent[i][j] - vpresent[i-1][j])/MESH.DeltasMU[i][j][0])*cos(PI + MESH.AngleMU[i][j]) + MESH.SupMP[i][j][1]*(0.50*(muBase[i][j]/(RhoPres[i][j]+1e-10) + muBase[i+1][j]/(RhoPres[i+1][j]+1e-10)) + 0.50*(vpresent[i][j] + vpresent[i+1][j]))*((vpresent[i+1][j] - vpresent[i][j])/MESH.DeltasMU[i+1][j][0])*cos(MESH.AngleMU[i+1][j]) + MESH.SupMP[i][j][2]*(0.50*(muBase[i][j]/(RhoPres[i][j]+1e-10) + muBase[i][j-1]/(RhoPres[i][j-1]+1e-10)) + 0.50*(vpresent[i][j] + vpresent[i][j-1]))*((vpresent[i][j] - vpresent[i][j-1])/MESH.DeltasMR[i][j][1])*cos(1.50*PI + MESH.AngleMR[i][j]) + MESH.SupMP[i][j][3]*(0.50*(muBase[i][j]/(RhoPres[i][j]+1e-10) + muBase[i][j+1]/(RhoPres[i][j+1]+1e-10)) + 0.50*(vpresent[i][j] + vpresent[i][j+1]))*((vpresent[i][j+1] - vpresent[i][j])/MESH.DeltasMR[i][j+1][1])*cos(0.50*PI + MESH.AngleMR[i][j+1]));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino2[0][j] = (1.0/(MESH.VolMP[0][j]*Sigma))*(MESH.SupMP[0][j][0]*(muWallsMU[0][j]/(RhoLeft[j]+1e-10) + vleft[j])*((vpresent[0][j] - vleft[j])/MESH.DeltasMU[0][j][0])*cos(PI + MESH.AngleMU[0][j]) + MESH.SupMP[0][j][1]*(0.50*(muBase[0][j]/(RhoPres[0][j]+1e-10) + muBase[1][j]/(RhoPres[1][j]+1e-10)) + 0.50*(vpresent[0][j] + vpresent[1][j]))*((vpresent[1][j] - vpresent[0][j])/MESH.DeltasMU[1][j][0])*cos(MESH.AngleMU[1][j]) + MESH.SupMP[0][j][2]*(0.50*(muBase[0][j]/(RhoPres[0][j]+1e-10) + muBase[0][j-1]/(RhoPres[0][j-1]+1e-10)) + 0.50*(vpresent[0][j] + vpresent[0][j-1]))*((vpresent[0][j] - vpresent[0][j-1])/MESH.DeltasMR[0][j][1])*cos(1.50*PI + MESH.AngleMR[0][j]) + MESH.SupMP[0][j][3]*(0.50*(muBase[0][j]/(RhoPres[0][j]+1e-10) + muBase[0][j+1]/(RhoPres[0][j+1]+1e-10)) + 0.50*(vpresent[0][j] + vpresent[0][j+1]))*((vpresent[0][j+1] - vpresent[0][j])/MESH.DeltasMR[0][j+1][1])*cos(0.50*PI + MESH.AngleMR[0][j+1]));

		//Parte derecha
		SA_Termino2[NA-1][j] = (1.0/(MESH.VolMP[NA-1][j]*Sigma))*(MESH.SupMP[NA-1][j][0]*(0.50*(muBase[NA-1][j]/(RhoPres[NA-1][j]+1e-10) + muBase[NA-2][j]/(RhoPres[NA-2][j]+1e-10)) + 0.50*(vpresent[NA-1][j] + vpresent[NA-2][j]))*((vpresent[NA-1][j] - vpresent[NA-2][j])/MESH.DeltasMU[NA-1][j][0])*cos(PI + MESH.AngleMU[NA-1][j]) + MESH.SupMP[NA-1][j][1]*(muWallsMU[NA][j]/(RhoRight[j]+1e-10) + vright[j])*((vright[j] - vpresent[NA-1][j])/MESH.DeltasMU[NA][j][0])*cos(MESH.AngleMU[NA][j]) + MESH.SupMP[NA-1][j][2]*(0.50*(muBase[NA-1][j]/(RhoPres[NA-1][j]+1e-10) + muBase[NA-1][j-1]/(RhoPres[NA-1][j-1]+1e-10)) + 0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1]))*((vpresent[NA-1][j] - vpresent[NA-1][j-1])/MESH.DeltasMR[NA-1][j][1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + MESH.SupMP[NA-1][j][3]*(0.50*(muBase[NA-1][j]/(RhoPres[NA-1][j]+1e-10) + muBase[NA-1][j+1]/(RhoPres[NA-1][j+1]+1e-10)) + 0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1]))*((vpresent[NA-1][j+1] - vpresent[NA-1][j])/MESH.DeltasMR[NA-1][j+1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]));

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino2[i][0] = (1.0/(MESH.VolMP[i][0]*Sigma))*(MESH.SupMP[i][0][0]*(0.50*(muBase[i][0]/(RhoPres[i][0]+1e-10) + muBase[i-1][0]/(RhoPres[i-1][0]+1e-10)) + 0.50*(vpresent[i][0] + vpresent[i-1][0]))*((vpresent[i][0] - vpresent[i-1][0])/MESH.DeltasMU[i][0][0])*cos(PI + MESH.AngleMU[i][0]) + MESH.SupMP[i][0][1]*(0.50*(muBase[i][0]/(RhoPres[i][0]+1e-10) + muBase[i+1][0]/(RhoPres[i+1][0]+1e-10)) + 0.50*(vpresent[i][0] + vpresent[i+1][0]))*((vpresent[i+1][0] - vpresent[i][0])/MESH.DeltasMU[i+1][0][0])*cos(MESH.AngleMU[i+1][0]) + MESH.SupMP[i][0][2]*(muWallsMR[i][0]/(RhoDown[i]+1e-10) + vdown[i])*((vpresent[i][0] - vdown[i])/MESH.DeltasMR[i][0][1])*cos(1.50*PI + MESH.AngleMR[i][0]) + MESH.SupMP[i][0][3]*(0.50*(muBase[i][0]/(RhoPres[i][0]+1e-10) + muBase[i][1]/(RhoPres[i][1]+1e-10)) + 0.50*(vpresent[i][0] + vpresent[i][1]))*((vpresent[i][1] - vpresent[i][0])/MESH.DeltasMR[i][1][1])*cos(0.50*PI + MESH.AngleMR[i][1]));

		//Parte arriba
		SA_Termino2[i][NRad-1] = (1.0/(MESH.VolMP[i][NRad-1]*Sigma))*(MESH.SupMP[i][NRad-1][0]*(0.50*(muBase[i][NRad-1]/(RhoPres[i][NRad-1]+1e-10) + muBase[i-1][NRad-1]/(RhoPres[i-1][NRad-1]+1e-10)) + 0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1]))*((vpresent[i][NRad-1] - vpresent[i-1][NRad-1])/MESH.DeltasMU[i][NRad-1][0])*cos(PI + MESH.AngleMU[i][NRad-1]) + MESH.SupMP[i][NRad-1][1]*(0.50*(muBase[i][NRad-1]/(RhoPres[i][NRad-1]+1e-10) + muBase[i+1][NRad-1]/(RhoPres[i+1][NRad-1]+1e-10)) + 0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1]))*((vpresent[i+1][NRad-1] - vpresent[i][NRad-1])/MESH.DeltasMU[i+1][NRad-1][0])*cos(MESH.AngleMU[i+1][NRad-1]) + MESH.SupMP[i][NRad-1][2]*(0.50*(muBase[i][NRad-1]/(RhoPres[i][NRad-1]+1e-10) + muBase[i][NRad-2]/(RhoPres[i][NRad-2]+1e-10)) + 0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2]))*((vpresent[i][NRad-1] - vpresent[i][NRad-2])/MESH.DeltasMR[i][NRad-1][1])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + MESH.SupMP[i][NRad-1][3]*(muWallsMR[i][NRad-1]/(RhoUp[i]+1e-10) + vup[i])*((vup[i] - vpresent[i][j])/MESH.DeltasMR[i][NRad][1])*cos(0.50*PI + MESH.AngleMR[i][NRad]));

	}

	//Esquina abajo izquierda
	SA_Termino2[0][0] = (1.0/(MESH.VolMP[0][0]*Sigma))*(MESH.SupMP[0][0][0]*(muWallsMU[0][0]/(RhoLeft[0]+1e-10) + vleft[0])*((vpresent[0][0] - vleft[0])/MESH.DeltasMU[0][0][0])*cos(PI + MESH.AngleMU[0][0]) + MESH.SupMP[0][0][1]*(0.50*(muBase[0][0]/(RhoPres[0][0]+1e-10) + muBase[1][0]/(RhoPres[1][0]+1e-10)) + 0.50*(vpresent[0][0] + vpresent[1][0]))*((vpresent[1][0] - vpresent[0][0])/MESH.DeltasMU[1][0][0])*cos(MESH.AngleMU[1][0]) + MESH.SupMP[0][0][2]*(muWallsMR[0][0]/(RhoDown[0]+1e-10) + vdown[0])*((vpresent[0][0] - vdown[0])/MESH.DeltasMR[0][0][1])*cos(1.50*PI + MESH.AngleMR[0][0]) + MESH.SupMP[0][0][3]*(0.50*(muBase[0][0]/(RhoPres[0][0]+1e-10) + muBase[0][1]/(RhoPres[0][1]+1e-10)) + 0.50*(vpresent[0][0] + vpresent[0][1]))*((vpresent[0][1] - vpresent[0][0])/MESH.DeltasMR[0][1][1])*cos(0.50*PI + MESH.AngleMR[0][1]));

	//Esquina arriba izquierda
	SA_Termino2[0][NRad-1] = (1.0/(MESH.VolMP[0][NRad-1]*Sigma))*(MESH.SupMP[0][NRad-1][0]*(muWallsMU[0][NRad-1]/(RhoLeft[NRad-1]+1e-10) + vleft[NRad-1])*((vpresent[0][j] - vleft[NRad-1])/MESH.DeltasMU[0][NRad-1][0])*cos(PI + MESH.AngleMU[0][NRad-1]) + MESH.SupMP[0][NRad-1][1]*(0.50*(muBase[0][NRad-1]/(RhoPres[0][NRad-1]+1e-10) + muBase[1][NRad-1]/(RhoPres[1][NRad-1]+1e-10)) + 0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1]))*((vpresent[1][NRad-1] - vpresent[0][NRad-1])/MESH.DeltasMU[1][NRad-1][0])*cos(MESH.AngleMU[1][NRad-1]) + MESH.SupMP[0][NRad-1][2]*(0.50*(muBase[0][NRad-1]/(RhoPres[0][NRad-1]+1e-10) + muBase[0][NRad-2]/(RhoPres[0][NRad-2]+1e-10)) + 0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2]))*((vpresent[0][NRad-1] - vpresent[0][NRad-2])/MESH.DeltasMR[0][NRad-1][1])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + MESH.SupMP[0][NRad-1][3]*(muWallsMR[0][NRad]/(RhoUp[0]+1e-10) + vup[0])*((vpresent[0][j+1] - vup[0])/MESH.DeltasMR[0][NRad][1])*cos(0.50*PI + MESH.AngleMR[0][NRad]));

	//Esquina abajo derecha
	SA_Termino2[NA-1][0] = (1.0/(MESH.VolMP[NA-1][0]*Sigma))*(MESH.SupMP[NA-1][0][0]*(0.50*(muBase[NA-1][0]/(RhoPres[NA-1][0]+1e-10) + muBase[NA-2][0]/(RhoPres[NA-2][0]+1e-10)) + 0.50*(vpresent[NA-1][0] + vpresent[NA-2][0]))*((vpresent[NA-1][0] - vpresent[NA-2][0])/MESH.DeltasMU[NA-1][0][0])*cos(PI + MESH.AngleMU[NA-1][0]) + MESH.SupMP[NA-1][0][1]*(muWallsMU[NA][0]/(RhoRight[0] + 1e-10) + vright[0])*((vright[0] - vpresent[NA-1][0])/MESH.DeltasMU[NA][0][0])*cos(MESH.AngleMU[NA][0]) + MESH.SupMP[NA-1][0][2]*(muWallsMR[NA-1][0]/(RhoDown[NA-1]+1e-10) + vdown[NA-1])*((vpresent[NA-1][0] - vdown[NA-1])/MESH.DeltasMR[NA-1][0][1])*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + MESH.SupMP[NA-1][0][3]*(0.50*(muBase[NA-1][0]/(RhoPres[NA-1][0]+1e-10) + muBase[NA-1][1]/(RhoPres[NA-1][1]+1e-10)) + 0.50*(vpresent[NA-1][0] + vpresent[NA-1][1]))*((vpresent[NA-1][1] - vpresent[NA-1][0])/MESH.DeltasMR[NA-1][1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]));
	
	//Esquina arriba derecha
	SA_Termino2[NA-1][NRad-1] = (1.0/(MESH.VolMP[NA-1][NRad-1]*Sigma))*(MESH.SupMP[NA-1][NRad-1][0]*(0.50*(muBase[NA-1][NRad-1]/(RhoPres[NA-1][NRad-1]+1e-10) + muBase[NA-2][NRad-1]/(RhoPres[NA-2][NRad-1]+1e-10)) + 0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1]))*((vpresent[NA-1][NRad-1] - vpresent[NA-2][NRad-1])/MESH.DeltasMU[NA-1][NRad-1][0])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + MESH.SupMP[NA-1][NRad-1][1]*(muWallsMU[NA][NRad-1]/(RhoRight[NRad-1]+1e-10) + vright[NRad-1])*((vright[NRad-1] - vpresent[NA-1][NRad-1])/MESH.DeltasMU[NA][NRad-1][0])*cos(MESH.AngleMU[NA][NRad-1]) + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(muBase[NA-1][NRad-1]/(RhoPres[NA-1][NRad-1]+1e-10) + muBase[NA-1][NRad-2]/(RhoPres[NA-1][NRad-2]+1e-10)) + 0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2]))*((vpresent[NA-1][NRad-1] - vpresent[NA-1][NRad-2])/MESH.DeltasMR[NA-1][NRad-1][1])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + MESH.SupMP[NA-1][NRad-1][3]*(muWallsMR[NA-1][NRad]/(RhoUp[NA-1]+1e-10) + vup[NA-1])*((vup[NA-1] - vpresent[NA-1][NRad-1])/MESH.DeltasMR[NA-1][NRad][1])*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]));

}


//Cálculo del término 3 de la ecuación de Spalart-Allmaras
void Solver::SA_Term3(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			SA_Termino3[i][j] = Cb1*(1.0 - ft2[i][j])*Smodel[i][j]*vpresent[i][j];
		}
	}

}


//Cálculo del término 4 de la ecuación de Spalart-Allmaras
void Solver::SA_Term4(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino4[i][j] = (Cb2/(Sigma*pow(MESH.VolMP[i][j],2.0)))*(pow(MESH.SupMP[i][j][0]*0.50*(vpresent[i][j] + vpresent[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + MESH.SupMP[i][j][1]*0.50*(vpresent[i][j] + vpresent[i+1][j])*cos(PI + MESH.AngleMU[i+1][j]),2.0) + pow(MESH.SupMP[i][j][2]*0.50*(vpresent[i][j] + vpresent[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + MESH.SupMP[i][j][3]*0.50*(vpresent[i][j] + vpresent[i][j+1])*cos(PI + MESH.AngleMR[i][j+1]),2.0));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino4[0][j] = (Cb2/(Sigma*pow(MESH.VolMP[0][j],2.0)))*(pow(MESH.SupMP[0][j][0]*vleft[j]*cos(PI + MESH.AngleMU[0][j]) + MESH.SupMP[0][j][1]*0.50*(vpresent[0][j] + vpresent[1][j])*cos(PI + MESH.AngleMU[1][j]),2.0) + pow(MESH.SupMP[0][j][2]*0.50*(vpresent[0][j] + vpresent[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + MESH.SupMP[0][j][3]*0.50*(vpresent[0][j] + vpresent[0][j+1])*cos(PI + MESH.AngleMR[i][1]),2.0));

		//Parte derecha
		SA_Termino4[NA-1][j] = (Cb2/(Sigma*pow(MESH.VolMP[NA-1][j],2.0)))*(pow(MESH.SupMP[NA-1][j][0]*0.50*(vpresent[NA-1][j] + vpresent[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + MESH.SupMP[NA-1][j][1]*vright[j]*cos(PI + MESH.AngleMU[NA][j]),2.0) + pow(MESH.SupMP[NA-1][j][2]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + MESH.SupMP[NA-1][j][3]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1])*cos(PI + MESH.AngleMR[NA-1][j+1]),2.0));
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino4[i][0] = (Cb2/(Sigma*pow(MESH.VolMP[i][0],2.0)))*(pow(MESH.SupMP[i][0][0]*0.50*(vpresent[i][0] + vpresent[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + MESH.SupMP[i][0][1]*0.50*(vpresent[i][0] + vpresent[i+1][0])*cos(PI + MESH.AngleMU[i+1][0]),2.0) + pow(MESH.SupMP[i][0][2]*vdown[i]*cos(1.50*PI + MESH.AngleMR[i][0]) + MESH.SupMP[i][0][3]*0.50*(vpresent[i][0] + vpresent[i][1])*cos(PI + MESH.AngleMR[i][1]),2.0));

		//Parte arriba
		SA_Termino4[i][NRad-1] = (Cb2/(Sigma*pow(MESH.VolMP[i][NRad-1],2.0)))*(pow(MESH.SupMP[i][NRad-1][0]*0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + MESH.SupMP[i][NRad-1][1]*0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1])*cos(PI + MESH.AngleMU[i+1][NRad-1]),2.0) + pow(MESH.SupMP[i][NRad-1][2]*0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + MESH.SupMP[i][NRad-1][3]*vup[i]*cos(PI + MESH.AngleMR[i][NRad]),2.0));
	}

	//Esquina abajo izquierda
	SA_Termino4[0][0] = (Cb2/(Sigma*pow(MESH.VolMP[0][0],2.0)))*(pow(MESH.SupMP[0][0][0]*vleft[0]*cos(PI + MESH.AngleMU[0][0]) + MESH.SupMP[0][0][1]*0.50*(vpresent[0][0] + vpresent[1][0])*cos(PI + MESH.AngleMU[1][0]),2.0) + pow(MESH.SupMP[0][0][2]*vdown[0]*cos(1.50*PI + MESH.AngleMR[0][0]) + MESH.SupMP[0][0][3]*0.50*(vpresent[0][0] + vpresent[0][1])*cos(PI + MESH.AngleMR[0][1]),2.0));

	//Esquina arriba izquierda
	SA_Termino4[0][NRad-1] = (Cb2/(Sigma*pow(MESH.VolMP[0][NRad-1],2.0)))*(pow(MESH.SupMP[0][NRad-1][0]*vleft[j]*cos(PI + MESH.AngleMU[0][NRad-1]) + MESH.SupMP[0][NRad-1][1]*0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1])*cos(PI + MESH.AngleMU[1][NRad-1]),2.0) + pow(MESH.SupMP[0][NRad-1][2]*0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + MESH.SupMP[0][NRad-1][3]*vup[0]*cos(PI + MESH.AngleMR[0][NRad]),2.0));

	//Esquina abajo derecha
	SA_Termino4[NA-1][0] = (Cb2/(Sigma*pow(MESH.VolMP[NA-1][0],2.0)))*(pow(MESH.SupMP[NA-1][0][0]*0.50*(vpresent[NA-1][0] + vpresent[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + MESH.SupMP[NA-1][0][1]*vright[0]*cos(PI + MESH.AngleMU[NA][0]),2.0) + pow(MESH.SupMP[NA-1][0][2]*vdown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + MESH.SupMP[NA-1][0][3]*0.50*(vpresent[NA-1][0] + vpresent[NA-1][1])*cos(PI + MESH.AngleMR[NA-1][1]),2.0));

	//Esquina arriba derecha
	SA_Termino4[NA-1][NRad-1] = (Cb2/(Sigma*pow(MESH.VolMP[NA-1][NRad-1],2.0)))*(pow(MESH.SupMP[NA-1][NRad-1][0]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + MESH.SupMP[NA-1][NRad-1][1]*vright[j]*cos(PI + MESH.AngleMU[NA][NRad-1]),2.0) + pow(MESH.SupMP[NA-1][NRad-1][2]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + MESH.SupMP[NA-1][NRad-1][3]*vup[NA-1]*cos(PI + MESH.AngleMR[NA-1][NRad]),2.0));
}

//Cálculo del término 5 de la ecuación de Spalart-Allmaras
void Solver::SA_Term5(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			SA_Termino5[i][j] = - (Cw1*fw[i][j] - (Cb1*ft2[i][j])/pow(K_VK,2.0))*pow(vpresent[i][j]/MESH.minDist[i][j],2.0);
		}
	}

}

//Cálculo de todas las contribuciones del modelo Spalart-Allmaras
void Solver::FmuSA(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Fmupresent[i][j] = SA_Termino1[i][j] + SA_Termino2[i][j] + SA_Termino3[i][j] + SA_Termino4[i][j] + SA_Termino5[i][j];
		}
	}

}

//Cálculo de variables y matrices necesarias para aplicar el modelo Spalart-Allmaras
void Solver::SpalartAllmarasPreparation(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			X[i][j] = vpresent[i][j]/(muBase[i][j]/(RhoPres[i][j]+1e-10) + 1e-10);

			fv1[i][j] = pow(X[i][j],3.0)/(pow(X[i][j],3.0) + pow(Cv1,3.0) + 1e-10);

			fv2[i][j] = 1.0 - X[i][j]/(1.0 + X[i][j]*fv1[i][j]);

			r[i][j] = vpresent[i][j]/(S[i][j]*pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0) + 1e-10);

			g[i][j] = r[i][j] + Cw2*(pow(r[i][j],6.0) - r[i][j]);

			fw[i][j] = g[i][j]*pow((1.0 + pow(Cw3,6.0))/(pow(g[i][j],6.0) + pow(Cw3,6.0)) + 1e-10,1.0/6.0);

			ft2[i][j] = Ct3*exp(-Ct4*pow(X[i][j],2.0));
 
			omega[i][j] = 0.50*(0.50*(GradV_DxMU[i][j] + GradV_DxMU[i+1][j]) - 0.50*(GradU_DyMR[i][j] + GradU_DyMR[i][j+1]));

			S[i][j] = sqrt(2.0*omega[i][j]*omega[i][j]);

			Smodel[i][j] = S[i][j] + vpresent[i][j]/(pow(K_VK,2.0)*pow(MESH.minDist[i][j],2.0)*fv2[i][j] + 1e-10);
		}
	}

}

//Cálculo de la viscosidad con el modelo Spalart-Allmaras
void Solver::SpalartAllmaras(){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			vfut[i][j] = (2.0*TimeBetta*vpresent[i][j] - (TimeBetta - 0.50)*vprevious[i][j])/(TimeBetta + 0.50) + DeltaT*((1.0 + TimeBetta)*Fmupresent[i][j] - TimeBetta*Fmuprevious[i][j]);
			
			muTurb[i][j] = RhoPres[i][j]*vfut[i][j]*fv1[i][j];
			
		}
	}

}