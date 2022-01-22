/Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Axial U)
void Solver::Get_MomemtumDifusiveU(Mesher MESH){
int i, j;
	
	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			MomentumDifusiveU[i][j] = (1.0/(sqrt(2)*MESH.VolMP[i][j]))*(
									  muWallsMU[i+1][j]*(sqrt(RhoPres[i+1][j])*Upres[i+1][j] - sqrt(RhoPres[i][j])*Upres[i][j])*(MESH.SupMP[i][j][1]/MESH.DeltasMU[i+1][j][0])*cos(MESH.AngleMU[i+1][j]) 
									+ muWallsMU[i][j]*(sqrt(RhoPres[i][j])*Upres[i][j] - sqrt(RhoPres[i-1][j])*Upres[i-1][j])*(MESH.SupMP[i][j][0]/MESH.DeltasMU[i][j][0])*cos(PI + MESH.AngleMU[i][j]) 
									+ muWallsMR[i][j+1]*(sqrt(RhoPres[i][j+1])*Upres[i][j+1] - sqrt(RhoPres[i][j])*Upres[i][j])*(MESH.SupMP[i][j][3]/MESH.DeltasMR[i][j+1][1])*sin(0.50*PI + MESH.AngleMR[i][j+1]) 
									+ muWallsMR[i][j]*(sqrt(RhoPres[i][j])*Upres[i][j] - sqrt(RhoPres[i][j-1])*Upres[i][j-1])*(MESH.SupMP[i][j][2]/MESH.DeltasMR[i][j][1])*sin(1.50*PI + MESH.AngleMR[i][j])
									);
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		MomentumDifusiveU[0][j] = (1.0/(sqrt(2)*MESH.VolMP[0][j]))*(
								  muWallsMU[1][j]*(sqrt(RhoPres[1][j])*Upres[1][j] - sqrt(RhoPres[0][j])*Upres[0][j])*(MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*cos(MESH.AngleMU[1][j]) 
							    + muWallsMU[0][j]*(sqrt(RhoPres[0][j])*Upres[0][j] - sqrt(RhoLeft[j])*Uleft[j])*(MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*cos(PI + MESH.AngleMU[0][j]) 
								+ muWallsMR[0][j+1]*(sqrt(RhoPres[0][j+1])*Upres[0][j+1] - sqrt(RhoPres[0][j])*Upres[0][j])*(MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*sin(0.50*PI + MESH.AngleMR[0][j+1]) 
								+ muWallsMR[0][j]*(sqrt(RhoPres[0][j])*Upres[0][j] - sqrt(RhoPres[0][j-1])*Upres[0][j-1])*(MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*sin(1.50*PI + MESH.AngleMR[0][j])
								);

		//Parte derecha
		MomentumDifusiveU[NA-1][j] = (1.0/(sqrt(2)*MESH.VolMP[NA-1][j]))*(
									 muWallsMU[NA][j]*(sqrt(RhoRight[j])*Uright[j] - sqrt(RhoPres[NA-1][j])*Upres[NA-1][j])*(MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*cos(MESH.AngleMU[NA][j]) 
								   + muWallsMU[NA-1][j]*(sqrt(RhoPres[NA-1][j])*Upres[NA-1][j] - sqrt(RhoPres[NA-1][j])*Upres[NA-2][j])*(MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*cos(PI + MESH.AngleMU[NA-1][j]) 
								   + muWallsMR[NA-1][j+1]*(sqrt(RhoPres[NA-1][j+1])*Upres[NA-1][j+1] - sqrt(RhoPres[NA-1][j])*Upres[NA-1][j])*(MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]) 
								   + muWallsMR[NA-1][j]*(sqrt(RhoPres[NA-1][j])*Upres[NA-1][j] - sqrt(RhoPres[NA-1][j-1])*Upres[NA-1][j-1])*(MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])
								   );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		MomentumDifusiveU[i][0] = (1.0/(sqrt(2)*MESH.VolMP[i][0]))*(
								  muWallsMU[i+1][0]*(sqrt(RhoPres[i+1][0])*Upres[i+1][0] - sqrt(RhoPres[i][0])*Upres[i][0])*(MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*cos(MESH.AngleMU[i+1][0]) 
								+ muWallsMU[i][0]*(sqrt(RhoPres[i][0])*Upres[i][0] - sqrt(RhoPres[i-1][0])*Upres[i-1][0])*(MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*cos(PI + MESH.AngleMU[i][0]) 
								+ muWallsMR[i][1]*(sqrt(RhoPres[i][1])*Upres[i][1] - sqrt(RhoPres[i][0])*Upres[i][0])*(MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*sin(0.50*PI + MESH.AngleMR[i][1]) 
								+ muWallsMR[i][0]*(sqrt(RhoPres[i][0])*Upres[i][0] - sqrt(RhoPres[i][0])*Upres[i][0])*(MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*sin(1.50*PI + MESH.AngleMR[i][0])
								);

		//Parte arriba
		MomentumDifusiveU[i][NRad-1] = (1.0/(sqrt(2)*MESH.VolMP[i][NRad-1]))*(
									   muWallsMU[i+1][NRad-1]*(sqrt(RhoPres[i+1][NRad-1])*Upres[i+1][NRad-1] - sqrt(RhoPres[i][NRad-1])*Upres[i][NRad-1])*(MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*cos(MESH.AngleMU[i+1][NRad-1]) 
									 + muWallsMU[i][NRad-1]*(sqrt(RhoPres[i][NRad-1])*Upres[i][NRad-1] - sqrt(RhoPres[i-1][NRad-1])*Upres[i-1][NRad-1])*(MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*cos(PI + MESH.AngleMU[i][NRad-1]) 
									 + muWallsMR[i][NRad]*(sqrt(RhoUp[i])*Uup[i] - sqrt(RhoPres[i][NRad-1])*Upres[i][NRad-1])*(MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*sin(0.50*PI + MESH.AngleMR[i][NRad]) 
									 + muWallsMR[i][NRad-1]*(sqrt(RhoPres[i][NRad-1])*Upres[i][NRad-1] - sqrt(RhoPres[i][NRad-2])*Upres[i][NRad-2])*(MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])
									 );
	}
	
	//Esquina abajo izquierda
	MomentumDifusiveU[0][0] = (1.0/(sqrt(2)*MESH.VolMP[0][0]))*(
							  muWallsMU[1][0]*(sqrt(RhoPres[1][0])*Upres[1][0] - sqrt(RhoPres[0][0])*Upres[0][0])*(MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*cos(MESH.AngleMU[1][0]) 
							+ muWallsMU[0][0]*(sqrt(RhoPres[0][0])*Upres[0][0] - sqrt(RhoLeft[0])*Uleft[0])*(MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*cos(PI + MESH.AngleMU[0][0]) 
							+ muWallsMR[0][1]*(sqrt(RhoPres[0][1])*Upres[0][1] - sqrt(RhoPres[0][0])*Upres[0][0])*(MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*sin(0.50*PI + MESH.AngleMR[0][1]) 
							+ muWallsMR[0][0]*(sqrt(RhoPres[0][0])*Upres[0][0] - sqrt(RhoPres[0][0])*Upres[0][0])*(MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*sin(1.50*PI + MESH.AngleMR[0][0])
							);

	//Esquina abajo derecha
	MomentumDifusiveU[NA-1][0] = (1.0/(sqrt(2)*MESH.VolMP[NA-1][0]))*(
								 muWallsMU[NA][0]*(sqrt(RhoRight[0])*Uright[0] - sqrt(RhoPres[NA-1][0])*Upres[NA-1][0])*(MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*cos(MESH.AngleMU[NA][0]) 
							   + muWallsMU[NA-1][0]*(sqrt(RhoPres[NA-1][0])*Upres[NA-1][0] - sqrt(RhoPres[NA-2][0])*Upres[NA-2][0])*(MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*cos(PI + MESH.AngleMU[NA-1][0]) 
							   + muWallsMR[NA-1][1]*(sqrt(RhoPres[NA-1][1])*Upres[NA-1][1] - sqrt(RhoPres[NA-1][0])*Upres[NA-1][0])*(MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]) 
							   + muWallsMR[NA-1][j]*(sqrt(RhoPres[NA-1][0])*Upres[NA-1][0] - sqrt(RhoPres[NA-1][0])*Upres[NA-1][0])*(MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*sin(1.50*PI + MESH.AngleMR[NA-1][0])
							   );
	
	//Esquina arriba izquierda
	MomentumDifusiveU[0][NRad-1] = (1.0/(sqrt(2)*MESH.VolMP[0][NRad-1]))*(
								   muWallsMU[1][NRad-1]*(sqrt(RhoPres[1][NRad-1])*Upres[1][NRad-1] - sqrt(RhoPres[0][NRad-1])*Upres[0][NRad-1])*(MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*cos(MESH.AngleMU[1][NRad-1]) 
								 + muWallsMU[0][NRad-1]*(sqrt(RhoPres[0][NRad-1])*Upres[0][NRad-1] - sqrt(RhoLeft[NRad-1])*Uleft[NRad-1])*(MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*cos(PI + MESH.AngleMU[0][NRad-1]) 
								 + muWallsMR[0][NRad]*(sqrt(RhoUp[0])*Uup[0] - sqrt(RhoPres[0][NRad-1])*Upres[0][NRad-1])*(MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*sin(0.50*PI + MESH.AngleMR[0][NRad]) 
								 + muWallsMR[0][NRad-1]*(sqrt(RhoPres[0][NRad-1])*Upres[0][NRad-1] - sqrt(RhoPres[0][NRad-2])*Upres[0][NRad-2])*(MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])
								 );

	//Esquina arriba derecha
	MomentumDifusiveU[NA-1][NRad-1] = (1.0/(sqrt(2)*MESH.VolMP[NA-1][NRad-1]))*(
									  muWallsMU[NA-1][NRad-1]*(sqrt(RhoRight[NRad-1])*Uright[NRad-1] - sqrt(RhoPres[NA-1][NRad-1])*Upres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*cos(MESH.AngleMU[NA][NRad-1]) 
									+ muWallsMU[NA-1][NRad-1]*(sqrt(RhoPres[NA-1][NRad-1])*Upres[NA-1][NRad-1] - sqrt(RhoPres[NA-2][NRad-1])*Upres[NA-2][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) 
									+ muWallsMR[NA-1][NRad]*(sqrt(RhoUp[NA-1])*Uup[NA-1] - sqrt(RhoPres[NA-1][NRad-1])*Upres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]) 
									+ muWallsMR[NA-1][NRad-1]*(sqrt(RhoPres[NA-1][NRad-1])*Upres[NA-1][NRad-1] - sqrt(RhoPres[NA-1][NRad-2])*Upres[NA-1][NRad-2])*(MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])
									);
}

//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveV(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			MomentumDifusiveV[i][j] = (1.0/(sqrt(2)*MESH.VolMP[i][j]))*(
									  muWallsMU[i+1][j]*(sqrt(RhoPres[i+1][j])*Vpres[i+1][j] - sqrt(RhoPres[i][j])*Vpres[i][j])*(MESH.SupMP[i][j][1]/MESH.DeltasMU[i+1][j][0])*cos(MESH.AngleMU[i+1][j]) 
									+ muWallsMU[i][j]*(sqrt(RhoPres[i][j])*Vpres[i][j] - sqrt(RhoPres[i-1][j])*Vpres[i-1][j])*(MESH.SupMP[i][j][0]/MESH.DeltasMU[i][j][0])*cos(PI + MESH.AngleMU[i][j]) 
									+ muWallsMR[i][j+1]*(sqrt(RhoPres[i][j+1])*Vpres[i][j+1] - sqrt(RhoPres[i][j])*Vpres[i][j])*(MESH.SupMP[i][j][3]/MESH.DeltasMR[i][j+1][1])*sin(0.50*PI + MESH.AngleMR[i][j+1]) 
									+ muWallsMR[i][j]*(sqrt(RhoPres[i][j])*Vpres[i][j] - sqrt(RhoPres[i][j-1])*Vpres[i][j-1])*(MESH.SupMP[i][j][2]/MESH.DeltasMR[i][j][1])*sin(1.50*PI + MESH.AngleMR[i][j])
									);
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		MomentumDifusiveV[0][j] = (1.0/(sqrt(2)*MESH.VolMP[0][j]))*(
								  muWallsMU[1][j]*(sqrt(RhoPres[1][j])*Vpres[1][j] - sqrt(RhoPres[0][j])*Vpres[0][j])*(MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*cos(MESH.AngleMU[1][j]) 
							    + muWallsMU[0][j]*(sqrt(RhoPres[0][j])*Vpres[0][j] - sqrt(RhoLeft[j])*Vleft[j])*(MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*cos(PI + MESH.AngleMU[0][j]) 
								+ muWallsMR[0][j+1]*(sqrt(RhoPres[0][j+1])*Vpres[0][j+1] - sqrt(RhoPres[0][j])*Vpres[0][j])*(MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*sin(0.50*PI + MESH.AngleMR[0][j+1]) 
								+ muWallsMR[0][j]*(sqrt(RhoPres[0][j])*Vpres[0][j] - sqrt(RhoPres[0][j-1])*Vpres[0][j-1])*(MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*sin(1.50*PI + MESH.AngleMR[0][j])
								);

		//Parte derecha
		MomentumDifusiveV[NA-1][j] = (1.0/(sqrt(2)*MESH.VolMP[NA-1][j]))*(
									 muWallsMU[NA][j]*(sqrt(RhoRight[j])*Vright[j] - sqrt(RhoPres[NA-1][j])*Vpres[NA-1][j])*(MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*cos(MESH.AngleMU[NA][j]) 
								   + muWallsMU[NA-1][j]*(sqrt(RhoPres[NA-1][j])*Vpres[NA-1][j] - sqrt(RhoPres[NA-1][j])*Vpres[NA-2][j])*(MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*cos(PI + MESH.AngleMU[NA-1][j]) 
								   + muWallsMR[NA-1][j+1]*(sqrt(RhoPres[NA-1][j+1])*Vpres[NA-1][j+1] - sqrt(RhoPres[NA-1][j])*Vpres[NA-1][j])*(MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]) 
								   + muWallsMR[NA-1][j]*(sqrt(RhoPres[NA-1][j])*Vpres[NA-1][j] - sqrt(RhoPres[NA-1][j-1])*Vpres[NA-1][j-1])*(MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])
								   );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		MomentumDifusiveV[i][0] = (1.0/(sqrt(2)*MESH.VolMP[i][0]))*(
								  muWallsMU[i+1][0]*(sqrt(RhoPres[i+1][0])*Vpres[i+1][0] - sqrt(RhoPres[i][0])*Vpres[i][0])*(MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*cos(MESH.AngleMU[i+1][0]) 
								+ muWallsMU[i][0]*(sqrt(RhoPres[i][0])*Vpres[i][0] - sqrt(RhoPres[i-1][0])*Vpres[i-1][0])*(MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*cos(PI + MESH.AngleMU[i][0]) 
								+ muWallsMR[i][1]*(sqrt(RhoPres[i][1])*Vpres[i][1] - sqrt(RhoPres[i][0])*Vpres[i][0])*(MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*sin(0.50*PI + MESH.AngleMR[i][1]) 
								+ muWallsMR[i][0]*(sqrt(RhoPres[i][0])*Vpres[i][0] - sqrt(RhoPres[i][0])*Vpres[i][0])*(MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*sin(1.50*PI + MESH.AngleMR[i][0])
								);

		//Parte arriba
		MomentumDifusiveV[i][NRad-1] = (1.0/(sqrt(2)*MESH.VolMP[i][NRad-1]))*(
									   muWallsMU[i+1][NRad-1]*(sqrt(RhoPres[i+1][NRad-1])*Vpres[i+1][NRad-1] - sqrt(RhoPres[i][NRad-1])*Vpres[i][NRad-1])*(MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*cos(MESH.AngleMU[i+1][NRad-1]) 
									 + muWallsMU[i][NRad-1]*(sqrt(RhoPres[i][NRad-1])*Vpres[i][NRad-1] - sqrt(RhoPres[i-1][NRad-1])*Vpres[i-1][NRad-1])*(MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*cos(PI + MESH.AngleMU[i][NRad-1]) 
									 + muWallsMR[i][NRad]*(sqrt(RhoUp[i])*Vup[i] - sqrt(RhoPres[i][NRad-1])*Vpres[i][NRad-1])*(MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*sin(0.50*PI + MESH.AngleMR[i][NRad]) 
									 + muWallsMR[i][NRad-1]*(sqrt(RhoPres[i][NRad-1])*Vpres[i][NRad-1] - sqrt(RhoPres[i][NRad-2])*Vpres[i][NRad-2])*(MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])
									 );
	}
	
	//Esquina abajo izquierda
	MomentumDifusiveV[0][0] = (1.0/(sqrt(2)*MESH.VolMP[0][0]))*(
							  muWallsMU[1][0]*(sqrt(RhoPres[1][0])*Vpres[1][0] - sqrt(RhoPres[0][0])*Vpres[0][0])*(MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*cos(MESH.AngleMU[1][0]) 
							+ muWallsMU[0][0]*(sqrt(RhoPres[0][0])*Vpres[0][0] - sqrt(RhoLeft[0])*Vleft[0])*(MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*cos(PI + MESH.AngleMU[0][0]) 
							+ muWallsMR[0][1]*(sqrt(RhoPres[0][1])*Vpres[0][1] - sqrt(RhoPres[0][0])*Vpres[0][0])*(MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*sin(0.50*PI + MESH.AngleMR[0][1]) 
							+ muWallsMR[0][0]*(sqrt(RhoPres[0][0])*Vpres[0][0] - sqrt(RhoPres[0][0])*Vpres[0][0])*(MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*sin(1.50*PI + MESH.AngleMR[0][0])
							);

	//Esquina abajo derecha
	MomentumDifusiveV[NA-1][0] = (1.0/(sqrt(2)*MESH.VolMP[NA-1][0]))*(
								 muWallsMU[NA][0]*(sqrt(RhoRight[0])*Vright[0] - sqrt(RhoPres[NA-1][0])*Vpres[NA-1][0])*(MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*cos(MESH.AngleMU[NA][0]) 
							   + muWallsMU[NA-1][0]*(sqrt(RhoPres[NA-1][0])*Vpres[NA-1][0] - sqrt(RhoPres[NA-2][0])*Vpres[NA-2][0])*(MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*cos(PI + MESH.AngleMU[NA-1][0]) 
							   + muWallsMR[NA-1][1]*(sqrt(RhoPres[NA-1][1])*Vpres[NA-1][1] - sqrt(RhoPres[NA-1][0])*Vpres[NA-1][0])*(MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]) 
							   + muWallsMR[NA-1][j]*(sqrt(RhoPres[NA-1][0])*Vpres[NA-1][0] - sqrt(RhoPres[NA-1][0])*Vpres[NA-1][0])*(MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*sin(1.50*PI + MESH.AngleMR[NA-1][0])
							   );
	
	//Esquina arriba izquierda
	MomentumDifusiveV[0][NRad-1] = (1.0/(sqrt(2)*MESH.VolMP[0][NRad-1]))*(
								   muWallsMU[1][NRad-1]*(sqrt(RhoPres[1][NRad-1])*Vpres[1][NRad-1] - sqrt(RhoPres[0][NRad-1])*Vpres[0][NRad-1])*(MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*cos(MESH.AngleMU[1][NRad-1]) 
								 + muWallsMU[0][NRad-1]*(sqrt(RhoPres[0][NRad-1])*Vpres[0][NRad-1] - sqrt(RhoLeft[NRad-1])*Vleft[NRad-1])*(MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*cos(PI + MESH.AngleMU[0][NRad-1]) 
								 + muWallsMR[0][NRad]*(sqrt(RhoUp[0])*Vup[0] - sqrt(RhoPres[0][NRad-1])*Vpres[0][NRad-1])*(MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*sin(0.50*PI + MESH.AngleMR[0][NRad]) 
								 + muWallsMR[0][NRad-1]*(sqrt(RhoPres[0][NRad-1])*Vpres[0][NRad-1] - sqrt(RhoPres[0][NRad-2])*Vpres[0][NRad-2])*(MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])
								 );

	//Esquina arriba derecha
	MomentumDifusiveV[NA-1][NRad-1] = (1.0/(sqrt(2)*MESH.VolMP[NA-1][NRad-1]))*(
									  muWallsMU[NA-1][NRad-1]*(sqrt(RhoRight[NRad-1])*Vright[NRad-1] - sqrt(RhoPres[NA-1][NRad-1])*Vpres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*cos(MESH.AngleMU[NA][NRad-1]) 
									+ muWallsMU[NA-1][NRad-1]*(sqrt(RhoPres[NA-1][NRad-1])*Vpres[NA-1][NRad-1] - sqrt(RhoPres[NA-2][NRad-1])*Vpres[NA-2][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) 
									+ muWallsMR[NA-1][NRad]*(sqrt(RhoUp[NA-1])*Vup[NA-1] - sqrt(RhoPres[NA-1][NRad-1])*Vpres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]) 
									+ muWallsMR[NA-1][NRad-1]*(sqrt(RhoPres[NA-1][NRad-1])*Vpres[NA-1][NRad-1] - sqrt(RhoPres[NA-1][NRad-2])*Vpres[NA-1][NRad-2])*(MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])
									);
}