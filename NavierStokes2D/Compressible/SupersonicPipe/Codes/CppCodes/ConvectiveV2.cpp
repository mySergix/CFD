//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Axial U)
void Solver::Get_MomentumConvectiveU(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			MomentumConvectiveU[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i-1][j])*0.50*(Upres[i][j] + Upres[i-1][j])
									  + MESH.SupMP[i][j][1]*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(MESH.AngleMU[i+1][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i+1][j])*0.50*(Upres[i][j] + Upres[i+1][j])
									  + MESH.SupMP[i][j][2]*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j-1])*0.50*(Upres[i][j] + Upres[i][j-1])
									  + MESH.SupMP[i][j][3]*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j+1])*0.50*(Upres[i][j] + Upres[i][j+1])
									  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		MomentumConvectiveU[0][j] = -(1.0/MESH.VolMP[0][j])*(
								    MESH.SupMP[0][j][0]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j]) + Vleft[j]*sin(PI + MESH.AngleMU[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoLeft[j])*0.50*(Upres[0][j] + Uleft[j])
								  + MESH.SupMP[0][j][1]*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(MESH.AngleMU[1][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[1][j])*0.50*(Upres[0][j] + Upres[1][j])
								  + MESH.SupMP[0][j][2]*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j-1])*0.50*(Upres[0][j] + Upres[0][j-1])
								  + MESH.SupMP[0][j][3]*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j+1])*0.50*(Upres[0][j] + Upres[0][j+1])
								  );
		//Parte derecha
		MomentumConvectiveU[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(
									   MESH.SupMP[NA-1][j][0]*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-2][j])*0.50*(Upres[NA-1][j] + Upres[NA-2][j])
								     + MESH.SupMP[NA-1][j][1]*(Uright[j]*cos(MESH.AngleMU[NA][j]) + Vright[j]*sin(MESH.AngleMU[NA][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoRight[j])*0.50*(Upres[NA-1][j] + Uright[j])
									 + MESH.SupMP[NA-1][j][2]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j-1])*0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])
									 + MESH.SupMP[NA-1][j][3]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j+1])*0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])
									 );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		MomentumConvectiveU[i][0] = -(1.0/MESH.VolMP[i][0])*(
									MESH.SupMP[i][0][0]*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i-1][0])*0.50*(Upres[i][0] + Upres[i-1][0])
								  + MESH.SupMP[i][0][1]*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(MESH.AngleMU[i+1][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i+1][0])*0.50*(Upres[i][0] + Upres[i+1][0])
								  + MESH.SupMP[i][0][2]*(Upres[i][0]*cos(1.50*PI + MESH.AngleMR[i][0]) + Vpres[i][0]*sin(1.50*PI + MESH.AngleMR[i][0]))*RhoPres[i][0]*Upres[i][0]
								  + MESH.SupMP[i][0][3]*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i][1])*0.50*(Upres[i][0] + Upres[i][1])
								  );
		//Parte arriba
		MomentumConvectiveU[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(
									     MESH.SupMP[i][NRad-1][0]*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i-1][NRad-1])*0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])
									   + MESH.SupMP[i][NRad-1][1]*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i+1][NRad-1])*0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])
									   + MESH.SupMP[i][NRad-1][2]*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i][NRad-2])*0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])
									   + MESH.SupMP[i][NRad-1][3]*(Uup[i]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + Vup[i]*sin(0.50*PI + MESH.AngleMR[i][NRad]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoUp[i])*0.50*(Upres[i][NRad-1] + Uup[i])
									   );
	}
	
	//Esquina abajo izquierda
	MomentumConvectiveU[0][0] = -(1.0/MESH.VolMP[0][0])*(
								MESH.SupMP[0][0][0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0]) + Vleft[0]*sin(PI + MESH.AngleMU[0][0]))*sqrt(RhoPres[0][0])*sqrt(RhoLeft[0])*0.50*(Upres[0][0] + Uleft[0])
							  + MESH.SupMP[0][0][1]*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(MESH.AngleMU[1][0]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[1][0])*0.50*(Upres[0][0] + Upres[1][0])
							  + MESH.SupMP[0][0][2]*(Upres[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vpres[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))*RhoPres[0][0]*Upres[0][0]
							  + MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[0][1])*0.50*(Upres[0][0] + Upres[0][1])
							  );

	//Esquina arriba izquierda
	MomentumConvectiveU[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
									 MESH.SupMP[0][NRad-1][0]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoLeft[NRad-1])*0.50*(Upres[0][NRad-1] + Uleft[NRad-1])
								   + MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[1][NRad-1])*0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])
								   + MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[0][NRad-2])*0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])
								   + MESH.SupMP[0][NRad-1][3]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoUp[0])*0.50*(Upres[0][NRad-1] + Uup[0])
								   );
	//Esquina abajo derecha
	MomentumConvectiveU[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								   MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-2][0])*0.50*(Upres[NA-1][0] + Upres[NA-2][0])
								 + MESH.SupMP[NA-1][0][1]*(Uright[0]*cos(MESH.AngleMU[NA][0]) + Vright[0]*sin(MESH.AngleMU[NA][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoRight[0])*0.50*(Upres[NA-1][0] + Uright[0])
								 + MESH.SupMP[NA-1][0][2]*(Upres[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vpres[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))*RhoPres[NA-1][0]*0.50*Upres[NA-1][0]
								 + MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-1][1])*0.50*(Upres[NA-1][0] + Upres[NA-1][1])
								 );
	//Esquina arriba derecha
	MomentumConvectiveU[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
									    MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-2][NRad-1])*0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoRight[NRad-1])*0.50*(Upres[NA-1][NRad-1] + Uright[NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-1][NRad-2])*0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])
									  + MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoUp[NA-1])*0.50*(Upres[NA-1][NRad-1] + Uup[NA-1])
									  );
									  
}

//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomentumConvectiveV(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			MomentumConvectiveV[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i-1][j])*0.50*(Vpres[i][j] + Vpres[i-1][j])
									  + MESH.SupMP[i][j][1]*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(MESH.AngleMU[i+1][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i+1][j])*0.50*(Vpres[i][j] + Vpres[i+1][j])
									  + MESH.SupMP[i][j][2]*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j-1])*0.50*(Vpres[i][j] + Vpres[i][j-1])
									  + MESH.SupMP[i][j][3]*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j+1])*0.50*(Vpres[i][j] + Vpres[i][j+1])
									  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		MomentumConvectiveV[0][j] = -(1.0/MESH.VolMP[0][j])*(
								    MESH.SupMP[0][j][0]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j]) + Vleft[j]*sin(PI + MESH.AngleMU[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoLeft[j])*0.50*(Vpres[0][j] + Vleft[j])
								  + MESH.SupMP[0][j][1]*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(MESH.AngleMU[1][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[1][j])*0.50*(Vpres[0][j] + Vpres[1][j])
								  + MESH.SupMP[0][j][2]*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j-1])*0.50*(Vpres[0][j] + Vpres[0][j-1])
								  + MESH.SupMP[0][j][3]*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j+1])*0.50*(Vpres[0][j] + Vpres[0][j+1])
								  );
		//Parte derecha
		MomentumConvectiveV[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(
									   MESH.SupMP[NA-1][j][0]*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-2][j])*0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])
								     + MESH.SupMP[NA-1][j][1]*(Uright[j]*cos(MESH.AngleMU[NA][j]) + Vright[j]*sin(MESH.AngleMU[NA][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoRight[j])*0.50*(Vpres[NA-1][j] + Vright[j])
									 + MESH.SupMP[NA-1][j][2]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j-1])*0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])
									 + MESH.SupMP[NA-1][j][3]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j+1])*0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])
									 );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		MomentumConvectiveV[i][0] = -(1.0/MESH.VolMP[i][0])*(
									MESH.SupMP[i][0][0]*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i-1][0])*0.50*(Vpres[i][0] + Vpres[i-1][0])
								  + MESH.SupMP[i][0][1]*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(MESH.AngleMU[i+1][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i+1][0])*0.50*(Vpres[i][0] + Vpres[i+1][0])
								  + MESH.SupMP[i][0][2]*(Upres[i][0]*cos(1.50*PI + MESH.AngleMR[i][0]) + Vpres[i][0]*sin(1.50*PI + MESH.AngleMR[i][0]))*RhoPres[i][0]*Vpres[i][0]
								  + MESH.SupMP[i][0][3]*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i][1])*0.50*(Vpres[i][0] + Vpres[i][1])
								  );
		//Parte arriba
		MomentumConvectiveV[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(
									     MESH.SupMP[i][NRad-1][0]*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i-1][NRad-1])*0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])
									   + MESH.SupMP[i][NRad-1][1]*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i+1][NRad-1])*0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])
									   + MESH.SupMP[i][NRad-1][2]*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i][NRad-2])*0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])
									   + MESH.SupMP[i][NRad-1][3]*(Uup[i]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + Vup[i]*sin(0.50*PI + MESH.AngleMR[i][NRad]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoUp[i])*0.50*(Vpres[i][NRad-1] + Vup[i])
									   );
	}
	
	//Esquina abajo izquierda
	MomentumConvectiveV[0][0] = -(1.0/MESH.VolMP[0][0])*(
								MESH.SupMP[0][0][0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0]) + Vleft[0]*sin(PI + MESH.AngleMU[0][0]))*sqrt(RhoPres[0][0])*sqrt(RhoLeft[0])*0.50*(Vpres[0][0] + Vleft[0])
							  + MESH.SupMP[0][0][1]*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(MESH.AngleMU[1][0]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[1][0])*0.50*(Vpres[0][0] + Vpres[1][0])
							  + MESH.SupMP[0][0][2]*(Upres[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vpres[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))*RhoPres[0][0]*Vpres[0][0]
							  + MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[0][1])*0.50*(Vpres[0][0] + Vpres[0][1])
							  );

	//Esquina arriba izquierda
	MomentumConvectiveV[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
									 MESH.SupMP[0][NRad-1][0]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoLeft[NRad-1])*0.50*(Upres[0][NRad-1] + Uleft[NRad-1])
								   + MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[1][NRad-1])*0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])
								   + MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[0][NRad-2])*0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])
								   + MESH.SupMP[0][NRad-1][3]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoUp[0])*0.50*(Vpres[0][NRad-1] + Vup[0])
								   );
	//Esquina abajo derecha
	MomentumConvectiveV[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								   MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-2][0])*0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])
								 + MESH.SupMP[NA-1][0][1]*(Uright[0]*cos(MESH.AngleMU[NA][0]) + Vright[0]*sin(MESH.AngleMU[NA][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoRight[0])*0.50*(Vpres[NA-1][0] + Vright[0])
								 + MESH.SupMP[NA-1][0][2]*(Upres[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vpres[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))*RhoPres[NA-1][0]*0.50*Vpres[NA-1][0]
								 + MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-1][1])*0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])
								 );
	//Esquina arriba derecha
	MomentumConvectiveV[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
									    MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-2][NRad-1])*0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoRight[NRad-1])*0.50*(Vpres[NA-1][NRad-1] + Vright[NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-1][NRad-2])*0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])
									  + MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoUp[NA-1])*0.50*(Vpres[NA-1][NRad-1] + Vup[NA-1])
									  );
}