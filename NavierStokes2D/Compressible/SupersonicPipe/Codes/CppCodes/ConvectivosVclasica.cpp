//Parte central
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			FRhoPresent[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i-1][j])
									  + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i+1][j])
									  + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j-1])
									  + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j+1])
									  );
		}
	}
	
	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		FRhoPresent[0][j] = -(1.0/MESH.VolMP[0][j])*(
									    MESH.SupMP[0][j][0]*(UwallsMU[0][j]*cos(PI + MESH.AngleMU[0][j]) + VwallsMU[0][j]*sin(PI + MESH.AngleMU[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoLeft[j])
									  + MESH.SupMP[0][j][1]*(UwallsMU[1][j]*cos(MESH.AngleMU[1][j]) + VwallsMU[1][j]*sin(MESH.AngleMU[1][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[1][j])
									  + MESH.SupMP[0][j][2]*(UwallsMR[0][j]*cos(1.50*PI + MESH.AngleMR[0][j]) + VwallsMR[0][j]*sin(1.50*PI + MESH.AngleMR[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j-1])
									  + MESH.SupMP[0][j][3]*(UwallsMR[0][j+1]*cos(0.50*PI + MESH.AngleMR[0][j+1]) + VwallsMR[0][j+1]*sin(0.50*PI + MESH.AngleMR[0][j+1]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j+1])
									  );

		//Parte derecha
		FRhoPresent[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(
									    MESH.SupMP[NA-1][j][0]*(UwallsMU[NA-1][j]*cos(PI + MESH.AngleMU[NA-1][j]) + VwallsMU[NA-1][j]*sin(PI + MESH.AngleMU[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-2][j])
									  + MESH.SupMP[NA-1][j][1]*(UwallsMU[NA][j]*cos(MESH.AngleMU[NA][j]) + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoRight[j])
									  + MESH.SupMP[NA-1][j][2]*(UwallsMR[NA-1][j]*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + VwallsMR[NA-1][j]*sin(1.50*PI + MESH.AngleMR[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j-1])
									  + MESH.SupMP[NA-1][j][3]*(UwallsMR[NA-1][j+1]*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + VwallsMR[NA-1][j+1]*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j+1])
									  );
		}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		FRhoPresent[i][0] = -(1.0/MESH.VolMP[i][0])*(
									MESH.SupMP[i][0][0]*(UwallsMU[i][0]*cos(PI + MESH.AngleMU[i][0]) + VwallsMU[i][0]*sin(PI + MESH.AngleMU[i][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i-1][0])
								  + MESH.SupMP[i][0][1]*(UwallsMU[i+1][0]*cos(MESH.AngleMU[i+1][0]) + VwallsMU[i+1][0]*sin(MESH.AngleMU[i+1][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i+1][0])
								  + MESH.SupMP[i][0][2]*(UwallsMR[i][0]*cos(1.50*PI + MESH.AngleMR[i][0]) + VwallsMR[i][0]*sin(1.50*PI + MESH.AngleMR[i][0]))*sqrt(RhoDown[i])*sqrt(RhoPres[i][0])
								  + MESH.SupMP[i][0][3]*(UwallsMR[i][1]*cos(0.50*PI + MESH.AngleMR[i][1]) + VwallsMR[i][1]*sin(0.50*PI + MESH.AngleMR[i][1]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i][1])
								  );
		//Parte arriba
		FRhoPresent[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(
									     MESH.SupMP[i][NRad-1][0]*(UwallsMU[i][NRad-1]*cos(PI + MESH.AngleMU[i][NRad-1]) + VwallsMU[i][NRad-1]*sin(PI + MESH.AngleMU[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i-1][NRad-1])
									   + MESH.SupMP[i][NRad-1][1]*(UwallsMU[i+1][NRad-1]*cos(MESH.AngleMU[i+1][NRad-1]) + VwallsMU[i+1][NRad-1]*sin(MESH.AngleMU[i+1][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i+1][NRad-1])
									   + MESH.SupMP[i][NRad-1][2]*(UwallsMR[i][NRad-1]*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + VwallsMR[i][NRad-1]*sin(1.50*PI + MESH.AngleMR[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i][NRad-2])
									   + MESH.SupMP[i][NRad-1][3]*(UwallsMR[i][NRad]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + VwallsMR[i][NRad]*sin(0.50*PI + MESH.AngleMR[i][NRad]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoUp[i])
									   );
	}

	//Esquina abajo izquierda
	FRhoPresent[0][0] = -(1.0/MESH.VolMP[0][0])*(
								MESH.SupMP[0][0][0]*(UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0]) + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0]))*sqrt(RhoPres[0][0])*sqrt(RhoLeft[0])
							  + MESH.SupMP[0][0][1]*(UwallsMU[1][0]*cos(MESH.AngleMU[1][0]) + VwallsMU[1][0]*sin(MESH.AngleMU[1][0]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[1][0])
							  + MESH.SupMP[0][0][2]*(UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))*sqrt(RhoDown[0])*sqrt(RhoPres[0][0])
							  + MESH.SupMP[0][0][3]*(UwallsMR[0][1]*cos(0.50*PI + MESH.AngleMR[0][1]) + VwallsMR[0][1]*sin(0.50*PI + MESH.AngleMR[0][1]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[0][1])
							  );

	//Esquina arriba izquierda
	FRhoPresent[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
									 MESH.SupMP[0][NRad-1][0]*(UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoLeft[NRad-1])
								   + MESH.SupMP[0][NRad-1][1]*(UwallsMU[1][NRad-1]*cos(MESH.AngleMU[1][NRad-1]) + VwallsMU[1][NRad-1]*sin(MESH.AngleMU[1][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[1][NRad-1])
								   + MESH.SupMP[0][NRad-1][2]*(UwallsMR[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + VwallsMR[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[0][NRad-2])
								   + MESH.SupMP[0][NRad-1][3]*(UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoUp[0])
								   );
	//Esquina abajo derecha
	FRhoPresent[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								   MESH.SupMP[NA-1][0][0]*(UwallsMU[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0]) + VwallsMU[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-2][0])
								 + MESH.SupMP[NA-1][0][1]*(UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0]) + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoRight[0])
								 + MESH.SupMP[NA-1][0][2]*(UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoDown[NA-1])
								 + MESH.SupMP[NA-1][0][3]*(UwallsMR[NA-1][1]*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + VwallsMR[NA-1][1]*sin(0.50*PI + MESH.AngleMR[NA-1][1]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-1][1])
								 );
	//Esquina arriba derecha
	FRhoPresent[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
									    MESH.SupMP[NA-1][NRad-1][0]*(UwallsMU[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + VwallsMU[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-2][NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][1]*(UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + VwallsMU[NA][NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoRight[NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][2]*(UwallsMR[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + VwallsMR[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-1][NRad-2])
									  + MESH.SupMP[NA-1][NRad-1][3]*(UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoUp[NA-1])
									  );





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
								  + MESH.SupMP[i][0][2]*(Udown[i]*cos(1.50*PI + MESH.AngleMR[i][0]) + Vdown[i]*sin(1.50*PI + MESH.AngleMR[i][0]))*sqrt(RhoDown[i])*sqrt(RhoPres[i][0])*0.50*(Upres[i][0] + Udown[i])
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
							  + MESH.SupMP[0][0][2]*(Udown[0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vdown[0]*sin(1.50*PI + MESH.AngleMR[0][0]))*sqrt(RhoDown[0])*sqrt(RhoPres[0][0])*0.50*(Upres[0][0] + Udown[0])
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
								 + MESH.SupMP[NA-1][0][2]*(Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vdown[NA-1]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoDown[NA-1])*0.50*(Upres[NA-1][0] + Udown[NA-1])
								 + MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-1][1])*0.50*(Upres[NA-1][0] + Upres[NA-1][1])
								 );
	//Esquina arriba derecha
	MomentumConvectiveU[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
									    MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-2][NRad-1])*0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoRight[NRad-1])*0.50*(Upres[NA-1][NRad-1] + Uright[NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-1][NRad-2])*0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])
									  + MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoUp[NA-1])*0.50*(Upres[NA-1][NRad-1] + Uup[NA-1])
									  );



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
								  + MESH.SupMP[i][0][2]*(Udown[i]*cos(1.50*PI + MESH.AngleMR[i][0]) + Vdown[i]*sin(1.50*PI + MESH.AngleMR[i][0]))*sqrt(RhoPres[i][0])*sqrt(RhoDown[i])*0.50*(Vpres[i][0] + Vdown[i])
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
							  + MESH.SupMP[0][0][2]*(Udown[0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vdown[0]*sin(1.50*PI + MESH.AngleMR[0][0]))*sqrt(RhoPres[0][0])*sqrt(RhoDown[0])*0.50*(Vpres[0][0] + Vdown[0])
							  + MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[0][1])*0.50*(Vpres[0][0] + Vpres[0][1])
							  );

	//Esquina arriba izquierda
	MomentumConvectiveV[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
									 MESH.SupMP[0][NRad-1][0]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoLeft[NRad-1])*0.50*(Vpres[0][NRad-1] + Vleft[NRad-1])
								   + MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[1][NRad-1])*0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])
								   + MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[0][NRad-2])*0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])
								   + MESH.SupMP[0][NRad-1][3]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoUp[0])*0.50*(Vpres[0][NRad-1] + Vup[0])
								   );
	//Esquina abajo derecha
	MomentumConvectiveV[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								   MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-2][0])*0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])
								 + MESH.SupMP[NA-1][0][1]*(Uright[0]*cos(MESH.AngleMU[NA][0]) + Vright[0]*sin(MESH.AngleMU[NA][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoRight[0])*0.50*(Vpres[NA-1][0] + Vright[0])
								 + MESH.SupMP[NA-1][0][2]*(Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vdown[NA-1]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoDown[NA-1])*0.50*(Vpres[NA-1][0] + Vdown[NA-1])
								 + MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-1][1])*0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])
								 );
	//Esquina arriba derecha
	MomentumConvectiveV[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
									    MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-2][NRad-1])*0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoRight[NRad-1])*0.50*(Vpres[NA-1][NRad-1] + Vright[NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-1][NRad-2])*0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])
									  + MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoUp[NA-1])*0.50*(Vpres[NA-1][NRad-1] + Vup[NA-1])
									  );



	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			EnergyConvective[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i-1][j])*sqrt(epres[i][j]*epres[i-1][j])
									  + MESH.SupMP[i][j][1]*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(MESH.AngleMU[i+1][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i+1][j])*sqrt(epres[i][j]*epres[i+1][j])
									  + MESH.SupMP[i][j][2]*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j-1])*sqrt(epres[i][j]*epres[i][j-1])
									  + MESH.SupMP[i][j][3]*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1]))*sqrt(RhoPres[i][j])*sqrt(RhoPres[i][j+1])*sqrt(epres[i][j]*epres[i][j+1])
									  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		EnergyConvective[0][j] = -(1.0/MESH.VolMP[0][j])*(
								    MESH.SupMP[0][j][0]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j]) + Vleft[j]*sin(PI + MESH.AngleMU[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoLeft[j])*sqrt(epres[0][j]*eleft[j])
								  + MESH.SupMP[0][j][1]*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(MESH.AngleMU[1][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[1][j])*sqrt(epres[0][j]*epres[1][j])
								  + MESH.SupMP[0][j][2]*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j-1])*sqrt(epres[0][j]*epres[0][j-1])
								  + MESH.SupMP[0][j][3]*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1]))*sqrt(RhoPres[0][j])*sqrt(RhoPres[0][j+1])*sqrt(epres[0][j]*epres[0][j+1])
								  );
		//Parte derecha
		EnergyConvective[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(
									   MESH.SupMP[NA-1][j][0]*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-2][j])*sqrt(epres[NA-1][j]*epres[NA-2][j])
								    + MESH.SupMP[NA-1][j][1]*(Uright[j]*cos(MESH.AngleMU[NA][j]) + Vright[j]*sin(MESH.AngleMU[NA][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoRight[j])*sqrt(epres[NA-1][j]*eright[j])
									 + MESH.SupMP[NA-1][j][2]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j-1])*sqrt(epres[NA-1][j]*epres[NA-1][j-1])
									 + MESH.SupMP[NA-1][j][3]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))*sqrt(RhoPres[NA-1][j])*sqrt(RhoPres[NA-1][j+1])*sqrt(epres[NA-1][j]*epres[NA-1][j+1])
									 );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		EnergyConvective[i][0] = -(1.0/MESH.VolMP[i][0])*(
									MESH.SupMP[i][0][0]*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i-1][0])*sqrt(epres[i][0]*epres[i-1][0])
								  + MESH.SupMP[i][0][1]*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(MESH.AngleMU[i+1][0]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i+1][0])*sqrt(epres[i][0]*epres[i+1][0])
								  + MESH.SupMP[i][0][2]*(Udown[i]*cos(1.50*PI + MESH.AngleMR[i][0]) + Vdown[i]*sin(1.50*PI + MESH.AngleMR[i][0]))*sqrt(RhoDown[i])*sqrt(RhoPres[i][0])*sqrt(epres[i][0]*edown[i])
								  + MESH.SupMP[i][0][3]*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1]))*sqrt(RhoPres[i][0])*sqrt(RhoPres[i][1])*sqrt(epres[i][0]*epres[i][1])
								  );
		//Parte arriba
		EnergyConvective[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(
									     MESH.SupMP[i][NRad-1][0]*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i-1][NRad-1])*sqrt(epres[i][NRad-1]*epres[i-1][NRad-1])
									   + MESH.SupMP[i][NRad-1][1]*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i+1][NRad-1])*sqrt(epres[i][NRad-1]*epres[i+1][NRad-1])
									   + MESH.SupMP[i][NRad-1][2]*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoPres[i][NRad-2])*sqrt(epres[i][NRad-1]*epres[i][NRad-2])
									   + MESH.SupMP[i][NRad-1][3]*(Uup[i]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + Vup[i]*sin(0.50*PI + MESH.AngleMR[i][NRad]))*sqrt(RhoPres[i][NRad-1])*sqrt(RhoUp[i])*sqrt(epres[i][NRad-1]*eup[i])
									   );
	}
	
	//Esquina abajo izquierda
	EnergyConvective[0][0] = -(1.0/MESH.VolMP[0][0])*(
								MESH.SupMP[0][0][0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0]) + Vleft[0]*sin(PI + MESH.AngleMU[0][0]))*sqrt(RhoPres[0][0])*sqrt(RhoLeft[0])*sqrt(epres[0][0]*eleft[0])
							  + MESH.SupMP[0][0][1]*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(MESH.AngleMU[1][0]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[1][0])*sqrt(epres[0][0]*epres[1][0])
							  + MESH.SupMP[0][0][2]*(Udown[0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vdown[0]*sin(1.50*PI + MESH.AngleMR[0][0]))*sqrt(RhoDown[0])*sqrt(RhoPres[0][0])*sqrt(epres[0][0]*edown[0])
							  + MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1]))*sqrt(RhoPres[0][0])*sqrt(RhoPres[0][1])*sqrt(epres[0][0]*epres[0][1])
							  );

	//Esquina arriba izquierda
	EnergyConvective[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(
									 MESH.SupMP[0][NRad-1][0]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoLeft[NRad-1])*sqrt(epres[0][NRad-1]*eleft[NRad-1])
								   + MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[1][NRad-1])*sqrt(epres[0][NRad-1]*epres[1][NRad-1])
								   + MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoPres[0][NRad-2])*sqrt(epres[0][NRad-1]*epres[0][NRad-2])
								   + MESH.SupMP[0][NRad-1][3]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad]))*sqrt(RhoPres[0][NRad-1])*sqrt(RhoUp[0])*sqrt(epres[0][NRad-1]*eup[0])
								   );
	//Esquina abajo derecha
	EnergyConvective[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(
								   MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-2][0])*sqrt(epres[NA-1][0]*epres[NA-2][0])
								 + MESH.SupMP[NA-1][0][1]*(Uright[0]*cos(MESH.AngleMU[NA][0]) + Vright[0]*sin(MESH.AngleMU[NA][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoRight[0])*sqrt(epres[NA-1][0]*eright[0])
								 + MESH.SupMP[NA-1][0][2]*(Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vdown[NA-1]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoDown[NA-1])*sqrt(epres[NA-1][0]*edown[NA-1])
								 + MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1]))*sqrt(RhoPres[NA-1][0])*sqrt(RhoPres[NA-1][1])*sqrt(epres[NA-1][0]*epres[NA-1][1])
								 );
	//Esquina arriba derecha
	EnergyConvective[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(
									    MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-2][NRad-1])*sqrt(epres[NA-1][NRad-1]*epres[NA-2][NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoRight[NRad-1])*sqrt(epres[NA-1][NRad-1]*eright[NRad-1])
									  + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoPres[NA-1][NRad-2])*sqrt(epres[NA-1][NRad-1]*epres[NA-1][NRad-2])
									  + MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))*sqrt(RhoPres[NA-1][NRad-1])*sqrt(RhoUp[NA-1])*sqrt(epres[NA-1][NRad-1]*eup[NA-1])
									  );