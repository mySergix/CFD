//Cálculo de la divergencia de la velocidad en cada nodo (Término viscoso ecuación de energía)
void Solver::Get_Divergence(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			Divergence[i][j] = (1.0/MESH.VolMP[i][j])*(
							   MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))
							 + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))
							 + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
							 + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
							 ) 
							 + Vpres[i][j]/MESH.MP[i][j][1]
			    			 ;
		}
	}

	//Nodos U

	//Centro
	for(j = 1; j < NRad; j++){
		for(i = 1; i < NA; i++){
			DivergenceMU[i][j] = 0.50*(Divergence[i][j] + Divergence[i-1][j]);
		}
	}

	for(j = 1; j < NRad; j++){

		//Parte izquierda
		DivergenceMU[0][j] = (1.0/(2.0*MESH.VolMP[0][j]))*(
						     MESH.SupMP[0][j][0]*(UwallsMU[0][j]*cos(PI + MESH.AngleMU[0][j]) + VwallsMU[0][j]*sin(PI + MESH.AngleMU[0][j]))
						   + MESH.SupMP[0][j][1]*(Upres[0][j]*cos(MESH.AngleMU[1][j]) + Vpres[0][j]*sin(MESH.AngleMU[1][j]))
						   + 0.50*MESH.SupMP[0][j][2]*(UwallsMR[0][j]*cos(1.50*PI + MESH.AngleMR[0][j]) + VwallsMR[0][j]*sin(1.50*PI + MESH.AngleMR[0][j]))
						   + 0.50*MESH.SupMP[0][j][3]*(UwallsMR[0][j+1]*cos(0.50*PI + MESH.AngleMR[0][j+1]) + VwallsMR[0][j+1]*sin(0.50*PI + MESH.AngleMR[0][j+1]))
						   )
						   + VwallsMU[0][j]/MESH.MU[0][j][1]	
						   ;
		//Parte derecha
		DivergenceMU[NA][j] = (1.0/(2.0*MESH.VolMP[NA-1][j]))*(
							  MESH.SupMP[NA-1][j][0]*(Upres[NA-1][j]*cos(PI + MESH.AngleMU[NA-1][j]) + Vpres[NA-1][j]*sin(PI + MESH.AngleMU[NA-1][j]))
							+ MESH.SupMP[NA-1][j][1]*(UwallsMU[NA][j]*cos(MESH.AngleMU[NA][j]) + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][j]))
							+ 0.50*MESH.SupMP[NA-1][j][2]*(UwallsMR[NA-1][j]*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + VwallsMR[NA-1][j]*sin(1.50*PI + MESH.AngleMR[NA-1][j]))
							+ 0.50*MESH.SupMP[NA-1][j][3]*(UwallsMR[NA-1][j+1]*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + VwallsMR[NA-1][j+1]*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]))
							)
							 + VwallsMU[NA][j]/MESH.MU[NA][j][1]	
							;

	}

	//Esquina abajo izquierda
	DivergenceMU[0][0] = (1.0/(2.0*MESH.VolMP[0][0]))*(
						  MESH.SupMP[0][0][0]*(UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0]) + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0]))
						   + MESH.SupMP[0][0][1]*(Upres[0][0]*cos(MESH.AngleMU[1][0]) + Vpres[0][0]*sin(MESH.AngleMU[1][0]))
						   + 0.50*MESH.SupMP[0][0][2]*(UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))
						   + 0.50*MESH.SupMP[0][0][3]*(UwallsMR[0][1]*cos(0.50*PI + MESH.AngleMR[0][1]) + VwallsMR[0][1]*sin(0.50*PI + MESH.AngleMR[0][1]))
						   )
							 + VwallsMU[0][0]/MESH.MU[0][0][1]	
							;

	//Esquina arriba izquierda
		DivergenceMU[0][NRad-1] = (1.0/(2.0*MESH.VolMP[0][NRad-1]))*(
						     MESH.SupMP[0][NRad-1][0]*(UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))
						   + MESH.SupMP[0][NRad-1][1]*(Upres[0][NRad-1]*cos(MESH.AngleMU[1][NRad-1]) + Vpres[0][NRad-1]*sin(MESH.AngleMU[1][NRad-1]))
						   + 0.50*MESH.SupMP[0][NRad-1][2]*(UwallsMR[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + VwallsMR[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad-1]))
						   + 0.50*MESH.SupMP[0][NRad-1][3]*(UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad]))
						   )
							 + VwallsMU[0][NRad-1]/MESH.MU[0][NRad-1][1]	
							;

	//Esquina abajo derecha
		DivergenceMU[NA][0] = (1.0/(2.0*MESH.VolMP[NA-1][0]))*(
							  MESH.SupMP[NA-1][0][0]*(Upres[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0]) + Vpres[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0]))
							+ MESH.SupMP[NA-1][0][1]*(UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0]) + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0]))
							+ 0.50*MESH.SupMP[NA-1][0][2]*(UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))
							+ 0.50*MESH.SupMP[NA-1][0][3]*(UwallsMR[NA-1][1]*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + VwallsMR[NA-1][1]*sin(0.50*PI + MESH.AngleMR[NA-1][1]))
							)
							 + VwallsMU[NA][0]/MESH.MU[NA][0][1]	
							;

	//Esquina arriba derecha
			DivergenceMU[NA][NRad-1] = (1.0/(2.0*MESH.VolMP[NA-1][NRad-1]))*(
							  MESH.SupMP[NA-1][NRad-1][0]*(Upres[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + Vpres[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1]))
							+ MESH.SupMP[NA-1][NRad-1][1]*(UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][NRad-1]))
							+ 0.50*MESH.SupMP[NA-1][NRad-1][2]*(UwallsMR[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + VwallsMR[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]))
							+ 0.50*MESH.SupMP[NA-1][NRad-1][3]*(UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))
							)
							 + VwallsMU[NA][NRad-1]/MESH.MU[NA][NRad-1][1]	
							;


	//Nodos R

	//Centro
	for(i = 0; i < NA; i++){
		for(j = 1; j < NRad-1; j++){
			DivergenceMR[i][j] =  0.50*(Divergence[i][j] + Divergence[i][j-1]);
		}
	}

	for(i = 0; i < NA; i++){
		
		//Parte abajo
		DivergenceMR[i][0] = (1.0/(2.0*MESH.VolMP[i][0]))*(
							   0.50*MESH.SupMP[i][0][0]*(UwallsMU[i][0]*cos(PI + MESH.AngleMU[i][0]) + VwallsMU[i][0]*sin(PI + MESH.AngleMU[i][0]))
							 + 0.50*MESH.SupMP[i][0][1]*(UwallsMU[i+1][0]*cos(MESH.AngleMU[i+1][0]) + VwallsMU[i+1][0]*sin(MESH.AngleMU[i+1][0]))
							 + MESH.SupMP[i][0][2]*(UwallsMR[i][0]*cos(1.50*PI + MESH.AngleMR[i][0]) + VwallsMR[i][0]*sin(1.50*PI + MESH.AngleMR[i][0]))
							 + MESH.SupMP[i][0][3]*(Upres[i][0]*cos(0.50*PI + MESH.AngleMR[i][1]) + Vpres[i][0]*sin(0.50*PI + MESH.AngleMR[i][1]))
							 )
							+ VwallsMR[i][0]/MESH.MR[i][0][1]
							;

		//Parte arriba
		DivergenceMR[i][NRad] = (1.0/(2.0*MESH.VolMP[i][NRad-1]))*(
							   0.50*MESH.SupMP[i][NRad-1][0]*(UwallsMU[i][NRad-1]*cos(PI + MESH.AngleMU[i][NRad-1]) + VwallsMU[i][NRad-1]*sin(PI + MESH.AngleMU[i][NRad-1]))
							 + 0.50*MESH.SupMP[i][NRad-1][1]*(UwallsMU[i+1][NRad-1]*cos(MESH.AngleMU[i+1][NRad-1]) + VwallsMU[i+1][NRad-1]*sin(MESH.AngleMU[i+1][NRad-1]))
							 + MESH.SupMP[i][NRad-1][2]*(Upres[i][NRad-1]*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + Vpres[i][NRad-1]*sin(1.50*PI + MESH.AngleMR[i][NRad]))
							 + MESH.SupMP[i][NRad-1][3]*(UwallsMR[i][NRad]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + VwallsMR[i][NRad]*sin(0.50*PI + MESH.AngleMR[i][NRad]))
							 )
							+ VwallsMR[i][NRad]/MESH.MR[i][NRad][1]
							;

	}

	//Esquina abajo izquierda
	DivergenceMR[0][0] = (1.0/(2.0*MESH.VolMP[0][0]))*(
							   0.50*MESH.SupMP[0][0][0]*(UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0]) + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0]))
							 + 0.50*MESH.SupMP[0][0][1]*(UwallsMU[1][0]*cos(MESH.AngleMU[1][0]) + VwallsMU[1][0]*sin(MESH.AngleMU[1][0]))
							 + MESH.SupMP[0][0][2]*(UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0]) + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0]))
							 + MESH.SupMP[0][0][3]*(Upres[0][0]*cos(0.50*PI + MESH.AngleMR[0][1]) + Vpres[0][0]*sin(0.50*PI + MESH.AngleMR[0][1]))
							 )
							+ VwallsMR[0][0]/MESH.MR[0][0][1]
							;

	//Esquina arriba izquierda
	DivergenceMR[0][NRad] = (1.0/(2.0*MESH.VolMP[0][NRad-1]))*(
							   0.50*MESH.SupMP[0][NRad-1][0]*(UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1]))
							 + 0.50*MESH.SupMP[0][NRad-1][1]*(UwallsMU[1][NRad-1]*cos(MESH.AngleMU[1][NRad-1]) + VwallsMU[1][NRad-1]*sin(MESH.AngleMU[1][NRad-1]))
							 + MESH.SupMP[0][NRad-1][2]*(Upres[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + Vpres[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad]))
							 + MESH.SupMP[0][NRad-1][3]*(UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad]))
							 )
							+ VwallsMR[0][NRad]/MESH.MR[0][NRad][1]
							;

	//Esquina abajo derecha
	DivergenceMR[NA-1][0] = (1.0/(2.0*MESH.VolMP[NA-1][0]))*(
							   0.50*MESH.SupMP[NA-1][0][0]*(UwallsMU[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0]) + VwallsMU[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0]))
							 + 0.50*MESH.SupMP[NA-1][0][1]*(UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0]) + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0]))
							 + MESH.SupMP[NA-1][0][2]*(UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0]))
							 + MESH.SupMP[NA-1][0][3]*(Upres[NA-1][0]*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + Vpres[NA-1][0]*sin(0.50*PI + MESH.AngleMR[NA-1][1]))
							 )
							+ VwallsMR[NA-1][0]/MESH.MR[NA-1][0][1]
							;

	//Esquina arriba derecha
	DivergenceMR[NA-1][NRad] = (1.0/(2.0*MESH.VolMP[NA-1][NRad-1]))*(
							   0.50*MESH.SupMP[NA-1][NRad-1][0]*(UwallsMU[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + VwallsMU[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1]))
							 + 0.50*MESH.SupMP[NA-1][NRad-1][1]*(UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + VwallsMU[NA][NRad-1]*sin(MESH.AngleMU[NA][NRad-1]))
							 + MESH.SupMP[NA-1][NRad-1][2]*(Upres[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + Vpres[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad]))
							 + MESH.SupMP[NA-1][NRad-1][3]*(UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]))
							 )
							+ VwallsMR[NA-1][NRad]/MESH.MR[NA-1][NRad][1]
							;

}