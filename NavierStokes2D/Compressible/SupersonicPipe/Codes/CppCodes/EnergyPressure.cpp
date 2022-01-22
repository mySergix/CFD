//Cálculo del término de presión en la ecuación de conservación de la energía
void Solver::Get_EnergyPressureTerm(Mesher MESH){
int i, j;
		
		//Centro
		for(i = 1; i < NA-1; i++){
			for(j = 1; j < NRad-1; j++){
				EnergyPressureTerm[i][j] = -(1.0/MESH.VolMP[i][j])*P[i][j]*(
										 + (MESH.SupMP[i][j][0]/MESH.DeltasMU[i][j][0])*(Upres[i][j] - Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) 
										 + (MESH.SupMP[i][j][1]/MESH.DeltasMU[i+1][j][0])*(Upres[i+1][j] - Upres[i][j])*cos(MESH.AngleMU[i+1][j]) 
										 + (MESH.SupMP[i][j][2]/MESH.DeltasMR[i][j][1])*(Upres[i][j] - Upres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j]) 
										 + (MESH.SupMP[i][j][3]/MESH.DeltasMR[i][j+1][1])*(Upres[i][j+1] - Upres[i][j])*sin(0.50*PI + MESH.AngleMR[i][j+1]) 
										 );	
			}
		}

		for(j = 1; j < NRad-1; j++){
			//Parte izquierda
			EnergyPressureTerm[0][j] = -(1.0/MESH.VolMP[0][j])*P[0][j]*(
										 + (MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*(Upres[0][j] - Uleft[j])*cos(PI + MESH.AngleMU[0][j]) 
										 + (MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*(Upres[1][j] - Upres[0][j])*cos(MESH.AngleMU[1][j]) 
										 + (MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*(Upres[0][j] - Upres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j]) 
										 + (MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*(Upres[0][j+1] - Upres[0][j])*sin(0.50*PI + MESH.AngleMR[0][j+1]) 
										 );	

			//Parte derecha
			EnergyPressureTerm[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*P[NA-1][j]*(
										 + (MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*(Upres[NA-1][j] - Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) 
										 + (MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*(Uright[j] - Upres[NA-1][j])*cos(MESH.AngleMU[NA][j]) 
										 + (MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*(Upres[NA-1][j] - Upres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j]) 
										 + (MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*(Upres[NA-1][j+1] - Upres[NA-1][j])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1]) 
										 );	
		}

		for(i = 1; i < NA-1; i++){
			//Parte abajo
			EnergyPressureTerm[i][0] = -(1.0/MESH.VolMP[i][0])*P[i][0]*(
										 + (MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*(Upres[i][0] - Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) 
										 + (MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*(Upres[i+1][0] - Upres[i][0])*cos(MESH.AngleMU[i+1][0]) 
										 + (MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*(Upres[i][0] - Udown[i])*sin(1.50*PI + MESH.AngleMR[i][0]) 
										 + (MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*(Upres[i][1] - Upres[i][0])*sin(0.50*PI + MESH.AngleMR[i][1]) 
										 );	

			//Parte arriba
			EnergyPressureTerm[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*P[i][NRad-1]*(
										 + (MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*(Upres[i][NRad-1] - Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) 
										 + (MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*(Upres[i+1][NRad-1] - Upres[i][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) 
										 + (MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*(Upres[i][NRad-1] - Upres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1]) 
										 + (MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*(Uup[i] - Upres[i][NRad-1])*sin(0.50*PI + MESH.AngleMR[i][NRad]) 
										 );	
		}

		//Esquina abajo izquierda
		EnergyPressureTerm[0][0] = -(1.0/MESH.VolMP[0][0])*P[0][0]*(
										 + (MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*(Upres[0][0] - Uleft[0])*cos(PI + MESH.AngleMU[0][0]) 
										 + (MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*(Upres[1][0] - Upres[0][0])*cos(MESH.AngleMU[1][0]) 
										 + (MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*(Upres[0][0] - Udown[0])*sin(1.50*PI + MESH.AngleMR[0][0]) 
										 + (MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*(Upres[0][1] - Upres[0][0])*sin(0.50*PI + MESH.AngleMR[0][1]) 
										 );	

		//Esquina arriba izquierda
		EnergyPressureTerm[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*P[0][NRad-1]*(
										 + (MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*(Upres[0][NRad-1] - Uleft[NRad-1])*cos(PI + MESH.AngleMU[0][NRad-1]) 
										 + (MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*(Upres[1][NRad-1] - Upres[0][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) 
										 + (MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*(Upres[0][NRad-1] - Upres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]) 
										 + (MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*(Uup[0] - Upres[0][NRad-1])*sin(0.50*PI + MESH.AngleMR[0][NRad]) 
										 );	

		//Esquina abajo derecha
		EnergyPressureTerm[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*P[NA-1][0]*(
									 + (MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*(Upres[NA-1][0] - Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) 
									 + (MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*(Uright[0] - Upres[NA-1][0])*cos(MESH.AngleMU[NA][0]) 
									 + (MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*(Upres[NA-1][0] - Udown[NA-1])*sin(1.50*PI + MESH.AngleMR[NA-1][0]) 
									 + (MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*(Upres[NA-1][1] - Upres[NA-1][0])*sin(0.50*PI + MESH.AngleMR[NA-1][1]) 
									 );	

		//Esquina arriba derecha
		EnergyPressureTerm[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*P[NA-1][NRad-1]*(
										 + (MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*(Upres[NA-1][NRad-1] - Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) 
										 + (MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*(Uright[NRad-1] - Upres[NA-1][NRad-1])*cos(MESH.AngleMU[NA][NRad-1]) 
										 + (MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*(Upres[NA-1][NRad-1] - Upres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) 
										 + (MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*(Uup[NA-1] - Upres[NA-1][NRad-1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad]) 
										 );	
}