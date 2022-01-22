//Cálculo de la viscosidad base de cada nodo con la ley de Sutherland
void Solver::Get_BaseViscosity(){
int i, j;
double muO2, muH2, muH2O;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			muO2 = mu0O2*((ToO2 + TsO2)/(Tpres[i][j] + TsO2))*pow(Tpres[i][j]/ToO2,1.50);

			muH2 = mu0H2*((ToH2 + TsH2)/(Tpres[i][j] + TsH2))*pow(Tpres[i][j]/ToH2,1.50);

			muH2O = mu0H2O*((ToH2O + TsH2O)/(Tpres[i][j] + TsH2O))*pow(Tpres[i][j]/ToH2O,1.50);

			muBase[i][j] = FracO2*muO2 + FracH2*muH2 + FracH2O*muH2O;
		}
	}
	
}

//Cálculo de la entalpía de los gases para una determinada temperatura
void Solver::Get_Enthalpy(){
int i, j, k;
double EnthalpyO2, EnthalpyH2, EnthalpyH2O;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			//Cálculo Cp O2
			if(Tpres[i][j] >= 100.0 && Tpres[i][j] <= 700.0){ k = 0; }//Rango 1 O2 (K)
			else if(Tpres[i][j] > 700.0 && Tpres[i][j] <= 2000.0){ k = 1; } //Rango 2 O2 (K)
			else{ k = 2; } //Rango 3 O2 (K)
				
			EnthalpyO2 = 1000.0*(1000.0/MmO2)*(HfO2 + JanafO2[k][0]*(Tpres[i][j]/1000.0) + 0.50*JanafO2[k][1]*pow(Tpres[i][j]/1000.0,2.0) + (1.0/3.0)*JanafO2[k][2]*pow(Tpres[i][j]/1000.0,3.0) + 0.25*JanafO2[k][3]*pow(Tpres[i][j]/1000.0,4.0) - JanafO2[k][4]/(Tpres[i][j]/1000.0) + JanafO2[k][5] - JanafO2[k][7]); //(J/Kg)


			//Cálculo Cp H2
			if(Tpres[i][j] >= 298.0 && Tpres[i][j] <= 1000.0){ k = 0; } //Rango 1 H2 (K)
			else if(Tpres[i][j] > 1000.0 && Tpres[i][j] <= 2500.0){ k = 1; } //Rango 2 H2 (K)
			else{ k = 2; } //Rango 3 H2 (K)
				
			EnthalpyH2 = 1000.0*(1000.0/MmH2)*(HfH2 + JanafH2[k][0]*(Tpres[i][j]/1000.0) + 0.50*JanafH2[k][1]*pow(Tpres[i][j]/1000.0,2.0) + (1.0/3.0)*JanafH2[k][2]*pow(Tpres[i][j]/1000.0,3.0) + 0.25*JanafH2[k][3]*pow(Tpres[i][j]/1000.0,4.0) - JanafH2[k][4]/(Tpres[i][j]/1000.0) + JanafH2[k][5] - JanafH2[k][7]); //(J/Kg)
			 

			//Cálculo Cp H2O
			if(Tpres[i][j] >= 500.0 && Tpres[i][j] <= 1700.0){ k = 0; } //Rango 1 H2O (K)
			else{ k = 1; } //Rango 2 H2O (K)
				
			EnthalpyH2O = 1000.0*(1000.0/MmH2O)*(HfH2O + JanafH2O[k][0]*(Tpres[i][j]/1000.0) + 0.50*JanafH2O[k][1]*pow(Tpres[i][j]/1000.0,2.0) + (1.0/3.0)*JanafH2O[k][2]*pow(Tpres[i][j]/1000.0,3.0) + 0.25*JanafH2O[k][3]*pow(Tpres[i][j]/1000.0,4.0) - JanafH2O[k][4]/(Tpres[i][j]/1000.0) + JanafH2O[k][5] - JanafH2O[k][7]); //(J/Kg)
			
			Hpres[i][j] = FracO2*EnthalpyO2 + FracH2*EnthalpyH2 + FracH2O*EnthalpyH2O;
		}
	}
}

//Calculo del Cp de cada uno de los nodos del dominio
void Solver::Get_CpCombustion(){
int i, j;
double CpO2, CpH2, CpH2O;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			//Cálculo Cp O2
			if(Tpres[i][j] >= 100.0 && Tpres[i][j] <= 700.0){ //K
				//Rango 1 O2
				CpO2 = (1000.0/MmO2)*(JanafO2[0][0] + JanafO2[0][1]*(Tpres[i][j]/1000.0) + JanafO2[0][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafO2[0][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafO2[0][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}
			else if(Tpres[i][j] > 700.0 && Tpres[i][j] <= 2000.0){ //K
				//Rango 2 O2
				CpO2 = (1000.0/MmO2)*(JanafO2[1][0] + JanafO2[1][1]*(Tpres[i][j]/1000.0) + JanafO2[1][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafO2[1][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafO2[1][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}
			else{
				//Rango 3 O2
				CpO2 = (1000.0/MmO2)*(JanafO2[2][0] + JanafO2[2][1]*(Tpres[i][j]/1000.0) + JanafO2[2][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafO2[2][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafO2[2][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}

			//Cálculo Cp H2
			if(Tpres[i][j] >= 298.0 && Tpres[i][j] <= 1000.0){ //K
				//Rango 1 H2
				CpH2 = (1000.0/MmH2)*(JanafH2[0][0] + JanafH2[0][1]*(Tpres[i][j]/1000.0) + JanafH2[0][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafH2[0][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafH2[0][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}
			else if(Tpres[i][j] > 1000.0 && Tpres[i][j] <= 2500.0){ //K
				//Rango 2 H2
				CpH2 = (1000.0/MmH2)*(JanafH2[1][0] + JanafH2[1][1]*(Tpres[i][j]/1000.0) + JanafH2[1][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafH2[1][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafH2[1][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}
			else{
				//Rango 3 H2
				CpH2 = (1000.0/MmH2)*(JanafH2[2][0] + JanafH2[2][1]*(Tpres[i][j]/1000.0) + JanafH2[2][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafH2[2][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafH2[2][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}

			//Cálculo Cp H2O
			if(Tpres[i][j] >= 500.0 && Tpres[i][j] <= 1700.0){ //K
				//Rango 1 H2O
				CpH2O = (1000.0/MmH2O)*(JanafH2O[0][0] + JanafH2O[0][1]*(Tpres[i][j]/1000.0) + JanafH2O[0][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafH2O[0][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafH2O[0][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}
			else{
				//Rango 2 H2O
				CpH2O = (1000.0/MmH2O)*(JanafH2O[1][0] + JanafH2O[1][1]*(Tpres[i][j]/1000.0) + JanafH2O[1][2]*pow(Tpres[i][j]/1000.0,2.0) + JanafH2O[1][3]*pow(Tpres[i][j]/1000.0,3.0) + JanafH2O[1][4]/pow(Tpres[i][j]/1000.0,2.0)); //(J/KgK)
			}
			
			Cp[i][j] = FracO2*CpO2 + FracH2*CpH2 + FracH2O*CpH2O;

		}
	}

}

//Pasar todos los coeficientes termoqímicos de las especies a sus matrices
void Solver::GetJanaf(ReadData R1){
int i;

	//Coeficientes de O2 Y H2
	for(i = 0; i < 24; i++){
		JanafO2[i/8][i - 8*(i/8)] = R1.JannafO2[i];
		JanafH2[i/8][i - 8*(i/8)] = R1.JannafH2[i];
	}

	//Coeficientes de H2O
	for(i = 0; i < 16; i++){
		JanafH2O[i/8][i - 8*(i/8)] = R1.JannafH2O[i];
	}

}


//Cálculo de la conductividad térmica de los gases en cada nodo y paredes
void Solver::Get_K(Mesher MESH){
int i, j;
double Predict;

	//Cálculo de K en los nodos
	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			K[i][j] = (Cp[i][j]*muTotal[i][j])/Pr;
		}
	}
//K[i][j] = (Cp[i][j]*muTotal[i][j])/Pr;
	//Cálculo de K en las paredes de los volúmenes de control
	
	//Nodos U
	for(j = 0; j < NRad; j++){
		KwallsMU[0][j] = KLeft[j]; //Parte izquierda
		KwallsMU[NA][j] = KRight[j]; //Parte derecha
		for(i = 1; i < NA; i++){
			KwallsMU[i][j] = 0.50*(K[i][j] + K[i-1][j]); //Centro
		}
	}
	
	//Nodos V
	for(i = 0; i < NA; i++){
		KwallsMV[i][0] = K[i][0]; //Parte abajo
		KwallsMV[i][NRad] = KUp[i]; //Parte arriba
		for(j = 1; j < NRad; j++){
			KwallsMV[i][j] = 0.50*(K[i][j] + K[i][j-1]); //Centro
		}
	}
}

//Cálculo del término difusivo de la ecuación de conservación de la energía
void Solver::Get_EnergyDifusive(Mesher MESH){
int i, j;
	
	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			EnergyDifusive[i][j] = (1.0/MESH.VolMP[i][j])*(
								   (KwallsMU[i][j]*(Tpres[i][j] - Tpres[i-1][j])*(MESH.SupMP[i][j][0])*cos(PI + MESH.AngleMU[i][j]))/MESH.DeltasMU[i][j][0]
								 + (KwallsMU[i+1][j]*(Tpres[i+1][j] - Tpres[i][j])*MESH.SupMP[i][j][1]*cos(MESH.AngleMU[i+1][j]))/MESH.DeltasMU[i+1][j][0] 
								 + (KwallsMV[i][j]*(Tpres[i][j] - Tpres[i][j-1])*MESH.SupMP[i][j][2]*sin(1.50*PI + MESH.AngleMR[i][j]))/MESH.DeltasMR[i][j][1] 
								 + (KwallsMV[i][j+1]*(Tpres[i][j+1] - Tpres[i][j])*MESH.SupMP[i][j][3]*sin(0.50*PI + MESH.AngleMR[i][j+1]))/MESH.DeltasMR[i][j+1][1]
								 );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		EnergyDifusive[0][j] = (1.0/MESH.VolMP[0][j])*(
							   KwallsMU[0][j]*(Tpres[0][j] - Tleft[j])*(MESH.SupMP[0][j][0]/MESH.DeltasMU[0][j][0])*cos(PI + MESH.AngleMU[0][j]) 
							 + KwallsMU[1][j]*(Tpres[1][j] - Tpres[0][j])*(MESH.SupMP[0][j][1]/MESH.DeltasMU[1][j][0])*cos(MESH.AngleMU[1][j]) 
							 + KwallsMV[0][j]*(Tpres[0][j] - Tpres[0][j-1])*(MESH.SupMP[0][j][2]/MESH.DeltasMR[0][j][1])*sin(1.50*PI + MESH.AngleMR[0][j]) 
							 + KwallsMV[0][j+1]*(Tpres[0][j+1] - Tpres[0][j])*(MESH.SupMP[0][j][3]/MESH.DeltasMR[0][j+1][1])*sin(0.50*PI + MESH.AngleMR[0][j+1])
							 );

		//Parte derecha
		EnergyDifusive[NA-1][j] = (1.0/MESH.VolMP[NA-1][j])*(
								  KwallsMU[NA-1][j]*(Tpres[NA-1][j] - Tpres[NA-2][j])*(MESH.SupMP[NA-1][j][0]/MESH.DeltasMU[NA-1][j][0])*cos(PI + MESH.AngleMU[NA-1][j]) 
								+ KwallsMU[NA][j]*(Tright[j] - Tpres[NA-1][j])*(MESH.SupMP[NA-1][j][1]/MESH.DeltasMU[NA][j][0])*cos(MESH.AngleMU[NA][j]) 
								+ KwallsMV[NA-1][j]*(Tpres[NA-1][j] - Tpres[NA-1][j-1])*(MESH.SupMP[NA-1][j][2]/MESH.DeltasMR[NA-1][j][1])*sin(1.50*PI + MESH.AngleMR[NA-1][j]) 
								+ KwallsMV[NA-1][j+1]*(Tpres[NA-1][j+1] - Tpres[NA-1][j])*(MESH.SupMP[NA-1][j][3]/MESH.DeltasMR[NA-1][j+1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])
								);
	}
	
	for(i = 1; i < NA-1; i++){
		//Parte abajo
		EnergyDifusive[i][0] = (1.0/MESH.VolMP[i][0])*(
							   KwallsMU[i][0]*(Tpres[i][0] - Tpres[i-1][0])*(MESH.SupMP[i][0][0]/MESH.DeltasMU[i][0][0])*cos(PI + MESH.AngleMU[i][0]) 
							 + KwallsMU[i+1][0]*(Tpres[i+1][0] - Tpres[i][0])*(MESH.SupMP[i][0][1]/MESH.DeltasMU[i+1][0][0])*cos(MESH.AngleMU[i+1][0]) 
							 + KwallsMV[i][0]*(Tpres[i][0] - Tdown[i])*(MESH.SupMP[i][0][2]/MESH.DeltasMR[i][0][1])*sin(1.50*PI + MESH.AngleMR[i][0]) 
							 + KwallsMV[i][1]*(Tpres[i][1] - Tpres[i][0])*(MESH.SupMP[i][0][3]/MESH.DeltasMR[i][1][1])*sin(0.50*PI + MESH.AngleMR[i][1])
							 );

		//Parte arriba
		EnergyDifusive[i][NRad-1] = (1.0/MESH.VolMP[i][NRad-1])*(
								    KwallsMU[i][NRad-1]*(Tpres[i][NRad-1] - Tpres[i-1][NRad-1])*(MESH.SupMP[i][NRad-1][0]/MESH.DeltasMU[i][NRad-1][0])*cos(PI + MESH.AngleMU[i][NRad-1]) 
								  + KwallsMU[i+1][NRad-1]*(Tpres[i+1][NRad-1] - Tpres[i][NRad-1])*(MESH.SupMP[i][NRad-1][1]/MESH.DeltasMU[i+1][NRad-1][0])*cos(MESH.AngleMU[i+1][NRad-1]) 
								  + KwallsMV[i][NRad-1]*(Tpres[i][NRad-1] - Tpres[i][NRad-2])*(MESH.SupMP[i][NRad-1][2]/MESH.DeltasMR[i][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[i][NRad-1]) 
								  + KwallsMV[i][NRad]*(Tup[i] - Tpres[i][NRad-1])*(MESH.SupMP[i][NRad-1][3]/MESH.DeltasMR[i][NRad][1])*sin(0.50*PI + MESH.AngleMR[i][NRad])
								  );
	}
	
	//Esquina abajo izquierda
	EnergyDifusive[0][0] = (1.0/MESH.VolMP[0][0])*(
						   KwallsMU[0][0]*(Tpres[0][0] - Tleft[0])*(MESH.SupMP[0][0][0]/MESH.DeltasMU[0][0][0])*cos(PI + MESH.AngleMU[0][0]) 
						 + KwallsMU[1][0]*(Tpres[1][0] - Tpres[0][0])*(MESH.SupMP[0][0][1]/MESH.DeltasMU[1][0][0])*cos(MESH.AngleMU[1][0]) 
						 + KwallsMV[0][0]*(Tpres[0][0] - Tdown[0])*(MESH.SupMP[0][0][2]/MESH.DeltasMR[0][0][1])*sin(1.50*PI + MESH.AngleMR[0][0]) 
						 + KwallsMV[0][1]*(Tpres[0][1] - Tpres[0][0])*(MESH.SupMP[0][0][3]/MESH.DeltasMR[0][1][1])*sin(0.50*PI + MESH.AngleMR[0][1])
						 );

	//Esquina arriba izquierda
	EnergyDifusive[0][NRad-1] = (1.0/MESH.VolMP[0][NRad-1])*(
						  	    KwallsMU[0][NRad-1]*(Tpres[0][NRad-1] - Tleft[NRad-1])*(MESH.SupMP[0][NRad-1][0]/MESH.DeltasMU[0][NRad-1][0])*cos(PI + MESH.AngleMU[0][NRad-1]) 
						 	  + KwallsMU[1][NRad-1]*(Tpres[1][NRad-1] - Tpres[0][NRad-1])*(MESH.SupMP[0][NRad-1][1]/MESH.DeltasMU[1][NRad-1][0])*cos(MESH.AngleMU[1][NRad-1]) 
						 	  + KwallsMV[0][NRad-1]*(Tpres[0][NRad-1] - Tpres[0][NRad-2])*(MESH.SupMP[0][NRad-1][2]/MESH.DeltasMR[0][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[0][NRad-1]) 
						 	  + KwallsMV[0][NRad]*(Tup[0] - Tpres[0][NRad-1])*(MESH.SupMP[0][NRad-1][3]/MESH.DeltasMR[0][NRad][1])*sin(0.50*PI + MESH.AngleMR[0][NRad])
						 	  );

	//Esquina abajo derecha
	EnergyDifusive[NA-1][0] = (1.0/MESH.VolMP[NA-1][0])*(
							  KwallsMU[NA-1][0]*(Tpres[NA-1][0] - Tpres[NA-2][0])*(MESH.SupMP[NA-1][0][0]/MESH.DeltasMU[NA-1][0][0])*cos(PI + MESH.AngleMU[NA-1][0]) 
							+ KwallsMU[NA][0]*(Tright[0] - Tpres[NA-1][0])*(MESH.SupMP[NA-1][0][1]/MESH.DeltasMU[NA][0][0])*cos(MESH.AngleMU[NA][0]) 
							+ KwallsMV[NA-1][0]*(Tpres[NA-1][0] - Tdown[NA-1])*(MESH.SupMP[NA-1][0][2]/MESH.DeltasMR[NA-1][0][1])*sin(1.50*PI + MESH.AngleMR[NA-1][0]) 
							+ KwallsMV[NA-1][1]*(Tpres[NA-1][1] - Tpres[NA-1][0])*(MESH.SupMP[NA-1][0][3]/MESH.DeltasMR[NA-1][1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1])
							);

	//Esquina arriba derecha
	EnergyDifusive[NA-1][NRad-1] = (1.0/MESH.VolMP[NA-1][NRad-1])*(
								   KwallsMU[NA-1][NRad-1]*(Tpres[NA-1][NRad-1] - Tpres[NA-2][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]/MESH.DeltasMU[NA-1][NRad-1][0])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) 
								 + KwallsMU[NA][NRad-1]*(Tright[NRad-1] - Tpres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][1]/MESH.DeltasMU[NA][NRad-1][0])*cos(MESH.AngleMU[NA][NRad-1]) 
								 + KwallsMV[NA-1][NRad-1]*(Tpres[NA-1][NRad-1] - Tpres[NA-1][NRad-2])*(MESH.SupMP[NA-1][NRad-1][2]/MESH.DeltasMR[NA-1][NRad-1][1])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) 
								 + KwallsMV[NA-1][NRad]*(Tup[NA-1] - Tpres[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][3]/MESH.DeltasMR[NA-1][NRad][1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])
								 );
}

void Solver::Get_eWalls(Mesher MESH){
int i, j;

		//Nodos U
		for(j = 0; j < NRad; j++){
			eWallsMU[0][j] = eleft[j];
			eWallsMU[NA][j] = eright[j];
			for(i = 1; i < NA; i++){
					eWallsMU[i][j] = sqrt(epres[i][j])*sqrt(epres[i-1][j]);
			}
		}

		//Nodos R
		for(i = 0; i < NA; i++){
			eWallsMR[i][0] = edown[i];
			eWallsMR[i][NRad] = eup[i];
			for(j = 1; j < NRad; j++){
				eWallsMR[i][j] = sqrt(epres[i][j])*sqrt(epres[i][j-1]);
			}
		}

}
//Cálculo del término convectivo de la energía
void Solver::Get_EnergyConvective(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			EnergyConvective[i][j] = -(1.0/MESH.VolMP[i][j])*(
									 MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*RhoWallsMU[i][j]*eWallsMU[i][j]
								   + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*RhoWallsMU[i+1][j]*eWallsMU[i+1][j]
								   + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*RhoWallsMR[i][j]*eWallsMR[i][j]
								   + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*RhoWallsMR[i][j+1]*eWallsMR[i][j+1]
								   );
		}
	}

}

//Cálculo del término de presión en la ecuación de conservación de la energía
void Solver::Get_EnergyPressureTerm(Mesher MESH){
int i, j;
		
		for(i = 0; i < NA; i++){
			for(j = 0; j < NRad; j++){
				EnergyPressureTerm[i][j] = -(1.0/MESH.VolMP[i][j])*(
										   MESH.SupMP[i][j][0]*PwallsU[i][j]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))
										 + MESH.SupMP[i][j][1]*PwallsU[i+1][j]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))
										 + MESH.SupMP[i][j][2]*PwallsR[i][j]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))	
										 + MESH.SupMP[i][j][3]*PwallsR[i][j+1]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
										 );

				
			}
		}
}

//Cálculo del término viscoso de la ecuaciń de conservación de la energía
void Solver::Get_EnergyViscousTerm(Mesher MESH){
int i, j;


	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			EnergyViscous[i][j] = (1.0/MESH.VolMP[i][j])*(
									MESH.SupMP[i][j][0]*(
										UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j])*TauZZ_mu[i][j] 
									  + UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j])*TauRZ_mu[i][j] 	
									  + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j])*TauRZ_mu[i][j] 
									  + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j])*TauRR_mu[i][j] 
									)
								  + MESH.SupMP[i][j][1]*(
										UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j])*TauZZ_mu[i+1][j] 
									  + UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j])*TauRZ_mu[i+1][j] 	
									  + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j])*TauRZ_mu[i+1][j] 
									  + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j])*TauRR_mu[i+1][j] 
									)	
								  + MESH.SupMP[i][j][2]*(
										UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j])*TauZZ_mr[i][j] 
									  + UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j])*TauRZ_mr[i][j] 	
									  + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j])*TauRZ_mr[i][j] 
									  + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j])*TauRR_mr[i][j] 
									)
								  + MESH.SupMP[i][j][3]*(
										UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1])*TauZZ_mr[i][j+1] 
									  + UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1])*TauRZ_mr[i][j+1] 	
									  + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1])*TauRZ_mr[i][j+1] 
									  + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1])*TauRR_mr[i][j+1] 
									)

								);
		}
	}
	
}

//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveW(Mesher MESH){
int i, j;
	
	for(i = 0; i <  NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumDifusiveW[i][j] = (1.0/MESH.VolMP[i][j])*(MESH.SupMP[i][j][1]*muWallsMU[i+1][j]*GradW_DxMU[i+1][j] - MESH.SupMP[i][j][0]*muWallsMU[i][j]*GradW_DxMU[i][j])
									+ (1.0/pow(MESH.MP[i][j][1],2.0))*(muWallsMR[i][j+1]*pow(MESH.MR[i][j+1][1],3.0)*GradW_rMR[i][j+1] - muWallsMR[i][j]*pow(MESH.MR[i][j][1],3.0)*GradW_rMR[i][j])
									- RhoPres[i][j]*Upres[i][j]*Wpres[i][j]/MESH.MP[i][j][1]
									;


		/*	MatrixNueva[i][j] =  (TauRR_mr[i][j+1] - TauRR_mr[i][j])/MESH.DeltasMP[i][j][1]
									+ (TauRZ_mu[i+1][j] - TauRZ_mu[i][j])/MESH.DeltasMP[i][j][0]
									+ (1.0/MESH.MP[i][j][1])*(RhoPres[i][j]*Vpres[i][j]*Vpres[i][j])
									+ (1.0/MESH.MP[i][j][1])*0.50*(TauRR_mr[i][j] + TauRR_mr[i][j+1])
									+ (1.0/MESH.MP[i][j][1])*((2.0*muTotal[i][j])/(3.0*ReynoldsM))*(-0.50*(GradU_DxMU[i+1][j] + GradU_DxMU[i][j]) -0.50*(GradV_DyMR[i][j+1] + GradV_DyMR[i][j]) + (2.0*Vpres[i][j])/MESH.MP[i][j][1])
									;*/
		}
	}


}

//Cálculo del término convectivo de la ecuaciones de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomentumConvectiveW(Mesher MESH){
int i, j;

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumConvectiveW[i][j] = -(1.0/MESH.VolMP[i][j])*(
									    MESH.SupMP[i][j][0]*(UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j]) + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j]))*RhoWallsMU[i][j]*WwallsMU[i][j]
									  + MESH.SupMP[i][j][1]*(UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j]) + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j]))*RhoWallsMU[i+1][j]*WwallsMU[i+1][j]
									  + MESH.SupMP[i][j][2]*(UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))*RhoWallsMR[i][j]*WwallsMR[i][j]
									  + MESH.SupMP[i][j][3]*(UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))*RhoWallsMR[i][j+1]*WwallsMR[i][j+1]
									  );
		}
	}
	
}