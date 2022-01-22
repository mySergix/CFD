//Cálculo de las presiones en las paredes de los volúmenes de control
void Solver::Get_Pwalls(Mesher MESH){
int i, j;
double Predict;

	//Nodos U
	for(j = 0; j < NRad; j++){
		PwallsU[0][j] = Pleft[j];
		PwallsU[NA][j] = Pright[j];
		
		for(i = 1; i < NA; i++){
			PwallsU[i][j] = 0.50*(P[i][j] + P[i-1][j]);
		}
	}
	
	//Nodos V
	for(i = 0; i < NA; i++){
		PwallsR[i][0] = P[i][0];
		PwallsR[i][NRad] = Pup[i];

		for(j = 1; j < NRad; j++){
				PwallsR[i][j] = 0.50*(P[i][j] + P[i][j-1]);
			}
	}

}



//Cálculo del término de presión en la ecuación de la energía
void Solver::Get_PressureTerm(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			EnergyPressTerm[i][j] = (1.0/MESH.VolMP[i][j])*(PwallsU[i][j]*MESH.SupMP[i][j][0]*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j])) + PwallsU[i+1][j]*MESH.SupMP[i][j][1]*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(MESH.AngleMU[i+1][j])) + PwallsR[i][j]*MESH.SupMP[i][j][2]*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j])) + PwallsR[i][j+1]*MESH.SupMP[i][j][3]*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1]) + (Vpres[i][j] + Vpres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1])));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		EnergyPressTerm[0][j] = (1.0/MESH.VolMP[0][j])*(PwallsU[0][j]*MESH.SupMP[0][j][0]*(0.50*(Upres[0][j] + Uleft[j])*cos(PI + MESH.AngleMU[0][j]) + 0.50*(Vpres[0][j] + Vleft[j])*sin(PI + MESH.AngleMU[0][j])) + PwallsU[1][j]*MESH.SupMP[0][j][1]*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(MESH.AngleMU[1][j])) + PwallsR[0][j]*MESH.SupMP[0][j][2]*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j])) + PwallsR[0][j+1]*MESH.SupMP[0][j][3]*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1]) + (Vpres[0][j] + Vpres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1])));

		//Parte derecha
		EnergyPressTerm[NA-1][j] = (1.0/MESH.VolMP[NA-1][j])*(PwallsU[NA-1][j]*MESH.SupMP[NA-1][j][0]*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j])) + PwallsU[NA][j]*MESH.SupMP[NA-1][j][1]*(0.50*(Upres[NA-1][j] + Uright[j])*cos(MESH.AngleMU[NA][j]) + 0.50*(Vpres[NA-1][j] + Vright[j])*sin(MESH.AngleMU[NA][j])) + PwallsR[NA-1][j]*MESH.SupMP[NA-1][j][2]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])) + PwallsR[NA-1][j+1]*MESH.SupMP[NA-1][j][3]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + (Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])));
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		EnergyPressTerm[i][0] = (1.0/MESH.VolMP[i][0])*(PwallsU[i][0]*MESH.SupMP[i][0][0]*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0])) + PwallsU[i+1][0]*MESH.SupMP[i][0][1]*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(MESH.AngleMU[i+1][0])) + PwallsR[i][0]*MESH.SupMP[i][0][2]*(0.50*(Upres[i][0] + Udown[i])*cos(1.50*PI + MESH.AngleMR[i][0]) + 0.50*(Vpres[i][0] + Vdown[i])*sin(1.50*PI + MESH.AngleMR[i][0])) + PwallsR[i][1]*MESH.SupMP[i][0][3]*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1]) + (Vpres[i][0] + Vpres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1])));

		//Parte arriba
		EnergyPressTerm[i][NRad-1] = (1.0/MESH.VolMP[i][NRad-1])*(PwallsU[i][NRad-1]*MESH.SupMP[i][NRad-1][0]*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1])) + PwallsU[i+1][NRad-1]*MESH.SupMP[i][NRad-1][1]*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1])) + PwallsR[i][NRad-1]*MESH.SupMP[i][NRad-1][2]*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])) + PwallsR[i][NRad]*MESH.SupMP[i][NRad-1][3]*(0.50*(Upres[i][NRad-1] + Uup[i])*cos(0.50*PI + MESH.AngleMR[i][NRad]) + (Vpres[i][NRad-1] + Vup[i])*sin(0.50*PI + MESH.AngleMR[i][NRad])));

	}

	//Esquina abajo izquierda
	EnergyPressTerm[0][0] = (1.0/MESH.VolMP[0][0])*(PwallsU[0][0]*MESH.SupMP[0][0][0]*(0.50*(Upres[0][0] + Uleft[0])*cos(PI + MESH.AngleMU[0][0]) + 0.50*(Vpres[0][0] + Vleft[0])*sin(PI + MESH.AngleMU[0][0])) + PwallsU[1][0]*MESH.SupMP[0][0][1]*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(MESH.AngleMU[1][0])) + PwallsR[0][0]*MESH.SupMP[0][0][2]*(0.50*(Upres[0][0] + Udown[0])*cos(1.50*PI + MESH.AngleMR[0][0]) + 0.50*(Vpres[0][0] + Vdown[0])*sin(1.50*PI + MESH.AngleMR[0][0])) + PwallsR[0][1]*MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + (Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1])));

	//Esquina arriba izquierda
	EnergyPressTerm[0][NRad-1] = (1.0/MESH.VolMP[0][NRad-1])*(PwallsU[0][NRad-1]*MESH.SupMP[0][NRad-1][0]*(0.50*(Upres[0][NRad-1] + Uleft[NRad-1])*cos(PI + MESH.AngleMU[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vleft[NRad-1])*sin(PI + MESH.AngleMU[0][NRad-1])) + PwallsU[1][NRad-1]*MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1])) + PwallsR[0][NRad-1]*MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])) + PwallsR[0][NRad]*MESH.SupMP[0][NRad-1][3]*(0.50*(Upres[0][NRad-1] + Uup[0])*cos(0.50*PI + MESH.AngleMR[0][NRad]) + (Vpres[0][NRad-1] + Vup[0])*sin(0.50*PI + MESH.AngleMR[0][NRad])));

	//Esquina abajo derecha
	EnergyPressTerm[NA-1][0] = (1.0/MESH.VolMP[NA-1][0])*(PwallsU[NA-1][0]*MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0])) + PwallsU[NA][0]*MESH.SupMP[NA-1][0][1]*(0.50*(Upres[NA-1][0] + Uright[0])*cos(MESH.AngleMU[NA][0]) + 0.50*(Vpres[NA-1][0] + Vright[0])*sin(MESH.AngleMU[NA][0])) + PwallsR[NA-1][0]*MESH.SupMP[NA-1][0][2]*(0.50*(Upres[NA-1][0] + Udown[NA-1])*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vdown[NA-1])*sin(1.50*PI + MESH.AngleMR[NA-1][0])) + PwallsR[NA-1][1]*MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + (Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1])));

	//Esquina arriba derecha
	EnergyPressTerm[NA-1][NRad-1] = (1.0/MESH.VolMP[NA-1][NRad-1])*(PwallsU[NA-1][NRad-1]*MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1])) + PwallsU[NA][NRad-1]*MESH.SupMP[NA-1][NRad-1][1]*(0.50*(Upres[NA-1][NRad-1] + Uright[NRad-1])*cos(MESH.AngleMU[NA][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vright[NRad-1])*sin(MESH.AngleMU[i+1][NRad-1])) + PwallsR[NA-1][NRad-1]*MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])) + PwallsR[NA-1][NRad]*MESH.SupMP[NA-1][NRad-1][3]*(0.50*(Upres[NA-1][NRad-1] + Uup[NA-1])*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + (Vpres[NA-1][NRad-1] + Vup[NA-1])*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])));

}



//Cálculo del término convectivo de la ecuación de energía
void Solver::Get_EnergyConvective(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			EnergyConvective[i][j] = -(1.0/MESH.VolMP[i][j])*(MESH.SupMP[i][j][0]*0.50*(RhoPres[i][j] + RhoPres[i-1][j])*0.50*(Hpres[i][j] + Hpres[i-1][j])*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j]))	+ MESH.SupMP[i][j][1]*0.50*(RhoPres[i][j] + RhoPres[i+1][j])*0.50*(Hpres[i][j] + Hpres[i+1][j])*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(MESH.AngleMU[i+1][j])) + MESH.SupMP[i][j][2]*0.50*(RhoPres[i][j] + RhoPres[i][j-1])*0.50*(Hpres[i][j] + Hpres[i][j-1])*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j])) + MESH.SupMP[i][j][3]*0.50*(RhoPres[i][j] + RhoPres[i][j+1])*0.50*(Hpres[i][j] + Hpres[i][j+1])*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1])));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		EnergyConvective[0][j] = -(1.0/MESH.VolMP[0][j])*(MESH.SupMP[0][j][0]*RhoLeft[j]*Hleft[j]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j]) + Vleft[j]*sin(PI + MESH.AngleMU[0][j])) + MESH.SupMP[0][j][1]*0.50*(RhoPres[0][j] + RhoPres[1][j])*0.50*(Hpres[0][j] + Hpres[1][j])*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(MESH.AngleMU[1][j])) + MESH.SupMP[0][j][2]*0.50*(RhoPres[0][j] + RhoPres[0][j-1])*0.50*(Hpres[0][j] + Hpres[0][j-1])*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j])) + MESH.SupMP[0][j][3]*0.50*(RhoPres[0][j] + RhoPres[0][j+1])*0.50*(Hpres[0][j] + Hpres[0][j+1])*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1])));

		//Parte derecha
		EnergyConvective[NA-1][j] = -(1.0/MESH.VolMP[NA-1][j])*(MESH.SupMP[NA-1][j][0]*0.50*(RhoPres[NA-1][j] + RhoPres[NA-2][j])*0.50*(Hpres[NA-1][j] + Hpres[NA-2][j])*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j])) + MESH.SupMP[NA-1][j][1]*RhoRight[j]*Hright[j]*(Uright[j]*cos(MESH.AngleMU[NA][j]) + Vright[j]*sin(MESH.AngleMU[NA][j])) + MESH.SupMP[NA-1][j][2]*0.50*(RhoPres[NA-1][j] + RhoPres[NA-1][j-1])*0.50*(Hpres[NA-1][j] + Hpres[NA-1][j-1])*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])) + MESH.SupMP[NA-1][j][3]*0.50*(RhoPres[NA-1][j] + RhoPres[NA-1][j+1])*0.50*(Hpres[NA-1][j] + Hpres[NA-1][j+1])*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])));
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		EnergyConvective[i][0] = -(1.0/MESH.VolMP[i][0])*(MESH.SupMP[i][0][0]*0.50*(RhoPres[i][0] + RhoPres[i-1][0])*0.50*(Hpres[i][0] + Hpres[i-1][0])*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0])) + MESH.SupMP[i][0][1]*0.50*(RhoPres[i][0] + RhoPres[i+1][0])*0.50*(Hpres[i][0] + Hpres[i+1][0])*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(MESH.AngleMU[i+1][0])) + MESH.SupMP[i][0][2]*RhoDown[i]*Hdown[i]*(Udown[i]*cos(1.50*PI + MESH.AngleMR[i][0]) + Vdown[i]*sin(1.50*PI + MESH.AngleMR[i][0])) + MESH.SupMP[i][0][3]*0.50*(RhoPres[i][0] + RhoPres[i][1])*0.50*(Hpres[i][0] + Hpres[i][1])*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1])));

		//Parte arriba
		EnergyConvective[i][NRad-1] = -(1.0/MESH.VolMP[i][NRad-1])*(MESH.SupMP[i][NRad-1][0]*0.50*(RhoPres[i][NRad-1] + RhoPres[i-1][NRad-1])*0.50*(Hpres[i][NRad-1] + Hpres[i-1][NRad-1])*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1]))	+ MESH.SupMP[i][NRad-1][1]*0.50*(RhoPres[i][NRad-1] + RhoPres[i+1][NRad-1])*0.50*(Hpres[i][NRad-1] + Hpres[i+1][NRad-1])*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1])) + MESH.SupMP[i][NRad-1][2]*0.50*(RhoPres[i][NRad-1] + RhoPres[i][NRad-2])*0.50*(Hpres[i][NRad-1] + Hpres[i][NRad-2])*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])) + MESH.SupMP[i][NRad-1][3]*RhoUp[i]*Hup[i]*(Uup[i]*cos(0.50*PI + MESH.AngleMR[i][NRad]) + Vup[i]*sin(0.50*PI + MESH.AngleMR[i][NRad])));
	}

	//Esquina abajo izquierda
	EnergyConvective[0][0] = -(1.0/MESH.VolMP[0][0])*(MESH.SupMP[0][0][0]*RhoLeft[0]*Hleft[0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0]) + Vleft[0]*sin(PI + MESH.AngleMU[0][0]))	+ MESH.SupMP[0][0][1]*0.50*(RhoPres[0][0] + RhoPres[1][0])*0.50*(Hpres[0][0] + Hpres[1][0])*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(MESH.AngleMU[1][0])) + MESH.SupMP[0][0][2]*RhoDown[0]*Hdown[0]*(Udown[0]*cos(1.50*PI + MESH.AngleMR[0][0]) + Vdown[0]*sin(1.50*PI + MESH.AngleMR[0][0])) + MESH.SupMP[0][0][3]*0.50*(RhoPres[0][0] + RhoPres[0][1])*0.50*(Hpres[0][0] + Hpres[0][1])*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1])));
	
	//Esquina arriba izquierda
	EnergyConvective[0][NRad-1] = -(1.0/MESH.VolMP[0][NRad-1])*(MESH.SupMP[0][NRad-1][0]*RhoLeft[NRad-1]*Hleft[NRad-1]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1])) + MESH.SupMP[0][NRad-1][1]*0.50*(RhoPres[0][NRad-1] + RhoPres[1][NRad-1])*0.50*(Hpres[0][NRad-1] + Hpres[1][NRad-1])*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1])) + MESH.SupMP[0][NRad-1][2]*0.50*(RhoPres[0][NRad-1] + RhoPres[0][NRad-2])*0.50*(Hpres[0][NRad-1] + Hpres[0][NRad-2])*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])) + MESH.SupMP[0][NRad-1][3]*RhoUp[0]*Hup[0]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad])));
	
	//Esquina abajo derecha
	EnergyConvective[NA-1][0] = -(1.0/MESH.VolMP[NA-1][0])*(MESH.SupMP[NA-1][0][0]*0.50*(RhoPres[NA-1][0] + RhoPres[NA-2][0])*0.50*(Hpres[NA-1][0] + Hpres[NA-2][0])*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0]))	+ MESH.SupMP[NA-1][0][1]*RhoRight[0]*Hright[0]*(Uright[0]*cos(MESH.AngleMU[NA][0]) + Vright[0]*sin(MESH.AngleMU[NA][0])) + MESH.SupMP[NA-1][0][2]*RhoDown[NA-1]*Hdown[NA-1]*(Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0]) + Vdown[NA-1]*sin(1.50*PI + MESH.AngleMR[NA-1][0])) + MESH.SupMP[NA-1][0][3]*0.50*(RhoPres[NA-1][0] + RhoPres[NA-1][1])*0.50*(Hpres[NA-1][0] + Hpres[NA-1][1])*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1])));

	//Esquina arriba derecha
	EnergyConvective[NA-1][NRad-1] = -(1.0/MESH.VolMP[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]*0.50*(RhoPres[NA-1][NRad-1] + RhoPres[NA-2][NRad-1])*0.50*(Hpres[NA-1][NRad-1] + Hpres[NA-2][NRad-1])*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1]))	+ MESH.SupMP[NA-1][NRad-1][1]*RhoRight[NRad-1]*Hright[NRad-1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1])) + MESH.SupMP[NA-1][NRad-1][2]*0.50*(RhoPres[NA-1][NRad-1] + RhoPres[NA-1][NRad-2])*0.50*(Hpres[NA-1][NRad-1] + Hpres[NA-1][NRad-2])*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])) + MESH.SupMP[NA-1][NRad-1][3]*RhoUp[NA-1]*Hup[NA-1]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])));

}



//Cálculo de la divergencia de la velocidad en cada nodo (Término viscoso eq energía)
void Solver::Get_Divergence(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			Divergence[i][j] = (1.0/MESH.VolMP[i][j])*(MESH.SupMP[i][j][0]*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*sin(PI + MESH.AngleMU[i][j])) + MESH.SupMP[i][j][1]*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(PI + MESH.AngleMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*sin(PI + MESH.AngleMU[i+1][j])) + MESH.SupMP[i][j][2]*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(PI + MESH.AngleMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*sin(PI + MESH.AngleMR[i][j])) + MESH.SupMP[i][j][3]*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(PI + MESH.AngleMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*sin(PI + MESH.AngleMR[i][j+1])));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		Divergence[0][j] = (1.0/MESH.VolMP[0][j])*(MESH.SupMP[0][j][0]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j]) + Vleft[j]*sin(PI + MESH.AngleMU[0][j])) + MESH.SupMP[0][j][1]*(0.50*(Upres[0][j] + Upres[1][j])*cos(PI + MESH.AngleMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*sin(PI + MESH.AngleMU[1][j])) + MESH.SupMP[0][j][2]*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(PI + MESH.AngleMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*sin(PI + MESH.AngleMR[0][j])) + MESH.SupMP[0][j][3]*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(PI + MESH.AngleMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*sin(PI + MESH.AngleMR[0][j+1])));

		//Parte derecha
		Divergence[NA-1][j] = (1.0/MESH.VolMP[NA-1][j])*(MESH.SupMP[NA-1][j][0]*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j])) + MESH.SupMP[NA-1][j][1]*(Uright[j]*cos(PI + MESH.AngleMU[NA][j]) + Vright[i]*sin(PI + MESH.AngleMU[NA][j])) + MESH.SupMP[NA-1][j][2]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(PI + MESH.AngleMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*sin(PI + MESH.AngleMR[NA-1][j])) + MESH.SupMP[NA-1][j][3]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(PI + MESH.AngleMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*sin(PI + MESH.AngleMR[NA-1][j+1])));

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		Divergence[i][0] = (1.0/MESH.VolMP[i][0])*(MESH.SupMP[i][0][0]*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*sin(PI + MESH.AngleMU[i][0])) + MESH.SupMP[i][0][1]*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(PI + MESH.AngleMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*sin(PI + MESH.AngleMU[i+1][0])) + MESH.SupMP[i][0][2]*(Udown[i]*cos(PI + MESH.AngleMR[i][0]) + Vdown[i]*sin(PI + MESH.AngleMR[i][0])) + MESH.SupMP[i][0][3]*(0.50*(Upres[i][0] + Upres[i][1])*cos(PI + MESH.AngleMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*sin(PI + MESH.AngleMR[i][1])));

		//Parte arriba
		Divergence[i][NRad-1] = (1.0/MESH.VolMP[i][NRad-1])*(MESH.SupMP[i][NRad-1][0]*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1])) + MESH.SupMP[i][NRad-1][1]*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(PI + MESH.AngleMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*sin(PI + MESH.AngleMU[i+1][NRad-1])) + MESH.SupMP[i][NRad-1][2]*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(PI + MESH.AngleMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*sin(PI + MESH.AngleMR[i][NRad-1])) + MESH.SupMP[i][NRad-1][3]*(Uup[i]*cos(PI + MESH.AngleMR[i][NRad]) + Vup[i]*sin(PI + MESH.AngleMR[i][NRad])));

	}

	//Esquina abajo izquierda
	Divergence[0][0] = (1.0/MESH.VolMP[0][0])*(MESH.SupMP[0][0][0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0]) + Vleft[0]*sin(PI + MESH.AngleMU[0][0])) + MESH.SupMP[0][0][1]*(0.50*(Upres[0][0] + Upres[1][0])*cos(PI + MESH.AngleMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*sin(PI + MESH.AngleMU[1][0])) + MESH.SupMP[0][0][2]*(Udown[0]*cos(PI + MESH.AngleMR[0][j]) + Vdown[0]*sin(PI + MESH.AngleMR[0][0])) + MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(PI + MESH.AngleMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*sin(PI + MESH.AngleMR[0][1])));

	//Esquina arriba izquierda
	Divergence[0][NRad-1] = (1.0/MESH.VolMP[0][NRad-1])*(MESH.SupMP[0][NRad-1][0]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1]) + Vleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1])) + MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(PI + MESH.AngleMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*sin(PI + MESH.AngleMU[1][NRad-1])) + MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(PI + MESH.AngleMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*sin(PI + MESH.AngleMR[0][NRad-1])) + MESH.SupMP[0][NRad-1][3]*(Uup[0]*cos(PI + MESH.AngleMR[0][NRad]) + Vup[0]*sin(PI + MESH.AngleMR[0][NRad])));

	//Esquina abajo derecha
	Divergence[NA-1][0] = (1.0/MESH.VolMP[NA-1][0])*(MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0])) + MESH.SupMP[NA-1][0][1]*(Uright[0]*cos(PI + MESH.AngleMU[i+1][0]) + Vright[0]*sin(PI + MESH.AngleMU[NA][0])) + MESH.SupMP[NA-1][0][2]*(Udown[NA-1]*cos(PI + MESH.AngleMR[NA-1][0]) + Vdown[NA-1]*sin(PI + MESH.AngleMR[NA-1][0])) + MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(PI + MESH.AngleMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*sin(PI + MESH.AngleMR[NA-1][1])));

	//Esquina arriba derecha
	Divergence[NA-1][NRad-1] = (1.0/MESH.VolMP[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1])) + MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(PI + MESH.AngleMU[NA][NRad-1]) + Vright[NRad-1]*sin(PI + MESH.AngleMU[NA][NRad-1])) + MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(PI + MESH.AngleMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*sin(PI + MESH.AngleMR[NA-1][NRad-1])) + MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(PI + MESH.AngleMR[NA-1][NRad]) + Vup[NA-1]*sin(PI + MESH.AngleMR[NA-1][NRad])));

}


//Cálculo de los gradientes de velocidades en cada una de las direcciones
void Solver::Get_SpeedGradients(Mesher MESH){
int i, j;

	//Matriz MU
	
	//Centro
	for(i = 1; i < NA; i++){
		for(j = 1; j < NRad-1; j++){
			//Gradiente de U respecto X
			GradU_DxMU[i][j] = (Upres[i][j] - Upres[i-1][j])/MESH.DeltasMU[i][j][0];

			//Gradiente de U respecto Y
			GradU_DyMU[i][j] = (0.50*(Upres[i][j+1] + Upres[i-1][j+1]) - 0.50*(Upres[i][j-1] + Upres[i-1][j-1]))/(0.50*(MESH.DeltasMU[i][j-1][1] + MESH.DeltasMU[i][j+1][1]) + MESH.DeltasMU[i][j][1]);

			//Gradiente de V respecto X
			GradV_DxMU[i][j] = (0.50*(Vpres[i][j+1] + Vpres[i-1][j+1]) - 0.50*(Vpres[i][j-1] + Vpres[i-1][j-1]))/(0.50*(MESH.DeltasMU[i][j-1][1] + MESH.DeltasMU[i][j+1][1]) + MESH.DeltasMU[i][j][1]);

			//Gradiente de V respecto Y
			GradV_DxMU[i][j] = (Vpres[i][j] - Vpres[i-1][j])/MESH.DeltasMU[i][j][0];
		}
	}

	
	for(j = 1; j < NRad-1; j++){
		//Parte izquierda

			//Gradiente de U respecto X
			GradU_DxMU[0][j] = (Upres[0][j] - Uleft[j])/MESH.DeltasMU[0][j][0];

			//Gradiente de U respecto Y
			GradU_DyMU[0][j] = (0.50*(Upres[0][j+1] + Uleft[j+1]) - 0.50*(Upres[0][j-1] + Uleft[j-1]))/(0.50*(MESH.DeltasMU[0][j-1][1] + MESH.DeltasMU[0][j+1][1]) + MESH.DeltasMU[0][j][1]);

			//Gradiente de V respecto X
			GradV_DxMU[0][j] = (0.50*(Vpres[0][j+1] + Vleft[j+1]) - 0.50*(Vpres[0][j-1] + Vleft[j-1]))/(0.50*(MESH.DeltasMU[0][j-1][1] + MESH.DeltasMU[0][j+1][1]) + MESH.DeltasMU[0][j][1]);

			//Gradiente de V respecto Y
			GradV_DxMU[0][j] = (Vpres[0][j] - Vleft[j])/MESH.DeltasMU[0][j][0];

		//Parte derecha

			//Gradiente de U respecto X
			GradU_DxMU[NA][j] = (Uright[j] - Upres[NA-1][j])/MESH.DeltasMU[NA][j][0];

			//Gradiente de U respecto Y
			GradU_DyMU[NA][j] = (0.50*(Uright[j+1] + Upres[NA-1][j+1]) - 0.50*(Upres[NA-1][j-1] + Upres[NA-1][j-1]))/(0.50*(MESH.DeltasMU[NA][j-1][1] + MESH.DeltasMU[NA][j+1][1]) + MESH.DeltasMU[NA][j][1]);

			//Gradiente de V respecto X
			GradV_DxMU[NA][j] = (0.50*(Vright[j+1] + Vpres[NA-1][j+1]) - 0.50*(Vright[j-1] + Vpres[NA-1][j-1]))/(0.50*(MESH.DeltasMU[NA][j-1][1] + MESH.DeltasMU[NA][j+1][1]) + MESH.DeltasMU[NA][j][1]);

			//Gradiente de V respecto Y
			GradV_DxMU[NA][j] = (Vright[j] - Vpres[NA-1][j])/MESH.DeltasMU[NA][j][0];

	}
	
	for(i = 1; i < NA-1; i++){
		//Parte abajo

			//Gradiente de U respecto X
			GradU_DxMU[i][0] = (Upres[i][0] - Upres[i-1][0])/MESH.DeltasMU[i][0][0];

			//Gradiente de U respecto Y
			GradU_DyMU[i][0] = (0.50*(Upres[i][1] + Upres[i-1][1]) - 0.50*(Udown[i] + Udown[i-1]))/(0.50*(MESH.DeltasMU[i][1][1]) + MESH.DeltasMU[i][0][1]);

			//Gradiente de V respecto X
			GradV_DxMU[i][0] = (0.50*(Vpres[i][1] + Vpres[i-1][1]) - 0.50*(Vdown[i] + Vdown[i-1]))/(0.50*(MESH.DeltasMU[i][1][1]) + MESH.DeltasMU[i][0][1]);

			//Gradiente de V respecto Y
			GradV_DxMU[i][0] = (Vpres[i][0] - Vpres[i-1][0])/MESH.DeltasMU[i][0][0];

		//Parte arriba

			//Gradiente de U respecto X
			GradU_DxMU[i][NRad-1] = (Upres[i][NRad-1] - Upres[i-1][NRad-1])/MESH.DeltasMU[i][NRad-1][0];

			//Gradiente de U respecto Y
			GradU_DyMU[i][NRad-1] = (0.50*(Uup[i] + Uup[i-1]) - 0.50*(Upres[i][NRad-2] + Upres[i-1][NRad-2]))/(0.50*(MESH.DeltasMU[i][NRad-2][1]) + MESH.DeltasMU[i][NRad-1][1]);

			//Gradiente de V respecto X
			GradV_DxMU[i][NRad-1] = (0.50*(Vup[i] + Vup[i-1]) - 0.50*(Vpres[i][NRad-2] + Vpres[i-1][NRad-2]))/(0.50*(MESH.DeltasMU[i][NRad-2][1]) + MESH.DeltasMU[i][NRad-1][1]);

			//Gradiente de V respecto Y
			GradV_DxMU[i][NRad-1] = (Vpres[i][NRad-1] - Vpres[i-1][NRad-1])/MESH.DeltasMU[i][NRad-1][0];

	}

	//Esquina abajo izquierda

		//Gradiente de U respecto X
		GradU_DxMU[0][0] = (Upres[0][0] - Uleft[0])/MESH.DeltasMU[0][0][0];

		//Gradiente de U respecto Y
		GradU_DyMU[0][0] = (0.50*(Upres[0][1] + Uleft[1]) - 0.50*(Udown[0]))/(0.50*(MESH.DeltasMU[0][1][1]) + MESH.DeltasMU[0][0][1]);

		//Gradiente de V respecto X
		GradV_DxMU[0][0] = (0.50*(Vpres[0][1] + Vleft[1]) - 0.50*(Vdown[0]))/(0.50*(MESH.DeltasMU[0][1][1]) + MESH.DeltasMU[0][0][1]);

		//Gradiente de V respecto Y
		GradV_DxMU[0][0] = (Vpres[0][0] - Vleft[0])/MESH.DeltasMU[0][0][0];


	//Esquina arriba izquierda

		//Gradiente de U respecto X
		GradU_DxMU[0][NRad-1] = (Upres[0][NRad-1] - Uleft[NRad-1])/MESH.DeltasMU[0][NRad-1][0];

		//Gradiente de U respecto Y
		GradU_DyMU[0][NRad-1] = (0.50*(Uup[0]) - 0.50*(Upres[0][NRad-2] + Uleft[NRad-2]))/(0.50*(MESH.DeltasMU[0][NRad-2][1]) + MESH.DeltasMU[0][NRad-1][1]);

		//Gradiente de V respecto X
		GradV_DxMU[0][NRad-1] = (0.50*(Vup[0]) - 0.50*(Vpres[0][NRad-2] + Vleft[NRad-2]))/(0.50*(MESH.DeltasMU[0][NRad-2][1] ) + MESH.DeltasMU[0][NRad-1][1]);

		//Gradiente de V respecto Y
		GradV_DxMU[0][NRad-1] = (Vpres[0][NRad-1] - Vleft[NRad-1])/MESH.DeltasMU[0][NRad-1][0];

	//Esquina abajo derecha

		//Gradiente de U respecto X
		GradU_DxMU[NA][0] = (Uright[0] - Upres[NA-1][0])/MESH.DeltasMU[NA][0][0];

		//Gradiente de U respecto Y
		GradU_DyMU[NA][0] = (0.50*(Uright[1] + Upres[NA-1][1]) - 0.50*(Udown[NA-1]))/(0.50*(MESH.DeltasMU[NA][1][1]) + MESH.DeltasMU[NA][0][1]);

		//Gradiente de V respecto X
		GradV_DxMU[NA][0] = (0.50*(Vright[1] + Vpres[NA-1][1]) - 0.50*(Vdown[NA-1]))/(0.50*(MESH.DeltasMU[NA][1][1]) + MESH.DeltasMU[NA][0][1]);

		//Gradiente de V respecto Y
		GradV_DxMU[NA][0] = (Vright[0] - Vpres[NA-1][0])/MESH.DeltasMU[NA][0][0];

	//Esquina arriba derecha
		//Gradiente de U respecto X
		GradU_DxMU[NA][NRad-1] = (Uright[NRad-1] - Upres[NA-1][NRad-1])/MESH.DeltasMU[NA][NRad-1][0];

		//Gradiente de U respecto Y
		GradU_DyMU[NA][NRad-1] = (0.50*(Uup[NA-1]) - 0.50*(Uright[NRad-2] + Upres[NA-1][NRad-2]))/(0.50*(MESH.DeltasMU[NA][NRad-2][1] + MESH.DeltasMU[NA][NRad-1][1]) + MESH.DeltasMU[NA][NRad-1][1]);


		//Gradiente de V respecto X
		GradV_DxMU[NA][NRad-1] = (0.50*(Vup[NA-1]) - 0.50*(Vright[NRad-2] + Vpres[NA-1][NRad-2]))/(0.50*(MESH.DeltasMU[NA][NRad-2][1] + MESH.DeltasMU[NA][NRad-1][1]) + MESH.DeltasMU[NA-1][NRad-1][1]);

		//Gradiente de V respecto Y
		GradV_DxMU[NA][NRad-1] = (Vright[NRad-1] - Vpres[NA-1][NRad-1])/MESH.DeltasMU[NA][NRad-1][0];

//------------------------------------------------------------------------------------------------------------------------------------

	//Matriz MR

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad; j++){
			//Gradiente de U respecto X
			GradU_DxMR[i][j] = (0.50*(Upres[i+1][j] + Upres[i+1][j-1]) - 0.50*(Upres[i-1][j] + Upres[i-1][j-1]))/(0.50*(MESH.DeltasMR[i-1][j][0] + MESH.DeltasMR[i+1][j][0]) + MESH.DeltasMR[i][j][0]);

			//Gradiente de U respecto Y
			GradU_DyMR[i][j] = (Upres[i][j] - Upres[i][j-1])/MESH.DeltasMR[i][j][1];

			//Gradiente de V respecto X
			GradV_DxMR[i][j] = (0.50*(Vpres[i+1][j] + Vpres[i+1][j-1]) - 0.50*(Vpres[i-1][j] + Vpres[i-1][j-1]))/(0.50*(MESH.DeltasMR[i-1][j][0] + MESH.DeltasMR[i+1][j][0]) + MESH.DeltasMR[i][j][0]); 

			//Gradiente de V respecto Y
			GradV_DxMR[i][j] = (Vpres[i][j] - Vpres[i][j-1])/MESH.DeltasMR[i][j][1];
		}
	}

	for(j = 1; j < NRad; j++){
		//Parte izquierda

			//Gradiente de U respecto X
			GradU_DxMR[0][j] = (0.50*(Upres[1][j] + Upres[1][j-1]) - 0.50*(Uleft[j] + Uleft[j-1]))/(0.50*(MESH.DeltasMR[1][j][0]) + MESH.DeltasMR[0][j][0]);

			//Gradiente de U respecto Y
			GradU_DyMR[0][j] = (Upres[0][j] - Upres[0][j-1])/MESH.DeltasMR[0][j][1];

			//Gradiente de V respecto X
			GradV_DxMR[0][j] = (0.50*(Vpres[1][j] + Vpres[1][j-1]) - 0.50*(Vleft[j] + Vleft[j-1]))/(0.50*(MESH.DeltasMR[1][j][0]) + MESH.DeltasMR[0][j][0]); 

			//Gradiente de V respecto Y
			GradV_DxMR[0][j] = (Vpres[0][j] - Vpres[0][j-1])/MESH.DeltasMR[0][j][1];

		//Parte derecha

			//Gradiente de U respecto X
			GradU_DxMR[NA-1][j] = (0.50*(Uright[j] + Uright[j-1]) - 0.50*(Upres[NA-2][j] + Upres[NA-2][j-1]))/(0.50*(MESH.DeltasMR[NA-2][j][0]) + MESH.DeltasMR[NA-1][j][0]);

			//Gradiente de U respecto Y
			GradU_DyMR[NA-1][j] = (Upres[NA-1][j] - Upres[NA-1][j-1])/MESH.DeltasMR[NA-1][j][1];

			//Gradiente de V respecto X
			GradV_DxMR[NA-1][j] = (0.50*(Vright[j] + Vright[j-1]) - 0.50*(Vpres[NA-2][j] + Vpres[NA-2][j-1]))/(0.50*(MESH.DeltasMR[NA-2][j][0]) + MESH.DeltasMR[NA-1][j][0]); 

			//Gradiente de V respecto Y
			GradV_DxMR[NA-1][j] = (Vpres[NA-1][j] - Vpres[NA-1][j-1])/MESH.DeltasMR[NA-1][j][1];

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo

			//Gradiente de U respecto X
			GradU_DxMR[i][0] = (0.50*(Upres[i+1][0] + Udown[i+1]) - 0.50*(Upres[i-1][0] + Udown[i-1]))/(0.50*(MESH.DeltasMR[i-1][0][0] + MESH.DeltasMR[i+1][0][0]) + MESH.DeltasMR[i][0][0]);

			//Gradiente de U respecto Y
			GradU_DyMR[i][0] = (Upres[i][0] - Udown[i])/MESH.DeltasMR[i][0][1];

			//Gradiente de V respecto X
			GradV_DxMR[i][0] = (0.50*(Vpres[i+1][0] + Vdown[i+1]) - 0.50*(Vpres[i-1][0] + Vdown[i-1]))/(0.50*(MESH.DeltasMR[i-1][0][0] + MESH.DeltasMR[i+1][0][0]) + MESH.DeltasMR[i][0][0]); 

			//Gradiente de V respecto Y
			GradV_DxMR[i][0] = (Vpres[i][0] - Vdown[i])/MESH.DeltasMR[i][0][1];

		//Parte arriba
		
			//Gradiente de U respecto X
			GradU_DxMR[i][NRad] = (0.50*(Uup[i+1] + Upres[i+1][NRad-1]) - 0.50*(Uup[i-1] + Upres[i-1][NRad-1]))/(0.50*(MESH.DeltasMR[i-1][NRad][0] + MESH.DeltasMR[i+1][NRad][0]) + MESH.DeltasMR[i][NRad][0]);

			//Gradiente de U respecto Y
			GradU_DyMR[i][NRad] = (Uup[i] - Upres[i][NRad-1])/MESH.DeltasMR[i][NRad][1];

			//Gradiente de V respecto X
			GradV_DxMR[i][NRad] = (0.50*(Vup[i+1] + Vpres[i+1][NRad-1]) - 0.50*(Vup[i-1] + Vpres[i-1][NRad-1]))/(0.50*(MESH.DeltasMR[i-1][NRad][0] + MESH.DeltasMR[i+1][NRad][0]) + MESH.DeltasMR[i][NRad][0]); 

			//Gradiente de V respecto Y
			GradV_DxMR[i][NRad] = (Vup[i] - Vpres[i][NRad-1])/MESH.DeltasMR[i][NRad][1];
	}	
		
	//Esquina abajo izquierda

		//Gradiente de U respecto X
		GradU_DxMR[0][0] = (0.50*(Upres[1][0] + Udown[1]) - 0.50*(Uleft[0]))/(0.50*(MESH.DeltasMR[1][0][0]) + MESH.DeltasMR[0][0][0]);

		//Gradiente de U respecto Y
		GradU_DyMR[0][0] = (Upres[0][0] - Udown[0])/MESH.DeltasMR[0][0][1];

		//Gradiente de V respecto X
		GradV_DxMR[0][0] = (0.50*(Vpres[1][0] + Vdown[1]) - 0.50*(Vleft[0]))/(0.50*(MESH.DeltasMR[1][0][0]) + MESH.DeltasMR[0][0][0]); 

		//Gradiente de V respecto Y
		GradV_DxMR[0][0] = (Vpres[0][0] - Vdown[0])/MESH.DeltasMR[0][0][1];

	//Esquina arriba izquierda

		//Gradiente de U respecto X
		GradU_DxMR[0][NRad] = (0.50*(Uup[1] + Upres[1][NRad-1]) - 0.50*(Uleft[NRad] + Uleft[NRad-1]))/(0.50*(MESH.DeltasMR[1][NRad][0]) + MESH.DeltasMR[0][NRad][0]);

		//Gradiente de U respecto Y
		GradU_DyMR[0][NRad] = (Uup[0] - Upres[0][NRad-1])/MESH.DeltasMR[0][NRad][1];

		//Gradiente de V respecto X
		GradV_DxMR[0][NRad] = (0.50*(Vup[1] + Vpres[1][NRad-1]) - 0.50*(Vleft[NRad] + Vleft[NRad-1]))/(0.50*(MESH.DeltasMR[1][NRad][0]) + MESH.DeltasMR[0][NRad][0]); 

		//Gradiente de V respecto Y
		GradV_DxMR[0][NRad] = (Vup[0] - Vpres[0][NRad-1])/MESH.DeltasMR[0][NRad][1];

	//Esquina abajo derecha

		//Gradiente de U respecto X
		GradU_DxMR[NA-1][0] = (0.50*(Uright[0]) - 0.50*(Upres[NA-2][0] + Udown[NA-2]))/(0.50*(MESH.DeltasMR[NA-2][0][0]) + MESH.DeltasMR[NA-1][0][0]);

		//Gradiente de U respecto Y
		GradU_DyMR[NA-1][0] = (Upres[NA-1][0] - Udown[NA-1])/MESH.DeltasMR[NA-1][0][1];

		//Gradiente de V respecto X
		GradV_DxMR[NA-1][0] = (0.50*(Vright[0]) - 0.50*(Vpres[NA-2][0] + Vdown[NA-2]))/(0.50*(MESH.DeltasMR[NA-2][0][0]) + MESH.DeltasMR[NA-1][0][0]); 

		//Gradiente de V respecto Y
		GradV_DxMR[NA-1][0] = (Vpres[NA-1][0] - Vdown[NA-1])/MESH.DeltasMR[NA-1][0][1];

	//Esquina arriba derecha	

		//Gradiente de U respecto X
		GradU_DxMR[NA-1][NRad] = (0.50*(Uright[NRad-1]) - 0.50*(Uup[NA-2] + Upres[NA-2][NRad-1]))/(0.50*(MESH.DeltasMR[NA-2][NRad][0]) + MESH.DeltasMR[NA-1][NRad][0]);

		//Gradiente de U respecto Y
		GradU_DyMR[NA-1][NRad] = (Uup[NA-1] - Upres[NA-1][NRad-1])/MESH.DeltasMR[NA-1][NRad][1];

		//Gradiente de V respecto X
		GradV_DxMR[NA-1][NRad] = (0.50*(Vright[NRad-1]) - 0.50*(Vup[NA-2] + Vpres[NA-2][NRad-1]))/(0.50*(MESH.DeltasMR[NA-2][NRad][0]) + MESH.DeltasMR[NA-1][NRad][0]); 

		//Gradiente de V respecto Y
		GradV_DxMR[NA-1][NRad] = (Vup[NA-1] - Vpres[NA-1][NRad-1])/MESH.DeltasMR[NA-1][NRad][1];

}




//Cálculo del término de disipación viscosa en la ecuación de energía
void Solver::Get_ViscousTerm(Mesher MESH){
int i, j;
	
	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			EnergyViscous[i][j] = (1.0/MESH.VolMP[i][j])*(MESH.SupMP[i][j][0]*(0.50*(Upres[i][j] + Upres[i-1][j])*cos(PI + MESH.AngleMU[i][j])*(2.0*muWallsMU[i][j]*GradU_DxMU[i][j] - (2.0/3.0)*muWallsMU[i][j]*0.50*(Divergence[i][j] + Divergence[i-1][j])) + 0.50*(Upres[i][j] + Upres[i-1][j])*sin(PI + MESH.AngleMU[i][j])*muWallsMU[i][j]*(GradU_DyMU[i][j] + GradV_DxMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*cos(PI + MESH.AngleMU[i][j])*muWallsMU[i][j]*(GradU_DyMU[i][j] + GradV_DxMU[i][j]) + 0.50*(Vpres[i][j] + Vpres[i-1][j])*(2.0*muWallsMU[i][j]*GradV_DyMU[i][j] - (2.0/3.0)*muWallsMU[i][j]*0.50*(Divergence[i][j] + Divergence[i-1][j]))) +

MESH.SupMP[i][j][1]*(0.50*(Upres[i][j] + Upres[i+1][j])*cos(MESH.AngleMU[i+1][j])*(2.0*muWallsMU[i+1][j]*GradU_DxMU[i+1][j] - (2.0/3.0)*muWallsMU[i+1][j]*0.50*(Divergence[i][j] + Divergence[i+1][j])) + 0.50*(Upres[i][j] + Upres[i+1][j])*sin(MESH.AngleMU[i+1][j])*muWallsMU[i+1][j]*(GradU_DyMU[i+1][j] + GradV_DxMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*cos(MESH.AngleMU[i+1][j])*muWallsMU[i+1][j]*(GradU_DyMU[i+1][j] + GradV_DxMU[i+1][j]) + 0.50*(Vpres[i][j] + Vpres[i+1][j])*(2.0*muWallsMU[i+1][j]*GradV_DyMU[i+1][j] - (2.0/3.0)*muWallsMU[i+1][j]*0.50*(Divergence[i][j] + Divergence[i+1][j]))) + 

MESH.SupMP[i][j][2]*(0.50*(Upres[i][j] + Upres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j])*(2.0*muWallsMR[i][j]*GradU_DxMR[i][j] - (2.0/3.0)*muWallsMR[i][j]*0.50*(Divergence[i][j] + Divergence[i][j-1])) + 0.50*(Upres[i][j] + Upres[i][j-1])*sin(1.50*PI + MESH.AngleMR[i][j])*muWallsMR[i][j]*(GradU_DyMR[i][j] + GradV_DxMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*cos(1.50*PI + MESH.AngleMR[i][j])*muWallsMR[i][j]*(GradU_DyMR[i][j] + GradV_DxMR[i][j]) + 0.50*(Vpres[i][j] + Vpres[i][j-1])*(2.0*muWallsMR[i][j]*GradV_DyMR[i][j] - (2.0/3.0)*muWallsMR[i][j]*0.50*(Divergence[i][j] + Divergence[i][j-1])))	+

MESH.SupMP[i][j][3]*(0.50*(Upres[i][j] + Upres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1])*(2.0*muWallsMR[i][j+1]*GradU_DxMR[i][j+1] - (2.0/3.0)*muWallsMR[i][j+1]*0.50*(Divergence[i][j] + Divergence[i][j+1])) + 0.50*(Upres[i][j] + Upres[i][j+1])*sin(0.50*PI + MESH.AngleMR[i][j+1])*muWallsMR[i][j+1]*(GradU_DyMR[i][j+1] + GradV_DxMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*cos(0.50*PI + MESH.AngleMR[i][j+1])*muWallsMR[i][j+1]*(GradU_DyMR[i][j+1] + GradV_DxMR[i][j+1]) + 0.50*(Vpres[i][j] + Vpres[i][j+1])*(2.0*muWallsMR[i][j+1]*GradV_DyMR[i][j+1] - (2.0/3.0)*muWallsMR[i][j+1]*0.50*(Divergence[i][j] + Divergence[i][j+1]))));
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		EnergyViscous[0][j] = (1.0/MESH.VolMP[0][j])*(MESH.SupMP[0][j][0]*(Uleft[j]*cos(PI + MESH.AngleMU[0][j])*(2.0*muWallsMU[0][j]*GradU_DxMU[0][j] - (2.0/3.0)*muWallsMU[0][j]*0.50*(Divergence[0][j])) + Uleft[j]*sin(PI + MESH.AngleMU[0][j])*muWallsMU[0][j]*(GradU_DyMU[0][j] + GradV_DxMU[0][j]) + Vleft[j]*cos(PI + MESH.AngleMU[0][j])*muWallsMU[0][j]*(GradU_DyMU[0][j] + GradV_DxMU[0][j]) + Vleft[j]*(2.0*muWallsMU[0][j]*GradV_DyMU[0][j] - (2.0/3.0)*muWallsMU[0][j]*0.50*(Divergence[i][j]))) +

MESH.SupMP[0][j][1]*(0.50*(Upres[0][j] + Upres[1][j])*cos(MESH.AngleMU[1][j])*(2.0*muWallsMU[1][j]*GradU_DxMU[1][j] - (2.0/3.0)*muWallsMU[1][j]*0.50*(Divergence[0][j] + Divergence[1][j])) + 0.50*(Upres[0][j] + Upres[1][j])*sin(MESH.AngleMU[1][j])*muWallsMU[1][j]*(GradU_DyMU[1][j] + GradV_DxMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*cos(MESH.AngleMU[1][j])*muWallsMU[1][j]*(GradU_DyMU[1][j] + GradV_DxMU[1][j]) + 0.50*(Vpres[0][j] + Vpres[1][j])*(2.0*muWallsMU[1][j]*GradV_DyMU[1][j] - (2.0/3.0)*muWallsMU[1][j]*0.50*(Divergence[0][j] + Divergence[1][j]))) + 

MESH.SupMP[0][j][2]*(0.50*(Upres[0][j] + Upres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j])*(2.0*muWallsMR[0][j]*GradU_DxMR[0][j] - (2.0/3.0)*muWallsMR[0][j]*0.50*(Divergence[0][j] + Divergence[0][j-1])) + 0.50*(Upres[0][j] + Upres[0][j-1])*sin(1.50*PI + MESH.AngleMR[0][j])*muWallsMR[0][j]*(GradU_DyMR[0][j] + GradV_DxMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*cos(1.50*PI + MESH.AngleMR[0][j])*muWallsMR[0][j]*(GradU_DyMR[0][j] + GradV_DxMR[0][j]) + 0.50*(Vpres[0][j] + Vpres[0][j-1])*(2.0*muWallsMR[0][j]*GradV_DyMR[0][j] - (2.0/3.0)*muWallsMR[0][j]*0.50*(Divergence[0][j] + Divergence[0][j-1])))	+

MESH.SupMP[0][j][3]*(0.50*(Upres[0][j] + Upres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1])*(2.0*muWallsMR[0][j+1]*GradU_DxMR[0][j+1] - (2.0/3.0)*muWallsMR[0][j+1]*0.50*(Divergence[0][j] + Divergence[0][j+1])) + 0.50*(Upres[0][j] + Upres[0][j+1])*sin(0.50*PI + MESH.AngleMR[0][j+1])*muWallsMR[0][j+1]*(GradU_DyMR[0][j+1] + GradV_DxMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*cos(0.50*PI + MESH.AngleMR[0][j+1])*muWallsMR[0][j+1]*(GradU_DyMR[0][j+1] + GradV_DxMR[0][j+1]) + 0.50*(Vpres[0][j] + Vpres[0][j+1])*(2.0*muWallsMR[0][j+1]*GradV_DyMR[i][j+1] - (2.0/3.0)*muWallsMR[0][j+1]*0.50*(Divergence[0][j] + Divergence[0][j+1]))));

		//Parte derecha
		EnergyViscous[NA-1][j] = (1.0/MESH.VolMP[NA-1][j])*(MESH.SupMP[NA-1][j][0]*(0.50*(Upres[NA-1][j] + Upres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j])*(2.0*muWallsMU[NA-1][j]*GradU_DxMU[NA-1][j] - (2.0/3.0)*muWallsMU[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-2][j])) + 0.50*(Upres[NA-1][j] + Upres[NA-2][j])*sin(PI + MESH.AngleMU[NA-1][j])*muWallsMU[NA-1][j]*(GradU_DyMU[NA-1][j] + GradV_DxMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*cos(PI + MESH.AngleMU[NA-1][j])*muWallsMU[NA-1][j]*(GradU_DyMU[NA-1][j] + GradV_DxMU[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-2][j])*(2.0*muWallsMU[NA-1][j]*GradV_DyMU[NA-1][j] - (2.0/3.0)*muWallsMU[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-2][j]))) +

MESH.SupMP[NA-1][j][1]*(Uright[j]*cos(MESH.AngleMU[NA][j])*(2.0*muWallsMU[NA][j]*GradU_DxMU[NA][j] - (2.0/3.0)*muWallsMU[NA][j]*Divergence[NA-1][j]) + Uright[j]*sin(MESH.AngleMU[NA][j])*muWallsMU[NA][j]*(GradU_DyMU[NA][j] + GradV_DxMU[NA][j]) + Vright[j]*cos(MESH.AngleMU[NA][j])*muWallsMU[NA][j]*(GradU_DyMU[NA][j] + GradV_DxMU[NA][j]) + Vright[j]*(2.0*muWallsMU[NA][j]*GradV_DyMU[NA][j] - (2.0/3.0)*muWallsMU[NA][j]*Divergence[NA-1][j])) + 

MESH.SupMP[NA-1][j][2]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j])*(2.0*muWallsMR[NA-1][j]*GradU_DxMR[NA-1][j] - (2.0/3.0)*muWallsMR[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j-1])) + 0.50*(Upres[NA-1][j] + Upres[NA-1][j-1])*sin(1.50*PI + MESH.AngleMR[NA-1][j])*muWallsMR[NA-1][j]*(GradU_DyMR[NA-1][j] + GradV_DxMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*cos(1.50*PI + MESH.AngleMR[NA-1][j])*muWallsMR[NA-1][j]*(GradU_DyMR[NA-1][j] + GradV_DxMR[NA-1][j]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j-1])*(2.0*muWallsMR[NA-1][j]*GradV_DyMR[NA-1][j] - (2.0/3.0)*muWallsMR[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j-1])))	+

MESH.SupMP[NA-1][j][3]*(0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1])*(2.0*muWallsMR[NA-1][j+1]*GradU_DxMR[NA-1][j+1] - (2.0/3.0)*muWallsMR[NA-1][j+1]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j+1])) + 0.50*(Upres[NA-1][j] + Upres[NA-1][j+1])*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])*muWallsMR[NA-1][j+1]*(GradU_DyMR[NA-1][j+1] + GradV_DxMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*cos(0.50*PI + MESH.AngleMR[NA-1][j+1])*muWallsMR[NA-1][j+1]*(GradU_DyMR[NA-1][j+1] + GradV_DxMR[NA-1][j+1]) + 0.50*(Vpres[NA-1][j] + Vpres[NA-1][j+1])*(2.0*muWallsMR[NA-1][j+1]*GradV_DyMR[NA-1][j+1] - (2.0/3.0)*muWallsMR[NA-1][j+1]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j+1]))));

	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		EnergyViscous[i][0] = (1.0/MESH.VolMP[i][0])*(MESH.SupMP[i][0][0]*(0.50*(Upres[i][0] + Upres[i-1][0])*cos(PI + MESH.AngleMU[i][0])*(2.0*muWallsMU[i][0]*GradU_DxMU[i][0] - (2.0/3.0)*muWallsMU[i][0]*0.50*(Divergence[i][0] + Divergence[i-1][0])) + 0.50*(Upres[i][0] + Upres[i-1][0])*sin(PI + MESH.AngleMU[i][0])*muWallsMU[i][0]*(GradU_DyMU[i][0] + GradV_DxMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*cos(PI + MESH.AngleMU[i][0])*muWallsMU[i][0]*(GradU_DyMU[i][0] + GradV_DxMU[i][0]) + 0.50*(Vpres[i][0] + Vpres[i-1][0])*(2.0*muWallsMU[i][0]*GradV_DyMU[i][0] - (2.0/3.0)*muWallsMU[i][0]*0.50*(Divergence[i][0] + Divergence[i-1][0]))) +

MESH.SupMP[i][0][1]*(0.50*(Upres[i][0] + Upres[i+1][0])*cos(MESH.AngleMU[i+1][0])*(2.0*muWallsMU[i+1][0]*GradU_DxMU[i+1][0] - (2.0/3.0)*muWallsMU[i+1][0]*0.50*(Divergence[i][0] + Divergence[i+1][0])) + 0.50*(Upres[i][0] + Upres[i+1][0])*sin(MESH.AngleMU[i+1][0])*muWallsMU[i+1][0]*(GradU_DyMU[i+1][0] + GradV_DxMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*cos(MESH.AngleMU[i+1][0])*muWallsMU[i+1][0]*(GradU_DyMU[i+1][0] + GradV_DxMU[i+1][0]) + 0.50*(Vpres[i][0] + Vpres[i+1][0])*(2.0*muWallsMU[i+1][0]*GradV_DyMU[i+1][0] - (2.0/3.0)*muWallsMU[i+1][0]*0.50*(Divergence[i][0] + Divergence[i+1][0]))) + 

MESH.SupMP[i][0][2]*(Udown[i]*cos(1.50*PI + MESH.AngleMR[i][0])*(2.0*muWallsMR[i][0]*GradU_DxMR[i][0] - (2.0/3.0)*muWallsMR[i][0]*Divergence[i][0]) + Udown[i]*sin(1.50*PI + MESH.AngleMR[i][0])*muWallsMR[i][0]*(GradU_DyMR[i][0] + GradV_DxMR[i][0]) + Vdown[i]*cos(1.50*PI + MESH.AngleMR[i][0])*muWallsMR[i][0]*(GradU_DyMR[i][0] + GradV_DxMR[i][0]) + Vdown[i]*(2.0*muWallsMR[i][0]*GradV_DyMR[i][0] - (2.0/3.0)*muWallsMR[i][0]*Divergence[i][0]))	+

MESH.SupMP[i][0][3]*(0.50*(Upres[i][0] + Upres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1])*(2.0*muWallsMR[i][1]*GradU_DxMR[i][1] - (2.0/3.0)*muWallsMR[i][1]*0.50*(Divergence[i][0] + Divergence[i][1])) + 0.50*(Upres[i][0] + Upres[i][1])*sin(0.50*PI + MESH.AngleMR[i][1])*muWallsMR[i][1]*(GradU_DyMR[i][1] + GradV_DxMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*cos(0.50*PI + MESH.AngleMR[i][1])*muWallsMR[i][1]*(GradU_DyMR[i][1] + GradV_DxMR[i][1]) + 0.50*(Vpres[i][0] + Vpres[i][1])*(2.0*muWallsMR[i][1]*GradV_DyMR[i][1] - (2.0/3.0)*muWallsMR[i][1]*0.50*(Divergence[i][0] + Divergence[i][1]))));

		//Parte arriba
		EnergyViscous[i][NRad-1] = (1.0/MESH.VolMP[i][NRad-1])*(MESH.SupMP[i][NRad-1][0]*(0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1])*(2.0*muWallsMU[i][NRad-1]*GradU_DxMU[i][NRad-1] - (2.0/3.0)*muWallsMU[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i-1][NRad-1])) + 0.50*(Upres[i][NRad-1] + Upres[i-1][NRad-1])*sin(PI + MESH.AngleMU[i][NRad-1])*muWallsMU[i][NRad-1]*(GradU_DyMU[i][NRad-1] + GradV_DxMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*cos(PI + MESH.AngleMU[i][NRad-1])*muWallsMU[i][NRad-1]*(GradU_DyMU[i][NRad-1] + GradV_DxMU[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i-1][NRad-1])*(2.0*muWallsMU[i][NRad-1]*GradV_DyMU[i][NRad-1] - (2.0/3.0)*muWallsMU[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i-1][NRad-1]))) +

MESH.SupMP[i][NRad-1][1]*(0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1])*(2.0*muWallsMU[i+1][NRad-1]*GradU_DxMU[i+1][NRad-1] - (2.0/3.0)*muWallsMU[i+1][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i+1][NRad-1])) + 0.50*(Upres[i][NRad-1] + Upres[i+1][NRad-1])*sin(MESH.AngleMU[i+1][NRad-1])*muWallsMU[i+1][NRad-1]*(GradU_DyMU[i+1][NRad-1] + GradV_DxMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*cos(MESH.AngleMU[i+1][NRad-1])*muWallsMU[i+1][NRad-1]*(GradU_DyMU[i+1][NRad-1] + GradV_DxMU[i+1][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i+1][NRad-1])*(2.0*muWallsMU[i+1][NRad-1]*GradV_DyMU[i+1][NRad-1] - (2.0/3.0)*muWallsMU[i+1][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i+1][NRad-1]))) + 

MESH.SupMP[i][NRad-1][2]*(0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1])*(2.0*muWallsMR[i][NRad-1]*GradU_DxMR[i][NRad-1] - (2.0/3.0)*muWallsMR[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i][NRad-2])) + 0.50*(Upres[i][NRad-1] + Upres[i][NRad-2])*sin(1.50*PI + MESH.AngleMR[i][NRad-1])*muWallsMR[i][NRad-1]*(GradU_DyMR[i][NRad-1] + GradV_DxMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*cos(1.50*PI + MESH.AngleMR[i][NRad-1])*muWallsMR[i][NRad-1]*(GradU_DyMR[i][NRad-1] + GradV_DxMR[i][NRad-1]) + 0.50*(Vpres[i][NRad-1] + Vpres[i][NRad-2])*(2.0*muWallsMR[i][NRad-1]*GradV_DyMR[i][NRad-1] - (2.0/3.0)*muWallsMR[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i][NRad-2])))	+

MESH.SupMP[i][NRad-1][3]*(Uup[i]*cos(0.50*PI + MESH.AngleMR[i][NRad])*(2.0*muWallsMR[i][NRad]*GradU_DxMR[i][NRad] - (2.0/3.0)*muWallsMR[i][NRad]*0.50*(Divergence[i][NRad-1])) + Uup[i]*sin(0.50*PI + MESH.AngleMR[i][NRad])*muWallsMR[i][NRad]*(GradU_DyMR[i][NRad] + GradV_DxMR[i][NRad]) + Vup[i]*cos(0.50*PI + MESH.AngleMR[i][NRad])*muWallsMR[i][NRad]*(GradU_DyMR[i][NRad] + GradV_DxMR[i][NRad]) + Vup[i]*(2.0*muWallsMR[i][NRad]*GradV_DyMR[i][NRad] - (2.0/3.0)*muWallsMR[i][NRad]*0.50*(Divergence[i][NRad-1] ))));

	}

	//Esquina abajo izquierda
	EnergyViscous[0][0] = (1.0/MESH.VolMP[0][0])*(MESH.SupMP[0][0][0]*(Uleft[0]*cos(PI + MESH.AngleMU[0][0])*(2.0*muWallsMU[0][0]*GradU_DxMU[0][0] - (2.0/3.0)*muWallsMU[0][0]*0.50*(Divergence[0][0])) + Uleft[0]*sin(PI + MESH.AngleMU[0][0])*muWallsMU[0][0]*(GradU_DyMU[0][0] + GradV_DxMU[0][0]) + Vleft[0]*cos(PI + MESH.AngleMU[0][0])*muWallsMU[0][0]*(GradU_DyMU[0][0] + GradV_DxMU[0][0]) + Vleft[0]*(2.0*muWallsMU[0][0]*GradV_DyMU[0][0] - (2.0/3.0)*muWallsMU[0][0]*0.50*(Divergence[0][0]))) +

MESH.SupMP[0][0][1]*(0.50*(Upres[0][0] + Upres[1][0])*cos(MESH.AngleMU[1][0])*(2.0*muWallsMU[1][0]*GradU_DxMU[1][0] - (2.0/3.0)*muWallsMU[1][0]*0.50*(Divergence[0][0] + Divergence[1][0])) + 0.50*(Upres[0][0] + Upres[1][0])*sin(MESH.AngleMU[1][0])*muWallsMU[1][0]*(GradU_DyMU[1][0] + GradV_DxMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*cos(MESH.AngleMU[1][0])*muWallsMU[1][0]*(GradU_DyMU[1][0] + GradV_DxMU[1][0]) + 0.50*(Vpres[0][0] + Vpres[1][0])*(2.0*muWallsMU[1][0]*GradV_DyMU[1][0] - (2.0/3.0)*muWallsMU[1][0]*0.50*(Divergence[0][0] + Divergence[1][0]))) + 

MESH.SupMP[0][0][2]*(Udown[0]*cos(1.50*PI + MESH.AngleMR[0][0])*(2.0*muWallsMR[0][0]*GradU_DxMR[0][0] - (2.0/3.0)*muWallsMR[0][0]*Divergence[0][0]) + Udown[0]*sin(1.50*PI + MESH.AngleMR[0][0])*muWallsMR[0][0]*(GradU_DyMR[0][0] + GradV_DxMR[0][0]) + Vdown[0]*cos(1.50*PI + MESH.AngleMR[0][0])*muWallsMR[0][0]*(GradU_DyMR[0][0] + GradV_DxMR[0][0]) + Vdown[0]*(2.0*muWallsMR[0][0]*GradV_DyMR[0][0] - (2.0/3.0)*muWallsMR[0][0]*Divergence[0][0]))	+

MESH.SupMP[0][0][3]*(0.50*(Upres[0][0] + Upres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1])*(2.0*muWallsMR[0][1]*GradU_DxMR[0][1] - (2.0/3.0)*muWallsMR[0][1]*0.50*(Divergence[0][0] + Divergence[0][1])) + 0.50*(Upres[0][0] + Upres[0][1])*sin(0.50*PI + MESH.AngleMR[0][1])*muWallsMR[0][1]*(GradU_DyMR[0][1] + GradV_DxMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*cos(0.50*PI + MESH.AngleMR[0][1])*muWallsMR[0][1]*(GradU_DyMR[0][1] + GradV_DxMR[0][1]) + 0.50*(Vpres[0][0] + Vpres[0][1])*(2.0*muWallsMR[0][1]*GradV_DyMR[0][1] - (2.0/3.0)*muWallsMR[0][1]*0.50*(Divergence[0][0] + Divergence[0][1]))));

	//Esquina arriba izquierda
	EnergyViscous[0][NRad-1] = (1.0/MESH.VolMP[0][NRad-1])*(MESH.SupMP[0][NRad-1][0]*(Uleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1])*(2.0*muWallsMU[0][NRad-1]*GradU_DxMU[0][NRad-1] - (2.0/3.0)*muWallsMU[0][NRad-1]*0.50*(Divergence[0][NRad-1])) + Uleft[NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1])*muWallsMU[0][NRad-1]*(GradU_DyMU[0][NRad-1] + GradV_DxMU[0][NRad-1]) + Vleft[NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1])*muWallsMU[0][NRad-1]*(GradU_DyMU[0][NRad-1] + GradV_DxMU[0][NRad-1]) + Vleft[NRad-1]*(2.0*muWallsMU[0][NRad-1]*GradV_DyMU[0][NRad-1] - (2.0/3.0)*muWallsMU[0][NRad-1]*0.50*(Divergence[0][NRad-1]))) +

MESH.SupMP[0][NRad-1][1]*(0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1])*(2.0*muWallsMU[1][NRad-1]*GradU_DxMU[1][NRad-1] - (2.0/3.0)*muWallsMU[1][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[1][NRad-1])) + 0.50*(Upres[0][NRad-1] + Upres[1][NRad-1])*sin(MESH.AngleMU[1][NRad-1])*muWallsMU[1][NRad-1]*(GradU_DyMU[1][NRad-1] + GradV_DxMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*cos(MESH.AngleMU[1][NRad-1])*muWallsMU[1][NRad-1]*(GradU_DyMU[1][NRad-1] + GradV_DxMU[1][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[1][NRad-1])*(2.0*muWallsMU[1][NRad-1]*GradV_DyMU[1][NRad-1] - (2.0/3.0)*muWallsMU[1][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[1][NRad-1]))) + 

MESH.SupMP[0][NRad-1][2]*(0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1])*(2.0*muWallsMR[0][NRad-1]*GradU_DxMR[0][NRad-1] - (2.0/3.0)*muWallsMR[0][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[0][NRad-2])) + 0.50*(Upres[0][NRad-1] + Upres[0][NRad-2])*sin(1.50*PI + MESH.AngleMR[0][NRad-1])*muWallsMR[0][NRad-1]*(GradU_DyMR[0][NRad-1] + GradV_DxMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*cos(1.50*PI + MESH.AngleMR[0][NRad-1])*muWallsMR[0][NRad-1]*(GradU_DyMR[0][NRad-1] + GradV_DxMR[0][NRad-1]) + 0.50*(Vpres[0][NRad-1] + Vpres[0][NRad-2])*(2.0*muWallsMR[0][NRad-1]*GradV_DyMR[0][NRad-1] - (2.0/3.0)*muWallsMR[0][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[0][NRad-2])))	+

MESH.SupMP[0][NRad-1][3]*(Uup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad])*(2.0*muWallsMR[0][NRad]*GradU_DxMR[0][NRad] - (2.0/3.0)*muWallsMR[0][NRad]*0.50*(Divergence[0][NRad-1])) + Uup[0]*sin(0.50*PI + MESH.AngleMR[0][NRad])*muWallsMR[0][NRad]*(GradU_DyMR[0][NRad] + GradV_DxMR[0][NRad]) + Vup[0]*cos(0.50*PI + MESH.AngleMR[0][NRad])*muWallsMR[0][NRad]*(GradU_DyMR[0][NRad] + GradV_DxMR[0][NRad]) + Vup[0]*(2.0*muWallsMR[0][NRad]*GradV_DyMR[0][NRad] - (2.0/3.0)*muWallsMR[0][NRad]*0.50*(Divergence[0][NRad-1] ))));
	
	//Esquina abajo derecha
	EnergyViscous[NA-1][0] = (1.0/MESH.VolMP[NA-1][0])*(MESH.SupMP[NA-1][0][0]*(0.50*(Upres[NA-1][0] + Upres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0])*(2.0*muWallsMU[NA-1][0]*GradU_DxMU[NA-1][0] - (2.0/3.0)*muWallsMU[NA-1][0]*0.50*(Divergence[NA-1][0] + Divergence[NA-2][0])) + 0.50*(Upres[NA-1][0] + Upres[NA-2][0])*sin(PI + MESH.AngleMU[NA-1][0])*muWallsMU[NA-1][0]*(GradU_DyMU[NA-1][0] + GradV_DxMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*cos(PI + MESH.AngleMU[NA-1][0])*muWallsMU[NA-1][0]*(GradU_DyMU[NA-1][0] + GradV_DxMU[NA-1][0]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-2][0])*(2.0*muWallsMU[NA-1][0]*GradV_DyMU[NA-1][0] - (2.0/3.0)*muWallsMU[NA-1][0]*0.50*(Divergence[NA-1][0] + Divergence[NA-2][0]))) +

MESH.SupMP[NA-1][0][1]*(Uright[0]*cos(MESH.AngleMU[NA][0])*(2.0*muWallsMU[NA][0]*GradU_DxMU[NA][0] - (2.0/3.0)*muWallsMU[NA][0]*Divergence[NA-1][0]) + Uright[0]*sin(MESH.AngleMU[NA][0])*muWallsMU[NA][0]*(GradU_DyMU[NA][0] + GradV_DxMU[NA][0]) + Vright[0]*cos(MESH.AngleMU[NA][0])*muWallsMU[NA][0]*(GradU_DyMU[NA][0] + GradV_DxMU[NA][0]) + Vright[0]*(2.0*muWallsMU[NA][0]*GradV_DyMU[NA][0] - (2.0/3.0)*muWallsMU[NA][0]*Divergence[NA-1][0])) + 

MESH.SupMP[NA-1][0][2]*(Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0])*(2.0*muWallsMR[NA-1][0]*GradU_DxMR[NA-1][0] - (2.0/3.0)*muWallsMR[NA-1][0]*Divergence[NA-1][0]) + Udown[NA-1]*sin(1.50*PI + MESH.AngleMR[NA-1][0])*muWallsMR[NA-1][0]*(GradU_DyMR[NA-1][0] + GradV_DxMR[NA-1][0]) + Udown[NA-1]*cos(1.50*PI + MESH.AngleMR[NA-1][0])*muWallsMR[NA-1][0]*(GradU_DyMR[NA-1][0] + GradV_DxMR[NA-1][0]) + Vdown[NA-1]*(2.0*muWallsMR[NA-1][0]*GradV_DyMR[NA-1][0] - (2.0/3.0)*muWallsMR[NA-1][0]*Divergence[NA-1][0]))	+

MESH.SupMP[NA-1][0][3]*(0.50*(Upres[NA-1][0] + Upres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1])*(2.0*muWallsMR[NA-1][1]*GradU_DxMR[NA-1][1] - (2.0/3.0)*muWallsMR[NA-1][1]*0.50*(Divergence[NA-1][0] + Divergence[NA-1][1])) + 0.50*(Upres[NA-1][0] + Upres[NA-1][1])*sin(0.50*PI + MESH.AngleMR[NA-1][1])*muWallsMR[NA-1][1]*(GradU_DyMR[NA-1][1] + GradV_DxMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*cos(0.50*PI + MESH.AngleMR[NA-1][1])*muWallsMR[NA-1][1]*(GradU_DyMR[NA-1][1] + GradV_DxMR[NA-1][1]) + 0.50*(Vpres[NA-1][0] + Vpres[NA-1][1])*(2.0*muWallsMR[NA-1][1]*GradV_DyMR[NA-1][1] - (2.0/3.0)*muWallsMR[NA-1][1]*0.50*(Divergence[NA-1][0] + Divergence[NA-1][1]))));

	//Esquina arriba derecha
	EnergyViscous[NA-1][NRad-1] = (1.0/MESH.VolMP[NA-1][NRad-1])*(MESH.SupMP[NA-1][NRad-1][0]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1])*(2.0*muWallsMU[NA-1][NRad-1]*GradU_DxMU[NA-1][NRad-1] - (2.0/3.0)*muWallsMU[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-2][NRad-1])) + 0.50*(Upres[NA-1][NRad-1] + Upres[NA-2][NRad-1])*sin(PI + MESH.AngleMU[NA-1][NRad-1])*muWallsMU[NA-1][NRad-1]*(GradU_DyMU[NA-1][NRad-1] + GradV_DxMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*cos(PI + MESH.AngleMU[NA-1][NRad-1])*muWallsMU[NA-1][NRad-1]*(GradU_DyMU[NA-1][NRad-1] + GradV_DxMU[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-2][NRad-1])*(2.0*muWallsMU[NA-1][NRad-1]*GradV_DyMU[NA-1][NRad-1] - (2.0/3.0)*muWallsMU[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-2][NRad-1]))) +

MESH.SupMP[NA-1][NRad-1][1]*(Uright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1])*(2.0*muWallsMU[NA][NRad-1]*GradU_DxMU[NA][NRad-1] - (2.0/3.0)*muWallsMU[NA][NRad-1]*Divergence[NA-1][NRad-1]) + Uright[NRad-1]*sin(MESH.AngleMU[NA][NRad-1])*muWallsMU[NA][NRad-1]*(GradU_DyMU[NA][NRad-1] + GradV_DxMU[NA][NRad-1]) + Vright[NRad-1]*cos(MESH.AngleMU[NA][NRad-1])*muWallsMU[NA][NRad-1]*(GradU_DyMU[NA][NRad-1] + GradV_DxMU[NA][NRad-1]) + Vright[NRad-1]*(2.0*muWallsMU[NA][NRad-1]*GradV_DyMU[NA][NRad-1] - (2.0/3.0)*muWallsMU[NA][NRad-1]*Divergence[NA-1][NRad-1])) + 

MESH.SupMP[NA-1][NRad-1][2]*(0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*(2.0*muWallsMR[NA-1][NRad-1]*GradU_DxMR[NA-1][NRad-1] - (2.0/3.0)*muWallsMR[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-1][NRad-2])) + 0.50*(Upres[NA-1][NRad-1] + Upres[NA-1][NRad-2])*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*muWallsMR[NA-1][NRad-1]*(GradU_DyMR[NA-1][NRad-1] + GradV_DxMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*muWallsMR[NA-1][NRad-1]*(GradU_DyMR[NA-1][NRad-1] + GradV_DxMR[NA-1][NRad-1]) + 0.50*(Vpres[NA-1][NRad-1] + Vpres[NA-1][NRad-2])*(2.0*muWallsMR[NA-1][NRad-1]*GradV_DyMR[NA-1][NRad-1] - (2.0/3.0)*muWallsMR[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-1][NRad-2])))	+

MESH.SupMP[NA-1][NRad-1][3]*(Uup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad])*(2.0*muWallsMR[NA-1][NRad]*GradU_DxMR[NA-1][NRad] - (2.0/3.0)*muWallsMR[NA-1][NRad]*0.50*(Divergence[NA-1][NRad-1])) + Uup[NA-1]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])*muWallsMR[NA-1][NRad]*(GradU_DyMR[NA-1][NRad] + GradV_DxMR[NA-1][NRad]) + Vup[NA-1]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad])*muWallsMR[NA-1][NRad]*(GradU_DyMR[NA-1][NRad] + GradV_DxMR[NA-1][NRad]) + Vup[NA-1]*(2.0*muWallsMR[NA-1][NRad]*GradV_DyMR[NA-1][NRad] - (2.0/3.0)*muWallsMR[NA-1][NRad]*Divergence[NA-1][NRad-1])));

}