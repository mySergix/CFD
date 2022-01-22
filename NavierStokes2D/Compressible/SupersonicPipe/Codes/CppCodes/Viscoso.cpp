for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			//Centro
			EnergyViscous[i][j] = (1.0/MESH.VolMP[i][j])*(
								  	MESH.SupMP[i][j][0]*(
								  		UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j])*(2.0*muWallsMU[i][j]*GradU_DxMU[i][j] - (2.0/3.0)*muWallsMU[i][j]*0.50*(Divergence[i][j] + Divergence[i-1][j]))
								  	  + UwallsMU[i][j]*cos(PI + MESH.AngleMU[i][j])*muWallsMU[i][j]*(GradU_DyMU[i][j] + GradV_DxMU[i][j])
								  	  + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j])*muWallsMU[i][j]*(GradU_DyMU[i][j] + GradV_DxMU[i][j])
								  	  + VwallsMU[i][j]*sin(PI + MESH.AngleMU[i][j])*(2.0*muWallsMU[i][j]*GradV_DyMU[i][j] - (2.0/3.0)*muWallsMU[i][j]*0.50*(Divergence[i][j] + Divergence[i-1][j]))
								    )
								  + MESH.SupMP[i][j][1]*(
								  		UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j])*(2.0*muWallsMU[i+1][j]*GradU_DxMU[i+1][j] - (2.0/3.0)*muWallsMU[i+1][j]*0.50*(Divergence[i][j] + Divergence[i+1][j]))
								  	  + UwallsMU[i+1][j]*cos(MESH.AngleMU[i+1][j])*muWallsMU[i+1][j]*(GradU_DyMU[i+1][j] + GradV_DxMU[i+1][j])
								  	  + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j])*muWallsMU[i+1][j]*(GradU_DyMU[i+1][j] + GradV_DxMU[i+1][j])
								  	  + VwallsMU[i+1][j]*sin(MESH.AngleMU[i+1][j])*(2.0*muWallsMU[i+1][j]*GradV_DyMU[i+1][j] - (2.0/3.0)*muWallsMU[i+1][j]*0.50*(Divergence[i][j] + Divergence[i+1][j]))
								    )
								  + MESH.SupMP[i][j][2]*(
								  		UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j])*(2.0*muWallsMR[i][j]*GradU_DxMR[i][j] - (2.0/3.0)*muWallsMR[i][j]*0.50*(Divergence[i][j] + Divergence[i][j-1]))
								  	  + UwallsMR[i][j]*cos(1.50*PI + MESH.AngleMR[i][j])*muWallsMR[i][j]*(GradU_DyMR[i][j] + GradV_DxMR[i][j])
								  	  + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j])*muWallsMR[i][j]*(GradU_DyMR[i][j] + GradV_DxMR[i][j])
								  	  + VwallsMR[i][j]*sin(1.50*PI + MESH.AngleMR[i][j])*(2.0*muWallsMR[i][j]*GradV_DyMR[i][j] - (2.0/3.0)*muWallsMR[i][j]*0.50*(Divergence[i][j] + Divergence[i][j-1]))
								    )
								  + MESH.SupMP[i][j][3]*(
								  		UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1])*(2.0*muWallsMR[i][j+1]*GradU_DxMR[i][j+1] - (2.0/3.0)*muWallsMR[i][j+1]*0.50*(Divergence[i][j] + Divergence[i][j+1]))
								  	  + UwallsMR[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1])*muWallsMR[i][j+1]*(GradU_DyMR[i][j+1] + GradV_DxMR[i][j+1])
								  	  + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1])*muWallsMR[i][j+1]*(GradU_DyMR[i][j+1] + GradV_DxMR[i][j+1])
								  	  + VwallsMR[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1])*(2.0*muWallsMR[i][j+1]*GradV_DyMR[i][j+1] - (2.0/3.0)*muWallsMR[i][j+1]*0.50*(Divergence[i][j] + Divergence[i][j+1]))
								    )
								);
		}
	}

	for(j = 1; j < NRad-1; j++){

		//Parte izquierda
		EnergyViscous[0][j] = (1.0/MESH.VolMP[0][j])*(
								  	MESH.SupMP[0][j][0]*(
								  		UwallsMU[0][j]*cos(PI + MESH.AngleMU[0][j])*(2.0*muWallsMU[0][j]*GradU_DxMU[0][j] - (2.0/3.0)*muWallsMU[0][j]*0.50*(Divergence[0][j]))
								  	  + UwallsMU[0][j]*cos(PI + MESH.AngleMU[0][j])*muWallsMU[0][j]*(GradU_DyMU[0][j] + GradV_DxMU[0][j])
								  	  + VwallsMU[0][j]*sin(PI + MESH.AngleMU[0][j])*muWallsMU[0][j]*(GradU_DyMU[0][j] + GradV_DxMU[0][j])
								  	  + VwallsMU[0][j]*sin(PI + MESH.AngleMU[0][j])*(2.0*muWallsMU[0][j]*GradV_DyMU[0][j] - (2.0/3.0)*muWallsMU[0][j]*0.50*(Divergence[0][j]))
								    )
								  + MESH.SupMP[0][j][1]*(
								  		UwallsMU[1][j]*cos(MESH.AngleMU[1][j])*(2.0*muWallsMU[1][j]*GradU_DxMU[1][j] - (2.0/3.0)*muWallsMU[1][j]*0.50*(Divergence[0][j] + Divergence[1][j]))
								  	  + UwallsMU[1][j]*cos(MESH.AngleMU[1][j])*muWallsMU[1][j]*(GradU_DyMU[1][j] + GradV_DxMU[1][j])
								  	  + VwallsMU[1][j]*sin(MESH.AngleMU[1][j])*muWallsMU[1][j]*(GradU_DyMU[1][j] + GradV_DxMU[1][j])
								  	  + VwallsMU[1][j]*sin(MESH.AngleMU[1][j])*(2.0*muWallsMU[1][j]*GradV_DyMU[1][j] - (2.0/3.0)*muWallsMU[1][j]*0.50*(Divergence[0][j] + Divergence[1][j]))
								    )
								  + MESH.SupMP[0][j][2]*(
								  		UwallsMR[0][j]*cos(1.50*PI + MESH.AngleMR[0][j])*(2.0*muWallsMR[0][j]*GradU_DxMR[0][j] - (2.0/3.0)*muWallsMR[0][j]*0.50*(Divergence[0][j] + Divergence[0][j-1]))
								  	  + UwallsMR[0][j]*cos(1.50*PI + MESH.AngleMR[0][j])*muWallsMR[0][j]*(GradU_DyMR[0][j] + GradV_DxMR[0][j])
								  	  + VwallsMR[0][j]*sin(1.50*PI + MESH.AngleMR[0][j])*muWallsMR[0][j]*(GradU_DyMR[0][j] + GradV_DxMR[0][j])
								  	  + VwallsMR[0][j]*sin(1.50*PI + MESH.AngleMR[0][j])*(2.0*muWallsMR[0][j]*GradV_DyMR[0][j] - (2.0/3.0)*muWallsMR[0][j]*0.50*(Divergence[0][j] + Divergence[0][j-1]))
								    )
								  + MESH.SupMP[0][j][3]*(
								  		UwallsMR[0][j+1]*cos(0.50*PI + MESH.AngleMR[0][j+1])*(2.0*muWallsMR[0][j+1]*GradU_DxMR[0][j+1] - (2.0/3.0)*muWallsMR[0][j+1]*0.50*(Divergence[0][j] + Divergence[0][j+1]))
								  	  + UwallsMR[0][j+1]*cos(0.50*PI + MESH.AngleMR[0][j+1])*muWallsMR[0][j+1]*(GradU_DyMR[0][j+1] + GradV_DxMR[0][j+1])
								  	  + VwallsMR[0][j+1]*sin(0.50*PI + MESH.AngleMR[0][j+1])*muWallsMR[0][j+1]*(GradU_DyMR[0][j+1] + GradV_DxMR[0][j+1])
								  	  + VwallsMR[0][j+1]*sin(0.50*PI + MESH.AngleMR[0][j+1])*(2.0*muWallsMR[0][j+1]*GradV_DyMR[0][j+1] - (2.0/3.0)*muWallsMR[0][j+1]*0.50*(Divergence[0][j] + Divergence[0][j+1]))
								    )
								);
		//Parte derecha
		EnergyViscous[NA-1][j] = (1.0/MESH.VolMP[NA-1][j])*(
								  	MESH.SupMP[NA-1][j][0]*(
								  		UwallsMU[NA-1][j]*cos(PI + MESH.AngleMU[NA-1][j])*(2.0*muWallsMU[NA-1][j]*GradU_DxMU[NA-1][j] - (2.0/3.0)*muWallsMU[i][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-2][j]))
								  	  + UwallsMU[NA-1][j]*cos(PI + MESH.AngleMU[NA-1][j])*muWallsMU[NA-1][j]*(GradU_DyMU[NA-1][j] + GradV_DxMU[NA-1][j])
								  	  + VwallsMU[NA-1][j]*sin(PI + MESH.AngleMU[NA-1][j])*muWallsMU[NA-1][j]*(GradU_DyMU[NA-1][j] + GradV_DxMU[NA-1][j])
								  	  + VwallsMU[NA-1][j]*sin(PI + MESH.AngleMU[NA-1][j])*(2.0*muWallsMU[NA-1][j]*GradV_DyMU[NA-1][j] - (2.0/3.0)*muWallsMU[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-2][j]))
								    )
								  + MESH.SupMP[NA-1][j][1]*(
								  		UwallsMU[NA][j]*cos(MESH.AngleMU[NA][j])*(2.0*muWallsMU[NA][j]*GradU_DxMU[NA][j] - (2.0/3.0)*muWallsMU[NA][j]*0.50*(Divergence[NA-1][j]))
								  	  + UwallsMU[NA][j]*cos(MESH.AngleMU[NA][j])*muWallsMU[NA][j]*(GradU_DyMU[NA][j] + GradV_DxMU[NA][j])
								  	  + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][j])*muWallsMU[NA][j]*(GradU_DyMU[NA][j] + GradV_DxMU[NA][j])
								  	  + VwallsMU[NA][j]*sin(MESH.AngleMU[NA][j])*(2.0*muWallsMU[NA][j]*GradV_DyMU[NA][j] - (2.0/3.0)*muWallsMU[NA][j]*0.50*(Divergence[NA-1][j]))
								    )
								  + MESH.SupMP[NA-1][j][2]*(
								  		UwallsMR[NA-1][j]*cos(1.50*PI + MESH.AngleMR[NA-1][j])*(2.0*muWallsMR[NA-1][j]*GradU_DxMR[NA-1][j] - (2.0/3.0)*muWallsMR[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j-1]))
								  	  + UwallsMR[NA-1][j]*cos(1.50*PI + MESH.AngleMR[NA-1][j])*muWallsMR[NA-1][j]*(GradU_DyMR[NA-1][j] + GradV_DxMR[NA-1][j])
								  	  + VwallsMR[NA-1][j]*sin(1.50*PI + MESH.AngleMR[NA-1][j])*muWallsMR[NA-1][j]*(GradU_DyMR[NA-1][j] + GradV_DxMR[NA-1][j])
								  	  + VwallsMR[NA-1][j]*sin(1.50*PI + MESH.AngleMR[NA-1][j])*(2.0*muWallsMR[NA-1][j]*GradV_DyMR[NA-1][j] - (2.0/3.0)*muWallsMR[NA-1][j]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j-1]))
								    )
								  + MESH.SupMP[NA-1][j][3]*(
								  		UwallsMR[NA-1][j+1]*cos(0.50*PI + MESH.AngleMR[NA-1][j+1])*(2.0*muWallsMR[NA-1][j+1]*GradU_DxMR[NA-1][j+1] - (2.0/3.0)*muWallsMR[NA-1][j+1]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j+1]))
								  	  + UwallsMR[NA-1][j+1]*cos(0.50*PI + MESH.AngleMR[NA-1][j+1])*muWallsMR[NA-1][j+1]*(GradU_DyMR[NA-1][j+1] + GradV_DxMR[NA-1][j+1])
								  	  + VwallsMR[NA-1][j+1]*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])*muWallsMR[NA-1][j+1]*(GradU_DyMR[NA-1][j+1] + GradV_DxMR[NA-1][j+1])
								  	  + VwallsMR[NA-1][j+1]*sin(0.50*PI + MESH.AngleMR[NA-1][j+1])*(2.0*muWallsMR[NA-1][j+1]*GradV_DyMR[NA-1][j+1] - (2.0/3.0)*muWallsMR[NA-1][j+1]*0.50*(Divergence[NA-1][j] + Divergence[NA-1][j+1]))
								    )
								);

	}

	for(i = 1; i < NA-1; i++){

		//Parte abajo
			EnergyViscous[i][0] = (1.0/MESH.VolMP[i][0])*(
								  	MESH.SupMP[i][0][0]*(
								  		UwallsMU[i][0]*cos(PI + MESH.AngleMU[i][0])*(2.0*muWallsMU[i][0]*GradU_DxMU[i][0] - (2.0/3.0)*muWallsMU[i][0]*0.50*(Divergence[i][0] + Divergence[i-1][0]))
								  	  + UwallsMU[i][0]*cos(PI + MESH.AngleMU[i][0])*muWallsMU[i][0]*(GradU_DyMU[i][0] + GradV_DxMU[i][0])
								  	  + VwallsMU[i][0]*sin(PI + MESH.AngleMU[i][0])*muWallsMU[i][0]*(GradU_DyMU[i][0] + GradV_DxMU[i][0])
								  	  + VwallsMU[i][0]*sin(PI + MESH.AngleMU[i][0])*(2.0*muWallsMU[i][0]*GradV_DyMU[i][0] - (2.0/3.0)*muWallsMU[i][0]*0.50*(Divergence[i][0] + Divergence[i-1][0]))
								    )
								  + MESH.SupMP[i][0][1]*(
								  		UwallsMU[i+1][0]*cos(MESH.AngleMU[i+1][0])*(2.0*muWallsMU[i+1][0]*GradU_DxMU[i+1][0] - (2.0/3.0)*muWallsMU[i+1][0]*0.50*(Divergence[i][0] + Divergence[i+1][0]))
								  	  + UwallsMU[i+1][0]*cos(MESH.AngleMU[i+1][0])*muWallsMU[i+1][0]*(GradU_DyMU[i+1][0] + GradV_DxMU[i+1][0])
								  	  + VwallsMU[i+1][0]*sin(MESH.AngleMU[i+1][0])*muWallsMU[i+1][0]*(GradU_DyMU[i+1][0] + GradV_DxMU[i+1][0])
								  	  + VwallsMU[i+1][0]*sin(MESH.AngleMU[i+1][0])*(2.0*muWallsMU[i+1][0]*GradV_DyMU[i+1][0] - (2.0/3.0)*muWallsMU[i+1][0]*0.50*(Divergence[i][0] + Divergence[i+1][0]))
								    )
								  + MESH.SupMP[i][0][2]*(
								  		UwallsMR[i][0]*cos(1.50*PI + MESH.AngleMR[i][0])*(2.0*muWallsMR[i][0]*GradU_DxMR[i][0] - (2.0/3.0)*muWallsMR[i][0]*0.50*(Divergence[i][0]))
								  	  + UwallsMR[i][0]*cos(1.50*PI + MESH.AngleMR[i][0])*muWallsMR[i][0]*(GradU_DyMR[i][0] + GradV_DxMR[i][0])
								  	  + VwallsMR[i][0]*sin(1.50*PI + MESH.AngleMR[i][0])*muWallsMR[i][0]*(GradU_DyMR[i][0] + GradV_DxMR[i][0])
								  	  + VwallsMR[i][0]*sin(1.50*PI + MESH.AngleMR[i][0])*(2.0*muWallsMR[i][0]*GradV_DyMR[i][0] - (2.0/3.0)*muWallsMR[i][0]*0.50*(Divergence[i][0]))
								    )
								  + MESH.SupMP[i][0][3]*(
								  		UwallsMR[i][1]*cos(0.50*PI + MESH.AngleMR[i][1])*(2.0*muWallsMR[i][1]*GradU_DxMR[i][1] - (2.0/3.0)*muWallsMR[i][1]*0.50*(Divergence[i][0] + Divergence[i][1]))
								  	  + UwallsMR[i][1]*cos(0.50*PI + MESH.AngleMR[i][1])*muWallsMR[i][1]*(GradU_DyMR[i][1] + GradV_DxMR[i][1])
								  	  + VwallsMR[i][1]*sin(0.50*PI + MESH.AngleMR[i][1])*muWallsMR[i][1]*(GradU_DyMR[i][1] + GradV_DxMR[i][1])
								  	  + VwallsMR[i][1]*sin(0.50*PI + MESH.AngleMR[i][1])*(2.0*muWallsMR[i][1]*GradV_DyMR[i][1] - (2.0/3.0)*muWallsMR[i][1]*0.50*(Divergence[i][0] + Divergence[i][1]))
								    )
								);
		//Parte arriba
		EnergyViscous[i][NRad-1] = (1.0/MESH.VolMP[i][NRad-1])*(
								  	MESH.SupMP[i][NRad-1][0]*(
								  		UwallsMU[i][NRad-1]*cos(PI + MESH.AngleMU[i][NRad-1])*(2.0*muWallsMU[i][NRad-1]*GradU_DxMU[i][NRad-1] - (2.0/3.0)*muWallsMU[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i-1][NRad-1]))
								  	  + UwallsMU[i][NRad-1]*cos(PI + MESH.AngleMU[i][NRad-1])*muWallsMU[i][NRad-1]*(GradU_DyMU[i][NRad-1] + GradV_DxMU[i][NRad-1])
								  	  + VwallsMU[i][NRad-1]*sin(PI + MESH.AngleMU[i][NRad-1])*muWallsMU[i][NRad-1]*(GradU_DyMU[i][NRad-1] + GradV_DxMU[i][NRad-1])
								  	  + VwallsMU[i][NRad-1]*sin(PI + MESH.AngleMU[i][NRad-1])*(2.0*muWallsMU[i][NRad-1]*GradV_DyMU[i][NRad-1] - (2.0/3.0)*muWallsMU[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i-1][NRad-1]))
								    )
								  + MESH.SupMP[i][NRad-1][1]*(
								  		UwallsMU[i+1][NRad-1]*cos(MESH.AngleMU[i+1][NRad-1])*(2.0*muWallsMU[i+1][NRad-1]*GradU_DxMU[i+1][NRad-1] - (2.0/3.0)*muWallsMU[i+1][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i+1][NRad-1]))
								  	  + UwallsMU[i+1][NRad-1]*cos(MESH.AngleMU[i+1][NRad-1])*muWallsMU[i+1][NRad-1]*(GradU_DyMU[i+1][NRad-1] + GradV_DxMU[i+1][NRad-1])
								  	  + VwallsMU[i+1][NRad-1]*sin(MESH.AngleMU[i+1][NRad-1])*muWallsMU[i+1][NRad-1]*(GradU_DyMU[i+1][NRad-1] + GradV_DxMU[i+1][NRad-1])
								  	  + VwallsMU[i+1][NRad-1]*sin(MESH.AngleMU[i+1][NRad-1])*(2.0*muWallsMU[i+1][NRad-1]*GradV_DyMU[i+1][NRad-1] - (2.0/3.0)*muWallsMU[i+1][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i+1][NRad-1]))
								    )
								  + MESH.SupMP[i][NRad-1][2]*(
								  		UwallsMR[i][NRad-1]*cos(1.50*PI + MESH.AngleMR[i][NRad-1])*(2.0*muWallsMR[i][NRad-1]*GradU_DxMR[i][NRad-1] - (2.0/3.0)*muWallsMR[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i][NRad-2]))
								  	  + UwallsMR[i][NRad-1]*cos(1.50*PI + MESH.AngleMR[i][NRad-1])*muWallsMR[i][NRad-1]*(GradU_DyMR[i][NRad-1] + GradV_DxMR[i][NRad-1])
								  	  + VwallsMR[i][NRad-1]*sin(1.50*PI + MESH.AngleMR[i][NRad-1])*muWallsMR[i][NRad-1]*(GradU_DyMR[i][NRad-1] + GradV_DxMR[i][NRad-1])
								  	  + VwallsMR[i][NRad-1]*sin(1.50*PI + MESH.AngleMR[i][NRad-1])*(2.0*muWallsMR[i][NRad-1]*GradV_DyMR[i][NRad-1] - (2.0/3.0)*muWallsMR[i][NRad-1]*0.50*(Divergence[i][NRad-1] + Divergence[i][NRad-2]))
								    )
								  + MESH.SupMP[i][NRad-1][3]*(
								  		UwallsMR[i][NRad]*cos(0.50*PI + MESH.AngleMR[i][NRad])*(2.0*muWallsMR[i][NRad]*GradU_DxMR[i][NRad] - (2.0/3.0)*muWallsMR[i][NRad]*0.50*(Divergence[i][NRad-1]))
								  	  + UwallsMR[i][NRad]*cos(0.50*PI + MESH.AngleMR[i][NRad])*muWallsMR[i][NRad]*(GradU_DyMR[i][NRad] + GradV_DxMR[i][NRad])
								  	  + VwallsMR[i][NRad]*sin(0.50*PI + MESH.AngleMR[i][NRad])*muWallsMR[i][NRad]*(GradU_DyMR[i][NRad] + GradV_DxMR[i][NRad])
								  	  + VwallsMR[i][NRad]*sin(0.50*PI + MESH.AngleMR[i][NRad])*(2.0*muWallsMR[i][NRad]*GradV_DyMR[i][NRad] - (2.0/3.0)*muWallsMR[i][NRad]*0.50*(Divergence[i][NRad-1]))
								    )
								);
									
	}

	//Esquina abajo izquierda
	EnergyViscous[0][0] = (1.0/MESH.VolMP[0][0])*(
								  	MESH.SupMP[0][0][0]*(
								  		UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0])*(2.0*muWallsMU[0][0]*GradU_DxMU[0][0] - (2.0/3.0)*muWallsMU[0][0]*0.50*(Divergence[0][0]))
								  	  + UwallsMU[0][0]*cos(PI + MESH.AngleMU[0][0])*muWallsMU[0][0]*(GradU_DyMU[0][0] + GradV_DxMU[0][0])
								  	  + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0])*muWallsMU[0][0]*(GradU_DyMU[0][0] + GradV_DxMU[0][0])
								  	  + VwallsMU[0][0]*sin(PI + MESH.AngleMU[0][0])*(2.0*muWallsMU[0][0]*GradV_DyMU[0][0] - (2.0/3.0)*muWallsMU[0][0]*0.50*(Divergence[0][0]))
								    )
								  + MESH.SupMP[0][0][1]*(
								  		UwallsMU[1][0]*cos(MESH.AngleMU[1][0])*(2.0*muWallsMU[1][0]*GradU_DxMU[1][0] - (2.0/3.0)*muWallsMU[1][0]*0.50*(Divergence[0][0] + Divergence[1][0]))
								  	  + UwallsMU[1][0]*cos(MESH.AngleMU[1][0])*muWallsMU[1][0]*(GradU_DyMU[1][0] + GradV_DxMU[1][0])
								  	  + VwallsMU[1][0]*sin(MESH.AngleMU[1][0])*muWallsMU[1][0]*(GradU_DyMU[1][0] + GradV_DxMU[1][0])
								  	  + VwallsMU[1][0]*sin(MESH.AngleMU[1][0])*(2.0*muWallsMU[1][0]*GradV_DyMU[1][0] - (2.0/3.0)*muWallsMU[1][0]*0.50*(Divergence[0][0] + Divergence[1][0]))
								    )
								  + MESH.SupMP[0][0][2]*(
								  		UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0])*(2.0*muWallsMR[0][0]*GradU_DxMR[0][0] - (2.0/3.0)*muWallsMR[0][0]*0.50*(Divergence[0][0]))
								  	  + UwallsMR[0][0]*cos(1.50*PI + MESH.AngleMR[0][0])*muWallsMR[0][0]*(GradU_DyMR[0][0] + GradV_DxMR[0][0])
								  	  + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0])*muWallsMR[0][0]*(GradU_DyMR[0][0] + GradV_DxMR[0][0])
								  	  + VwallsMR[0][0]*sin(1.50*PI + MESH.AngleMR[0][0])*(2.0*muWallsMR[0][0]*GradV_DyMR[i][j] - (2.0/3.0)*muWallsMR[0][0]*0.50*(Divergence[0][0]))
								    )
								  + MESH.SupMP[0][0][3]*(
								  		UwallsMR[0][1]*cos(0.50*PI + MESH.AngleMR[0][1])*(2.0*muWallsMR[0][1]*GradU_DxMR[0][1] - (2.0/3.0)*muWallsMR[0][1]*0.50*(Divergence[0][0] + Divergence[0][1]))
								  	  + UwallsMR[0][1]*cos(0.50*PI + MESH.AngleMR[0][1])*muWallsMR[0][1]*(GradU_DyMR[0][1] + GradV_DxMR[0][1])
								  	  + VwallsMR[0][1]*sin(0.50*PI + MESH.AngleMR[0][1])*muWallsMR[0][1]*(GradU_DyMR[0][1] + GradV_DxMR[0][1])
								  	  + VwallsMR[0][1]*sin(0.50*PI + MESH.AngleMR[0][1])*(2.0*muWallsMR[0][1]*GradV_DyMR[0][1] - (2.0/3.0)*muWallsMR[0][1]*0.50*(Divergence[0][0] + Divergence[0][1]))
								    )
								);

	//Esquina arriba izquierda
	EnergyViscous[0][NRad-1] = (1.0/MESH.VolMP[0][NRad-1])*(
								  	MESH.SupMP[0][NRad-1][0]*(
								  		UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1])*(2.0*muWallsMU[0][NRad-1]*GradU_DxMU[0][NRad-1] - (2.0/3.0)*muWallsMU[0][NRad-1]*0.50*(Divergence[0][NRad-1]))
								  	  + UwallsMU[0][NRad-1]*cos(PI + MESH.AngleMU[0][NRad-1])*muWallsMU[0][NRad-1]*(GradU_DyMU[0][NRad-1] + GradV_DxMU[0][NRad-1])
								  	  + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1])*muWallsMU[0][NRad-1]*(GradU_DyMU[0][NRad-1] + GradV_DxMU[0][NRad-1])
								  	  + VwallsMU[0][NRad-1]*sin(PI + MESH.AngleMU[0][NRad-1])*(2.0*muWallsMU[0][NRad-1]*GradV_DyMU[0][NRad-1] - (2.0/3.0)*muWallsMU[0][NRad-1]*0.50*(Divergence[0][NRad-1]))
								    )
								  + MESH.SupMP[0][NRad-1][1]*(
								  		UwallsMU[1][NRad-1]*cos(MESH.AngleMU[1][NRad-1])*(2.0*muWallsMU[1][NRad-1]*GradU_DxMU[1][NRad-1] - (2.0/3.0)*muWallsMU[1][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[1][NRad-1]))
								  	  + UwallsMU[1][NRad-1]*cos(MESH.AngleMU[1][NRad-1])*muWallsMU[1][NRad-1]*(GradU_DyMU[1][NRad-1] + GradV_DxMU[1][NRad-1])
								  	  + VwallsMU[1][NRad-1]*sin(MESH.AngleMU[1][NRad-1])*muWallsMU[1][NRad-1]*(GradU_DyMU[1][NRad-1] + GradV_DxMU[1][NRad-1])
								  	  + VwallsMU[1][NRad-1]*sin(MESH.AngleMU[1][NRad-1])*(2.0*muWallsMU[1][NRad-1]*GradV_DyMU[1][NRad-1] - (2.0/3.0)*muWallsMU[1][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[1][NRad-1]))
								    )
								  + MESH.SupMP[0][NRad-1][2]*(
								  		UwallsMR[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1])*(2.0*muWallsMR[0][NRad-1]*GradU_DxMR[0][NRad-1] - (2.0/3.0)*muWallsMR[0][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[0][NRad-2]))
								  	  + UwallsMR[0][NRad-1]*cos(1.50*PI + MESH.AngleMR[0][NRad-1])*muWallsMR[0][NRad-1]*(GradU_DyMR[0][NRad-1] + GradV_DxMR[0][NRad-1])
								  	  + VwallsMR[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad-1])*muWallsMR[0][NRad-1]*(GradU_DyMR[0][NRad-1] + GradV_DxMR[0][NRad-1])
								  	  + VwallsMR[0][NRad-1]*sin(1.50*PI + MESH.AngleMR[0][NRad-1])*(2.0*muWallsMR[0][NRad-1]*GradV_DyMR[0][NRad-1] - (2.0/3.0)*muWallsMR[0][NRad-1]*0.50*(Divergence[0][NRad-1] + Divergence[0][NRad-2]))
								    )
								  + MESH.SupMP[0][NRad-1][3]*(
								  		UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad])*(2.0*muWallsMR[0][NRad]*GradU_DxMR[0][NRad] - (2.0/3.0)*muWallsMR[0][NRad]*0.50*(Divergence[0][NRad-1]))
								  	  + UwallsMR[0][NRad]*cos(0.50*PI + MESH.AngleMR[0][NRad])*muWallsMR[0][NRad]*(GradU_DyMR[0][NRad] + GradV_DxMR[0][NRad])
								  	  + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad])*muWallsMR[0][NRad]*(GradU_DyMR[0][NRad] + GradV_DxMR[0][NRad])
								  	  + VwallsMR[0][NRad]*sin(0.50*PI + MESH.AngleMR[0][NRad])*(2.0*muWallsMR[0][NRad]*GradV_DyMR[0][NRad] - (2.0/3.0)*muWallsMR[0][NRad]*0.50*(Divergence[0][NRad-1]))
								    )
								);

	//Esquina abajo derecha
	EnergyViscous[NA-1][0] = (1.0/MESH.VolMP[NA-1][0])*(
								  	MESH.SupMP[NA-1][0][0]*(
								  		UwallsMU[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0])*(2.0*muWallsMU[NA-1][0]*GradU_DxMU[NA-1][0] - (2.0/3.0)*muWallsMU[NA-1][0]*0.50*(Divergence[NA-1][0] + Divergence[NA-2][0]))
								  	  + UwallsMU[NA-1][0]*cos(PI + MESH.AngleMU[NA-1][0])*muWallsMU[NA-1][0]*(GradU_DyMU[NA-1][0] + GradV_DxMU[NA-1][0])
								  	  + VwallsMU[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0])*muWallsMU[NA-1][0]*(GradU_DyMU[NA-1][0] + GradV_DxMU[NA-1][0])
								  	  + VwallsMU[NA-1][0]*sin(PI + MESH.AngleMU[NA-1][0])*(2.0*muWallsMU[NA-1][0]*GradV_DyMU[NA-1][0] - (2.0/3.0)*muWallsMU[NA-1][0]*0.50*(Divergence[NA-1][0] + Divergence[NA-2][0]))
								    )
								  + MESH.SupMP[NA-1][0][1]*(
								  		UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0])*(2.0*muWallsMU[NA][0]*GradU_DxMU[NA][0] - (2.0/3.0)*muWallsMU[NA][0]*0.50*(Divergence[NA-1][0]))
								  	  + UwallsMU[NA][0]*cos(MESH.AngleMU[NA][0])*muWallsMU[NA][0]*(GradU_DyMU[NA][0] + GradV_DxMU[NA][0])
								  	  + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0])*muWallsMU[NA][0]*(GradU_DyMU[NA][0] + GradV_DxMU[NA][0])
								  	  + VwallsMU[NA][0]*sin(MESH.AngleMU[NA][0])*(2.0*muWallsMU[NA][0]*GradV_DyMU[NA][0] - (2.0/3.0)*muWallsMU[NA][0]*0.50*(Divergence[NA-1][0]))
								    )
								  + MESH.SupMP[NA-1][0][2]*(
								  		UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0])*(2.0*muWallsMR[NA-1][0]*GradU_DxMR[NA-1][0] - (2.0/3.0)*muWallsMR[NA-1][0]*0.50*(Divergence[NA-1][0]))
								  	  + UwallsMR[NA-1][0]*cos(1.50*PI + MESH.AngleMR[NA-1][0])*muWallsMR[NA-1][0]*(GradU_DyMR[NA-1][0] + GradV_DxMR[NA-1][0])
								  	  + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0])*muWallsMR[NA-1][0]*(GradU_DyMR[NA-1][0] + GradV_DxMR[NA-1][0])
								  	  + VwallsMR[NA-1][0]*sin(1.50*PI + MESH.AngleMR[NA-1][0])*(2.0*muWallsMR[NA-1][0]*GradV_DyMR[NA-1][0] - (2.0/3.0)*muWallsMR[NA-1][0]*0.50*(Divergence[NA-1][0]))
								    )
								  + MESH.SupMP[NA-1][0][3]*(
								  		UwallsMR[NA-1][1]*cos(0.50*PI + MESH.AngleMR[NA-1][1])*(2.0*muWallsMR[NA-1][1]*GradU_DxMR[NA-1][1] - (2.0/3.0)*muWallsMR[i][1]*0.50*(Divergence[NA-1][0] + Divergence[NA-1][1]))
								  	  + UwallsMR[NA-1][1]*cos(0.50*PI + MESH.AngleMR[NA-1][1])*muWallsMR[NA-1][1]*(GradU_DyMR[NA-1][1] + GradV_DxMR[NA-1][1])
								  	  + VwallsMR[NA-1][1]*sin(0.50*PI + MESH.AngleMR[NA-1][1])*muWallsMR[NA-1][1]*(GradU_DyMR[NA-1][1] + GradV_DxMR[NA-1][1])
								  	  + VwallsMR[NA-1][1]*sin(0.50*PI + MESH.AngleMR[NA-1][1])*(2.0*muWallsMR[NA-1][1]*GradV_DyMR[NA-1][1] - (2.0/3.0)*muWallsMR[NA-1][1]*0.50*(Divergence[NA-1][0] + Divergence[NA-1][1]))
								    )
								);

	//Esquina arriba derecha
	EnergyViscous[NA-1][NRad-1] = (1.0/MESH.VolMP[NA-1][NRad-1])*(
								  	MESH.SupMP[NA-1][NRad-1][0]*(
								  		UwallsMU[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1])*(2.0*muWallsMU[NA-1][NRad-1]*GradU_DxMU[NA-1][NRad-1] - (2.0/3.0)*muWallsMU[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-2][NRad-1]))
								  	  + UwallsMU[NA-1][NRad-1]*cos(PI + MESH.AngleMU[NA-1][NRad-1])*muWallsMU[NA-1][NRad-1]*(GradU_DyMU[NA-1][NRad-1] + GradV_DxMU[NA-1][NRad-1])
								  	  + VwallsMU[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1])*muWallsMU[NA-1][NRad-1]*(GradU_DyMU[NA-1][NRad-1] + GradV_DxMU[NA-1][NRad-1])
								  	  + VwallsMU[NA-1][NRad-1]*sin(PI + MESH.AngleMU[NA-1][NRad-1])*(2.0*muWallsMU[NA-1][NRad-1]*GradV_DyMU[NA-1][NRad-1] - (2.0/3.0)*muWallsMU[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-2][NRad-1]))
								    )
								  + MESH.SupMP[NA-1][NRad-1][1]*(
								  		UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1])*(2.0*muWallsMU[NA][NRad-1]*GradU_DxMU[NA][NRad-1] - (2.0/3.0)*muWallsMU[NA][NRad-1]*0.50*(Divergence[NA-1][NRad-1]))
								  	  + UwallsMU[NA][NRad-1]*cos(MESH.AngleMU[NA][NRad-1])*muWallsMU[NA][NRad-1]*(GradU_DyMU[NA][NRad-1] + GradV_DxMU[NA][NRad-1])
								  	  + VwallsMU[NA][NRad-1]*sin(MESH.AngleMU[NA][NRad-1])*muWallsMU[NA][NRad-1]*(GradU_DyMU[NA][NRad-1] + GradV_DxMU[NA][NRad-1])
								  	  + VwallsMU[NA][NRad-1]*sin(MESH.AngleMU[NA][NRad-1])*(2.0*muWallsMU[NA][NRad-1]*GradV_DyMU[NA][NRad-1] - (2.0/3.0)*muWallsMU[NA][NRad-1]*0.50*(Divergence[NA-1][NRad-1]))
								    )
								  + MESH.SupMP[NA-1][j][2]*(
								  		UwallsMR[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*(2.0*muWallsMR[NA-1][NRad-1]*GradU_DxMR[NA-1][NRad-1] - (2.0/3.0)*muWallsMR[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-1][NRad-2]))
								  	  + UwallsMR[NA-1][NRad-1]*cos(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*muWallsMR[NA-1][NRad-1]*(GradU_DyMR[NA-1][NRad-1] + GradV_DxMR[NA-1][NRad-1])
								  	  + VwallsMR[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*muWallsMR[NA-1][NRad-1]*(GradU_DyMR[NA-1][NRad-1] + GradV_DxMR[NA-1][NRad-1])
								  	  + VwallsMR[NA-1][NRad-1]*sin(1.50*PI + MESH.AngleMR[NA-1][NRad-1])*(2.0*muWallsMR[NA-1][NRad-1]*GradV_DyMR[NA-1][NRad-1] - (2.0/3.0)*muWallsMR[NA-1][NRad-1]*0.50*(Divergence[NA-1][NRad-1] + Divergence[NA-1][NRad-2]))
								    )
								  + MESH.SupMP[NA-1][NRad-1][3]*(
								  		UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad])*(2.0*muWallsMR[NA-1][NRad]*GradU_DxMR[NA-1][NRad] - (2.0/3.0)*muWallsMR[NA-1][NRad]*0.50*(Divergence[NA-1][NRad-1]))
								  	  + UwallsMR[NA-1][NRad]*cos(0.50*PI + MESH.AngleMR[NA-1][NRad])*muWallsMR[NA-1][NRad]*(GradU_DyMR[NA-1][NRad] + GradV_DxMR[NA-1][NRad])
								  	  + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])*muWallsMR[NA-1][NRad]*(GradU_DyMR[NA-1][NRad] + GradV_DxMR[NA-1][NRad])
								  	  + VwallsMR[NA-1][NRad]*sin(0.50*PI + MESH.AngleMR[NA-1][NRad])*(2.0*muWallsMR[NA-1][NRad]*GradV_DyMR[NA-1][NRad] - (2.0/3.0)*muWallsMR[NA-1][NRad]*0.50*(Divergence[NA-1][NRad-1]))
								    )
								);
