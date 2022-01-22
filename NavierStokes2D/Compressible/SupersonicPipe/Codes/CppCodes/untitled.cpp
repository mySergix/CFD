void Solver::Get_Stresses(Mesher MESH){
int i, j;

	for(i = 0; i < NA+1; i++){
		for(j = 0; j < NRad; j++){
			TauRR_mu[i][j] = -muWallsMU[i][j]*(2.0*GradV_DyMU[i][j] - (2.0/3.0)*DivergenceMU[i][j]);

			TauRZ_mu[i][j] = -muWallsMU[i][j]*(GradV_DxMU[i][j] + GradU_DyMU[i][j]);

			TauZZ_mu[i][j] = -muWallsMU[i][j]*(2.0*GradU_DxMU[i][j] - (2.0/3.0)*DivergenceMU[i][j]);
		}
	}

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad+1; j++){
			TauRR_mr[i][j] = -muWallsMR[i][j]*(2.0*GradV_DyMR[i][j] - (2.0/3.0)*DivergenceMR[i][j]);

			TauRZ_mr[i][j] = -muWallsMR[i][j]*(GradV_DxMR[i][j] + GradU_DyMR[i][j]);

			TauZZ_mr[i][j] = -muWallsMR[i][j]*(2.0*GradU_DxMR[i][j] - (2.0/3.0)*DivergenceMR[i][j]);
		}
	}
}


for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumDifusiveU[i][j] = 
									  (1.0/MESH.VolMP[i][j])*(
								      +	MESH.SupMP[i][j][0]*(TauZZ_mu[i][j]*cos(PI + MESH.AngleMU[i][j]) + TauRZ_mu[i][j]*sin(PI + MESH.AngleMU[i][j]))
								      +	MESH.SupMP[i][j][1]*(TauZZ_mu[i+1][j]*cos(MESH.AngleMU[i+1][j]) + TauRZ_mu[i+1][j]*sin(MESH.AngleMU[i+1][j]))
								      +	MESH.SupMP[i][j][2]*(TauZZ_mr[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + TauRZ_mr[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
								      +	MESH.SupMP[i][j][3]*(TauZZ_mr[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + TauRZ_mr[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
									  )
									   (1.0/MESH.MP[i][j][1])*(muTotal[i][j]/ReynoldsM)*(0.50*(GradU_DyMU[i][j] + GradU_DyMU[i][j]) + 0.50*(GradV_DxMU[i][j] + GradV_DxMU[i][j]))
									  ;
		}
	}








	void Solver::Get_Stresses(Mesher MESH){
int i, j;

	for(i = 0; i < NA+1; i++){
		for(j = 0; j < NRad; j++){
			TauRR_mu[i][j] = muWallsMU[i][j]*(2.0*GradV_DyMU[i][j] - (2.0/3.0)*DivergenceMU[i][j]);

			TauRZ_mu[i][j] = muWallsMU[i][j]*(GradV_DxMU[i][j] + GradU_DyMU[i][j]);

			TauZZ_mu[i][j] = muWallsMU[i][j]*(2.0*GradU_DxMU[i][j] - (2.0/3.0)*DivergenceMU[i][j]);
		}
	}

	for(i = 0; i < NA; i++){
		for(j = 0; j < NRad+1; j++){
			TauRR_mr[i][j] = muWallsMR[i][j]*(2.0*GradV_DyMR[i][j] - (2.0/3.0)*DivergenceMR[i][j]);

			TauRZ_mr[i][j] = muWallsMR[i][j]*(GradV_DxMR[i][j] + GradU_DyMR[i][j]);

			TauZZ_mr[i][j] = muWallsMR[i][j]*(2.0*GradU_DxMR[i][j] - (2.0/3.0)*DivergenceMR[i][j]);
		}
	}
}

//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveU(Mesher MESH){
int i, j;

for(i = 0; i < NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumDifusiveU[i][j] = 
									  (1.0/MESH.VolMP[i][j])*(
								      +	MESH.SupMP[i][j][0]*(TauZZ_mu[i][j]*cos(PI + MESH.AngleMU[i][j]) + TauRZ_mu[i][j]*sin(PI + MESH.AngleMU[i][j]))
								      +	MESH.SupMP[i][j][1]*(TauZZ_mu[i+1][j]*cos(MESH.AngleMU[i+1][j]) + TauRZ_mu[i+1][j]*sin(MESH.AngleMU[i+1][j]))
								      +	MESH.SupMP[i][j][2]*(TauZZ_mr[i][j]*cos(1.50*PI + MESH.AngleMR[i][j]) + TauRZ_mr[i][j]*sin(1.50*PI + MESH.AngleMR[i][j]))
								      +	MESH.SupMP[i][j][3]*(TauZZ_mr[i][j+1]*cos(0.50*PI + MESH.AngleMR[i][j+1]) + TauRZ_mr[i][j+1]*sin(0.50*PI + MESH.AngleMR[i][j+1]))
									  )
									  ;
		}
	}

}



//Cálculo del término difusivo de la ecuación de cantidad de movimiento (Velocidad Radial V)
void Solver::Get_MomemtumDifusiveV(Mesher MESH){
int i, j;
	
	for(i = 0; i <  NA; i++){
		for(j = 0; j < NRad; j++){
			MomentumDifusiveV[i][j] = 
									+ (TauRR_mr[i][j+1] - TauRR_mr[i][j])/MESH.DeltasMP[i][j][1]
									+ (TauRZ_mu[i+1][j] - TauRZ_mu[i][j])/MESH.DeltasMP[i][j][0]
									+ ((2.0*muTotal[i][j])/MESH.MP[i][j][1])*((VwallsMR[i][j+1] - VwallsMR[i][j])/MESH.DeltasMP[i][j][1] - Vpres[i][j]/MESH.MP[i][j][1]) 
									;
		}
	}
}