/*	if(Problema == 2){ //Problema Differentially Heated

			BoussinesqU = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
			BoussinesqV = M1.AllocateDouble(Fx - Ix + 2*Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
			BoussinesqW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ+1 + 2*HP, 1); 

			//Temperatura
			TLPAS = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TLPRES = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TLFUT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			//Matrices de las contribuciones de cada una de las ecuaciones
			ConvectiveT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			DiffusiveT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			TcontributionPast = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TcontributionPres = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

		}*/


        /*	if(Problema == 2){ //Problema Differentially Heated

			BoussinesqU = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
			BoussinesqV = M1.AllocateDouble(Fx - Ix + 2*Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
			BoussinesqW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ+1 + 2*HP, 1); 


			//Temperatura globales
			TGPRES = M1.AllocateDouble(NX + 2*HP,NY + 2*HP,NZ + 2*HP,1); //Temperatura T Global Presente
			TGFUT = M1.AllocateDouble(NX + 2*HP,NY + 2*HP,NZ + 2*HP,1); //Temperatura T Global Futuro

			//Temperatura
			TLPAS = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TLPRES = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TLFUT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			//Matrices de contribuciones
			ConvectiveT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			DiffusiveT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			TcontributionPast = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TcontributionPres = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

		}*/


        /*	if(Problema == 2){ //Problema Dfferentially Heated

			BoussinesqU = M1.AllocateDouble(Fx - Ix + 1 + 2*Halo, NY + 2*HP, NZ + 2*HP, 1);
			BoussinesqV = M1.AllocateDouble(Fx - Ix + 2*Halo, NY+1 + 2*HP, NZ + 2*HP, 1);
			BoussinesqW = M1.AllocateDouble(Fx - Ix + 2*Halo, NY + 2*HP, NZ+1 + 2*HP, 1); 

			//Temperatura
			TLPAS = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TLPRES = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TLFUT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			//Matrices de contribuciones
			ConvectiveT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			DiffusiveT = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);

			TcontributionPast = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			TcontributionPres = M1.AllocateDouble(Fx - Ix + Halo + HP, NY + 2*HP, NZ + 2*HP, 1);
			
		}*/

        /*	//Temperatura
	if(Problema == 2){ //Problema Differentially Heated

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){

					TLPAS[LPC(i,j,k,0)] = 0.0;
					TLPRES[LPC(i,j,k,0)] = 0.0;
					TLFUT[LPC(i,j,k,0)] = 0.0;

					//Matrices de contribuciones
					ConvectiveT[LPC(i,j,k,0)] = 0.0;

					DiffusiveT[LPC(i,j,k,0)] = 0.0;

					TcontributionPast[LPC(i,j,k,0)] = 0.0;
					TcontributionPres[LPC(i,j,k,0)] = 0.0;

				}
			}
		}
		
	}*/
}

/*
//Cálculo del término difusivo de la temperatura T
void Solver::Get_DiffusiveT(Mesher MESH, ParPro MPI1){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix; i < Fx; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}
		

	}
	else if(Rank == 0){

		//Centro
		for(i = Ix + 1; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix + 1; i < Fx; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}


		//Parte Izquierda
		i = 0;

		//Centro
		for(j = 1; j < NY - 1; j++){
			for(k = 1; k < NZ - 1; k++){
				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}
		}

		//Partes Superior e Inferior
		for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

		//Partes Here y There
		for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}


			//Esquina Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLEFT[PLEFT(0,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

	}

	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx - 1; i++){
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix; i < Fx - 1; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
		}

		//Parte Derecha
		i = NX - 1;

			//Centro
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}



			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte There
				k = NZ - 1;

				DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TLPRES[LPC(i,j+1,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TBOT[PBOT(i,0,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(TLPRES[LPC(i,j,k+1,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - There[PHERE(i,j,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveT[LPC(i,j,k,0)] = (K/(Rho*Cp*MESH.VolMP[GP(i,j,k,0)]))*(
											+ MESH.SupMP[GP(i,j,k,1)]*(TLPRES[LPC(i+1,j,k,0)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMP[GP(i,j,k,0)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMP[GP(i,j,k,3)]*(TTOP[PTOP(i,NY,k)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMP[GP(i,j,k,2)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMP[GP(i,j,k,5)]*(Tthere[PTHERE(i,j,NZ)] - TLPRES[LPC(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMP[GP(i,j,k,4)]*(TLPRES[LPC(i,j,k,0)] - TLPRES[LPC(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);


	}

}
*/

//Función con los diferentes esquemas convectivos utilizados
inline double Solver::ConvectiveScheme(double CoordObjetivo, double Velocity, double Coord1, double Phi1, double Coord2, double Phi2, double Coord3, double Phi3, double Coord4, double Phi4, string Esquema){

double PhiObjetivo;

double CoordD;
double PhiD;

double CoordC;
double PhiC;

double CoordU;
double PhiU;

	if (Velocity <= 0.0 || (Phi1 == 0.0 && Coord1 == 0.0)){

		CoordD = Coord2;
		PhiD = Phi2;
		CoordC = Coord3;
		PhiC = Phi3;
		CoordU = Coord4;
		PhiU = Phi4;

	}
	else if(Velocity > 0.0 || (Phi4 == 0.0 && Coord4 == 0.0)){

		CoordD = Coord3;
		PhiD = Phi3;
		CoordC = Coord2;
		PhiC = Phi2;
		CoordU = Coord1;
		PhiU = Phi1;

	}

	//Adimensionalizacion
	double PhiAdimC;
	double AdimCoordC;
	double AdimCoordE;

	PhiAdimC = (PhiC - PhiU)/(PhiD - PhiU);

	AdimCoordC = (CoordC - CoordU)/(CoordD - CoordU);

	AdimCoordE = (CoordObjetivo - CoordU)/(CoordD - CoordU);

	//Evaluacion
	double PhiF;

	if (PhiD == PhiU){
		PhiObjetivo = PhiD;
	}
	else{
		if(Esquema == "CDS"){
			PhiF = ((AdimCoordE - AdimCoordC)/(1.0 - AdimCoordC)) + ((AdimCoordE - 1.0)/(AdimCoordC - 1.0))*PhiAdimC;	
		}
		else if(Esquema == "UDS"){
			PhiF = PhiAdimC;	
		}
		else if(Esquema == "SUDS"){
			PhiF = (AdimCoordE/AdimCoordC)*PhiAdimC;
		}
		else if(Esquema == "QUICK"){
			PhiF = AdimCoordE + (((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0))))*(PhiAdimC - AdimCoordC);
		}
		else if(Esquema == "SMART"){
			if(PhiAdimC > 0 && PhiAdimC < AdimCoordC/3.0){
				PhiF = -((AdimCoordE*(1.0 - 3.0*AdimCoordC + 2.0*AdimCoordE))/(AdimCoordC*(AdimCoordC - 1.0)))*PhiAdimC;
			}
			else if(PhiAdimC > AdimCoordC/6.0 && PhiAdimC < (AdimCoordC/AdimCoordE)*(1.0 + AdimCoordE - AdimCoordC)){
				PhiF = ((AdimCoordE*(AdimCoordE - AdimCoordC))/(1.0 - AdimCoordC)) + ((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0)))*PhiAdimC;
			}
			else if(PhiAdimC > (AdimCoordC/AdimCoordE)*(1.0 + AdimCoordE - AdimCoordC) && PhiAdimC < 1.0){
				PhiF = 1.0;
			}
			else{
				PhiF = PhiAdimC;
			}
		}

		//Dimensionalizacion
		PhiObjetivo = PhiU + (PhiD - PhiU)*PhiF;
	}

	return PhiObjetivo;

}

/*
//Cálculo del término convectivo de la temperatura T
void Solver::Get_ConvectiveT(Mesher MESH, ParPro MPI1){
int i, j, k;
double Tw, Te, Ts, Tn, Th, Tt;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro 
		for(i = Ix; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){

					Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], UFIELD[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], UFIELD[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VFIELD[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WFIELD[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WFIELD[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VFIELD[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WFIELD[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );
				}
			}
		}

		for(i = Ix; i < Fx; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VFIELD[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VFIELD[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte Superior
				j = NY - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte There
				k = NZ - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

		}


	}
	else if(Rank == 0){

		//Centro 
		for(i = Ix + 1; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){

					Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );
				}
			}
		}

		for(i = Ix + 1; i < Fx; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte Superior
				j = NY - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte There
				k = NZ - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

		}

		//Parte Izquierda
		i = 0;

		//Centro
		for(j = 1; j < NY - 1; j++){
			for(k = 1; k < NZ - 1; k++){

					Tw = TLEFT[PLEFT(0,j,k)];
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );
			}
		}

		//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte Superior
				j = NY - 1;

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte There
				k = NZ - 1;

				Tw = TLEFT[PLEFT(0,j,k)];
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = TLEFT[PLEFT(0,j,k)];
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*Uleft[ULEFT(0,j,k)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[WTHERE(i,j,NZ)]
											  );


	}
	else if(Rank == Procesos - 1){

		//Centro 
		for(i = Ix; i < Fx - 1; i++){
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){

					Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
					Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

					Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
					Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

					Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
					Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

					ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );
				}
			}
		}

		for(i = Ix; i < Fx - 1; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte Superior
				j = NY - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte There
				k = NZ - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], ULPRES[LUC(i+1,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], MESH.MP[GP(i+2,j,k,0)], TLPRES[LPC(i+2,j,k,0)], EsquemaLargo);

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*ULPRES[LUC(i+1,j,k,0)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

		}


		//Parte Derecha
		i = NX - 1;

		//Centro 
		for(j = 1; j < NY - 1; j++){
			for(k = 1; k < NZ - 1; k++){

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );
			}
		}

		for(i = Ix; i < Fx - 1; i++){

			//Partes Superior e Inferior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = TBOT[PBOT(i,0,k)];
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte Superior
				j = NY - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = TTOP[PTOP(i,NY,k)];

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = There[PHERE(i,j,0)];
				Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

				//Parte There
				k = NZ - 1;

				Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
				Te = TRIGHT[PRIGHT(NX,j,k)];

				Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
				Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

				Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
				Tt = Tthere[PTHERE(i,j,NZ)];

				ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			}


			//Fila Abajo Here
			j = 0;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = TBOT[PBOT(i,0,k)];
			Tn = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], VLPRES[LVC(i,j+1,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], MESH.MP[GP(i,j+2,k,1)], TLPRES[LPC(i,j+2,k,0)], EsquemaLargo);

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*Vbot[VBOT(i,0,k)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*VLPRES[LVC(i,j+1,k,0)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = There[PHERE(i,j,0)];
			Tt = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], WLPRES[LWC(i,j,k+1,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,0)], MESH.MP[GP(i,j,k+2,2)], TLPRES[LPC(i,j,k+2,2)], EsquemaLargo);

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*Where[WHERE(i,j,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*WLPRES[LWC(i,j,k+1,0)]
											  );

			//Fla Arriba There
			j = NY - 1;
			k = NZ - 1;

			Tw = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
			Te = TRIGHT[PRIGHT(NX,j,k)];

			Ts = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
			Tn = TTOP[PTOP(i,NY,k)];

			Th = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
			Tt = Tthere[PTHERE(i,j,NZ)];

			ConvectiveT[LPC(i,j,k,0)] = (1.0/MESH.VolMP[GP(i,j,k,0)])*(
											  - MESH.SupMP[GP(i,j,k,0)]*Tw*ULPRES[LUC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,1)]*Te*Uright[URIGHT(NX,j,k)]
											  - MESH.SupMP[GP(i,j,k,2)]*Ts*VLPRES[LVC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,3)]*Tn*Vtop[VTOP(i,NY,k)]
											  - MESH.SupMP[GP(i,j,k,4)]*Th*WLPRES[LWC(i,j,k,0)]
											  + MESH.SupMP[GP(i,j,k,5)]*Tt*Wthere[PTHERE(i,j,NZ)]
											  );

		}


	}

}
*/
/*
//Cálculo de término de Boussinesq de la ecuación de Cantidad de Movimiento
void Solver::Get_BoussinesqTerm(Mesher MESH){
int i, j, k;
double Tu, Tv, Tw;

	if(gx != 0.0){
		//Velocidades U
		if(Rank != 0 && Rank != Procesos - 1){
	
			for(i = Ix; i < Fx; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){

						Tu = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
					
						BoussinesqU[LUC(i,j,k,0)] = gx*(1.0 - Beta*(Tu - To));

					}
				}
			}	

		}
		else if(Rank == 0){

			for(i = Ix + 1; i < Fx; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){

						Tu = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);
					
						BoussinesqU[LUC(i,j,k,0)] = gx*(1.0 - Beta*(Tu - To));

					}
				}
			}

		}
		else if(Rank == Procesos - 1){

			for(i = Ix; i < Fx; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){

						Tu = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], ULPRES[LUC(i,j,k,0)], MESH.MP[GP(i-2,j,k,0)], TLPRES[LPC(i-2,j,k,0)], MESH.MP[GP(i-1,j,k,0)], TLPRES[LPC(i-1,j,k,0)], MESH.MP[GP(i,j,k,0)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i+1,j,k,0)], TLPRES[LPC(i+1,j,k,0)], EsquemaLargo);

						BoussinesqU[LUC(i,j,k,0)] = gx*(1.0 - Beta*(Tu - To));

					}
				}
			}

		}

	}
	

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(j = 1; j < NY; j++){
			for(k = 0; k < NZ; k++){

				Tv = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LVC(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LPC(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LPC(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LPC(i,j+1,k,0)], EsquemaLargo);
					
				BoussinesqV[LVC(i,j,k,0)] = gy*(1.0 - Beta*(Tv - To));

			}
		}
	}

	if(gz != 0.0){

		//Velocidades W
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 1; k < NZ; k++){

					Tw = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], WLPRES[LWC(i,j,k,0)], MESH.MP[GP(i,j,k-2,2)], TLPRES[LPC(i,j,k-2,0)], MESH.MP[GP(i,j,k-1,2)], TLPRES[LPC(i,j,k-1,0)], MESH.MP[GP(i,j,k,2)], TLPRES[LPC(i,j,k,0)], MESH.MP[GP(i,j,k+1,2)], TLPRES[LPC(i,j,k+1,2)], EsquemaLargo);
					
					BoussinesqW[LWC(i,j,k,0)] = gz*(1.0 - Beta*(Tw - To));

				}
			}
		}

	}
	
}
*/

	/*	else if(Problema == 2){

			for(i = Ix; i < Fx + 1; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){
						UcontributionPres[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)] + BoussinesqU[LUC(i,j,k,0)];
					}
				}
			}

		}*/
		
/*	else if(Problema == 2){

			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){

					//Parte Izquierda
					UcontributionPres[LUC(0,j,k,0)] = 0.0;

					for(i = Ix + 1; i < Fx + 1; i++){
						UcontributionPres[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)] + BoussinesqU[LUC(i,j,k,0)];
					}
				}
			}

		}*/
		
        /*	else if(Problema == 2){

			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){

					//Parte Derecha
					UcontributionPres[LUC(NX,j,k,0)] = 0.0;

					for(i = Ix; i < Fx; i++){
						UcontributionPres[LUC(i,j,k,0)] = - ConvectiveU[LUC(i,j,k,0)] + DiffusiveU[LUC(i,j,k,0)] + BoussinesqU[LUC(i,j,k,0)];
					}
				}
			}

		}*/
		
        /*	else if(Problema == 2){

		//Velocidades V
		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){

				//Parte Inferior
				VcontributionPres[LVC(i,0,k,0)] = 0.0;

				//Parte Superior
				VcontributionPres[LVC(i,NY,k,0)] = 0.0;

				for(j = 1; j < NY; j++){
					VcontributionPres[LVC(i,j,k,0)] = - ConvectiveV[LVC(i,j,k,0)] + DiffusiveV[LVC(i,j,k,0)] + BoussinesqV[LVC(i,j,k,0)];
				}
			}
		}

		//Velocidades W
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){

				//Parte Here
				WcontributionPres[LWC(i,j,0,0)] = 0.0;

				//Parte There
				WcontributionPres[LWC(i,j,NZ,0)] = 0.0;

				for(k = 1; k < NZ; k++){
					WcontributionPres[LWC(i,j,k,0)] = - ConvectiveW[LWC(i,j,k,0)] + DiffusiveW[LWC(i,j,k,0)] + BoussinesqW[LWC(i,j,k,0)];
				}
			}
		}

	}*/


	/*
//Cálculo de las contribuciones de Temperatura
void Solver::Get_TemperatureContributions(){
int i, j, k;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				TcontributionPres[LPC(i,j,k,0)] = - ConvectiveT[LPC(i,j,k,0)] + DiffusiveT[LPC(i,j,k,0)];
			}
		}
	}

}
*/


/*
//Cálculo máxima diferecia Gauss Seidel
void Solver::Get_MaxDifGS(ParPro MPI1){
int i, j, k;
MaxDiffGS = 0.0;
MPI_Status ST;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				if(abs(PLFUT[LPC(i,j,k,0)] - PLFUTsup[LPC(i,j,k,0)]) >= MaxDiffGS){
					MaxDiffGS = abs(PLFUT[LPC(i,j,k,0)] - PLFUTsup[LPC(i,j,k,0)]);	
				}
			}
			
		}
	}

	MPI1.SendDataToZero(MaxDiffGS, PDiff);

	double Diff;
	if(Rank == 0){
		Diff = PDiff[0];
		for(i = 1; i < Procesos; i++){
			if(PDiff[i] >= Diff){
				Diff = PDiff[i];
			}
		}
	}

	MPI1.SendDataToAll(Diff, MaxDiffGS);

}
*/

/*
//Resolución del Sistema con Gauss-Seidel
void Solver::Get_Pressure(ParPro MPI1){
int i, j, k;
double MaxDiff;
MaxDiffGS = 2.0*ConvergenciaGS;

	while(MaxDiffGS > ConvergenciaGS){

		MPI1.CommunicateDataLP(PLFUT, PLFUT, Ix, Fx); //Enviar Presiones P

		//Parte Central
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)]  + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}
		}
		// 
		Get_MaxDifGS(MPI1);

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					PLFUTsup[LPC(i,j,k,0)] = PLFUT[LPC(i,j,k,0)];
				}
			}
		}

		
	}
	
}
*/

/*
//Cálculo de las temperaturas del sistema
void Solver::Get_Temperatures(){
int i, j, k;

	for(i = Ix; i < Fx; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				TLFUT[LPC(i,j,k,0)] = TLPRES[LPC(i,j,k,0)] + DeltaT*(1.50*TcontributionPres[LPC(i,j,k,0)] - 0.50*TcontributionPast[LPC(i,j,k,0)]);
			}
		}
	}

}
*/

/*
		if(Problema == 2){

			//Comprobación Temperatura T
			for(i = 0; i < NX; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){
						if(abs((TGFUT[GP(i,j,k,0)] - TGPRES[GP(i,j,k,0)])/(TGPRES[GP(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal){
							MaxDiffGlobal = abs((TGFUT[GP(i,j,k,0)] - TGPRES[GP(i,j,k,0)])/(TGPRES[GP(i,j,k,0)] + 1e-10));
						}
					}
				}
			}

		}*/


		/*	if(Problema == 2){

		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					TLPAS[LPC(i,j,k,0)] = TLPRES[LPC(i,j,k,0)];
					TLPRES[LPC(i,j,k,0)] = TLFUT[LPC(i,j,k,0)];

					TcontributionPast[LPC(i,j,k,0)] = TcontributionPres[LPC(i,j,k,0)];
				}
			}
		}

	}*/

/*
		if(Problema == 2){

			MPI1.CommunicateDataLP(TLPRES, TLPRES, Ix, Fx); //Enviar Temperatura T
			Get_DiffusiveT(MESH, MPI1); //Cálculo Difusivo Temperatura T
			Get_ConvectiveT(MESH, MPI1); //Cálculo Convectivo T

		}
		*/
	
		
		//if(Problema == 2){
		//	Get_BoussinesqTerm(MESH); //Cálculo del Término de Boussinesq
		//}
		
	

		//if(Problema == 2){
		//	Get_TemperatureContributions(); //Cálculo contribuciones Temperatura
		//}
		

	//if(Problema == 2){
		//	Get_Temperatures();
		//}
		