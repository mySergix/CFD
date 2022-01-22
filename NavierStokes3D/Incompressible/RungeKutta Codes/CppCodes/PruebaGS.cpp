if(Rank != 0 && Rank != Procesos - 1){

			//Parte Central
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}
			}

			for(i = Ix; i < Fx; i++){
				
				//Partes Inferior y Superior
				for(k = 1; k < NZ - 1; k++){

					//Parte Inferior
					j = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
					//Parte Superior
					j = NY - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/5.0;
				
				}

				//Partes Here y There
				for(j = 1; j < NY - 1; j++){

					//Parte Here
					k = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

					//Parte There
					k = NZ - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				}

				//Esquinas

				//Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

			}

		}
        else if(Rank == 0){

			//Parte Central
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}
			}

			for(i = Ix + 1; i < Fx; i++){
				
				//Partes Inferior y Superior
				for(k = 1; k < NZ - 1; k++){

					//Parte Inferior
					j = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
					//Parte Superior
					j = NY - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				}

				//Partes Here y There
				for(j = 1; j < NY - 1; j++){

					//Parte Here
					k = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

					//Parte There
					k = NZ - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				}

				//Esquinas

				//Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

			}

			//Parte Izquierda
			i = 0;
				//Parte Central
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}

				//Partes Inferior y Superior
				for(k = 1; k < NZ - 1; k++){

					//Parte Inferior
					j = 0;
					PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
					//Parte Superior
					j = NY - 1;
					PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				}

				//Partes Here y There
				for(j = 1; j < NY - 1; j++){

					//Parte Here
					k = 0;
					PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

					//Parte There
					k = NZ - 1;
					PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				}

				//Esquinas

				//Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

		}
		else if(Rank == Procesos - 1){

			//Parte Central
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}
			}

			for(i = Ix; i < Fx - 1; i++){
				
				//Partes Inferior y Superior
				for(k = 1; k < NZ - 1; k++){

					//Parte Inferior
					j = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
					//Parte Superior
					j = NY - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				}

				//Partes Here y There
				for(j = 1; j < NY - 1; j++){

					//Parte Here
					k = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

					//Parte There
					k = NZ - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				}

				//Esquinas

				//Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LPC(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

			}

			//Parte Derecha
			i = NX - 1;

			//Parte Central
			for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}

			//Partes Inferior y Superior
				for(k = 1; k < NZ - 1; k++){

					//Parte Inferior
					j = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
					//Parte Superior
					j = NY - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				}

				//Partes Here y There
				for(j = 1; j < NY - 1; j++){

					//Parte Here
					k = 0;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

					//Parte There
					k = NZ - 1;
					PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				}

				//Esquinas

				//Abajo Here
				j = 0;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LPC(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + at[LA(i,j,k,0)]*PLFUT[LPC(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

				//Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LPC(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LPC(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LPC(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LPC(i,j,k-1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];

		}
