void Solver::ComputeParameters(Mesher MESH, int Step){
int i, j;
double ReynoldsTau, MachTau;
int LimInf, LimSup;
double SumVel = 0.0;
double MeanVel, MeanRho;
double ReynoldsMean;

double RhoM = 0.0;

RhoW = 0.0;
TauW = 0.0;

LimInf = 0.80*NA;
LimSup = 0.95*NA;
	
	TiempoEstadistico += DeltaT;

	for(j = 0; j < NRad; j++){
		DensidadMedia[j] = 0.0;
		TemperaturaMedia[j] = 0.0;
		PresionMedia[j] = 0.0;
	}


	

	for(j = 0; j < NRad; j++){
		for(i = LimInf; i <= LimSup; i++){
			DensidadMedia[j] = DensidadMedia[j] + (RhoPres[i][j]*MESH.DeltasMP[i][j][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
			TemperaturaMedia[j] = TemperaturaMedia[j] + (Tpres[i][j]*MESH.DeltasMP[i][j][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
			PresionMedia[j] = PresionMedia[j] + (P[i][j]*MESH.DeltasMP[i][j][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
		}
	}

	for(j = 0; j < NRad; j++){
			DensidadMediaTemporal[j] = DensidadMediaTemporal[j] + (DensidadMedia[j]*DeltaT);
			TemperaturaMediaTemporal[j] = TemperaturaMediaTemporal[j] + (TemperaturaMedia[j]*DeltaT);
			PresionMediaTemporal[j] = PresionMediaTemporal[j] + (PresionMedia[j]*DeltaT);
	}

	
	if(Step%1000 == 0){

	for(i = LimInf; i <= LimSup; i++){
		RhoW = RhoW + (RhoPres[i][NRad-1]*MESH.DeltasMP[i][NRad-1][0])/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);
		TauW = TauW + (muTotal[i][NRad-1]*((Upres[i][NRad-1]/MESH.DeltasMR[i][NRad][1])*MESH.DeltasMP[i][NRad-1][0]))/(MESH.AxialCoord[LimSup+1] - MESH.AxialCoord[LimInf]);

	}

	for(j = 0; j < NRad; j++){
			SumVel = SumVel + (Upres[LimSup+1][j]*MESH.DeltasMP[LimSup+1][j][1])/VelocidadSonido;
	}

	for(j = 0; j < NRad; j++){
		RhoM = RhoM + (RhoPres[LimSup+1][j]*MESH.DeltasMP[LimSup+1][j][1]);
	}
	MeanRho = RhoM/(0.50*PipeDiameter);
	MeanVel = SumVel/(PipeDiameter/2.0);

	ReynoldsMean = (MeanRho*MeanVel*VelocidadSonido*0.50*PipeDiameter)/muW;
	uT = sqrt(abs(TauW)/RhoW);
	
	ReynoldsTau = (RhoW*uT*PipeDiameter)/(2.0*muW);
	MachTau = uT/sqrt(1.40*Rideal*Twall);


		for(j = 0; j < NRad; j++){
			DensidadMediaTemporalDato[j] = DensidadMediaTemporal[j]/TiempoEstadistico;
			TemperaturaMediaTemporalDato[j] = TemperaturaMediaTemporal[j]/TiempoEstadistico;
			PresionMediaTemporalDato[j] = PresionMediaTemporal[j]/TiempoEstadistico;
	}

	FILE *fp1;
		fp1 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/ResultadosParametros.txt","w");
			fprintf(fp1,"ReynoldsMean: %f \n", ReynoldsMean);
			fprintf(fp1,"ReynoldsTau: %f \n", ReynoldsTau);
			fprintf(fp1,"MachTau: %f \n", MachTau);
			fprintf(fp1,"Limite Inferior: %d \n", LimInf);
			fprintf(fp1,"Limite Superior: %d \n", LimSup);
			fprintf(fp1,"RhoW: %f \n", RhoW);
			fprintf(fp1,"Mean cross velocity: %f \n", MeanVel);
			fprintf(fp1,"TauW: %f \n", TauW);
			fprintf(fp1,"uT: %f \n", uT);
			fprintf(fp1,"TiempoEstadistico: %f \n", TiempoEstadistico);
					
	fclose(fp1);

	FILE *fp2;
		fp2 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/MeanDensity.txt","w");
			for(j = 0; j < NRad; j++){
				fprintf(fp2,"%f \t %f \n", MESH.MP[LimInf][NRad-1-j][1]/(PipeDiameter/2.0), DensidadMediaTemporalDato[j]/RhoW);
			}
	fclose(fp2);

	FILE *fp3;
		fp3 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/MeanTemperature.txt","w");
			for(j = 0; j < NRad; j++){
				fprintf(fp3,"%f \t %f \n", MESH.MP[LimInf][NRad-1-j][1]/(PipeDiameter/2.0), TemperaturaMediaTemporalDato[j]/Twall);
			}
	fclose(fp3);

	FILE *fp4;
		fp4 = fopen("/home/sergio/Desktop/TFG/RegimenTurbulento/MeshResults/Pipe/MeanPressure.txt","w");
			for(j = 0; j < NRad; j++){
				fprintf(fp4,"%f \t %f \n", MESH.MP[LimInf][NRad-1-j][1]/(PipeDiameter/2.0), PresionMediaTemporalDato[j]/PresionMediaTemporalDato[NRad-1]);
			}
	fclose(fp4);
	}
}