//Cálculo del término 4 de la ecuación de Spalart-Allmaras
void Solver::SA_Term4(Mesher MESH){
int i, j;

	//Centro
	for(i = 1; i < NA-1; i++){
		for(j = 1; j < NRad-1; j++){
			SA_Termino4[i][j] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[i][j])*(
									 MESH.SupMP[i][j][0]*0.50*(vpresent[i][j] + vpresent[i-1][j]) 
								   + MESH.SupMP[i][j][1]*0.50*(vpresent[i][j] + vpresent[i+1][j])
								   ),2.0)
							  + pow((1.0/MESH.VolMP[i][j])*(
							  		 MESH.SupMP[i][j][2]*0.50*(vpresent[i][j] + vpresent[i][j-1])
							  	   + MESH.SupMP[i][j][3]*0.50*(vpresent[i][j] + vpresent[i][j+1])
							  	   ),2.0)
							  );
		}
	}

	for(j = 1; j < NRad-1; j++){
		//Parte izquierda
		SA_Termino4[0][j] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[0][j])*(
									 MESH.SupMP[0][j][0]*vleft[j] 
								   + MESH.SupMP[0][j][1]*0.50*(vpresent[0][j] + vpresent[1][j])
								   ),2.0)
							  + pow((1.0/MESH.VolMP[0][j])*(
							  		 MESH.SupMP[0][j][2]*0.50*(vpresent[0][j] + vpresent[0][j-1])
							  	   + MESH.SupMP[0][j][3]*0.50*(vpresent[0][j] + vpresent[0][j+1])
							  	   ),2.0)
							  );

		//Parte derecha
			SA_Termino4[NA-1][j] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[NA-1][j])*(
									 MESH.SupMP[NA-1][j][0]*0.50*(vpresent[NA-1][j] + vpresent[NA-2][j]) 
								   + MESH.SupMP[NA-1][j][1]*vright[j]
								   ),2.0)
							  + pow((1.0/MESH.VolMP[NA-1][j])*(
							  		 MESH.SupMP[NA-1][j][2]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j-1])
							  	   + MESH.SupMP[NA-1][j][3]*0.50*(vpresent[NA-1][j] + vpresent[NA-1][j+1])
							  	   ),2.0)
							  );
	}

	for(i = 1; i < NA-1; i++){
		//Parte abajo
		SA_Termino4[i][0] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[i][0])*(
									 MESH.SupMP[i][0][0]*0.50*(vpresent[i][0] + vpresent[i-1][0]) 
								   + MESH.SupMP[i][0][1]*0.50*(vpresent[i][0] + vpresent[i+1][0])
								   ),2.0)
							  + pow((1.0/MESH.VolMP[i][0])*(
							  		 MESH.SupMP[i][0][2]*vdown[i]
							  	   + MESH.SupMP[i][0][3]*0.50*(vpresent[i][0] + vpresent[i][1])
							  	   ),2.0)
							  );

		//Parte arriba
		SA_Termino4[i][NRad-1] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[i][NRad-1])*(
									 MESH.SupMP[i][NRad-1][0]*0.50*(vpresent[i][NRad-1] + vpresent[i-1][NRad-1]) 
								   + MESH.SupMP[i][NRad-1][1]*0.50*(vpresent[i][NRad-1] + vpresent[i+1][NRad-1])
								   ),2.0)
							  + pow((1.0/MESH.VolMP[i][NRad-1])*(
							  		 MESH.SupMP[i][NRad-1][2]*0.50*(vpresent[i][NRad-1] + vpresent[i][NRad-2])
							  	   + MESH.SupMP[i][NRad-1][3]*vup[i]
							  	   ),2.0)
							  );
	}

	//Esquina abajo izquierda
	SA_Termino4[0][0] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[0][0])*(
									 MESH.SupMP[0][0][0]*vleft[0]
								   + MESH.SupMP[0][0][1]*0.50*(vpresent[0][0] + vpresent[1][0])
								   ),2.0)
							  + pow((1.0/MESH.VolMP[0][0])*(
							  		 MESH.SupMP[0][0][2]*vdown[0]
							  	   + MESH.SupMP[0][0][3]*0.50*(vpresent[0][0] + vpresent[0][1])
							  	   ),2.0)
							  );

	//Esquina arriba izquierda
	SA_Termino4[0][NRad-1] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[0][NRad-1])*(
									 MESH.SupMP[0][NRad-1][0]*vleft[NRad-1]
								   + MESH.SupMP[0][NRad-1][1]*0.50*(vpresent[0][NRad-1] + vpresent[1][NRad-1])
								   ),2.0)
							  + pow((1.0/MESH.VolMP[0][NRad-1])*(
							  		 MESH.SupMP[0][NRad-1][2]*0.50*(vpresent[0][NRad-1] + vpresent[0][NRad-2])
							  	   + MESH.SupMP[0][NRad-1][3]*vup[0]
							  	   ),2.0)
							  );

	//Esquina abajo derecha
	SA_Termino4[NA-1][0] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[NA-1][0])*(
									 MESH.SupMP[NA-1][0][0]*0.50*(vpresent[NA-1][0] + vpresent[NA-2][0]) 
								   + MESH.SupMP[NA-1][0][1]*vright[0]
								   ),2.0)
							  + pow((1.0/MESH.VolMP[NA-1][0])*(
							  		 MESH.SupMP[NA-1][0][2]*vdown[NA-1]
							  	   + MESH.SupMP[NA-1][0][3]*0.50*(vpresent[NA-1][0] + vpresent[NA-1][1])
							  	   ),2.0)
							  );

	//Esquina arriba derecha
	SA_Termino4[NA-1][NRad-1] = (Cb2/Sigma)*(
								pow((1.0/MESH.VolMP[NA-1][j])*(
									 MESH.SupMP[NA-1][j][0]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-2][NRad-1]) 
								   + MESH.SupMP[NA-1][j][1]*vright[NRad-1]
								   ),2.0)
							  + pow((1.0/MESH.VolMP[NA-1][j])*(
							  		 MESH.SupMP[NA-1][j][2]*0.50*(vpresent[NA-1][NRad-1] + vpresent[NA-1][NRad-2])
							  	   + MESH.SupMP[NA-1][j][3]*vup[NA-1]
							  	   ),2.0)
							  );
}