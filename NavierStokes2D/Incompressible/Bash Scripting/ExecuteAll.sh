#!/bin/bash

Directorio1=/home/sergiogus/Desktop/CTTC/NavierStokes

cd $Directorio1/Codigos

InputData=$Directorio1/DatosIniciales/ProblemasTotal.txt

Counter1=1;

while IFS= read -r Datos
do
	if [ "$Datos" == "Problema Driven Cavity" ]
	then
		N1=$Counter1
	elif [ "$Datos" == "Problema Differentially Heated" ]
	then 
		N2=$Counter1
	elif [ "$Datos" == "Problema Square Cylinder" ]
	then
		N3=$Counter1
	elif [ "$Datos" == "FINAL" ]
	then 
		N4=$Counter1
	fi

	((Counter1++))
done < "$InputData"

CasosDriven=$(($N2-$N1-2))
CasosDifferentially=$(($N3-$N2-2))
CasosSquare=$(($N4-$N3-2))

DrivenCavity=1
DifferentiallyHeated=2
SquareCylinder=3

#make clean
#make

Counter2=1
while IFS= read -r Datos
do
	if [ $Counter2 -gt $N1 ] && [ $Counter2 -lt $(($N2-1)) ]
	then	
		Carpeta=$(printf 'RE_%d' $Datos)
		mkdir $Directorio1/ResultadosScripting/DrivenCavity/$Carpeta
		sed -e 's/RE_NUM/'"$Datos"'/g' $Directorio1/DatosIniciales/ProblemDataBase.txt > $Directorio1/DatosIniciales/ProblemData.txt

		sed -i 's/NUM_PRO/'"$DrivenCavity"'/g' $Directorio1/DatosIniciales/ProblemData.txt
		sed -i 's/RA_NUM/'"$DrivenCavity"'/g' $Directorio1/DatosIniciales/ProblemData.txt
	
		#./Main
	elif [ $Counter2 -gt $N2 ]  && [ $Counter2 -lt $(($N3-1)) ]
	then
		N=`echo "l($Datos)/l(10)" | bc -l`
		Carpeta=$(printf 'RA_1E%d' $N)
		mkdir $Directorio1/ResultadosScripting/DifferentiallyHeated/$Carpeta
		sed -e 's/RA_NUM/'"$Datos"'/g' $Directorio1/DatosIniciales/ProblemDataBase.txt > $Directorio1/DatosIniciales/ProblemData.txt

		sed -i 's/NUM_PRO/'"$DifferentiallyHeated"'/g' $Directorio1/DatosIniciales/ProblemData.txt
		sed -i 's/RE_NUM/'"$DifferentiallyHeated"'/g' $Directorio1/DatosIniciales/ProblemData.txt
		
		#./Main

	elif [ $Counter2 -gt $N3 ]  && [ $Counter2 -lt $(($N4-1)) ]
	then	
		Carpeta=$(printf 'RE_%d'  $Datos)
		mkdir $Directorio1/ResultadosScripting/SquareCylinder/$Carpeta
		sed -e 's/RE_NUM/'"$Datos"'/g' $Directorio1/DatosIniciales/ProblemDataBase.txt > $Directorio1/DatosIniciales/ProblemData.txt

		sed -i 's/NUM_PRO/'"$SquareCylinder"'/g' $Directorio1/DatosIniciales/ProblemData.txt
		sed -i 's/RA_NUM/'"$SquareCylinder"'/g' $Directorio1/DatosIniciales/ProblemData.txt
		
		#./Main
	fi
	((Counter2++))
done < "$InputData"



