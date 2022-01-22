#!/bin/bash

Directorio1=/home/sergiogus/Desktop/CTTC/BashScripting
#echo $Directorio1
cd $Directorio1

g++ EscribirRE.cpp -o Ejecutable

input1=$Directorio1/NumerosReynolds.txt

while IFS= read -r line
do
	Carpeta=$(printf 'RE_%d' "$line")
	mkdir $Directorio1/$Carpeta
	
	sed -e 's/RE_NUM/'$line'/g' $Directorio1/ProblemDataBase.txt > $Directorio1/ProblemData.txt
	$Directorio1/Ejecutable
  echo "$line"
done < "$input1"
