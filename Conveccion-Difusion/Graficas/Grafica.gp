
set title "Ecuación de convección difusión, Problema 4" font "Arial-Bold,30"
set xlabel "Coordenada X" font "Arial-Bold,20" 
set ylabel "Coordenada Y" font "Arial-Bold,20"
set zlabel "Valor de la propiedad" font "Arial-Bold,20"
set pm3d
set hidden3d
set dgrid3d 
splot 'Datos.txt' with pm3d
