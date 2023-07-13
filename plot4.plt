#PLOT 4 - Instruções para plotagem do grafíco do hamilitoniano para a rede Union Jack autovalor 2
set encoding utf8
set title font "Arial, 12"
set title "Rede Union Jack - Autovalor 2"
set xrange [0:pi]
set yrange [0:pi]
set zrange [0.5:1.5]
set dgrid3d
set nokey
set xtics 0.5
set ytics 0.5
set ztics 0.5
set xtics font "Arial, 10"
set ytics font "Arial, 10"
set ztics font "Arial, 10"
set xlabel font "Arial, 12"
set xlabel font "Arial, 12"
set xtics nomirror rotate by -45
set ytics nomirror rotate by 45
set xlabel "k_x"
set ylabel "k_y"
set surface
set pm3d
set cntrparam levels 10
set palette defined (0 "green", 1 "blue", 2 "purple")
set hidden3d
set isosamples 80
#Constantes usadas no cálculo
j1= -0.2
j2= -0.2
D=0.2
S= sqrt(3)/2
#
splot (4.0*j2*S*cos(x)*cos(y)- 8.0*j1*S-4.0*j2*S)-sqrt((2.0*j1*S*(cos(x)+cos(y)))**2+((D/2.0)*S*(cos(x)+cos(y)))**2+(4.0*j2*(cos(x)*cos(y)- 1.0))**2)
