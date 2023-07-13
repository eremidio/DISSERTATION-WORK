#PLOT 10 - Instruções para plotagem do grafíco do hamilitoniano para a rede de chevron checando se há gap para um determiando #valor das constantes de acoplamento 
set encoding utf8
set title font "Arial, 12"
set title "Rede Union Jack - Gap"
set xrange [0:pi]
set yrange [0:pi]
set zrange [-5:10]
set dgrid3d
set nokey
set xtics 0.5
set ytics 0.5
set ztics 2.0
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
j1= -0.4
j2= -0.4
D=0.5
S= sqrt(3)/2
#
splot 2*sqrt((2.0*j1*S*(cos(x)+cos(y)))**2+((D/2.0)*S*(cos(x)+cos(y)))**2+(4.0*j2*(cos(x)*cos(y)- 1.0))**2)
