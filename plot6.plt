#PLOT 6 - Instruções para plotagem do grafíco da curvatura de Berry para a rede Union Jack
set encoding utf8
set title font "Arial, 12"
set title "Rede Union Jack - curvatura de Berry"
set xrange [0:pi]
set yrange [0:pi]
set zrange [-1:2]
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
set palette defined (0 "red", 1 "orange", 2 "yellow")
set hidden3d
set isosamples 80
#Constantes usadas no cálculo
j1= -0.2
j2= -0.2
D = 0.5
S= sqrt(3)/2
#
splot 0,5*((((D/2.0) + (2.0*j1))**2*(cos(x)+cos(y))**2)+(16.0*j2**2*(cos(x)*cos(y)-1)**2))**(-1,5)*((4.0*j1*j2*D*(cos(x)+cos(y)))*(cos(y)*(sin(x))**2 + cos(x)*(sin(y))**2-sin(2.0*x)*sin(y)))
