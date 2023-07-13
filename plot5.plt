#PLOT 5 - Instruções para plotagem do grafíco da curvatura de Berry para a rede de chevron
set encoding utf8
set title font "Arial, 12"
set title "Rede de Chevron - curvatura de Berry"
set xrange [0:pi]
set yrange [0:pi]
set zrange [-1:2.5]
set dgrid3d
set nokey
set xtics 0.5
set ytics 0.5
set ztics 2
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
j1= -0.4
j2= -0.3
j3= -0.3
D=0.4
r=0.7
S= sqrt(3)/2
#
splot (1/16.0)*((D**2*(sin(x)+sin(y))**2)+(2.0*j2*D*cos(y)*(sin(x)+sin(y)))+(j2**2*cos(y)**2)+(4.0*j3**2*cos(y)**2*sin(x)**2)+(j1**2*cos(x)**2*(r-1.0)**2))**(-1.5)*((16.0*j1*j3*D*(sin(x)+sin(y))*sin(x)**2*sin(y)*(1.0-r))+(16.0*j1*j3*D*cos(y)**2*sin(x)**2*(1.0-r))+(8.0*j1*j3*D*cos(x)*sin(y)*sin(2.0*x))+( 16.0*j1*j3*D*cos(x)**2*cos(y)**2*(1.0-r))+(8.0*j1*j2*j3*cos(x)**2*sin(2.0*y)))
