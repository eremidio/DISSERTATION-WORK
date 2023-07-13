#PLOT 8 - Instruções para plotagem do gráfico coeficiente de spin Nerst x Temperatura
set encoding utf8
set title font "Arial, 12"
set title "Rede Union Jack"
set xrange [0:28]
set yrange [0:14]
set nokey
set nogrid
set xtics 4.0
set ytics 2.0
set xtics font "Arial, 10"
set ytics font "Arial, 10"
set xlabel font "Arial, 12"
set xlabel font "Arial, 12"
set xlabel "T(K)"
set ylabel "k_N"
plot "simulação7.txt" with lines lt rgb "red"

