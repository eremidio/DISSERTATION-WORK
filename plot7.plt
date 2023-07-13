#PLOT 7 - Instruções para plotagem do gráfico condutividade Hall de spin x Temperatura
set encoding utf8
set title font "Arial, 12"
set title "Rede Union Jack"
set xrange [0:28]
set yrange [0:2.4]
set nokey
set nogrid
set xtics 4.0
set ytics 0.4
set xtics font "Arial, 10"
set ytics font "Arial, 10"
set xlabel font "Arial, 12"
set xlabel font "Arial, 12"
set xlabel "T(K)"
set ylabel "{/Symbol s}^{xy}"
plot "simulação6.txt" with lines lt rgb "blue"

