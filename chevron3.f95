program chevron3
use dissertacao
implicit none

!Variaveis
real(kind=kind(1.0d0))::autovalor1, autovalor2, T1, berry, soma1, soma2, hallspin, halltermica
real(kind=kind(1.0d0)), dimension(1000000)::kx, ky, x1, x2, x3, x4, x5, x6
integer::i, fu, fi, ios

!Procedimentos externos
interface
function autovalor1(kx, ky)
real(kind=kind(1.0d0))::energia
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S= sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.4d0
real(kind=kind(1.0d0)), parameter::j1=-0.4d0
real(kind=kind(1.0d0)), parameter::j2=-0.3d0
real(kind=kind(1.0d0)), parameter::j3=-0.3d0
real(kind=kind(1.0d0)), parameter::r=0.5d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
end function autovalor1

function autovalor2(kx, ky) 
real(kind=kind(1.0d0))::energia
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.4d0
real(kind=kind(1.0d0)), parameter::j1=-0.4d0
real(kind=kind(1.0d0)), parameter::j2=-0.3d0
real(kind=kind(1.0d0)), parameter::j3=-0.3d0
real(kind=kind(1.0d0)), parameter::r=0.5d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
end function autovalor2

function berry(kx, ky)
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.4d0
real(kind=kind(1.0d0)), parameter::j1=-0.4d0
real(kind=kind(1.0d0)), parameter::j2=-0.3d0
real(kind=kind(1.0d0)), parameter::j3=-0.3d0
real(kind=kind(1.0d0)), parameter::r=0.5d0
real(kind=kind(1.0d0))::p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
end function berry
end interface




!Definição


!1. Lendo dados de arquivos e abrindo arquivo para exportar dados se necessário
open(newunit=fu, file='meshredequadrada.txt', status='old', action='read', iostat=ios)
open(newunit=fi, file='chevron3.txt', status='replace', action='readwrite', iostat=ios)

!2.Entrando dados no programa manualmente (temperatura e número de pontos na rede quadrada)
T1=25.0d0

!3.Computando as funções correspodentes na malha quadrada

!Determinando os valores de kx e ky para calculo da energia da banda na malha quadrada
do i=1, 1000000, 1
read(unit=fu, fmt=*) kx(i), ky(i)
end do

!Calculando o autovalores do Hamiltoniano
do i=1, 1000000, 1
x1(i)= autovalor1(kx(i), ky(i))
end do

do i=1, 1000000, 1
x2(i)= autovalor2(kx(i), ky(i))
end do

!Cálculo do número de estados de ocupação
do i=1, 1000000, 1
if (x1(i)/=0.0d0) then
x3(i)= boseeinstein(T1, x1(i))
else  
x3(i)=0.0d0
endif
end do

do i=1, 1000000, 1
if (x2(i)/=0.0d0) then
x4(i)= boseeinstein(T1, x2(i))
else 
x4(i)=0.0d0
endif
end do

!Cálculo das funções características que aparecem nas integrais dos coeficientes de transporte
do i=1, 1000000, 1
x5(i)= constant1(x3(i))
end do

do i=1, 1000000, 1
x6(i)= constant1(x4(i))
end do

!4.Cálculo da integral que dá a coeficiente de Spin Nerst e a condutividade Hall térmica
!4.1 Coeficiente Nerst de spin
soma1=(((x5(1)+x6(1))*berry(kx(1),ky(1))+(x5(1000000)+x6(1000000))*berry(kx(1000000),ky(1000000)))/2.0d0)*(pi**2/1000000d0)
do i=2, 999999, 1
soma1=soma1+(x5(i)+x6(i))*berry(kx(i),ky(i))*(pi**2/1000000d0)
end do
hallspin=soma1*T1/(pi**2)
write(*,*) "O valor do coeficiente Nerst de spin usando uma malha quadrada de 1000x1000 pontos a", T1, "K é", hallspin

!5.Exportando os dados para um arquivo externo
do i=1, 1000000, 1
write(unit=fi, fmt=*) kx(i), ky(i), autovalor1(kx(i), ky(i)), autovalor2(kx(i), ky(i)), x3(i), x4(i), x5(i), x6(i)
end do

close(unit=fu, status="keep", iostat=ios)
close(unit=fu, status="keep", iostat=ios)

end program chevron3



function autovalor1(kx, ky) result(energia)
use dissertacao
!Esta função nos dá os autovalores do hamiltoniano em termos do parametro kx e ky no espaço recíproco e deve ser Calculado analiticamente para cada caso
!Por esta razão é recomendável usá-la como um procedimento externo para cada  rede cristalina em particular
 implicit none
!Resultado da função
real(kind=kind(1.0d0))::energia
!Variáveis
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.4d0
real(kind=kind(1.0d0)), parameter::j1=-0.4d0
real(kind=kind(1.0d0)), parameter::j2=-0.3d0
real(kind=kind(1.0d0)), parameter::j3=-0.3d0
real(kind=kind(1.0d0)), parameter::r=0.5d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
!Definição
h0=j1*cos(kx)*(r+1.0d0)- 2.0d0*j1*S - 2.0d0*j2*S- 4.0d0*j3*S
hx=2.0d0*D*S*(sin(kx)+sin(ky)) + 2.0d0*j2*S*cos(ky)
hy= 4.0*j3*S*cos(ky)*sin(kx)
hz=2.0d0*j1*S*cos(kx)*(r-1.0d0)
h=sqrt(hx*hx+hy*hy*hz*hz)
energia=h0+h
end function autovalor1



function autovalor2(kx, ky) result(energia)

use dissertacao
!Esta função nos dá os autovalores do hamiltoniano em termos do parametro kx e ky no espaço recíproco e deve ser Calculado analiticamente para cada caso
!Por esta razão é recomendável usá-la como um procedimento externo para cada  rede cristalina em particular
 implicit none
!Resultado da função
real(kind=kind(1.0d0))::energia
!Variáveis
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.4d0
real(kind=kind(1.0d0)), parameter::j1=-0.4d0
real(kind=kind(1.0d0)), parameter::j2=-0.3d0
real(kind=kind(1.0d0)), parameter::j3=-0.3d0
real(kind=kind(1.0d0)), parameter::r=0.5d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
!Definição
h0=j1*cos(kx)*(r+1.0d0)- 2.0d0*j1*S - 2.0d0*j2*S- 4.0d0*j3*S
hx=2.0d0*D*S*(sin(kx)+sin(ky)) + 2.0d0*j2*S*cos(ky)
hy= 4.0*j3*S*cos(ky)*sin(kx)
hz=2.0d0*j1*S*cos(kx)*(r-1.0d0)
h=sqrt(hx*hx+hy*hy*hz*hz)
energia=h0-h
end function autovalor2


double precision function berry(kx, ky)
use dissertacao
!Esta função nos dá os valores da curvatura de Berry no espaço recíproco e deve ser Calculado analiticamente para cada ponto
!Por esta razão é recomendável usá-la como um procedimento externo para cada  rede cristalina em particular
 implicit none
!Variáveis
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.4d0
real(kind=kind(1.0d0)), parameter::j1=-0.4d0
real(kind=kind(1.0d0)), parameter::j2=-0.3d0
real(kind=kind(1.0d0)), parameter::j3=-0.3d0
real(kind=kind(1.0d0)), parameter::r=0.5d0
real(kind=kind(1.0d0))::p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
!Definição
p1=D**2*(sin(kx)+sin(ky))**2
p2=2.0d0*j2*D*cos(ky)*(sin(kx)+sin(ky))
p3=j2**2*cos(ky)**2
p4=4.0d0*j3**2*cos(ky)**2*sin(kx)**2
p5=j1**2*cos(kx)**2*(r-1.0d0)**2
p6=16.0d0*j1*j3*D*(sin(kx)+sin(ky))*sin(kx)**2*sin(ky)*(1.0d0-r)
p7=16.0d0*j1*j3*D*cos(ky)**2*sin(kx)**2*(1.0d0-r)
p8=8.0d0*j1*j3*D*cos(kx)*sin(ky)*sin(2.0d0*kx)
p9= 16.0d0*j1*j3*D*cos(kx)**2*cos(ky)**2*(1.0d0-r)
p10=8.0d0*j1*j2*j3*cos(kx)**2*sin(2.0d0*ky)

if ((p1+p2+p3+p4+p5)/=0.0d0) then
berry=(1/16.0d0)*(p1+p2+p3+p4+p5)**(-1.5)*(p6+p7+p8+p9+p10)
else
berry=0.0d0
endif


end function berry
