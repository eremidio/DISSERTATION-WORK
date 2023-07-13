program unionjack1
use dissertacao
implicit none
!Variaveis
real(kind=kind(1.0d0))::autovalor1, autovalor2, T1, berry, soma, hall
real(kind=kind(1.0d0)), dimension(1000000)::kx, ky, x1, x2, x3, x4
integer::i, fu, fi, ios

!Procedimentos externos
interface
function autovalor1(kx, ky)
real(kind=kind(1.0d0))::energia
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.5d0
real(kind=kind(1.0d0)), parameter::j1= -0.2d0
real(kind=kind(1.0d0)), parameter::j2= -0.2d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
end function autovalor1

function autovalor2(kx, ky)
real(kind=kind(1.0d0))::energia
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.5d0
real(kind=kind(1.0d0)), parameter::j1=-0.2d0
real(kind=kind(1.0d0)), parameter::j2=-0.2d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
end function autovalor2

function berry(kx, ky)
real(kind=kind(1.0d0)),intent(in)::kx, ky
real(kind=kind(1.0d0)), parameter::S=sqrt(3.0d0)/2.0d0
real(kind=kind(1.0d0)), parameter::D=0.5d0
real(kind=kind(1.0d0)), parameter::j1=-0.2d0
real(kind=kind(1.0d0)), parameter::j2=-0.2d0
real(kind=kind(1.0d0))::p1, p2, p3, p4
end function berry
end interface



!Definição

!1. Lendo dados de arquivos e abrindo arquivo para exportar dados se necessário
open(newunit=fu, file='meshredequadrada.txt', status='old', action='read', iostat=ios)
open(newunit=fi, file='unionjack1.txt', status='replace', action='readwrite', iostat=ios)

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

!4.Cálculo da integral que dá a condutividade Hall de Spin
 soma= (((x3(1)+x4(1))*berry(kx(1), ky(1)) +(x3(1000000)+x4(1000000))*berry(kx(1000000), ky(1000000)))/2.0d0)*(pi**2/1000000d0)
do i= 2, 999999, 1
soma= soma+(x3(i)+x4(i))*berry(kx(i), ky(i))*(pi**2/1000000d0)
end do
hall= soma/(pi**2)

write(*,*) "O valor da condutividade Hall usando uma malha quadrada de 1000x1000 pontos a", T1, "K é", hall

!5.Exportando os dados para um arquivo externo
do i=1, 1000000, 1
write(unit=fi, fmt=*) kx(i), ky(i), autovalor1(kx(i), ky(i)), autovalor2(kx(i), ky(i)), x3(i), x4(i)
end do

close(unit=fu, status="keep", iostat=ios)
close(unit=fi, status="keep", iostat=ios)

end program unionjack1




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
real(kind=kind(1.0d0)), parameter::D=0.5d0
real(kind=kind(1.0d0)), parameter::j1=-0.2d0
real(kind=kind(1.0d0)), parameter::j2=-0.2d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
!Definição
h0= 4.0d0*j2*S*cos(kx)*cos(ky)- 8.0d0*j1*S-4.0d0*j2*S
hx=2.0d0*j1*S*(cos(kx)+cos(ky))
hy=(D/2.0d0)*S*(cos(kx)+cos(ky))
hz=4.0d0*j2*(cos(kx)*cos(ky)- 1.0d0)
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
real(kind=kind(1.0d0)), parameter::D=0.5d0
real(kind=kind(1.0d0)), parameter::j1=-0.2d0
real(kind=kind(1.0d0)), parameter::j2=-0.2d0
real(kind=kind(1.0d0))::h0, hx, hy, hz, h
!Definição
h0= 4.0d0*j2*S*cos(kx)*cos(ky)- 8.0d0*j1*S-4.0d0*j2*S
hx=2.0d0*j1*S*(cos(kx)+cos(ky))
hy=(D/2.0d0)*S*(cos(kx)+cos(ky))
hz=4.0d0*j2*(cos(kx)*cos(ky)- 1.0d0)
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
real(kind=kind(1.0d0)), parameter::D=0.5d0
real(kind=kind(1.0d0)), parameter::j1=-0.2d0
real(kind=kind(1.0d0)), parameter::j2=-0.2d0
real(kind=kind(1.0d0))::p1, p2, p3, p4

!Definição
p1=((D/2.0d0) + (2.0d0*j1))**2*(cos(kx)+cos(ky))**2
p2=16.0d0*j2**2*(cos(kx)*cos(ky)-1)**2
p3=4.0d0*j1*j2*D*(cos(kx)+cos(ky))
p4= cos(ky)*(sin(kx))**2 + cos(kx)*(sin(ky))**2-sin(2.0d0*kx)*sin(ky)
if ((p1+p2)/=0.0d0) then
berry=0.5d0*(p1+p2)**(-1.5)*p3*p4
else
berry=0.0d0
endif
end function berry
