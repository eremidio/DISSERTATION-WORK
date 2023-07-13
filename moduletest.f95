
program moduletest
use dissertacao

implicit none

!variaveis de entrada para cálculo da função
real(kind=kind(1.0d0)), dimension(10)::kx, ky, x1
real(kind=kind(1.0d0))::autovalor, T
integer::i

!Priocedimentos externos
interface
function autovalor(kx, ky)
real(kind=kind(1.0d0))::energia
real(kind=kind(1.0d0)),intent(in)::kx, ky
end function  autovalor
end interface

!Procedimentos executáveis
!1. Exibindo as constantes usadas
call show()

!1.Construido as entradas para o programa em questão
do i=1,10
kx(i)=2.0*pi*i/10.0
ky(i)=2.0*pi*i/10.0
end do

!2.Calculando as funções os respectivos pontos
!A temperatura é um parâmetro a ser declarado pelo progrmador
 T=10.0
write(*,*)"kx(i)                  ky(i)                  energia               número de estados ocupados" 

do i=1,10
write(*,*) kx(i), ky(i), autovalor(kx(i), ky(i)), boseeinstein(T, autovalor(kx(i), ky(i)))
enddo

!Alocando o valor do número de estados ocupados em em outro arranjo
do i=1,10
x1(i) = boseeinstein(T, autovalor(kx(i), ky(i)))
enddo

write(*,*)"número de estados ocupados         Li_2(x)              c1(x)                   c2(x)"

do i=1,10
write(*,*) x1(i), polilog2(x1(i)), constant1(x1(i)), constant2(x1(i))
enddo

end program moduletest


function autovalor(kx, ky) result(energia)
!Esta função nos dá os autovalores do hamiltoniano em termos do parametro K no espaço recíproco e deve ser Calculado analiticamente para cada caso
!Por esta razão é recomendável usá-la como um procedimento externo para cada situação
!Este é um modelo de brinquedo 
 implicit none
!Resultado da função
real(kind=kind(1.0d0))::energia
!Variáveis
real(kind=kind(1.0d0)),intent(in)::kx, ky
!Definição
energia=1.6e-34*(sin(kx))**2*cos(ky) + 0.4e-34*sin(kx)*sin(ky)**2+0.5e-34*cos(kx)*cos(ky)**3
end function autovalor
