
module dissertacao
implicit none
save
public::show, polilog2, constant1, constant2, boseeinstein
!Constantes usadas nos cálculos
!O valor de n define o número de pontos usados no cálculo da função polilogarítmica
real(kind=kind(1.0d0)), parameter, public::pi=4.0*atan(1.0)
real(kind=kind(1.0d0)), parameter, public::euler=exp(1.0)
real(kind=kind(1.0d0)), parameter, public::Boltzmann= 1.0d0
real(kind=kind(1.0d0)), parameter, public::Planck=1.0d0
integer, parameter, public::n=10000

real(kind=kind(1.0d0)):: polilog2, constant1, constant2, boseeinstein
contains
!Procedimentos

!Checagem das constantes usadas no cálculo dos coeficientes de transporte
subroutine show()
write(*,*)"n=", n
write(*,*) "pi=", pi
write(*,*)"e=",euler
write(*,*)"Constante de Planck=", Planck
write(*,*) "Constante de Boltzmann=", Boltzmann
end subroutine

!Lista defunções usadas

function constant1(x)
!Variáveis
real(kind=kind(1.0d0)), intent(in)::x
!Definição
constant1=(1.0d0+x)*log(1.0d0+x)-x*log(x)
end function constant1


function constant2(x)
!Variáveis
real(kind=kind(1.0d0)), intent(in)::x
constant2=(1.0d0+x)*(log((1.0d0+x)/x))**2-(log(x))**2-2.0d0*polilog2(-x)
end function constant2

function polilog2(x)
!Variaveis
real(kind=kind(1.0d0)), intent(in)::x
real(kind=kind(1.0d0))::soma
integer::i

if (abs(x)<1) then
soma=0.0
do i=1, n
soma= soma+((x**2)/(i**2))
enddo
elseif (abs(x)>1) then 
soma=((pi**2)/16)-(1/2)*(log(x))**2
do i=1, n
soma=soma+((-1)**(i-1))/(i**2*x**i)
enddo
endif

polilog2=soma
end function polilog2

function boseeinstein(T, x)
!Variaveis
real(kind=kind(1.0d0)), intent(in)::x, T
! Definição
boseeinstein=1.0d0/(exp(x/(Boltzmann*T))- 1.0d0)
endfunction boseeinstein

end module dissertacao



