program redequadrada
use dissertacao
implicit none
!Variáveis
real(kind=kind(1.0d0)), dimension(:), allocatable::kx, ky
real(kind=kind(1.0d0))::q
character(len=*), parameter::filename="meshredequadrada.txt"
integer::i, j, n1, ios, fi

!Entrada de dados
write(*,*) "Quantos pontos você deseja usar na construção do mesh uniforme para a rede quadrada?"
read(*,*)n1

!Abrindo arquivo para registrar os valores de kx e ky
open(newunit=fi, file=filename, action="write", status="replace", iostat=ios)
!Calculando os pontos da malha quadrada para se fazer a integração

allocate(kx(n1), ky(n1))
q=pi/n1

do i =1, n1
 do j= 1, n1
kx(i)=q*(i-1)
ky(j)=q*(j-1)

write(unit=fi, fmt=*)kx(i), ky(j)
 end do
enddo
deallocate(kx, ky)

close(unit=fi, status="keep", iostat=ios)
end program redequadrada
