module parameters
implicit none
real,dimension(257)::arr_x,arr_v,arr_k
complex:: iota=(0,1)
real::pi=acos(-1.0)
real:: alpha=20.0, m=14500, xmin=-2.0,dx=0.02,dt=0.1,p0=20.0,x0=-0.5,h=1.0
integer::t=5000,N=257
end module

program prog
use parameters
implicit none
complex,dimension(257)::psi_t0,psi_t1,phi,psi_prime,psi0p,psi0dp,psi1part
integer::i,j

open(unit=10, file='500.txt')
do i=1,N
arr_x(i)= xmin+ (i-1)*dx
if(arr_x(i)<0) then
arr_v(i)=0.0
endif
if(arr_x(i) .ge. 0) then
arr_v(i)=1.0
endif
enddo

do i=0,(N/2)-1
arr_k(i+1) = (2*pi*i)/(N*dx)
enddo
do i=N/2,N-1
arr_k(i+1) = (2*pi*(i-N))/(N*dx)
enddo

do i=1,N
psi_t0(i)= (exp(iota*p0*(arr_x(i)-x0)))*(exp(-alpha*((arr_x(i)-x0)**2)))*(((2*alpha)/(pi))**0.25)
enddo

do j=1,5000

do i=1,N
psi0p(i)=exp(-iota*(arr_v(i)/2.0)*dt)*psi_t0(i)
enddo

call fwd_ft(psi0p,psi0dp)
do i=1,N
psi0dp(i)=psi0dp(i)*exp(-iota*arr_k(i)*arr_k(i)*dt/(2*m))
enddo
call rev_ft(psi0dp,psi1part)

do i=1,N
psi_t1(i) = exp(-iota*dt*arr_v(i)/2)*psi1part(i)
enddo

if(j .eq. 4500) then
do i=1,N
write(10,*) arr_x(i), (abs(psi_t1(i)))**2
enddo
endif

do i=1,N
psi_t0(i)=psi_t1(i)
enddo

enddo


contains
subroutine fwd_ft(psi,phi)
use parameters
complex,dimension(N)::psi,phi
integer::i,j
do j=1,N
phi(j)=0
do i=1,N
phi(j) = phi(j) + psi(i)*exp(-iota*arr_k(j)*arr_x(i))
enddo
phi(j) = phi(j) * (1/sqrt(257.0))
enddo 
end subroutine

subroutine rev_ft(phi,psi_prime)
use parameters
integer::i,j
complex,dimension(N)::phi,psi_prime
do j=1,N
psi_prime(j) = 0
do i=1,N
psi_prime(j) = psi_prime(j) + phi(i)*exp(iota*arr_k(i)*arr_x(j))
enddo
psi_prime(j) = psi_prime(j)*(1/sqrt(257.0))
enddo 

end subroutine
end
