module chem
real(kind=8),dimension(:,:),allocatable::q
real(kind=8)::kab=1,kba=3,kbc=4.2,kcb=7.3,kcd=0.4,dt=0.0001,t0=0,tmax=60,a0=5,b0=0,c0=0,d0=0,fa,fb,fc,fd
integer(kind=8)::n,i,j

end module chem

program kinetics
use chem
implicit none
n = tmax/dt
allocate(q(0:n,5))

open(1,file='output.txt')

call euler
do i=0,n
write(1,*) (q(i,j),j=1,5)
enddo

end program kinetics



subroutine euler()
    use chem 
    q(0,:) = (/t0,a0,b0,c0,d0/)
    do i=1,n
    q(i,1) = i*dt
    q(i,2) = q(i-1,2) + dt*(-kab*q(i-1,2) + kba*q(i-1,3))
    q(i,3) = q(i-1,3) + dt*(kab*q(i-1,2) - kba*q(i-1,3)-kbc*q(i-1,3) + kcb*(q(i-1,4)**2))
    q(i,4) = q(i-1,4) + dt*(2*kbc*q(i-1,3) - 2*kcb*(q(i-1,4)**2) - kcd*q(i-1,4))
    q(i,5) = q(i-1,5) + dt*(kcd*q(i-1,4))
	
	enddo
end subroutine euler

