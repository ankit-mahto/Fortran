module q1
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx=0.02, dk, dt=0.1, alpha=20, x0=-0.5, p0=40, m=14500, xmin=-2,l,pmax,tprob=0.0,rprob=0.0,total=0.0
    real(kind=8),dimension(:),allocatable::x,v,k
    complex(kind=8),dimension(:,:),allocatable::psi
    complex(kind=8),dimension(:),allocatable::phi,d2psi,psip,psidp
    complex(kind=8),parameter::iota=(0,1)
    complex(kind=8)::sum
    integer::i,j,t,nt=5000,nx=256
    character(len= 1024)::number
end module q1

program hello
    use q1
    implicit none
    l=nx*dx
    dk = 2*pi/l
    pmax = pi/dx
    open(1,file='out.txt')
    allocate(x(1:nx),k(1:nx),v(1:nx),phi(1:nx),psip(1:nx),psidp(1:nx),psi(1:nx,0:nt),d2psi(1:nx))
    do i =1,nx
        !setting x grid
        x(i) = xmin + i*dx
        
        !setting potential grid
        if(x(i)>0 .and. x(i)<0.5) then
            v(i)=0.1
        else
            v(i)=0
        endif
        !initialising psi at t=0
        psi(i,0) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*p0*(x(i)-x0))
    enddo
     !setting momentum grid
    do i=1,nx/2
            k(i) =  (i-1)*dk
    enddo
    do i=nx/2+1 , nx
        k(i) =  (i-1-nx)*dk
    enddo
	
    do t=1,nt

        !psip contains the first part of split operator
        do i=1,nx
        psip(i)=exp(-iota*(v(i)/2.0)*dt)*psi(i,t-1)
        enddo
        !phi contains the forward fourier transformation
        do j=1,nx
        phi(j)=0
        do i=1,nx
        phi(j) = phi(j) + psip(i)*exp(-iota*k(j)*x(i))
        enddo
        phi(j) = phi(j) * (1/sqrt(real(nx)))
        enddo 
        !psidp contains the momentum array multiplied by -k**2/2m
        do i=1,nx
        psidp(i)=phi(i)*exp(-iota*k(i)*k(i)*dt/(2*m))
        enddo
        !d2psi contains the backard fourier transformation
        do j=1,nx
        d2psi(j) = 0
        do i=1,nx
        d2psi(j) = d2psi(j) + psidp(i)*exp(iota*k(i)*x(j))
        enddo
        d2psi(j) = d2psi(j)*(1/sqrt(real(nx)))
        enddo 
        !last part of multiplication by v/2 of split operator
        do i=1,nx
        psi(i,t) = exp(-iota*dt*v(i)/2)*d2psi(i)
        enddo

    enddo

    do i=1,nx
    total = total + abs(psi(i,5000)**2)
    enddo
    do i=125,nx
    tprob = tprob + abs(psi(i,5000)**2)
    enddo
    do i=1,100
    rprob = rprob + abs(psi(i,5000)**2)
    enddo
        print* , "transmission probability =", tprob/total , "reflection probability =", rprob/total
    
end program hello
