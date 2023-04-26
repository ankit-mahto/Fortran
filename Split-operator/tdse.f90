module q1
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx=0.02, dk, dt=0.1, alpha=20, x0=-0.5, p0=20, m=14500, xmin=-2,l,pmax
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
    dk = 2*pi/nx
    pmax = pi/dx
    open(1,file='out1.txt')
    allocate(x(0:nx),k(0:nx),v(0:nx),phi(0:nx),psip(0:nx),psidp(0:nx),psi(0:nx,0:nt),d2psi(0:nx))
    do i =0,nx
        x(i) = xmin + i*dx
        k(i) = -pmax + i*dk
        if(x(i)<0) then
            v(i)=0
        else
            v(i)=1
        endif
        psi(i,0) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*p0*(x(i)-x0))
    enddo

    do t=1,nt
        do i=1,nx
        psip(i)=exp(-iota*(v(i)/2.0)*dt)*psi(i,t-1)
        enddo

        do j=1,nx
        phi(j)=0
        do i=1,nx
        phi(j) = phi(j) + psip(i)*exp(-iota*k(j)*x(i))
        enddo
        phi(j) = phi(j) * (1/sqrt(real(nx)))
        enddo 

        do i=1,nx
        psidp(i)=phi(i)*exp(-iota*k(i)*k(i)*dt/(2*m))
        enddo
        
        do j=1,nx
        d2psi(j) = 0
        do i=1,nx
        d2psi(j) = d2psi(j) + psidp(i)*exp(iota*k(i)*x(j))
        enddo
        d2psi(j) = d2psi(j)*(1/sqrt(real(nx)))
        enddo 

        do i=1,nx
        psi(i,t) = exp(-iota*dt*v(i)/2)*d2psi(i)
        enddo

    enddo
    do i=0,nx
        write(1,*) x(i), (abs(psi(i,j)**2) , j=0,nt,100)
    enddo
end program hello
