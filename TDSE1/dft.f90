module ft
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx=0.02, dp, dt=0.1, alpha=20, x0=-0.5, p0=20, m=14500, xmin=-2,l,pmin
    real(kind=8),dimension(:),allocatable::x,v,p
    complex(kind=8),dimension(:,:),allocatable::psi
    complex(kind=8),dimension(:),allocatable::phi,psi2
    complex(kind=8),parameter::iota=(0,1)
    complex(kind=8)::sum
    integer(kind=8)::i,j,t,nt=5000,nx=256
    character(len= 1024)::number
end module ft
!this program uses finite difference method for calculating kinetic part
program ha
    use ft
    implicit none
    l=nx*dx
    dp = 2*pi/l
    pmin = -pi/dx
    open(1,file='out.txt')
    allocate(x(0:nx),p(0:nx),v(0:nx),phi(0:nx),psi(0:nx,0:nt),psi2(0:nx))
    do i =0,nx
        x(i) = xmin + i*dx
        if(x(i)<0) then
            v(i)=0
        else
            v(i)=1
        endif
        p(i) = pmin + i*dp
        psi(i,0) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*p0*(x(i)-x0))
    enddo
    do i=0,nx
        sum = 0
        do j =0,nx
            sum=sum+ psi(j,0)*exp(-iota*p(i)*x(j))
        enddo
        phi(i) = -(p(i)**2)*sum/sqrt(real(nx))
    enddo
    do i=0,nx
        sum=0
        do j =0,nx
            sum = sum + phi(j)*exp(iota*x(i)*p(j))
        enddo
        psi2(i)=sum/sqrt(real(nx))
    enddo
    do i=0,nx
        psi(i,1) = psi(i,0) + iota*dt*(psi2(i)/(2*m) - v(i)*psi(i,0))
    enddo
    
    do t=1,nt-1
        do i=0,nx
            sum = 0
            do j =0,nx
                sum=sum+ psi(j,t)*exp(-iota*p(i)*x(j))
            enddo
            phi(i) = -(p(i)**2)*sum/sqrt(real(nx))
        enddo
        do i=0,nx
            sum=0
            do j =0,nx
                sum = sum + phi(j)*exp(iota*x(i)*p(j))
            enddo
            psi2(i)=sum/sqrt(real(nx))
        enddo
        do i=0,nx
            psi(i,t+1) = psi(i,t-1) + 2*iota*dt*(psi2(i)/(2*m) - v(i)*psi(i,t))
        enddo
    enddo

    open(unit=10,file='write1.txt')
    do i=0,256
    write(10,*) x(i),abs(psi(i,5000)**2)
    enddo

end program 