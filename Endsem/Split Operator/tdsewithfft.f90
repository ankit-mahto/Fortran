module q1
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx=0.02, dk, dt=0.1, alpha=20, x0=-0.5, p0=20, m=14500, xmin=-2,l,pmax,tprob,sum1,sum2,rprob
    real(kind=8),dimension(:),allocatable::x,v,k
    complex(kind=8),dimension(:,:),allocatable::psi
    complex(kind=8),dimension(:),allocatable::phi,d2psi,psip,psidp
    complex(kind=8),parameter::iota=(0,1)
    complex(kind=8)::sum
    integer::i,j,t,nt=5000,nx=256,np
end module q1

program hello
    use q1
    implicit none
    l=nx*dx
    dk = 2*pi/l
    pmax = pi/dx
    open(1,file='out.txt')
    allocate(x(nx),k(nx),v(nx),phi(nx),psip(nx),psidp(nx),psi(nx,0:nt),d2psi(nx))
    do i=1,nx/2
        k(i)=2.0*pi*(i-1)/(nx*dx)
    enddo
    do i=nx/2+1,nx
        k(i)=2.0*pi*(i-1-nx)/(nx*dx)
    enddo
    do i =1,nx
        x(i) = xmin + (i-1)*dx
        ! k(i) = -pmax + (i-1)*dk
        if(x(i)>=0 .and. x(i)<=0.5) then
            v(i)=0.1
        else
            v(i)=0
        endif
        psi(i,0) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*p0*(x(i)-x0))
    enddo

    do t=1,nt
    !for v/2 ooperator
        do i=1,nx
        psip(i)=exp(-iota*(v(i)/2.0)*dt)*psi(i,t-1)
        enddo
        !for T operator
	! fourier transformation
        call fft(psip,nx,+1)
       
        do i=1,nx
        psidp(i)=psip(i)*exp(-iota*k(i)*k(i)*dt/(2*m))
        enddo
	 ! inverse fourier transformation
       call fft(psidp,nx,-1)
       psidp=psidp/nx
	! for v/2 operator 
        do i=1,nx
        psi(i,t) = exp(-iota*dt*v(i)/2)*psidp(i)	!calculating final psi at time t
        enddo

    enddo
    !do i=0,nx
    !    write(1,*) x(i), (abs(psi(i,j)**2) , j=0,nt,100)
    !enddo
    do i=1,nx
    write(1,*) x(i), abs(psi(i,nt)**2), v(i)
    enddo
    
    np=2.5/dx +1	!end point of barrier
    sum1=0
    do i=np,nx
    sum1=sum1+abs(psi(i,nt)**2)	!numerator for transmission prob
    enddo
    sum2=0
    do i=1,nx
    sum2=sum2+abs(psi(i,nt)**2)	!denominator for transmission prob
    enddo
    tprob = sum1/sum2
    rprob = 1-tprob
    print*, "transmission probability=",tprob
    print*, "reflection probability=",rprob
   

end program hello
