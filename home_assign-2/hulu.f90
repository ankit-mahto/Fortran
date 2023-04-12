module q2
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx=0.01,dt=0.0001
    real(kind=8),dimension(:),allocatable::x
    complex(kind=8),dimension(:,:),allocatable::psi
    complex(kind=8)::z=(0,1),s1,s2,s3,s4
    integer(kind=8)::i,j
end module q2
program wave
    use q2
    implicit none
    allocate(psi(0:1000,0:300),x(0:1000))
    open(1,file='psi.txt')
    do i=0,1000
        x(i) = i*dx
        psi(i,0)= ((2*20/pi)**0.25) * exp(-((x(i)-5.0)**2)/(2*0.1*0.1))*exp(z*20*x(i))
    enddo
    do j=0,299
        psi(0,j+1)=psi(0,0)
        psi(1000,j+1)=psi(1000,0)
        do i=1,999
            s1 = z*(psi(i+1,j)+psi(i-1,j)-2*psi(i,j))/(2*dx*dx)
            s2 = z*(psi(i+1,j)+psi(i-1,j)-2*(psi(i,j)+dt*s1/2))/(2*dx*dx)
            s3 = z*(psi(i+1,j)+psi(i-1,j)-2*(psi(i,j)+dt*s2/2))/(2*dx*dx)
            s4 = z*(psi(i+1,j)+psi(i-1,j)-2*(psi(i,j)+dt*s3))/(2*dx*dx)
            ! s2 = z*((psi(i+1,j)+dt*s1/2)+(psi(i-1,j)+dt*s1/2)-2*(psi(i,j)+dt*s1/2))/(2*dx*dx)
            ! s3 = z*((psi(i+1,j)+dt*s2/2)+(psi(i-1,j)+dt*s2/2)-2*(psi(i,j)+dt*s2/2))/(2*dx*dx)
            ! s4 = z*((psi(i+1,j)+dt*s3)+(psi(i-1,j)+dt*s3)-2*(psi(i,j)+dt*s3))/(2*dx*dx)
            psi(i,j+1) = psi(i,j) + dt*(s1+2*s2+2*s3+s4)/6
        enddo
    enddo
    do i=0,1000
        write(1,*) x(i),abs(psi(i,0)**2),abs(psi(i,1)**2), abs(psi(i,150)**2),abs(psi(i,300)**2)
    enddo
end program wave


