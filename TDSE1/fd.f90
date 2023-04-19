module ft2
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx=0.02, dp, dt=0.1, alpha=20, x0=-0.5, p0=20, m=14500, xmin=-2,l
    real(kind=8),dimension(:),allocatable::x,v,p
    complex(kind=8),dimension(:,:),allocatable::psi
    complex(kind=8),dimension(:),allocatable::phi,psi2
    complex(kind=8),parameter::iota=(0,1)
    complex(kind=8)::sum
    integer(kind=8)::i,j,t,nt=5000,nx=256
    character(len= 1024)::number
end module ft2
program ha
    use ft2
    implicit none
    l=nx*dx
    dp = 2*pi/l
    open(1,file='out2.txt')
    allocate(x(0:nx),p(0:nx),v(0:nx),phi(0:nx),psi(0:nx,0:nt),psi2(0:nx))
    do i =0,nx
        x(i) = xmin + i*dx
        if(x(i)<0) then
            v(i)=0
        else
            v(i)=1
        endif
        p(i) = dp*(i-nx/2)
        psi(i,0) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*p0*(x(i)-x0))
    enddo
    psi(0,1)=0
    psi(nx,1)=0
    do i=1,nx-1
        psi(i,1) = psi(i,0) + iota*dt*(((psi(i+1,0)+psi(i-1,0)-2*psi(i,0))/(1*(dx**2)))/(2*m) - v(i)*psi(i,0))
    enddo
    do t=1,nt-1
        psi(0,t)=0
        psi(nx,t)=0
        do i=1,nx-1
            psi(i,t+1) = psi(i,t-1) + 2*iota*dt*(((psi(i+1,t)+psi(i-1,t)-2*psi(i,t))/(1*(dx**2)))/(2*m) - v(i)*psi(i,t))
        enddo
    enddo
    ! do i=0,nx
    !     do j=0,nt,100
    !     write(number,*) j
    !     open(j,file=trim(adjustl(number))//"dft.txt")
    !     write(j,*) x(i), abs(psi(i,j)**2)
    !     enddo
    ! enddo
    open(unit=10,file='write.txt')
    do i=0,256
    write(10,*) x(i),abs(psi(i,5000)**2)
    enddo

end program 