! this program is for solving tdse- (finite difference for time + finite difference for position)
program fd
    implicit none
    real(kind=8),parameter::pi=acos(-1.0)
    complex(kind=8),parameter::iota=(0,1)
    real(kind=8),dimension(:),allocatable::x,v
    real(kind=8)::dx=0.02,dt=0.1,m=14500.0
    complex(kind=8),dimension(:,:),allocatable::psi
    integer::i,j

    allocate(psi(0:255,0:5000))
    allocate(x(0:255))
    allocate(v(0:255))


    do i=0,255
        x(i)=-2.0+i*dx
        if(x(i)<0.0) then
            v(i)=0.0
        else
            v(i)=1.0
        endif
    psi(i,0) = ((2*20.0/pi)**0.25)* exp(iota*20.0*(x(i)+0.5))*exp(-20.0* ((x(i)+0.5)**2))
    enddo

    do i=0,255
        psi(i,1)=psi(i,0)+iota*dt*( ((1.0/2*m)*diff(i,0)) -v(i)*psi(i,0)) 
    enddo
    
    do j=1,4999
        do i=0,255
        psi(i,j+1)= psi(i,j-1)-2*iota*0.1*((-1.0/2.0*m)*diff(i,j)+v(i)*psi(i,j))
    enddo
enddo

    open(unit=10,file='write.txt')
    do i=0,255
        ! write(10,*) x(i),abs(psi(i,3)**2)
        ! print* , abs((psi(i,0)+iota*dt*(((1.0/2*m)*diff(i,0))-v(i)*psi(i,0)) )**2)
        write(10,*)  abs(psi(i,1)**2)
    ! write(10,*) (abs(psi**2), j=0,5000)
    enddo

    contains
    complex(kind=8) function diff(i,j)
    integer::i,j
    if(i==0) then
        diff = (1/0.02**2)*(-psi(i+3,j)+4*psi(i+2,j)-5*psi(i+1,j)+2*psi(i,j)) !forward diff
    elseif(i==255) then
        diff = (1/0.02**2)*(2*psi(i,j)-5*psi(i-1,j)+4*psi(i-2,j)-psi(i-3,j))  !backward diff
    else
    diff = (1/0.02**2)*(psi(i+1,j)-2*psi(i,j)+psi(i-1,j)) !central diff
    endif
    end function diff
end program fd