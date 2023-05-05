module q1
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::dx,dt=8,p0,xmin,xmax,hbar=1,l,dk,mass,dxm
    real(kind=8),dimension(:),allocatable::x,v11,v22,v12,k,d
    complex(kind=8),dimension(:,:),allocatable::psi1,psi2
    complex(kind=8),dimension(:),allocatable::phi1,phi2,phi1p,phi2p
    complex(kind=8),parameter::iota=(0,1)
    integer(kind=8)::nx=2048,nt=1000,i

end module q1

program hello
    use q1
    implicit none
    open(1,file='potential.txt')
    open(2,file='wave1.txt')
    open(3,file='wave2.txt')
    allocate(x(nx),k(nx),d(nx),v11(nx),v22(nx),v12(nx),psi1(nx,0:nt),psi2(nx,0:nt),phi1(nx),phi2(nx),phi1p(nx),phi2p(nx))
    call grid
    call potential
    call init
    call propogation
    
    do i=1,nx
        write(1,*) x(i),k(i),d(i),v11(i),v22(i),v12(i)
        write(2,*) x(i),abs(psi1(i,0)**2),abs(psi1(i,200)**2),abs(psi1(i,400)**2),abs(psi1(i,600)**2),abs(psi1(i,1000)**2)
        write(3,*) x(i),abs(psi2(i,0)**2),abs(psi2(i,200)**2),abs(psi2(i,400)**2),abs(psi2(i,600)**2),abs(psi2(i,1000)**2)
    enddo
end program hello
subroutine grid
    use q1
    real(kind=8)::xmask
    xmin=-45
    xmax=45
    dx=(xmax-xmin)/real(nx-1)
    l=nx*dx
 dk = 2*pi/l
    do i=1,nx
        x(i)=xmin+(i-1)*dx
        if(i<=nx/2) then
        k(i)=(i-1)*dk
        else 
        k(i)=(i-1-nx)*dk
        endif
        
        if(x(i)<=-30) then
            xmask=-30
            dxm=xmin-xmask
            d(i)=sin(pi*(xmask+dxm-x(i))/(2*dxm))
        else if(x(i)>=10) then
            xmask=10
            dxm=xmax-xmask
            d(i)=sin(pi*(xmask+dxm-x(i))/(2*dxm))
        else
            d(i)=1
        endif
    enddo
end subroutine grid

subroutine init
    use q1
    real(kind=8)::x0,sigma,beta
    integer(kind=8)::nx0
    x0=9
    sigma=0.3
    beta=1./(4*(sigma**2))
    mass=3474.057
    nx0=int((x0-xmin)/dx) + 1
    p0=-sqrt(2*mass*(0.029-v11(nx0)))

    psi1(:,0)= ((1/(2*pi*(sigma**2))**(0.25)))*exp(-beta*((x-x0)**2))*exp(iota*p0*(x-x0))*d
    psi2(:,0)=0
end subroutine init

subroutine potential
    use q1
    real(kind=8)::v1,v2,v3,v4,v5,vasym,vlower,x1,x2,x3,x4,beta1,beta2,beta3,beta4
    real(kind=8),dimension(nx)::v1ad,v2ad,f

    v1=4.0167971782296E-2
    v2=4.79833373E-3
    v3=9.8998917754E-1
    v4=1.122019E-2
    v5=7.9781762366E-1

    vasym=3.61196179E-1
    vlower=0.0

    x1=-4.364721325998E-2
    x2=5.0012635420E-2
    x3=-7.6042693477E-1
    x4=8.1790045179E-1

    beta1=5.5
    beta2=4.9818195151
    beta3=2.3471780470
    beta4=1.0487590725

    v1ad=(v1*exp(beta1*(x-x1)))/((1+exp(beta1*(x-x1)))**2) + (v2*exp(beta1*(x-x1)))/(1+exp(beta1*(x-x1)))
    v2ad=vasym-(v3*exp(beta2*(x-x2)))/((1+exp(beta2*(x-x2)))**2)-(v4*exp(beta2*(x-x2)))/(1+exp(beta2*(x-x2)))&
    -(v5*exp(beta3*(x-x3)))/((1+exp(beta3*(x-x3)))**2) - vlower

    f=(1-tanh(beta4*(x-x4)))/2.0

    v11=(1-f)*v1ad+f*v2ad
    v22=f*v1ad+(1-f)*v2ad
    v12=-sqrt(f*(1-f))*(v2ad-v1ad)
    
end subroutine potential

subroutine propogation
    use q1
    do i=1,nt
    phi1=psi1(:,i-1)
    phi2=psi2(:,i-1)
    
    call split

    call fft(phi1,nx,1)
    phi1=phi1*exp(-iota*k*k*dt/(2*mass))
    call fft(phi1,nx,-1)
    phi1=phi1/real(nx)

    call fft(phi2,nx,1)
    phi2=phi2*exp(-iota*k*k*dt/(2*mass))
    call fft(phi2,nx,-1)
    phi2=phi2/real(nx)
   
    call split
    
    psi1(:,i)=phi1*d
    psi2(:,i)=phi2*d
    
    enddo
    
end subroutine propogation

subroutine split
    use q1
    
phi1=exp(-iota*dt*v11/(4*hbar))*phi1
phi2=exp(-iota*dt*v22/(4*hbar))*phi2

phi1p=cos(v12*dt/(2*hbar))*phi1 - iota*sin(v12*dt/(2*hbar))*phi2
phi2p=-iota*sin(v12*dt/(2*hbar))*phi1 + cos(v12*dt/(2*hbar))*phi2

phi1=exp(-iota*dt*v11/(4*hbar))*phi1p
phi2=exp(-iota*dt*v22/(4*hbar))*phi2p
    
    
end subroutine split

