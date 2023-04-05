program Rk4

    integer::i,n=200
    real(kind=8),parameter::pi=acos(-1.0)
    real(kind=8)::x,s1,s2,E,r1,r2,h 
    real(kind=8), dimension(:), allocatable::z,r,y
    ! real(kind=8), dimension(:,:), allocatable::y
    
    allocate(y(0:n-1))
    allocate(z(0:n-1))
    allocate(r(0:n-1))

    h=0.02*2*pi
    
        y(0) = 1.0
        z(0) = 0.0
        x=0.0

    do i=0,n-2
        s1 = h*Rprime(x,y(i),z(i))
        r1 = h*pprime(x,y(i),z(i))
        s2 = h*Rprime(x+h,y(i)+s1,z(i)+r1)
        r2 = h*pprime(x+h,y(i)+s1,z(i)+r1)   
    y(i+1)= y(i) + (s1+s2)/2
    z(i+1)= z(i) + (r1+r2)/2

    x=x+h 
    enddo

    do i=0,n-1
    r(i)= (i)*0.02*2*pi
    enddo


open(10,file='XvsP1.txt')
open(20,file='Evst1.txt')
open(30,file='XvstT1.txt')
    do i=0,n-1
    write(10,*) y(i),z(i)
enddo


    do i=0,n-1
    write(20,*) r(i)/2*pi ,y(i)**2+z(i)**2 
enddo


    do i=0,n-1
    write(30,*) r(i)/2*pi,y(i)
enddo

    contains
    real(kind=8) function Rprime(x,y,z)
    real(kind=8)::y,x,z
    Rprime=z
    end function Rprime

    real(kind=8) function pprime(x,y,z)
    real(kind=8)::y,x,z
    pprime = -y
    end function pprime
    
    end program Rk4