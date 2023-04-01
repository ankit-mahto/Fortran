program huen

    integer::i,n=10000
    real(kind=8)::x,s1,s2,s3,s4,E,r1,r2,r3,r4 
    real(kind=8), dimension(:), allocatable::z,r
    real(kind=8), dimension(:,:), allocatable::y
    
    allocate(y(0:n-1,0:20))
    allocate(z(0:n-1))
    allocate(r(0:n-1))

    h=0.0005
    do j=0,20
        E = -0.6 +j*0.01
        y(0,j) = 0.000001
        z(0) = -1000.0
        x=0.0005
    do i=0,n-2
        s1 = h*Rprime(x,y(i,j),z(i))
        r1 = h*pprime(x,y(i,j),z(i))
        s2 = h*Rprime(x+h/2,y(i,j)+s1/2,z(i)+r1/2)
        r2 = h*pprime(x+h/2,y(i,j)+s1/2,z(i)+r1/2)
        s3 = h*Rprime(x+h/2,y(i,j)+s2/2,z(i)+r2/2)
        r3 = h*pprime(x+h/2,y(i,j)+s2/2,z(i)+r2/2)
        s4 = h*Rprime(x+h,y(i,j)+s3,z(i)+r3)
        r4 = h*pprime(x+h,y(i,j)+s3,z(i)+r3)   
    y(i+1,j)= y(i,j) + (s1+2*s2+2*s3+s4)/6
    z(i+1)= z(i) + (r1+2*r2+2*r3+r4)/6

    x=x+h 
    enddo

    do i=0,n-1
    r(i)= (i+1)*0.0005
    enddo
enddo

    open(10,file='write.txt')
    do i=0,n-1
    write(10,*) r(i),(abs(r(i)*y(i,j))**2, j=0,20)
enddo

    contains
    real(kind=8) function Rprime(x,y,z)
    real(kind=8)::y,x,z
    Rprime=z
    end function Rprime

    real(kind=8) function pprime(x,y,z)
    real(kind=8)::y,x,z
    pprime = -2*((1.0/x)*z + (E+1.0/x)*y)
    end function pprime
    
    end program huen