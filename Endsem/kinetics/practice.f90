module q
real(kind=8),dimension(:),allocatable::a,b,c,d
real(kind=8)::k1,k11,k2,kcb,kcd,kdc,dt,t
real(kind=8)::s1,s2,p1,p2,q1,q2,r1,r2
integer::i,j,n
end module q

program hello 
    use q
    implicit none

    n=6000
        k1=0.1
        k11=0.05
        k2 =0.21  
        dt=0.05
        allocate(a(n),b(n),c(n),d(n))
        a(1)=1.0
        b(1)=0.0
        c(1)=0.0
        ! d(1)=0.0
    
    t=0.0

    do i=1,n-1
        s1=dt*aprime(a(i),b(i),c(i))
        p1=dt*bprime(a(i),b(i),c(i))
        q1=dt*cprime(a(i),b(i),c(i))
        ! r1=dt*dprime(a(i),b(i),c(i))
        s2=dt*aprime(a(i)+s1,b(i)+p1,c(i)+q1)
        p2=dt*bprime(a(i)+s1,b(i)+p1,c(i)+q1)
        q2=dt*cprime(a(i)+s1,b(i)+p1,c(i)+q1)
        ! r2=dt*dprime(a(i)+s1,b(i)+p1,c(i)+q1)
        a(i+1)= a(i)+(s1+s2)/2.0
        b(i+1)= b(i)+(p1+p2)/2.0
        c(i+1)= c(i)+(q1+q2)/2.0
        ! d(i+1)= d(i)+(r1+r2)/2.0
        t=t+dt
        if(a(i)<0.5) then
            print* , b(i),i
            exit
        endif
    enddo
    open(10,file='out.txt')
    ! do i=1,n
    !     write(10,*) i*dt,a(i), b(i) , c(i) 
    ! enddo
    
contains
    real(kind=8) function aprime(a,b,c)
        real(kind=8)::a,b,c
    aprime=-k1*a+k11*b
    end function aprime
    real(kind=8) function bprime(a,b,c)
        real(kind=8)::a,b,c
        bprime=k1*a-k11*b-k2*b
    end function bprime
    real(kind=8) function cprime(a,b,c)
        real(kind=8)::a,b,c
        cprime= k2*b
    end function cprime
    ! real(kind=8) function dprime(a,b,c)
    !     real(kind=8)::a,b,c
    !     dprime = kcd*c
    ! end function dprime

end program hello

