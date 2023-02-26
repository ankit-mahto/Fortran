module modl
    real(kind=8),dimension(:),allocatable::f,l,x
    real(kind=8)::ans=0,y
    integer::i,j,k,n=22
    
end module modl

program lagr
    use modl
    implicit none
    allocate(x(0:n))
    allocate(l(0:n))
    allocate(f(0:n))

    print*,"enter value of point"
    read*,y

    do i=0,n
        x(i)=0+i*8.0/22.0
        ! read*, x(i)
        ! read*, f(i)
        f(i)= 0.1745*(1-(exp(-0.90*((x(i)-1.4)**2))))
    enddo

    do i=0,n
        l(i)=1
        do j=0,n
            if(j/=i) then
            l(i)= l(i)*(y-x(j))/(x(i)-x(j))
            endif
        enddo
    enddo
    
    do i =0,n
    ans=ans+f(i)*l(i)
    enddo

    print* ,ans
end program lagr
