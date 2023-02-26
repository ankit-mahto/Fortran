module modn
    real(kind=8),dimension(:),allocatable::a,x
    real(kind=8),dimension(:,:),allocatable::f
    real(kind=8)::ans=0,pro,y
    integer::i,j,k,n=2
    
end module modn

program new
    use modn
    implicit none
    allocate(x(0:n))
    allocate(f(0:n,0:n))

    call calc()
    call calf()
    call diff()

    print*,"enter value of point"
    read*,y

    ans=f(0,0)
    do i=1,n
        pro=1
        do j=0,i-1
            pro=pro*(y-x(j))
        enddo
        ans = ans+f(0,i)*pro
    enddo

    print*,f(0,0),f(0,1),f(0,2),ans
end program new

subroutine calc()
    use modn
    print*, "enter values of x"
    do i =0,n
        read* , x(i)
    enddo
end subroutine calc

subroutine calf()
    use modn
    print*, "enter values of f(x)"
    do i=0,n
        read* , f(i,0)
    enddo
end subroutine calf

subroutine diff()
    use modn
    integer::r,c
    do c=1,n
        do r=0,n-c
            f(r,c)=(f(r+1,c-1)-f(r,c-1))/(x(r+c)-x(r))
        enddo
    enddo
end subroutine diff