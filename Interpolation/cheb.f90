module modu
    real(kind=8),parameter::pi = acos(-1.0)
    real(kind=8),dimension(:),allocatable::t,f,c,x
    real(kind=8)::ans=0,sum
    integer::i,j,k,n=3
    
end module modu

program cheb
    use modu
    implicit none
    allocate(x(0:n))
    allocate(t(0:n))
    allocate(f(0:n))
    allocate(c(0:n))

    do i=0,n
        x(i)=cos((2*i+1)*pi/(2*n+2))
        f(i)=fx(x(i))
    enddo

    call calc()
    call tx(1.0)

    do i =0,n
        ans = ans+c(i)*t(i)
        print*,x(i),c(i)
    enddo
    print* , ans

    contains

    real(kind=8) function fx(y)
    use modu
    real(kind=8)::y
    fx=exp(y)
    return
    end function fx

    subroutine tx(y)
    use modu
    real::y
    t(0)=1
    t(1)=y
    do i=1,n-1
        t(i+1)=2*y*t(i)-t(i-1)
    enddo
    ! do i=0,n
    ! t(i)=cos(i*acos(y))
    ! enddo
    end subroutine tx

    subroutine calc()
    use modu
    do j=0,n
        sum=0
        do k=0,n
            sum = sum + f(k)*cos(j*pi*(k*2+1)/(2*n+2))
        enddo
        if (j==0) then
            c(j)=sum/(n+1)
        else 
            c(j)=2*sum/(n+1)
        endif
        enddo
    end subroutine calc

end program cheb

    