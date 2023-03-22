Program Hello
    implicit none
    integer::i,n=100
    real(kind=8)::h,a,b,sum1=0.0,ans,x
    a=2.0 
    b=3.0
    
    h=(b-a)/real(n)
    do i=1,n/2
    x = a+(2*i-1)*h 
    sum1 = sum1+ 4*fx(x)
    enddo
    do i=1,(n/2-1)
    x = a+(2*i)*h 
    sum1 = sum1+ 2*fx(x)
    enddo
    
    ans = ( fx(a)+ sum1+fx(b) ) *h/3
    print* , ans
    contains 
    real(kind=8) function fx(a)
    real(kind=8)::a
    fx=a
    end function fx
    End Program Hello
    