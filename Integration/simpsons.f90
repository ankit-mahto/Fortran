Program Hello
    implicit none
    integer::i,n=100
    real,dimension(:),allocatable::x
    real::h,a,b,sum1=0,ans
    allocate(x(n))
    a=2 
    b=3
    
    h=(b-a)/(n)
    do i=1,n
    x(i) = a+(i-1)*h 
    if(mod(i,2)==1) then 
    sum1 = sum1+ 4*fx(x(i))
    else
    sum1 = sum1+ 2*fx(x(i))
    endif
    enddo
    
    ans = (fx(a)+ sum1+ fx(b) )*h/3
    print* , ans
    contains 
    real function fx(a)
    real::a
    fx=a
    end function fx
    End Program Hello
    