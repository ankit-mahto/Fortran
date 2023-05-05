program newton 
    real,parameter::epsilon=1.0E-6
    INTEGER::i
    real,dimension(:),allocatable::x
    real::a,b

    a=4.0
    b=2.0
    allocate(x(0:1000))
    x(0)=a
    x(1)=b

    do i=2,1000
        x(i) = x(i-1)- ( f(x(i-1))* (x(i-1)-x(i-2)) /(f(x(i-1))-f(x(i-2))))
        if(f(x(i))==0) THEN
            EXIT
            endif
    
        if(abs(x(i)-x(i-1))<epsilon) THEN
            EXIT
        endif   
    enddo
    print* , "the root of the equation is ", x(i)

    contains
    real function f(x)
    real x
    f=x**2-4*x-10
    end function f
end program newton