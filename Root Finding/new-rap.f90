program newton 
    real,parameter::epsilon=1.0E-6
    INTEGER::i
    real,dimension(:),allocatable::x
    real::a

    a=0.0
    allocate(x(0:1000))
    x(0)=a

    do i=1,1000
        x(i) = x(i-1)-f(x(i-1))/df(x(i-1))
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
    f=x**2-3*x+2
    end function f

    real function df(x)
    real x
    df=x*2-3
    end function df
end program newton