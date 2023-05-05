PROGRAM bisection_root_finding

    IMPLICIT NONE
    real,parameter::epsilon = 1.0E-6
    INTEGER :: i
    REAL :: a, b
    real,dimension(:),allocatable::x

    ! Set the initial interval 
    a = 0.0
    b = 4.0
    allocate(x(100))
    open(10,file='write.txt')
    WRITE(10,*) "a,b, x(i),error"
    ! Loop until the interval is sufficiently small
    DO i = 1, 100
         ! Evaluate the function at the endpoints and midpoint of the interval
        x(i) = (a + b) / 2.0
        ! Output the result
        WRITE(10,*) a,b, x(i),(x(i)-a)/x(i) 
        !check if the midvalue is exactly a root
        IF (f(x(i)) ==0) THEN
            EXIT
        END IF
        ! Check if the interval is sufficiently small
        IF (ABS((x(i) -a )/x(i)) < epsilon) THEN
            EXIT
        END IF
        ! Check if the root is in the left or right half of the interval
        IF (f(a) * f(x(i)) < 0.0) THEN
            b = x(i)
        ELSE
            a = x(i)
        END IF        
    END DO

    contains
    real function f(x)
    real::x
    f = x**2-3*x-9
    end function f

END PROGRAM bisection_root_finding
