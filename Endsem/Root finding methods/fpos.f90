PROGRAM bisection_root_finding

    IMPLICIT NONE
    real,parameter::epsilon = 1.0E-6
    INTEGER :: i
    REAL :: a, b
    real,dimension(:),allocatable::x

    ! Set the initial interval 
    a = 0.0
    b = 3.0
    allocate(x(100))

    ! Loop until the interval is sufficiently small
    DO i = 1, 100
         ! Evaluate the function at the endpoints and midpoint of the interval
        x(i) = (a*f(b)-b*f(a) ) /(f(b)-f(a))
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

    ! Output the result
    WRITE(*,*) 'The root of the equation is: ', x(i)

    contains
    real function f(x)
    real::x
    f = x**3-2*x-5
    end function f

END PROGRAM bisection_root_finding
