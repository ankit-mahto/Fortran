module mod2

    real,parameter::pi=acos(-1.0)

    contains
    !_________________________________

    real function fact(n)
    integer::x
    if(n==0) then
        fact=1
ELSE
    do x=1,n
        fact = fact*n
    enddo
ENDIF
END FUNCTION FACT

    !_________________________________


    real function xcalc(i)
    real:: xmax =2*pi,delx,xmin=0
    delx=(xmax-xmin)/100
    xcalc = xmin + delx*i
    return
end function xcalc

    !_________________________________


real function sinex(x,n)
    real:: sum =0.0
    do i =1,4
        sum = sum + ((-1)**n)*(x**(2*n+1))/fact(2*n+1)
        enddo
        sinex=sum
        return
end function sinex

    !_________________________________


end module mod2 