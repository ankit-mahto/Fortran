module pibs
    implicit none
    real, parameter::pi = acos(-1.0)
    ! pi = acos(-1.0)

    contains
    
    real function xcalc(i)
    real ::xmin,xmax,deltax
    integer::i 
    xmin = 0
    xmax = 1
    deltax = (xmax-xmin)/200
    xcalc = xmin+deltax*i
    return
end function xcalc

real function psi(x,n)
    real::x
    integer::n
    psi = sin(n*pi*x) 

end function psi

end module pibs