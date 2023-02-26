program calc
    use mod2
    implicit none
  

    integer::n,i
    real::x
    open(unit=10,file='write2.txt')
    print*,"enter value of n"
    read*,n

    do i = 1,100
        x = xcalc(i)
        write(10,*) x, sinex(x,n)
    enddo
    end program calc