program pib
    use pibs
    implicit none
    ! main()
    integer :: n,i
    real::x
    ! real :: pi=3.14
    open(unit=10,file = '1write.txt')
    print*,"enter the quantum state"
    read*, n
    
    do i=0,200 
        x = xcalc(i)
    write(10,*) x, psi(x,n)
enddo     

end program pib



