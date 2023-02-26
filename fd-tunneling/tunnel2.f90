program tunnel
    implicit none
    
    real(kind=8),dimension(:),allocatable::x,V,En
    complex(kind=8),dimension(:),allocatable::psi
    complex(kind=8),parameter::iota = (0,1)
    real(kind=8)::k,E=9,m=1,h=1,max,min,pavg,t
    integer::i,j,n
    
    open (unit =10,file='write2.txt')  
    ! open (unit =10,file='write2.txt')   
    
    allocate(v(0:1000))
    allocate(x(0:1000))
    allocate(psi(0:1000))
    allocate(En(0:1000))
    
    do n=0,250           
    En(n) = 1+n*0.1           
    
    do i=0,1000
        x(i) = i*0.01
        V(i) = 0
        if(i>=400 .and. i<=500) then
            V(i) = 9
        endif
    enddo

    k=sqrt( 2*m*(En(n)-v(1))/(h**2))
    ! k=sqrt(2*m*( E-v(1) )/(h**2)) 

    psi(0)=1
    psi(1)=exp(-1*iota*k*x(1))
    
    do j=1,999
        psi(j+1) = (2-(2*m/(h**2)*(En(n)-V(j))*0.0001))*psi(j)-psi(j-1)    
        ! psi(j+1) = (2-(2*m/(h**2)*(E-V(j))*0.0001))*psi(j)-psi(j-1)      
    enddo
    
    max = abs(psi(601))**2
    min = abs(psi(601))**2
    do i=601,1000
        if(max<abs(psi(i))**2) then
            max=abs(psi(i))**2
        endif
        if(min>abs(psi(i))**2) then
            min=abs(psi(i))**2
        endif
    enddo
    
    pavg = (max+min)/2
    t=2/(1+pavg)
    
    write(10,*) En(n),t        
    enddo                      
    
    ! do i=0,1000                          
    !     write(10,*) x(i),abs(psi(i))**2  
    ! enddo                                
    
    end program tunnel