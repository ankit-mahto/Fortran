program hamiltonian
    implicit none

    real ,parameter::pi = acos(-1.0)
    real(kind=8),dimension(:),allocatable::x,v,evalue
    real(kind=8),dimension(:,:),allocatable::ham,evector
    real(kind=8)::length=5.0,hbar=1.0,mass=1.0,dx,t,ex_ev
    integer::i,j,n,state

    print*, "enter value of state"
    read*,state
    print*, "enter number of grid points"
    read*,n

    dx= length/dble(n-1)
    t = hbar**2/(2*mass*dx**2)

    allocate(x(1:n))
    allocate(v(1:n))
    allocate(ham(1:n,1:n))
    ! allocate(x(n,n))

    do i=1,n
        x(i)= (i-1)*dx
        v(i) = 0.1*x(i)
    enddo

    ex_ev = ((state**2)*(hbar**2)) / (8.0*mass*(length**2))

    do i=1,n
        do j=1,n
            if(i == j) then
                ham(i,j) = v(i) + 2*t
            else if(abs(i-j)==1) then
                ham(i,j) = -t
            else 
                ham(i,j)=0.d0
            endif
        enddo
    enddo

    open(unit=10,file='write.txt')

    ! d
    print*,ex_ev


end program hamiltonian




