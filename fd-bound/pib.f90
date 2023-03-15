program diagonalize_hamiltonian
    implicit none
    
    integer, parameter :: n = 10     ! number of grid points
    integer :: i, j, info
    real(kind=8), parameter :: length = 5.d0, hbar = 1.d0, mass = 1.d0
    real(kind=8), parameter :: dx = length/dble(n-1)
    real(kind=8), dimension(1:n) :: x, v
    real(kind=8), dimension(n,n) :: hamiltonian
    real(kind=8), dimension(n) :: eigenvalues
    real(kind=8), dimension(n,n) :: eigenvectors
    
    open(unit=10,file='write1.txt')
    ! Construct x grid and potential v(x)
    do i = 1, n
        x(i) = dble(i-1)*dx
        v(i) = 0.1d0*x(i)
    end do
    
    ! Construct Hamiltonian matrix
    do i = 1, n
        do j = 1, n
            if (i == j) then
                hamiltonian(i,j) = 2.d0*hbar**2/(2.d0*mass*dx**2) + v(i)
            else if (abs(i-j) == 1) then
                hamiltonian(i,j) = -hbar**2/(2.d0*mass*dx**2)
            else
                hamiltonian(i,j) = 0.d0
            end if
        end do
    end do
    
    ! Diagonalize Hamiltonian matrix
    ! call dsyev('V', 'U', n, hamiltonian, n, eigenvalues, eigenvectors, n, info)
    
    do i = 1, n
        write(10,*) (hamiltonian(i,j), j=1,n)
    end do

    ! Print out eigenvalues
    ! write(*,*) "Eigenvalues:"
    ! do i = 1, n
    !     write(*,*) eigenvalues(i)
    ! end do
end program diagonalize_hamiltonian
