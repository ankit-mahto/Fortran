Program Hello
implicit none
integer::i,n=100000
real,dimension(:),allocatable::x
real::h,a,b,sum=0 
allocate(x(n))
a=2 
b=3

h=(b-a)/(n-1)
do i=1,n
x(i) = a+(i-1)*h 
sum = sum+ fx(x(i))*h
enddo
print* , sum

contains 
real function fx(a)
    real::a
    fx=a
end function fx

End Program Hello
