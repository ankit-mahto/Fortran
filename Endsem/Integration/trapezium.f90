Program Hello
implicit none
integer::i,n=10
real(kind=8),dimension(:),allocatable::x
real(kind=8)::h,a,b,sum=0,ans
allocate(x(n))
a=2.0 
b=3.0

h=(b-a)/real(n)
do i=1,n-1
x(i) = a+i*h 
sum = sum + fx(x(i))
enddo

ans=h*( (fx(a) +fx(b))/2 + sum)
print* , ans

contains 
real(kind=8) function fx(a)
    real(kind=8)::a
    fx= exp(a)
end function fx

End Program Hello
