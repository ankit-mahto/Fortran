Program Hello
implicit none
integer::i,n=100
real,dimension(:),allocatable::x
real::h,a,b,sum=0,ans
allocate(x(n))
a=2 
b=3

h=(b-a)/(n)
do i=1,n-1
x(i) = a+(i)*h 
sum = sum + fx(x(i))
enddo

ans=h*( (fx(a) +fx(b))/2 + sum)
print* , ans

contains 
real function fx(a)
    real::a
    fx=a
end function fx

End Program Hello
