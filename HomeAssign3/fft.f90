!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Subroutines from Numerical Recipes
!     X: The Input data
!     N: SIZE of the data
!     ISIGN: +1 for FORWARD TRANSFORMATION
!     ISIGN: -1 for INVERSE TRANSFORMATION
!     A simple call for position to momentum transformation to this routine will be as follows
!     call fft(PSI,NX,1)
!     While compiling this routine with the main program, there may appear four "WARNINGS". Ignore
!     warnings.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fft(x,n,isign)
        implicit real*8(a-h,o-z)

        parameter(npts=10,nptx=2*2**npts)
        complex*16 s,v,w,x(n),cstore
!       dimension of cstore should be at the least twice the grid size
        dimension cstore(nptx)
        complex conjg
        data ntbl/0/
!       roots of unity in the first call
        if(n.gt.ntbl)then
        ntbl=n
        pi=4.0*atan(1.00)
        j=1
        icnt=0
10      s=pi*(0.0,1.0)/float(j)
        do 20 k=0,j-1
        icnt=icnt+1
20      cstore(icnt)=exp(s*float(k))
        j=j+j
        if(j.lt.n)go to 10
        end if
!       permutation of x(j)
        j=1
        do 30 i=1,n
        if(i.le.j)then
        v=x(j)
        x(j)=x(i)
        x(i)=v
        end if
        m=n/2
25      continue
        if(j.gt.m)then
        j=j-m
        m=m/2
        if(m.ge.1)go to 25
        else
        j=j+m
        end if
30       continue
!       multiply x(j) and the roots of unity
        j=1
        icnt=0
40       jj=j+j
         do 50 k=1,j
        icnt=icnt+1
        w=cstore(icnt)
        if(isign.lt.0)w=conjg(w)
        do 50 i=k,n,jj
        v=w*x(i+j)
        x(i+j)=x(i)-v
50      x(i)=x(i)+v
        j=jj
        if(j.lt.n) go to 40

        return
        end
