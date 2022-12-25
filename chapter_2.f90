! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes for chapter 2.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHAPTER 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.1.1
subroutine IPVR(v1,v2,n,c)
implicit none
integer::n,i
real*8::v1(0:n-1),v2(0:n-1),c
c=0.0d0
do i=0,n-1
c=c+v1(i)*v2(i)
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.1.2
subroutine IPVC(v1,v2,n,c)
implicit none
integer::n,i
complex*16::v1(0:n-1),v2(0:n-1),c
c=0.0d0
do i=0,n-1
c=c+dconjg(v1(i))*v2(i)
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.2.1
subroutine NORMVR(v,n,norm)
implicit none
integer::n,i
real*8::v(0:n-1),c,norm
c=0.0d0
do i=0,n-1
c=c+v(i)*v(i)
end do
norm=dsqrt(c)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.2.2
subroutine NORMVC(v,n,norm)
implicit none
integer::n,i
complex*16::v(0:n-1),c
real*8::norm
c=0.0d0
do i=0,n-1
c=c+conjg(v(i))*v(i)
end do
norm=dsqrt(abs(c))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.3.1
subroutine NORMALIZEVR(v,n)
implicit none
integer::n,i
real*8::v(0:n-1),c
c=0.0d0
do i=0,n-1
c=c+v(i)*v(i)
end do
v=v/dsqrt(c)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.3.2
subroutine NORMALIZEVC(v,n)
implicit none
integer::n,i
complex*16::v(0:n-1),c
c=0.0d0
do i=0,n-1
c=c+conjg(v(i))*v(i)
end do
v=v/dsqrt(abs(c))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.4.1
subroutine OPVR(v,w,n1,n2,A)
implicit none
integer::n1,n2,i,j
real*8::v(0:n1-1),w(0:n2-1),A(0:n1-1,0:n2-1)
do i=0,n1-1
do j=0,n2-1
A(i,j)=v(i)*w(j)
end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.4.2
subroutine OPVC(v,w,n1,n2,A)
implicit none
integer::n1,n2,i,j
complex*16::v(0:n1-1),w(0:n2-1),A(0:n1-1,0:n2-1)
do i=0,n1-1
do j=0,n2-1
A(i,j)=v(i)*conjg(w(j))
end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.5.1
subroutine MMULMR(A,B,m,n,p,C)
implicit none
integer::n,m,l,i,j,k,p
real*8::A(0:m-1,0:n-1),B(0:n-1,0:p-1),C(0:m-1,0:p-1),temp
do i=0,m-1,1
    do j=0,p-1,1
        temp = 0.0d0
        do k=0,n-1,1
            temp=temp+A(i,k)*B(k,j)
        enddo
        C(i,j) = temp
    enddo
enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.5.2
subroutine MMULMC(A,B,m,n,p,C)
implicit none
integer::n,m,l,i,j,k,p
complex*16::A(0:m-1,0:n-1),B(0:n-1,0:p-1),C(0:m-1,0:p-1),temp
do i=0,m-1,1
    do j=0,p-1,1
        temp = 0.0d0
        do k=0,n-1,1
            temp=temp+A(i,k)*B(k,j)
        enddo
        C(i,j) = temp
    enddo
enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.6.1
subroutine KPMR(m,n,A,p,q,B,C)
implicit none
integer::m,n,p,q,i,j,k,l
real*8,dimension(0:m-1,0:n-1)::A
real*8,dimension(0:p-1,0:q-1)::B
real*8,dimension(0:m*p-1,0:n*q-1)::C
do k=0,m-1,1
    do l=0,p-1,1
        do i=0,n-1,1
            do j=0,q-1,1
                C(p*k+l,q*i+j)=A(k,i)*B(l,j)
            end do
        end do
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.6.2
subroutine KPMC(m,n,A,p,q,B,C)
implicit none
integer::m,n,p,q,i,j,k,l
complex*16,dimension(0:m-1,0:n-1)::A
complex*16,dimension(0:p-1,0:q-1)::B
complex*16,dimension(0:m*p-1,0:n*q-1)::C
do k=0,m-1,1
    do l=0,p-1,1
        do i=0,n-1,1
            do j=0,q-1,1
                C(p*k+l,q*i+j)=A(k,i)*B(l,j)
            end do
        end do
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.7.1
subroutine TRMR(A,n,trace)
integer::i,n
real*8::A(0:n-1,0:n-1),trace
trace=0.0d0
do i=0,n-1
    trace=trace+A(i,i)
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.7.2
subroutine TRMC(A,n,trace)
integer::i,n
complex*16::A(0:n-1,0:n-1),trace
trace=0.0d0
do i=0,n-1
    trace=trace+A(i,i)
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.8.1
subroutine COMMR(n,A,B,C)
implicit none
integer::n
real*8::A(0:n-1,0:n-1),B(0:n-1,0:n-1),C(0:n-1,0:n-1)
C=matmul(A,B)-matmul(B,A)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.8.2
subroutine COMMC(n,A,B,C)
implicit none
integer::n
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),C(0:n-1,0:n-1)
C=matmul(A,B)-matmul(B,A)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.9.1
subroutine ANTICOMMR(n,A,B,C)
implicit none
integer::n
real*8::A(0:n-1,0:n-1),B(0:n-1,0:n-1),C(0:n-1,0:n-1)
C=matmul(A,B)+matmul(B,A)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.9.2
subroutine ANTICOMMC(n,A,B,C)
implicit none
integer::n
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),C(0:n-1,0:n-1)
C=matmul(A,B)+matmul(B,A)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.10
subroutine BTOD(m,n,s)
implicit none
integer*4::s
integer,dimension(0:s-1)::m
integer::n,i,k
n=0
do i=0,s-1
    k=(2**(i))*(m(s-1-i))
    n=n+k
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.11
subroutine DTOB(m,tt,s)
implicit none
integer::s
integer,dimension(0:s-1)::tt
integer::m,k,a2
tt=0
a2=m
do k = 0,s-1,1
    tt(s-k-1) = mod(a2,2)
    a2 = a2/2
    if (a2== 0) then
        exit
    end if
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.12
subroutine DTOBONEBIT(m,ia,i0,s)
implicit none
integer*4::s
integer,dimension(0:s-1)::tt
integer::m,ia,i0
call DTOB(m,tt,s)
ia=tt(i0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.13
subroutine LEFTSHIFT(a,s,k)
implicit none
integer::s,k
integer,dimension(0:s-1)::a,b
integer::i,j
do j=0,k-1
    do i=0,s-2,1
        b(i)=a(i+1)
    end do
    b(s-1)=a(0)
    a=b
    b=0
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.14
subroutine RIGHTSHIFT(a,s,k)
implicit none
integer::s,k
integer,dimension(0:s-1)::a,b
integer::i,j
do j=0,k-1
    do i=0,s-2,1
        b(i+1)=a(i)
    end do
    b(0)=a(s-1)
    a=b
    b=0
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!2.15
subroutine SWAP(s,i,j,a)
implicit none
integer::s,i,j,temp
integer,dimension(0:s-1)::a
temp=a(i)
a(i)=a(j)
a(j)=temp
end subroutine

