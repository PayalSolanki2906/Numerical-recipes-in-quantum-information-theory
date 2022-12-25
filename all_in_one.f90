! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes chapterwise.

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHAPTER 3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.1.1
subroutine MATINVMR(N,A,B)
implicit none
integer::N
real*8,dimension(0:N-1,0:N-1)::A,B
integer::LWORK,INFO,i
integer,dimension(0:N-1)::IPIV
double precision::WORK(0:N-1)
B=A
call DGETRF( N, N, B, N, IPIV, INFO )
call DGETRI( N, B, N, IPIV, WORK, N, INFO )
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.1.2
subroutine MATINVMC(N,A,B)
implicit none
integer::N
complex*16,dimension(0:N-1,0:N-1)::A,B
integer::LWORK,INFO
integer,dimension(0:N-1)::IPIV
complex*16::WORK(0:60*N-1)
B=A
call ZGETRF( N, N, B, N, IPIV, INFO )
call ZGETRI( N, B, N, IPIV, WORK, 60*N, INFO)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.2.1
subroutine LANCZOSSR(n,A,D,E)
implicit none
integer::n,INFO
real*8::A(0:n-1,0:n-1),D(0:n-1),E(0:n-2),TAU(0:n-1),WORK(0:n-1)
call DSYTRD( 'U', n, A, n, D, E, TAU, WORK,n, INFO )
call DORGTR( 'U', n, A, n, TAU, WORK, n, INFO )
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.2.2
subroutine LANCZOSHC(n,A,D,E)
implicit none
integer::n,INFO
complex*16::A(0:n-1,0:n-1)
complex*16::TAU(0:n-1),WORK(0:n-1)
real*8::D(0:n-1),E(0:n-2)
call ZHETRD( 'U', n, A, n, D, E, TAU, WORK,n , INFO )
call ZUNGTR( 'U', n, A, n, TAU, WORK, n, INFO )
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.3.1
subroutine QRMR(n,k,A,R)
implicit none
integer::n,k,INFO,i,j,temp
real*8::A(0:n-1,0:k-1),R(0:k-1,0:k-1)
real*8::WORK(0:k-1)
real*8,allocatable,dimension(:)::TAU
R=0.0d0
temp=min(n,k)
allocate(TAU(0:temp-1))
call DGEQRF( n, k, A, n, TAU, WORK, k, INFO )
do i=0,k-1
    do j=i,k-1
        R(i,j)=A(i,j)
    end do
end do
call DORGQR( n, k, k, A, n, TAU, WORK, n, INFO )
deallocate(TAU)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.3.2
subroutine QRMC(n,k,A,R)
implicit none
integer::n,k,INFO,i,j,temp
complex*16::A(0:n-1,0:k-1),R(0:k-1,0:k-1)
complex*16::WORK(0:k-1)
complex*16,allocatable,dimension(:)::TAU
R=0.0d0
temp=min(n,k)
allocate(TAU(0:temp-1))
call ZGEQRF( n, k, A, n, TAU, WORK, k, INFO )
do i=0,k-1
    do j=i,k-1
        R(i,j)=A(i,j)
    end do
end do
call ZUNGQR( n, k, k, A, n, TAU, WORK, k, INFO )
deallocate(TAU)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.4.1
subroutine FUNSR(A,n,F,B)
implicit none
real*8::A(0:n-1,0:n-1),B(0:n-1,0:n-1),D(0:n-1,0:n-1)
integer::n,i
character*3::JOBZ, UPLO ,F  !!!JOBZ='N' or 'V', UPLO='U' or 'L'
integer::INFO
real*8::W(0:N-1),WORK(0:3*N-2),pi
JOBZ='V'
UPLO='U'
pi=4.0d0*atan(1.0d0)
call DSYEV( JOBZ, UPLO, N, A, N, W, WORK, 3*N-1, INFO )
D=0.0d0
if (F=='exp') then
    do i=0,n-1
        D(i,i)=exp(W(i))
    end do
end if

if (F=='log') then
    do i=0,n-1
        if (W(i) .le. 0.0000000000001) then
            print*, ' Error occured' 
            print*, 'Logarithm of zero and & 
                     & negative number are not defined'
            exit
        end if
        D(i,i)=log(W(i))
    end do
end if

if (F=='sin') then
    do i=0,n-1
        D(i,i)=sin(W(i))
    end do
end if

if (F=='cos') then
    do i=0,n-1
        D(i,i)=cos(W(i))
    end do
end if

if (F=='tan') then
    do i=0,n-1
        D(i,i)=tan(W(i))
    end do
end if
B=matmul(matmul(A,D),transpose(A))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.4.2
subroutine FUNHC(A,n,F,B)
implicit none
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),D(0:n-1,0:n-1)
integer::n,i
character*3::JOBZ, UPLO,F   
integer::INFO
real*8::W(0:N-1),RWORK(0:3*N-3)
complex*16::WORK(0:2*N-2)
JOBZ='V'
UPLO='U'
call ZHEEV( JOBZ, UPLO, N, A, N, W, WORK, 2*N-1, RWORK,INFO )
D=0.0d0
if (F=='exp') then
    do i=0,n-1
        D(i,i)=exp(W(i))
    end do
end if

if (F=='log') then
    do i=0,n-1
        if (W(i) .le. 0.0000000000001) then
            print*, ' Error occured' 
            print*, 'Logarithm of zero and & 
                  &negative number are not defined'
            exit
        end if
        D(i,i)=log(W(i))
    end do
end if

if (F=='sin') then
    do i=0,n-1
        D(i,i)=sin(W(i))
    end do
end if

if (F=='cos') then
    do i=0,n-1
        D(i,i)=cos(W(i))
    end do
end if

if (F=='tan') then
    do i=0,n-1
        D(i,i)=tan(W(i))
    end do
end if
B=matmul(matmul(A,D),conjg(transpose(A)))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.4.3
subroutine FUNMR(A,n,F,E)
implicit none
integer::n,i,j
complex*16::B(0:n-1,0:n-1),D(0:n-1,0:n-1)
complex*16::E(0:n-1,0:n-1),C(0:n-1,0:n-1)
character*3::JOBVL, JOBVR ,F 
integer::INFO
complex*16::W(0:N-1),VL(0:N-1,0:N-1),VR(0:N-1,0:N-1)
complex*16::WORK(0:2*N-1)
real*8::RWORK(0:2*N-1),A(0:N-1,0:N-1)
JOBVL='N'
JOBVR='V'
do i=0,n-1
    do j=0,n-1
        C(i,j)=dcmplx(A(i,j),0)
    end do
end do
call ZGEEV( JOBVL, JOBVR, N, C, N, W, VL, N, VR, N, & 
            WORK, 2*N, RWORK, INFO )
D=0.0d0
if (F=='exp') then
    do i=0,n-1
        D(i,i)=exp(W(i))
    end do
end if

if (F=='log') then
    do i=0,n-1
        if (abs(W(i)) .le. 0.0000000000001) then
            print*, ' Error occured' 
            print*, 'Logarithm of zero and & 
                  &negative number are not defined'
            exit
        end if
        D(i,i)=log(W(i))
    end do
end if

if (F=='sin') then
    do i=0,n-1
        D(i,i)=sin(W(i))
    end do
end if

if (F=='cos') then
    do i=0,n-1
        D(i,i)=cos(W(i))
    end do
end if

if (F=='tan') then
    do i=0,n-1
        D(i,i)=tan(W(i))
    end do
end if
call MATINVMC(n,VR,B)
E=matmul(matmul(VR,D),B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.4.4
subroutine FUNMC(A,n,F,E)
implicit none
integer::n,i
complex*16::B(0:n-1,0:n-1),D(0:n-1,0:n-1)
complex*16::E(0:n-1,0:n-1)
character*3::JOBVL, JOBVR ,F 
integer::INFO
complex*16::W(0:N-1),VL(0:N-1,0:N-1),VR(0:N-1,0:N-1)
complex*16::WORK(0:2*N-1),A(0:N-1,0:N-1)
real*8::RWORK(0:2*N-1)
JOBVL='N'
JOBVR='V'
call ZGEEV( JOBVL, JOBVR, N, A, N, W, VL, N, VR, &
            N, WORK, 2*N, RWORK, INFO )
D=0.0d0
if (F=='exp') then
    do i=0,n-1
        D(i,i)=exp(W(i))
    end do
end if

if (F=='log') then
    do i=0,n-1
        if (abs(W(i)) .le. 0.0000000000001) then
            print*, ' Error occured' 
            print*, 'Logarithm of zero and & 
                  &negative number are not defined'
            exit
        end if
        D(i,i)=log(W(i))
    end do
end if

if (F=='sin') then
    do i=0,n-1
        D(i,i)=sin(W(i))
    end do
end if

if (F=='cos') then
    do i=0,n-1
        D(i,i)=cos(W(i))
    end do
end if

if (F=='tan') then
    do i=0,n-1
        D(i,i)=tan(W(i))
    end do
end if
call MATINVMC(n,VR,B)
 E=matmul(matmul(VR,D),B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.5.1
subroutine POWNSR(A,n,k,B)
implicit none
real*8::A(0:n-1,0:n-1)
complex*16::D(0:n-1,0:n-1),im,B(0:n-1,0:n-1)
integer::n,i
character::JOBZ, UPLO  
integer::INFO
real*8::W(0:N-1),WORK(0:3*N-2),k
im=dcmplx(0.d0,1.0d0)
JOBZ='V'
UPLO='U'
call DSYEV( JOBZ, UPLO, N, A, N, W, WORK, 3*N-1, INFO )
D=0.0d0
do i=0,n-1
    if (abs(W(i))<0.0000000000001) then
        W(i)=0.0d0
    end if
    if (W(i) .lt. 0.0d0) then
        D(i,i)=(abs(W(i))**(k))*(im**(2.0d0*k))
    else
        D(i,i)=W(i)**k
    end if
end do
B=matmul(matmul(A,D),transpose(A))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.5.2
subroutine POWNHC(A,n,k,B)
implicit none
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),D(0:n-1,0:n-1),im
integer::n,i
character::JOBZ, UPLO   
integer::INFO
real*8::W(0:N-1),RWORK(0:3*N-3),k
complex*16::WORK(0:2*N-2)
im=dcmplx(0.d0,1.0d0)
JOBZ='V'
UPLO='U'
call ZHEEV( JOBZ, UPLO, N, A, N, W, WORK, 2*N-1, RWORK,INFO )
D=0.0d0
do i=0,n-1
    if (abs(W(i))<0.0000000000001) then
        W(i)=0.0d0
    end if
    if (W(i) .lt. 0.0d0) then
        D(i,i)=(abs(W(i))**(k))*(im**(2.0d0*k))
    else
        D(i,i)=W(i)**k
    end if
end do
B=matmul(matmul(A,D),conjg(transpose(A)))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.5.3
subroutine POWNMR(A,n,k,E)
implicit none
integer::n,i,j
complex*16::E(0:n-1,0:n-1),D(0:n-1,0:n-1)
complex*16::B(0:n-1,0:n-1),C(0:n-1,0:n-1)
character*3::JOBVL, JOBVR 
integer::INFO
complex*16::W(0:N-1),VL(0:N-1,0:N-1),VR(0:N-1,0:N-1)
complex*16::WORK(0:2*N-1)
real*8::RWORK(0:2*N-1),k,A(0:n-1,0:n-1)
JOBVL='N'
JOBVR='V'
do i=0,n-1
    do j=0,n-1
        C(i,j)=dcmplx(A(i,j),0.0d0)
    end do
end do
call ZGEEV( JOBVL, JOBVR, N, C, N, W, VL, N, &
     VR, N, WORK, 2*N, RWORK, INFO )
D=0.0d0
do i=0,n-1
    if (abs(W(i))<0.0000000000001) then
        W(i)=0
    end if
    D(i,i)=W(i)**(k)
end do
call MATINVMC(n,VR,B)
E=matmul(matmul(VR,D),B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.5.4
subroutine POWNMC(A,n,k,E)
implicit none
integer::n,i
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1)
complex*16::D(0:n-1,0:n-1),E(0:n-1,0:n-1)
character*3::JOBVL, JOBVR 
integer::INFO
complex*16::W(0:N-1),VL(0:N-1,0:N-1),VR(0:N-1,0:N-1)
complex*16::WORK(0:2*N-1)
real*8::RWORK(0:2*N-1),k
JOBVL='N'
JOBVR='V'
call ZGEEV( JOBVL, JOBVR, N, A, N, W, VL, N, VR, &
         N, WORK, 2*N, RWORK, INFO )
D=0.0d0
do i=0,n-1
if (abs(W(i))<0.0000000000001) then
W(i)=0
end if
D(i,i)=W(i)**(k)
end do
call MATINVMC(n,VR,B)
E=matmul(matmul(VR,D),B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.6.1
subroutine TRPOWNSR(A,n,k,trace)
implicit none
real*8::A(0:n-1,0:n-1)
complex*16::trace,im
integer::n,i
character::JOBZ, UPLO  
integer::INFO
real*8::W(0:N-1),WORK(0:3*N-2),k
JOBZ='N'
UPLO='U'
im=dcmplx(0.0d0,1.0d0)
call DSYEV( JOBZ, UPLO, N, A, N, W, WORK, 3*N-1, INFO )
trace=0.0d0
do i=0,n-1
    if (abs(W(i))<0.0000000000001) then
        W(i)=0.0d0
    end if
    if (W(i) .lt. 0.0d0) then
        trace=trace+(abs(W(i))**(k))*(im**(2.0d0*k))
    else
        trace=trace+W(i)**k
    end if
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.6.2
subroutine TRPOWNHC(A,n,k,trace)
implicit none
complex*16::A(0:n-1,0:n-1),trace,im
integer::n,i
character::JOBZ, UPLO
integer::INFO
real*8::W(0:N-1),RWORK(0:3*N-3),k
complex*16::WORK(0:2*N-2)
JOBZ='N'
UPLO='U'
im=dcmplx(0.0d0,1.0d0)
call ZHEEV( JOBZ, UPLO, N, A, N, W, WORK, 2*N-1, RWORK,INFO )
trace=0.0d0
do i=0,n-1
    if (abs(W(i))<0.0000000000001) then
        W(i)=0.0d0
    end if
    if (W(i) .lt. 0.0d0) then
        trace=trace+(abs(W(i))**(k))*(im**(2.0d0*k))
    else
        trace=trace+W(i)**k
    end if
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.6.3
subroutine TRPOWNMR(A,n,k,trace)
implicit none
integer::n,i
real*8::k
complex*16::trace,o
character::JOBVL, JOBVR   
integer::INFO
real*8::WR(0:N-1),WI(0:N-1),VL(0:N-1,0:N-1)
real*8::VR(0:N-1,0:N-1),WORK(0:4*N-1),A(0:N-1,0:N-1)
JOBVL='N'
JOBVR='N'
call DGEEV( JOBVL, JOBVR, N, A, N, WR, WI,   &
                   VL, N, VR,N, WORK, 4*N, INFO )
trace=0.0d0
do i=0,n-1,1
    o=dcmplx(WR(i),WI(i))
    if (abs(o)<0.0000000000001) then
        o=0.0d0
    end if
trace=trace+o**k
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.6.4
subroutine TRPOWNMC(A,n,k,trace)
implicit none
integer::n,i
real*8::k
complex*16::trace,o
character::JOBVL, JOBVR  
integer::INFO
complex*16::W(0:N-1),VL(0:N-1,0:N-1)
complex*16::VR(0:N-1,0:N-1),WORK(0:2*N-1),A(0:N-1,0:N-1)
real*8::RWORK(0:2*N-1)
JOBVL='N'
JOBVR='N'
call ZGEEV( JOBVL, JOBVR, N, A, N, W, VL, N, &
                   VR, N, WORK, 2*N, RWORK, INFO )
trace=0.0d0
do i=0,n-1,1
    if (abs(W(i))<0.0000000000001) then
        W(i)=0.0d0
    end if
o=W(i)
trace=trace+o**k
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.7.1
subroutine DETSR(A,n,det)
implicit none
real*8::A(0:n-1,0:n-1),det
integer::n,i
character::JOBZ, UPLO  
integer::INFO
real*8::W(0:N-1),WORK(0:3*N-2)
JOBZ='N'
UPLO='U'
call DSYEV( JOBZ, UPLO, N, A, N, W, WORK, 3*N-1, INFO )
det=1.0d0
do i=0,n-1
    det=det*W(i) 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.7.2
subroutine DETHC(A,n,det)
implicit none
complex*16::A(0:n-1,0:n-1)
integer::n,i,INFO
character::JOBZ, UPLO   
real*8::W(0:N-1),RWORK(0:3*N-3),det
complex*16::WORK(0:2*N-2)
JOBZ='N'
UPLO='U'
call ZHEEV( JOBZ, UPLO, N, A, N, W, WORK, 2*N-1, RWORK,INFO )
det=1.0d0
do i=0,n-1
    det=det*W(i) 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.7.3
subroutine DETMR(A,n,det)
implicit none
integer::n,i
real*8::det
complex*16::det1,o
character::JOBVL, JOBVR  
integer::INFO
real*8::WR(0:N-1),WI(0:N-1),VL(0:N-1,0:N-1)
real*8::VR(0:N-1,0:N-1),WORK(0:4*N-1),A(0:N-1,0:N-1)
JOBVL='N'
JOBVR='N'
call DGEEV( JOBVL, JOBVR, N, A, N, WR, WI, VL, &
                 N, VR,N, WORK, 4*N, INFO )
det1=dcmplx(1.0d0,0.d0)
do i=0,n-1,1
    o=dcmplx(WR(i),WI(i))
    det1=det1*o
end do
det=real(det1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.7.4
subroutine DETMC(A,n,det)
implicit none
integer::n,i
complex*16::det
character::JOBVL, JOBVR  
integer::INFO
complex*16::W(0:N-1),VL(0:N-1,0:N-1)
complex*16::VR(0:N-1,0:N-1),WORK(0:2*N-1),A(0:N-1,0:N-1)
real*8::RWORK(0:2*N-1)
JOBVL='N'
JOBVR='N'
call ZGEEV( JOBVL, JOBVR, N, A, N, W, VL, &
                          N, VR, N, WORK, 2*N, RWORK, INFO )
det=dcmplx(1.0d0,0.d0)
do i=0,n-1,1
    det=det*W(i)
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.8.1
subroutine TRNORMMR(A,n,trace_norm)
implicit none
integer::n
real*8::A(0:N-1,0:N-1),B(0:N-1,0:N-1),trace_norm
complex*16::trace
B=matmul(transpose(A),A)
call TRPOWNSR(B,n,0.5d0,trace)
trace_norm=real(trace)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.8.2
subroutine TRNORMMC(A,n,trace_norm)
implicit none
integer::n
real*8::trace_norm
complex*16::A(0:N-1,0:N-1),B(0:N-1,0:N-1),trace
B=matmul(conjg(transpose(A)),A)
call TRPOWNHC(B,n,0.5d0,trace)
trace_norm=real(trace)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.9.1
subroutine HISCNORMMR(A,n,HS_norm)
implicit none
integer::n
real*8::trace,HS_norm,A(0:n-1,0:n-1),B(0:n-1,0:n-1)
B=matmul(transpose(A),A)
call TRMR(B,n,trace)
HS_norm=sqrt(trace)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.9.2
subroutine HISCNORMMC(A,n,HS_norm)
implicit none
integer::n
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),trace
real*8::HS_norm
B=matmul(conjg(transpose(A)),A)
call TRMC(B,n,trace)
HS_norm=sqrt(abs(trace))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.10.1
subroutine ABSMR(A,n,B)
implicit none
real*8::A(0:n-1,0:n-1),C(0:n-1,0:n-1),B(0:n-1,0:n-1)
complex*16::D(0:n-1,0:n-1)
integer::n
C=matmul(transpose(A),A)
call POWNSR(C,n,0.5d0,D)
B=real(D)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.10.2
subroutine ABSMC(A,n,B)
implicit none
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),C(0:N-1,0:N-1)
integer::n
C=matmul(conjg(transpose(A)),A)
call POWNHC(C,n,0.5d0,B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.11.1
subroutine HISCIPMR(n,A,B,HS_inp)
integer::n
real*8::A(0:n-1,0:n-1),HS_inp,A2(0:n-1,0:n-1),B(0:n-1,0:n-1)
A2=matmul(transpose(A),B)
call TRMR(A2,n,HS_inp)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.11.2
subroutine HISCIPMC(n,A,B,HS_inp)
integer::n
complex*16::A(0:n-1,0:n-1),HS_inp,A2(0:n-1,0:n-1),B(0:n-1,0:n-1)
A2=matmul(transpose(conjg(A)),B)
call TRMC(A2,n,HS_inp)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.12.1
subroutine GSOMR(n,k,wi,vi)
implicit none
integer::n,k,i,j
real*8::wi(0:n-1,0:k-1),vi(0:n-1,0:k-1),c
real*8,dimension(0:n-1)::u1,u2,u3,u4,u5
vi=0.0d0
u1=wi(:,0)
call NORMALIZEVR(u1,n)
vi(:,0)=u1
do i=1,k-1
    u2=wi(:,i)
    u3=0.0d0
    u4=0.0d0
    do j=0,i-1
        u5=vi(:,j)
        call IPVR(u5,u2,n,c)
        u3=u3+c*u5
    end do   
    u4=u2-u3
    call NORMALIZEVR(u4,n)
    vi(:,i)=u4
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.12.2
subroutine GSOMC(n,k,wi,vi)
implicit none
integer::n,k,i,j
complex*16::wi(0:n-1,0:k-1),vi(0:n-1,0:k-1),c
complex*16,dimension(0:n-1)::u1,u2,u3,u4,u5
vi=0.0d0
u1=wi(:,0)
call NORMALIZEVC(u1,n)
vi(:,0)=u1
do i=1,k-1
    u2=wi(:,i)
    u3=0.0d0
    u4=0.0d0
    do j=0,i-1
        u5=vi(:,j)
        call IPVC(u5,u2,n,c)
        u3=u3+c*u5
    end do   
    u4=u2-u3
    call NORMALIZEVC(u4,n)
    vi(:,i)=u4
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.13.1
subroutine SVDMR(A,M,N,U,Si,VT)
integer::M,N
integer::LWORK
real*8::A(0:m-1,0:n-1),S(0:min(M,N)-1),U(0:M-1,0:M-1)
real*8::VT(0:N-1,0:N-1),B(0:m-1,0:n-1),si(0:m-1,0:n-1)
character::JOBU, JOBVT
real*8,allocatable,dimension(:)::WORK
integer::INFO
JOBU='A'
JOBVT='A'
LWORK=MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
allocate(WORK(0:LWORK-1))
call  DGESVD( JOBU, JOBVT, M, N, A, M, S, &
              U, M, VT, N,WORK, LWORK, INFO )
si=0
do i=0,min(M,N)
si(i,i)=S(i)
end do
deallocate(WORK)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!3.13.2
subroutine SVDMC(A,M,N,U,Si,VT)
integer::m,n
integer::LWORK
real*8::S(0:min(M,N)-1),si(0:m-1,0:n-1),RWORK(0:5*min(M,N)-1)
complex*16::A(0:m-1,0:n-1),U(0:M-1,0:M-1),VT(0:N-1,0:N-1)
complex*16::B(0:m-1,0:n-1)
complex*16,allocatable,dimension(:)::WORK
character::JOBU, JOBVT
integer::INFO
JOBU='A'
JOBVT='A'
LWORK=MAX(1,2*MIN(M,N)+MAX(M,N))
allocate(WORK(0:LWORK-1))
call ZGESVD( JOBU, JOBVT, M, N, A, M, S, &
                      U, M, VT, N,WORK, LWORK, RWORK, INFO )
si=0
do i=0,min(M,N)
si(i,i)=S(i)
end do
deallocate(WORK)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHAPTER 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.1
subroutine PAULIX(sigmax)
real*8::sigmax(0:1,0:1)
sigmax(0,:)=(/0,1/)
sigmax(1,:)=(/1,0/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.2
subroutine PAULIY(sigmay)
complex*16::sigmay(0:1,0:1)
sigmay(0,0)=dcmplx(0,0)
sigmay(0,1)=dcmplx(0,-1)
sigmay(1,0)=dcmplx(0,1)
sigmay(1,1)=dcmplx(0,0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.3
subroutine PAULIZ(sigmaz)
real*8::sigmaz(0:1,0:1)
sigmaz(0,:)=(/1,0/)
sigmaz(1,:)=(/0,-1/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.4
subroutine HADAMARDG(H)
real*8::H(0:1,0:1),ele
ele=1/dsqrt(2.0d0)
H(0,:)=(/ele,ele/)
H(1,:)=(/ele,-ele/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.5
subroutine PHASEG(S)
complex*16::S(0:1,0:1)
S(0,0)=dcmplx(1,0)
S(0,1)=dcmplx(0,0)
S(1,0)=dcmplx(0,0)
S(1,1)=dcmplx(0,1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.6
subroutine ROTG(R,k)
integer::k
real*8::pi
complex*16::R(0:1,0:1),i1,ele
pi=4*atan(1.0d0)
i1=dcmplx(0,1)
ele=exp((2*pi*i1)/(2**k))
R(0,0)=dcmplx(1,0)
R(0,1)=dcmplx(0,0)
R(1,0)=dcmplx(0,0)
R(1,1)=ele
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.7
subroutine CONOTG(A)
real*8::A(0:3,0:3)
A(0,:)=(/1,0,0,0/)
A(1,:)=(/0,1,0,0/)
A(2,:)=(/0,0,0,1/)
A(3,:)=(/0,0,1,0/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.8
subroutine COZG(A)
real*8::A(0:3,0:3)
A(0,:)=(/1,0,0,0/)
A(1,:)=(/0,1,0,0/)
A(2,:)=(/0,0,1,0/)
A(3,:)=(/0,0,0,-1/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.9
subroutine SWAPG(A)
real*8::A(0:3,0:3)
A(0,:)=(/1,0,0,0/)
A(1,:)=(/0,0,1,0/)
A(2,:)=(/0,1,0,0/)
A(3,:)=(/0,0,0,1/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.10
subroutine TOFFOLIG(A)
real*8::A(0:7,0:7)
A(0,:)=(/1,0,0,0,0,0,0,0/)
A(1,:)=(/0,1,0,0,0,0,0,0/)
A(2,:)=(/0,0,1,0,0,0,0,0/)
A(3,:)=(/0,0,0,1,0,0,0,0/)
A(4,:)=(/0,0,0,0,1,0,0,0/)
A(5,:)=(/0,0,0,0,0,1,0,0/)
A(6,:)=(/0,0,0,0,0,0,0,1/)
A(7,:)=(/0,0,0,0,0,0,1,0/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.1.11
subroutine FREDKING(A)
real*8::A(0:7,0:7)
A(0,:)=(/1,0,0,0,0,0,0,0/)
A(1,:)=(/0,1,0,0,0,0,0,0/)
A(2,:)=(/0,0,1,0,0,0,0,0/)
A(3,:)=(/0,0,0,1,0,0,0,0/)
A(4,:)=(/0,0,0,0,1,0,0,0/)
A(5,:)=(/0,0,0,0,0,0,1,0/)
A(6,:)=(/0,0,0,0,0,1,0,0/)
A(7,:)=(/0,0,0,0,0,0,0,1/)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.2
subroutine NQBHADAMARDG(s,A)
implicit none
integer::i,j,s,t1(0:s-1),t2(0:s-1),k,xy
real*8::A(0:2**s-1,0:2**s-1)
A=0.0d0
do i=0,2**s-1
    do j=0,2**s-1
        call DTOB(i,t1,s)
        call DTOB(j,t2,s)
        xy=0
        do k=0,s-1
            xy=xy+t1(k)*t2(k)
        end do
        A(j,i)=(1/dsqrt(dfloat(2**s)))*((-1)**(xy))
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.3
subroutine NQBCB(n,II)
implicit none
integer::n,j
real*8::II(0:n-1,0:n-1)
do j=0,n-1,1
    II(j,j)=1.0d0
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.4
subroutine BELL(b1,b2,b3,b4)
implicit none
real*8::II(0:3,0:3),b1(0:3),b2(0:3),b3(0:3),b4(0:3)
II(0,:)=(/1,0,0,0/)
II(1,:)=(/0,1,0,0/)
II(2,:)=(/0,0,1,0/)
II(3,:)=(/0,0,0,1/)
b1=(II(:,0)+II(:,3))/dsqrt(2.0d0)
b2=(II(:,0)-II(:,3))/dsqrt(2.0d0)
b3=(II(:,1)+II(:,2))/dsqrt(2.0d0)
b4=(II(:,1)-II(:,2))/dsqrt(2.0d0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.5
subroutine NQBGHZ(s,psi)
implicit none
integer::s
real*8::psi(0:2**s-1)
psi=0.0d0
psi(0)=1.0d0/dsqrt(2.0d0)
psi(2**s-1)=1.0d0/dsqrt(2.0d0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.6
subroutine NQBWSTATE(s,psi)
integer::s,m
real*8::psi(0:2**s-1)
integer::state(0:s-1)
psi=0.0d0
state=0
state(s-1)=1
do i=0,s-1
    call RIGHTSHIFT(state,s,1)
    call BTOD(state,m,s)
    psi(m)=1.0/dsqrt(dfloat(s))
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.7
subroutine NQBWERNER(s,x,rho)
implicit none
integer::s,i,j
real*8,dimension(0:2**s-1,0:2**s-1)::rho,rho1,II
real*8::psi(0:2**s-1),x
call NQBGHZ(s,psi)
call OPVR(psi,psi,2**s,2**s,rho1)
call NQBCB(2**s,II)
rho=x*rho1+((1-x)/2**s)*II
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.8
subroutine SE(p,sp,n)
implicit none
integer::n,i
real*8::sp,p(0:n-1),sptemp
sp=0
do i=0,n-1,1
    sptemp=-p(i)*dlog(p(i))/dlog(2.0d0)
    sp=sp+sptemp 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.9.1
subroutine LEDMR(A,n,le)
implicit none
integer::i,n
real*8::A(0:n-1,0:n-1),le
complex*16::trace
call TRPOWNSR(A,n,2.0d0,trace)
le=1-real(trace)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.9.2
subroutine LEDMC(A,n,le)
implicit none
integer::i,n
real*8::le
complex*16::A(0:n-1,0:n-1),trace
call TRPOWNHC(A,n,2.0d0,trace)
le=1-real(trace)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.10.1
subroutine REDMR(A,B,n,re)
integer::n,i
real*8::A(0:n-1,0:n-1),B(0:n-1,0:n-1),re
real*8::logA(0:n-1,0:n-1),logB(0:n-1,0:n-1)
real*8::temp(0:n-1,0:n-1),C(0:n-1,0:n-1),D(0:n-1,0:n-1)
 C=A
D=B
call FUNSR(C,n,'log',logA)
call FUNSR(D,n,'log',logB)
temp=matmul(A,logA)-matmul(A,logB)
call TRMR(temp,n,re)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.10.2
subroutine REDMC(A,B,n,re)
integer::n,i
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),re1
complex*16::logA(0:n-1,0:n-1),logB(0:n-1,0:n-1)
complex*16::temp(0:n-1,0:n-1),C(0:n-1,0:n-1),D(0:n-1,0:n-1)
real*8::re
C=A
D=B
call FUNHC(C,n,'log',logA)
call FUNHC(D,n,'log',logB)
temp=matmul(A,logA)-matmul(A,logB)
call TRMC(temp,n,re1)
re=real(re1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.11.1
subroutine TRDISDMR(A,B,n,T)
implicit none
integer::n
real*8::A(0:n-1,0:n-1),B(0:n-1,0:n-1),T
real*8::temp(0:n-1,0:n-1),trace
temp=A-B
call TRNORMMR(temp,n,trace)
T=trace/2.0d0
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.11.2
subroutine TRDISDMC(A,B,n,T)
implicit none
integer::n,i
complex*16::A(0:n-1,0:n-1),B(0:n-1,0:n-1),T
complex*16::temp(0:n-1,0:n-1)
real*8::trace
temp=A-B
call TRNORMMC(temp,n,trace)
T=trace/2.0d0
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.12.1
subroutine FIDVR(psi1,psi2,n,F)
implicit none
integer::n,i
real*8::psi1(0:n-1),psi2(0:n-1),F
call IPVR(psi1,psi2,n,F)
F=F*F
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.12.2
subroutine FIDVC(psi1,psi2,n,F)
implicit none
integer::n,i
real*8::F
complex*16::psi1(0:n-1),psi2(0:n-1),F1
call IPVC(psi1,psi2,n,F1)
F1=F1*conjg(F1)
F=real(F1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.12.3
subroutine FIDDMR(rho1,rho2,n,F)
implicit none
integer::n
real*8::rho1(0:n-1,0:n-1),rho2(0:n-1,0:n-1),temp(0:n-1,0:n-1)
complex*16::rtrho1(0:n-1,0:n-1),rtrho2(0:n-1,0:n-1)
complex*16::temp1(0:n-1,0:n-1)
real*8::F,rho3(0:n-1,0:n-1),rho4(0:n-1,0:n-1)
real*8::pow
pow=0.5d0
rho3=rho1
rho4=rho2
call POWNSR(rho3,n,pow,rtrho1)
call POWNSR(rho4,n,pow,rtrho2)
temp1=matmul(rtrho1,rtrho2)
temp=real(temp1)
call TRNORMMR(temp,n,F)
F=F**2
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.12.4
subroutine FIDDMC(rho1,rho2,n,F)
implicit none
complex*16::rho1(0:n-1,0:n-1),rho2(0:n-1,0:n-1)
complex*16::rtrho1(0:n-1,0:n-1),rtrho2(0:n-1,0:n-1)
complex*16::temp(0:n-1,0:n-1)
complex*16::rho3(0:n-1,0:n-1),rho4(0:n-1,0:n-1)
integer::n
real*8::pow,F
rho3=rho1
rho4=rho2
call POWNHC(rho3,n,0.5d0,rtrho1)
call POWNHC(rho4,n,0.5d0,rtrho2)
temp=matmul(rtrho1,rtrho2)
call TRNORMMC(temp,n,F)
F=F**2
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.12.5
subroutine FIDVRDMR(n,rho,psi,F)
implicit none
integer::n,i
real*8::rho(0:n-1,0:n-1),psi(0:n-1),F,temp(0:n-1)
temp=matmul(rho,psi)
call IPVR(temp,psi,n,F)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.12.6
subroutine FIDVCDMC(n,rho,psi,F)
implicit none
integer::n,i
complex*16::rho(0:n-1,0:n-1),psi(0:n-1),F1,temp(0:n-1)
real*8::F
temp=matmul(rho,psi)
call IPVC(temp,psi,n,F1)
F=real(F1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.13.1
subroutine SUPFIDDMR(rho1,rho2,n,F)
implicit none
integer::n
real*8::rho1(0:n-1,0:n-1),rho2(0:n-1,0:n-1)
real*8::a1(0:n-1,0:n-1),a2(0:n-1,0:n-1),a3(0:n-1,0:n-1)
real*8::F,trace1,trace2,trace3,temp1,temp2
a1=matmul(rho1,rho2)
a2=matmul(rho1,rho1)
a3=matmul(rho2,rho2)
call TRMR(a1,n,trace1)
call TRMR(a2,n,trace2)
call TRMR(a3,n,trace3)
temp1=1-trace2
if (abs(temp1) .le. 0.0000000000001) then
    temp1=0.0d0
end if
temp2=1-trace3
if (abs(temp2) .le. 0.0000000000001) then
    temp2=0.0d0
end if
F=trace1+dsqrt(temp1)*dsqrt(temp2)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.13.2
subroutine SUPFIDDMC(rho1,rho2,n,F)
implicit none
integer::n
complex*16::rho1(0:n-1,0:n-1),rho2(0:n-1,0:n-1)
complex*16::a1(0:n-1,0:n-1),a2(0:n-1,0:n-1),a3(0:n-1,0:n-1)
complex*16::trace1,trace2,trace3
real*8::F,temp1,temp2
a1=matmul(rho1,rho2)
a2=matmul(rho1,rho1)
a3=matmul(rho2,rho2)
call TRMC(a1,n,trace1)
call TRMC(a2,n,trace2)
call TRMC(a3,n,trace3)
temp1=1-real(trace2)
if (abs(temp1) .le. 0.0000000000001) then
    temp1=0.0d0
end if
temp2=1-real(trace3)
if (abs(temp2) .le. 0.0000000000001) then
    temp2=0.0d0
end if
F=real(trace1)+dsqrt(temp1)*dsqrt(temp2)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.14.1
subroutine BURDISVR(psi1,psi2,n,D)
implicit none
integer::n
real*8::psi1(0:n-1),psi2(0:n-1),D
real*8::F,temp
call FIDVR(psi1,psi2,n,F)
temp=1-dsqrt(F)
If (temp .le. 0.0000000000001) then
temp=0.0d0
end if
D=dsqrt(2*(temp))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.14.2
subroutine BURDISVC(psi1,psi2,n,D)
implicit none
integer::n
complex*16::psi1(0:n-1),psi2(0:n-1)
real*8::D,F,temp
call FIDVC(psi1,psi2,n,F)
temp=1-dsqrt(F)
If (temp .le. 0.0000000000001) then
temp=0.0d0
end if
D=dsqrt(2*(temp))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.14.3
subroutine BURDISDMR(rho1,rho2,n,D)
implicit none
integer::n
real*8::rho1(0:n-1,0:n-1),rho2(0:n-1,0:n-1),D
real*8::F,temp
call FIDDMR(rho1,rho2,n,F)
temp=1-dsqrt(F)
If (temp .le. 0.0000000000001) then
temp=0.0d0
end if
D=dsqrt(2*(temp))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.14.4
subroutine BURDISDMC(rho1,rho2,n,D)
implicit none
integer::n
complex*16::rho1(0:n-1,0:n-1),rho2(0:n-1,0:n-1)
real*8::D,F,temp
call FIDDMC(rho1,rho2,n,F)
temp=1-dsqrt(F)
If (temp .le. 0.0000000000001) then
temp=0.0d0
end if
D=dsqrt(2*(temp))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.15.1
subroutine AVGVRMR(n,A,psi,val)
implicit none
integer::n,i,j
real*8::A(0:n-1,0:n-1),rho(0:n-1,0:n-1)
real*8::val,temp(0:n-1,0:n-1),psi(0:n-1)
call OPVR(psi,psi,n,n,rho)
temp=matmul(A,rho)
call TRMR(temp,n,val)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.15.2
subroutine AVGVCMC(n,A,psi,val)
implicit none
integer::n,i,j
complex*16::A(0:n-1,0:n-1),rho(0:n-1,0:n-1)
complex*16::val1,temp(0:n-1,0:n-1),psi(0:n-1)
real*8::val
call OPVC(psi,psi,n,n,rho)
temp=matmul(A,rho)
call TRMC(temp,n,val1)
val=real(val1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.15.3
subroutine AVGDMRMR(n,A,rho,val)
implicit none
integer::n,i
real*8::A(0:n-1,0:n-1),rho(0:n-1,0:n-1),val,temp(0:n-1,0:n-1)
temp=matmul(rho,A)
call TRMR(temp,n,val)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.15.4
subroutine AVGDMCMC(n,A,rho,val)
implicit none
integer::n,i
complex*16::A(0:n-1,0:n-1),rho(0:n-1,0:n-1),temp(0:n-1,0:n-1)
complex*16::val1
real*8::val
temp=matmul(rho,A)
call TRMC(temp,n,val1)
val=real(val1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.16.1
subroutine MEASUREVR(s,n,psi,o,psim,pm)
implicit none
integer::s,n,i,j
character::o
real*8::psi(0:2**s-1),psim(0:2**s-1),pm,id(0:2**s-1,0:2**s-1)
id=0.0d0
if (o=='0') then
    do i=0,2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
if (o=='1') then
    do i=2**(s-n),2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
psim=matmul(id,psi)
call IPVR(psim,psim,2**s,pm)
psim=psim/sqrt(pm)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.16.2
subroutine MEASUREVC(s,n,psi,o,psim,pm)
implicit none
integer::s,n,i,j
character::o
complex*16::psi(0:2**s-1),psim(0:2**s-1),pm1
real*8::pm,id(0:2**s-1,0:2**s-1)
id=0.0d0
if (o=='0') then
    do i=0,2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
if (o=='1') then
    do i=2**(s-n),2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
psim=matmul(id,psi)
call IPVC(psim,psim,2**s,pm1)
pm=real(pm1)
psim=psim/sqrt(pm)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.16.3
subroutine MEASUREDMR(s,n,rho,o,rhom,pm)
implicit none
integer::s,n,i,j
character::o
real*8::rho(0:2**s-1,0:2**s-1),rhom(0:2**s-1,0:2**s-1),pm
real*8::id(0:2**s-1,0:2**s-1)
id=0.0d0
if (o=='0') then
    do i=0,2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
if (o=='1') then
    do i=2**(s-n),2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
rhom=matmul(id,matmul(rho,transpose(id)))
call TRMR(rhom,2**s,pm)
rhom=rhom/pm
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!4.16.4
subroutine MEASUREDMC(s,n,rho,o,rhom,pm)
implicit none
integer::s,n,i,j
character::o
complex*16::rho(0:2**s-1,0:2**s-1),rhom(0:2**s-1,0:2**s-1),pm1
real*8::id(0:2**s-1,0:2**s-1),pm
id=0.0d0
if (o=='0') then
    do i=0,2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
if (o=='1') then
    do i=2**(s-n),2**s-1,2**(s-n+1)
        do j=i,i+2**(s-n)-1
            id(j,j)=1.0d0
        end do
    end do
end if
rhom=matmul(id,matmul(rho,transpose(id)))
call TRMC(rhom,2**s,pm1)
rhom=rhom/pm1
pm=real(pm1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHAPTER 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.1.2
subroutine  PTRVR(s,s1,site,vin,rdm)
implicit none
integer::s, s1,s2
real*8::trace
real*8,dimension(0:2**s-1)::vin
real*8,dimension(0:2**s1-1,0:2**s1-1)::rdm
integer::ii,i1,a,i0,ia,j1,oo,k1,x1,y1,j,t1,t2,z1,i
integer,dimension(0:s1-1)::site,site2
integer,dimension(0:s1-1)::bin
integer,dimension(0:(2**(s-s1))*(2**s1)-1)::ind
s2=s-s1
x1=0
y1=0
ind=0
do j1=2**s1-1,0,-1
        do ii=0,2**s-1
                   a=ii  
                   do i1=0,s1-1,1
                       i0=site(i1)
                       call DTOBONEBIT(a,ia,i0,s)
                       site2(i1)=ia
                   enddo
                   call DTOB(j1,bin,s1)
                   oo=1
                   do k1=0,s1-1,1
                   oo=oo*(bin(k1)-site2(k1))   
                   end do
                   if(abs(oo)==1) then
                      ind(x1)=ii
                        x1=x1+1
                   end if
         end do
end do

rdm=0.0d0
do t1=0,2**s1-1,1
   do t2=0,2**s1-1,1
      do z1=0,2**s2-1,1
         rdm(t1,t2)=rdm(t1,t2)+vin(ind((2**s2)*t1+z1)) & 
                            *vin(ind((2**s2)*t2+z1))
      enddo
   end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.1.3
subroutine  PTRVC(s,s1,site,vin,rdm)
implicit none
integer::s, s1,s2
real*8::trace
complex*16,dimension(0:2**s-1)::vin
complex*16,dimension(0:2**s1-1,0:2**s1-1)::rdm
integer::ii,i1,a,i0,ia,j1,oo,k1,x1,y1,j,t1,t2,z1,i
integer,dimension(0:s1-1)::site,site2
integer,dimension(0:s1-1)::bin
integer,dimension(0:(2**(s-s1))*(2**s1)-1)::ind
s2=s-s1
x1=0
y1=0
ind=0

do j1=2**s1-1,0,-1
        do ii=0,2**s-1
                   a=ii  
                   do i1=0,s1-1,1
                       i0=site(i1)
                       call DTOBONEBIT(a,ia,i0,s)
                       site2(i1)=ia
                   enddo
                   call DTOB(j1,bin,s1)
                   oo=1
                   do k1=0,s1-1,1
                   oo=oo*(bin(k1)-site2(k1))   
                   end do
                   if(abs(oo)==1) then
                      ind(x1)=ii
                        x1=x1+1
                   end if
         end do
end do
rdm=0.0d0
do t1=0,2**s1-1,1
   do t2=0,2**s1-1,1
      do z1=0,2**s2-1,1
         rdm(t1,t2)=rdm(t1,t2)+vin(ind((2**s2)*t1+z1))* &
                          conjg(vin(ind((2**s2)*t2+z1)))
      enddo
   end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.1.4
subroutine PTRDMR(s,s1,site,rho,rdm)
implicit none
integer::s, s1,s2
real*8::trace
real*8,dimension(0:2**s-1,0:2**s-1)::rho
real*8,dimension(0:2**s1-1,0:2**s1-1)::rdm
integer::ii,i1,a,i0,ia,j1,oo,k1,x1,y1,j,t1,t2,z1,i
integer,dimension(0:s1-1)::site,site2
integer,dimension(0:s1-1)::bin
integer,dimension(0:(2**(s-s1))*(2**s1)-1)::ind
s2=s-s1
x1=0
y1=0
ind=0
do j1=2**s1-1,0,-1
        do ii=0,2**s-1
                   a=ii  
                   do i1=0,s1-1,1
                       i0=site(i1)
                       call DTOBONEBIT(a,ia,i0,s)
                       site2(i1)=ia
                   enddo
                   call DTOB(j1,bin,s1)
                   oo=1
                   do k1=0,s1-1,1
                   oo=oo*(bin(k1)-site2(k1))   
                   end do
                   if(abs(oo)==1) then
                      ind(x1)=ii
                        x1=x1+1
                   end if
         end do
end do
rdm=0.0d0
do t1=0,2**s1-1,1
   do t2=0,2**s1-1,1
      do z1=0,2**s2-1,1
         rdm(t1,t2)=rdm(t1,t2)+ &
            rho(ind((2**s2)*t1+z1),ind((2**s2)*t2+z1))
      enddo
   end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.1.5
subroutine PTRDMC(s,s1,site,rho,rdm)
implicit none
integer::s, s1,s2
real*8::trace
complex*16,dimension(0:2**s-1,0:2**s-1)::rho
complex*16,dimension(0:2**s1-1,0:2**s1-1)::rdm
integer::ii,i1,a,i0,ia,j1,oo,k1,x1,y1,j,t1,t2,z1,i
integer,dimension(0:s1-1)::site,site2
integer,dimension(0:s1-1)::bin
integer,dimension(0:(2**(s-s1))*(2**s1)-1)::ind
s2=s-s1
x1=0
y1=0
ind=0

do j1=2**s1-1,0,-1
        do ii=0,2**s-1
                   a=ii  
                   do i1=0,s1-1,1
                       i0=site(i1)
                       call DTOBONEBIT(a,ia,i0,s)
                       site2(i1)=ia
                   enddo
                   call DTOB(j1,bin,s1)
                   oo=1
                   do k1=0,s1-1,1
                   oo=oo*(bin(k1)-site2(k1))   
                   end do
                   if(abs(oo)==1) then
                      ind(x1)=ii
                        x1=x1+1
                   end if
         end do
end do
rdm=0.0d0
do t1=0,2**s1-1,1
   do t2=0,2**s1-1,1
      do z1=0,2**s2-1,1
         rdm(t1,t2)=rdm(t1,t2)+ &
              rho(ind((2**s2)*t1+z1),ind((2**s2)*t2+z1))
      enddo
   end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.2.1
subroutine PTDMR(s,s1,site,rho,rhotp)
implicit none
integer::s,i1,j1,i,j,oo,f,flag,t2,t3,t1,x1,s1
real*8::t4,t5,trace
integer,dimension(0:s1-1)::site
real*8,dimension(0:2**s-1,0:2**s-1):: rho,rhotp,rhot
integer,dimension(0:2**s-1)::ind1,ind2
integer,dimension(0:s-1)::m1,m2,m4,m5
integer,dimension(0:2*s-1)::m3
rhotp=rho
ind1=0
ind2=0
oo=0
do i=0,2**s-1,1
do j=0,2**s-1,1
f=0
    do j1=0,2**s-1
         if ( i==ind1(j1) .and. j==ind2(j1) ) then
              flag=1
         else
              flag=0
         end if
      f=f+flag
    end do

if(f==0 .and. i .ne. j) then 
     call DTOB(i,m1,s)
     call DTOB(j,m2,s)

    do i1=0,2*s/2-1,1
       m3(i1)=m1(i1)
    end do

    do i1=0,2*s/2-1,1
      m3(2*s/2+i1)=m2(i1)
    end do

    do x1=0,s1-1
       t1=m3(site(x1))
       m3(site(x1))=m3(site(x1)+s)
       m3(site(x1)+s)=t1
    end do

    do i1=0,2*s/2-1,1
       m4(i1)=m3(i1)
    end do
    do i1=0,2*s/2-1,1
       m5(i1)=m3(2*s/2+i1)
    end do

call BTOD(m4,t2,s)
call BTOD(m5,t3,s)

  t4=rho(t2,t3)
  t5=rho(i,j)
  rhotp(i,j)=t4
  rhotp(t2,t3)=t5
  ind1(oo)=t2
  ind2(oo)=t3

end if

end do
oo=oo+1
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.2.2
subroutine PTDMC(s,s1,site,rho,rhotp)
implicit none
integer::s,i1,j1,i,j,oo,f,flag,t2,t3,t1,x1,s1
complex*16::t4,t5
real*8::trace
integer,dimension(0:s1-1)::site
complex*16,dimension(0:2**s-1,0:2**s-1):: rho,rhotp
integer,dimension(0:2**s-1)::ind1,ind2
integer,dimension(0:s-1)::m1,m2,m4,m5
integer,dimension(0:2*s-1)::m3
rhotp=rho
ind1=0
ind2=0
oo=0
do i=0,2**s-1,1
do j=0,2**s-1,1
       f=0
    do j1=0,2**s-1
         if ( i==ind1(j1) .and. j==ind2(j1) ) then
              flag=1
         else
              flag=0
         end if
      f=f+flag
    end do

if(f==0 .and. i .ne. j) then 
     call DTOB(i,m1,s)
     call DTOB(j,m2,s)

    do i1=0,2*s/2-1,1
       m3(i1)=m1(i1)
    end do

    do i1=0,2*s/2-1,1
      m3(2*s/2+i1)=m2(i1)
    end do

    do x1=0,s1-1
       t1=m3(site(x1))
       m3(site(x1))=m3(site(x1)+s)
       m3(site(x1)+s)=t1
    end do

    do i1=0,2*s/2-1,1
       m4(i1)=m3(i1)
    end do
    do i1=0,2*s/2-1,1
       m5(i1)=m3(2*s/2+i1)
    end do

call BTOD(m4,t2,s)
call BTOD(m5,t3,s)

  t4=rho(t2,t3)
  t5=rho(i,j)
  rhotp(i,j)=t4
  rhotp(t2,t3)=t5
  ind1(oo)=t2
  ind2(oo)=t3

end if

end do
oo=oo+1
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.2.3
subroutine PTVR(s,s1,site,vin,rhotp)
implicit none
integer::s,s1,i,j
integer,dimension(0:s1-1)::site
real*8,dimension(0:2**s-1)::vin
real*8,dimension(0:2**s-1,0:2**s-1):: rho,rhotp
call OPVR(vin,vin,2**s,2**s,rho)
call PTDMR(s,s1,site,rho,rhotp)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.2.4
subroutine PTVC(s,s1,site,vin,rhotp)
implicit none
integer::s,s1,i,j
integer,dimension(0:s1-1)::site
complex*16,dimension(0:2**s-1,0:2**s-1):: rho,rhotp
complex*16,dimension(0:2**s-1)::vin
call OPVC(vin,vin,2**s,2**s,rho)
call PTDMC(s,s1,site,rho,rhotp)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.3.1
subroutine CONVR(s,i0,j0,vin,con)
implicit none
integer::s,i0,j0 ,site(0:1),i
real*8::vin(0:2**s-1),evs
real*8,dimension(0:3,0:3)::rdm,MM1,MM2,MM3,MM4
INTEGER::INFO, LDA, LWORK, NN,LDVL,LDVR
parameter(NN=4,LDA=NN,LWORK=4*NN,LDVL=NN,LDVR=NN)
real*8::WORK(0:LWORK-1),AA2(0:LDA-1,0:NN-1)
real*8::VL(0:NN-1,0:NN-1),VR(0:NN-1,0:NN-1)
real*8::WR(0:NN-1),WI(0:NN-1)
real*8::ev(NN) ,con,c2,x11,const,sigmay(0:3,0:3)   
sigmay(0,:)=(/0.0,0.0,0.0,-1.0/)
sigmay(1,:)=(/0.0,0.0,1.0,0.0/)
sigmay(2,:)=(/0.0,1.0,0.0,0.0/)
sigmay(3,:)=(/-1.0,0.0,0.0,0.0/)
site(0)=i0
site(1)=j0
call PTRVR(s,2,site,vin,rdm)
MM1=matmul(rdm,sigmay)      
MM2=matmul(rdm,sigmay)
MM3=matmul(MM1,MM2)
call dgeev( 'N', 'N', NN, MM3, LDA, WR,    &
           WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
do i=1,NN,1
    if (WR(i-1) .lt. 0.0000000000001) then
        WR(i-1)=0.0d0
    end if
    ev(i)=WR(i-1)
end do

do i=2,4
    if(ev(i)>=ev(1))then
        evs=ev(i)
        ev(i)=ev(1)
        ev(1)=evs   
    end if
end do

do i=3,4
    if(ev(i)>=ev(2))then
        evs=ev(i)
        ev(i)=ev(2)
        ev(2)=evs   
    end if
end do

do i=4,4
    if(ev(i)>=ev(3))then
        evs=ev(i)
        ev(i)=ev(3)
        ev(3)=evs   
    end if
end do
c2=0          
con=0
c2=dsqrt(ev(1))-dsqrt(ev(2))-dsqrt(ev(3))-dsqrt(ev(4))
if(c2>0)then
    con=c2
else
    con=0
end if
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.3.2
subroutine CONVC(s,i0,j0,vin,con)
implicit none
integer::s,i0,j0 ,site(0:1),i
complex*16::vin(0:2**s-1),evs
complex*16,dimension(0:3,0:3)::rdm,MM1,MM2,MM3
INTEGER::INFO, LDA, LWORK, N,LDVL,LDVR
parameter(N=4,LDA=N,LWORK=2*N,LDVL=N,LDVR=N)
complex*16::WORK(0:LWORK-1),VL(0:N-1,0:N-1)
complex*16::VR(0:N-1,0:N-1),W(0:N-1)
real*8::RWORK(0:2*N-1),ev(N) ,con,c2,x11,sigmay(0:3,0:3)      
sigmay(0,:)=(/0.0,0.0,0.0,-1.0/)
sigmay(1,:)=(/0.0,0.0,1.0,0.0/)
sigmay(2,:)=(/0.0,1.0,0.0,0.0/)
sigmay(3,:)=(/-1.0,0.0,0.0,0.0/)
site(0)=i0
site(1)=j0
call  PTRVC(s,2,site,vin,rdm)
MM1=matmul(rdm,sigmay)
MM2=matmul(conjg(rdm),sigmay)
MM3=matmul(MM1,MM2)
call zgeev( 'N', 'N', N, MM3, LDA, W, VL, LDVL, &
	      VR, LDVR, WORK, LWORK, RWORK, INFO )
do i=1,N,1
    if (real(W(i-1))<0.0000000000001) then
        W(i-1)=0.0d0
    end if
    ev(i)=W(i-1)
end do

do i=2,4
    if(ev(i)>=ev(1))then
        evs=ev(i)
        ev(i)=ev(1)
        ev(1)=evs   
    end if
end do

do i=3,4
    if(ev(i)>=ev(2))then
        evs=ev(i)
        ev(i)=ev(2)
        ev(2)=evs   
    end if
end do

do i=4,4
    if(ev(i)>=ev(3))then
        evs=ev(i)
        ev(i)=ev(3)
        ev(3)=evs   
    end if
end do
c2=0          
con=0
c2=dsqrt(ev(1))-dsqrt(ev(2))-dsqrt(ev(3))-dsqrt(ev(4))
if(c2>0)then
    con=c2
else
    con=0
end if
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.3.3
subroutine CONDMR(s,i0,j0,rho,con)
implicit none
integer::s,i0,j0 ,site(0:1),i
real*8::rho(0:2**s-1,0:2**s-1)
real*8,dimension(0:3,0:3)::rdm,MM1,MM2,MM3,MM4
INTEGER::INFO, LDA, LWORK, N,LDVL,LDVR
parameter(N=4,LDA=N,LWORK=4*N,LDVL=N,LDVR=N)
real*8::WORK(0:LWORK-1),VL(0:N-1,0:N-1)
real*8::VR(0:N-1,0:N-1),WR(0:N-1),WI(0:N-1)
real*8::RWORK(0:2*N-1),ev(N),con,c2,sigmay(0:3,0:3),evs    
sigmay(0,:)=(/0.0,0.0,0.0,-1.0/)
sigmay(1,:)=(/0.0,0.0,1.0,0.0/)
sigmay(2,:)=(/0.0,1.0,0.0,0.0/)
sigmay(3,:)=(/-1.0,0.0,0.0,0.0/)
site(0)=i0
site(1)=j0
call PTRDMR(s,2,site,rho,rdm)
MM1=matmul(rdm,sigmay)
MM2=matmul(rdm,sigmay)
MM3=matmul(MM1,MM2)
call dgeev( 'N', 'N', N, MM3, LDA, WR, WI, VL, &
                  LDVL, VR, LDVR, WORK, LWORK, INFO )
do i=1,N,1
    if (WR(i-1)<0.0000000000001) then
        WR(i-1)=0.0d0
    end if
    ev(i)=WR(i-1)
end do

do i=2,4
    if(ev(i)>=ev(1))then
        evs=ev(i)
        ev(i)=ev(1)
        ev(1)=evs   
    end if
end do

do i=3,4
    if(ev(i)>=ev(2))then
        evs=ev(i)
        ev(i)=ev(2)
        ev(2)=evs   
    end if
end do

do i=4,4
    if(ev(i)>=ev(3))then
        evs=ev(i)
        ev(i)=ev(3)
        ev(3)=evs   
    end if
end do
c2=0          
con=0
c2=dsqrt(ev(1))-dsqrt(ev(2))-dsqrt(ev(3))-dsqrt(ev(4))
if(c2>0)then
    con=c2
else
    con=0
end if
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.3.4
subroutine CONDMC(s,i0,j0,rho,con) 
implicit none
integer::s,i0,j0 ,site(0:1),i
complex*16::rho(0:2**s-1,0:2**s-1)
complex*16,dimension(0:3,0:3)::rdm,MM1,MM2,MM3,MM4
INTEGER::INFO, LDA, LWORK, N,LDVL,LDVR
parameter(N=4,LDA=N,LWORK=2*N,LDVL=N,LDVR=N)
complex*16::WORK(0:0,0:LWORK-1),AA2(0:LDA-1,0:N-1)
complex*16::VL(0:N-1,0:N-1),VR(0:N-1,0:N-1),W(0:N-1)
real*8::RWORK(0:2*N-1),ev(N) ,con,c2,x11,const 
real*8::sigmay(0:3,0:3),evs  
sigmay(0,:)=(/0.0,0.0,0.0,-1.0/)
sigmay(1,:)=(/0.0,0.0,1.0,0.0/)
sigmay(2,:)=(/0.0,1.0,0.0,0.0/)
sigmay(3,:)=(/-1.0,0.0,0.0,0.0/)
site(0)=i0
site(1)=j0
call PTRDMC(s,2,site,rho,rdm)
MM1=matmul(rdm,sigmay)
MM2=matmul(conjg(rdm),sigmay)
MM3=matmul(MM1,MM2)
call zgeev( 'N', 'N', N, MM3, LDA, W, VL, LDVL, &
	         VR, LDVR, WORK, LWORK, RWORK, INFO )
do i=1,N,1
    if (real(W(i-1))<0.0000000000001) then
        W(i-1)=0.0d0
    end if
    ev(i)=W(i-1)
end do

do i=2,4
    if(ev(i)>=ev(1))then
        evs=ev(i)
        ev(i)=ev(1)
        ev(1)=evs   
    end if
end do

do i=3,4
    if(ev(i)>=ev(2))then
        evs=ev(i)
        ev(i)=ev(2)
        ev(2)=evs   
    end if
end do

do i=4,4
    if(ev(i)>=ev(3))then
        evs=ev(i)
        ev(i)=ev(3)
        ev(3)=evs   
    end if
end do
c2=0          
con=0
c2=dsqrt(ev(1))-dsqrt(ev(2))-dsqrt(ev(3))-dsqrt(ev(4))
if(c2>0)then
    con=c2
else
    con=0
end if
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.4.1
subroutine  BERPVR(s,s1,site,vin,be)
implicit none
integer::s, s1,site(0:s1-1)
real*8::be,betemp,vin(0:2**s-1),rdm(0:2**s1-1,0:2**s1-1)
integer::i,INFO
real*8::WORK(0:3*(2**s1)-1),W(0:2**s1-1)
call PTRVR(s,s1,site,vin,rdm)
call  DSYEV( 'N','U', 2**s1, rdm, 2**s1, &
               W, WORK,3*(2**s1)-1,INFO )
be=0
betemp=0
do i=0,2**s1-1,1 
    if (W(i)>0.0000000000001) then
        betemp=-(dlog(abs(W(i)))*abs(W(i)))/dlog(2.0d0)
        be=be+betemp 
    end if 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.4.2
subroutine  BERPVC(s,s1,site,vin,be)
implicit none
integer::s, s1,site(0:s1-1)
real*8::be,betemp
complex*16::vin(0:2**s-1),rdm(0:2**s1-1,0:2**s1-1)
integer::i,INFOO
complex*16::WORKK(0:2*(2**s1)-1)
real*8::RWORKK(0:3*(2**s1)-3),WW(0:2**s1-1)
call PTRVC(s,s1,site,vin,rdm)
call  zheev( 'N','U', 2**s1, rdm, 2**s1, & 
             WW, WORKK,2*(2**s1)-1,RWORKK,INFOO )
be=0
betemp=0
do i=0,2**s1-1,1 
    if (WW(i)>0.0000000000001) then
        betemp=-(dlog(abs(WW(i)))*abs(WW(i)))/dlog(2.0d0)
        be=be+betemp 
    end if 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.4.3
subroutine BERPDMR(s,s1,site,rho,be)
implicit none
integer::s, s1,i,INFOO,site(0:s1-1)
real*8::trace,be,betemp,rho(0:2**s-1,0:2**s-1)
real*8::rdm(0:2**s1-1,0:2**s1-1)
real*8::WORKK(0:3*(2**s1)-1),WW(0:2**s1-1)
call PTRDMR(s,s1,site,rho,rdm)
call  DSYEV( 'N','U', 2**s1, rdm, 2**s1, &
                     WW, WORKK,3*(2**s1)-1,INFOO )
be=0
betemp=0
do i=0,2**s1-1,1 
    if (WW(i)>0.0000000000001) then
        betemp=-(dlog(abs(WW(i)))*abs(WW(i)))/dlog(2.0d0)
        be=be+betemp 
    end if 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.4.4
subroutine BERPDMC(s,s1,site,rho,be)
implicit none
integer::s, s1,i,INFOO,site(0:s1-1)
real*8::trace,be,betemp
complex*16::rho(0:2**s-1,0:2**s-1),rdm(0:2**s1-1,0:2**s1-1)
complex*16::WORKK(0:2*(2**s1)-1)
real*8::RWORKK(0:3*(2**s1)-3),WW(0:2**s1-1)
call PTRDMC(s,s1,site,rho,rdm)
call  zheev( 'V','U', 2**s1, rdm, 2**s1, &
                    WW, WORKK,2*(2**s1)-1,RWORKK,INFOO )
be=0
betemp=0
do i=0,2**s1-1,1 
    if (WW(i)>0.0000000000001) then
        betemp=-(dlog(abs(WW(i)))*abs(WW(i)))/dlog(2.0d0)
        be=be+betemp 
    end if 
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.5.1
subroutine RERPDMR(s,rho,alpha,re)
implicit none
integer::s
real*8::rho(0:2**s-1,0:2**s-1)
real*8::alpha,re
complex*16::rdma(0:2**s-1,0:2**s-1),trace
call POWNSR(rho,2**s,alpha,rdma)
call TRMC(rdma,2**s,trace)
re=(1/(1-alpha))*dlog(abs(trace))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.5.2
subroutine RERPDMC(s,rho,alpha,re)
implicit none
integer::s
real*8::re,alpha
complex*16::rho(0:2**s-1,0:2**s-1),trace
complex*16::rdma(0:2**s-1,0:2**s-1)
call POWNHC(rho,2**s,alpha,rdma)
call TRMC(rdma,2**s,trace)
re=(1/(1-alpha))*log(abs(trace))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.6.1
subroutine NEGVR(s,s1,site,vin,neg,logneg)
implicit none
integer::s,s1
real*8::trace,neg,logneg
integer,dimension(0:s1-1)::site
real*8,dimension(0:2**s-1)::vin
real*8,dimension(0:2**s-1,0:2**s-1):: rho,rhotp
call PTVR(s,s1,site,vin,rhotp)
call TRNORMMR(rhotp,2**s,trace)
neg=(trace-1)/2.0d0
logneg=dlog(trace)/dlog(2.0d0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.6.2
subroutine NEGVC(s,s1,site,vin,neg,logneg)
implicit none
integer::s,s1
real*8::trace,neg,logneg
integer,dimension(0:s1-1)::site
complex*16,dimension(0:2**s-1,0:2**s-1)::rho,rhotp
complex*16,dimension(0:2**s-1)::vin
call PTVC(s,s1,site,vin,rhotp)
call TRNORMMC(rhotp,2**s,trace)
neg=(trace-1)/2.0d0
logneg=dlog(trace)/dlog(2.0d0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.6.3
subroutine NEGDMR(s,s1,site,rho,neg,logneg)
implicit none
integer::s,s1
real*8::trace,neg,logneg
integer,dimension(0:s1-1)::site
real*8,dimension(0:2**s-1,0:2**s-1):: rho,rhotp
call PTDMR(s,s1,site,rho,rhotp)
call TRNORMMR(rhotp,2**s,trace)
neg=(trace-1)/2.0d0
logneg=dlog(trace)/dlog(2.0d0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.6.4
subroutine NEGDMC(s,s1,site,rho,neg,logneg)
implicit none
integer::s,s1
real*8::trace,neg,logneg
integer,dimension(0:s1-1)::site
complex*16,dimension(0:2**s-1,0:2**s-1):: rho,rhotp
call PTDMC(s,s1,site,rho,rhotp)
call TRNORMMC(rhotp,2**s,trace)
neg=(trace-1)/2.0d0
logneg=dlog(trace)/dlog(2.0d0)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.7.1
subroutine QMVR(s,vin,Q)
implicit none
integer::s,site(0:0),i
real*8::vin(0:2**s-1),rdm(0:1,0:1)
real*8::rdm2(0:1,0:1),trace,summ,Q
summ=0.0d0
do i=0,s-1,1
    site(0)=i
    call PTRVR(s,1,site,vin,rdm)
    rdm2=matmul(rdm,rdm)
    call TRMR(rdm2,2,trace)
    summ=summ+trace
end do
Q=2*(1-summ/s)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.7.2
subroutine QMVC(s,vin,Q)
implicit none
integer::s,site(0:0),i
complex*16::vin(0:2**s-1),trace
complex*16,dimension(0:1,0:1)::rdm,rdm2
real*8::summ,Q
summ=0.0d0
do i=0,s-1,1
    site(0)=i
    call PTRVC(s,1,site,vin,rdm)
    rdm2=matmul(rdm,rdm)
    call TRMC(rdm2,2,trace)
    summ=summ+real(trace)
end do
Q=2*(1-summ/s)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.7.3
subroutine QMDMR(s,rho,Q)
implicit none
integer::s,site(0:0),i
real*8::rho(0:2**s-1,0:2**s-1),rdm(0:1,0:1)
real*8::rdm2(0:1,0:1),trace,summ,Q
summ=0.0d0
do i=0,s-1,1
    site(0)=i
    call PTRDMR(s,1,site,rho,rdm)
    rdm2=matmul(rdm,rdm)
    call TRMR(rdm2,2,trace)
    summ=summ+trace
end do
Q=2*(1-summ/s)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.7.4
subroutine QMDMC(s,rho,Q)
implicit none
integer::s,site(0:0),i
complex*16::rho(0:2**s-1,0:2**s-1),trace
complex*16,dimension(0:1,0:1)::rdm,rdm2
real*8::summ,Q
summ=0.0d0
do i=0,s-1,1
    site(0)=i
    call PTRDMC(s,1,site,rho,rdm)
    rdm2=matmul(rdm,rdm)
    call TRMC(rdm2,2,trace)
    summ=summ+real(trace)
end do
Q=2*(1-summ/s)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.8.1
subroutine ENTSPECDMR(s,rho,eval,lneval)
integer::s,i,INFO
real*8::rho(0:2**s-1,0:2**s-1),eval(0:2**s-1)
real*8::WORK(0:3*(2**s)-1),W(0:2**s-1),lneval(0:2**s-1)
call  DSYEV( 'N','U', 2**s, rho, 2**s, W, WORK,3*(2**s)-1,INFO )
do i=0,2**s-1
    eval(i)=W(i)
    lneval(i)=-log(W(i))
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.8.2
subroutine ENTSPECDMC(s,rho,eval,lneval)
integer::s,i,INFO
complex*16::rho(0:2**s-1,0:2**s-1)
real*8::eval(0:2**s-1),lneval(0:2**s-1)
complex*16::WORK(0:2*(2**s)-1)
real*8::RWORK(0:3*(2**s)-3),W(0:2**s-1)
call  zheev( 'N','U', 2**s, rho, 2**s, W, WORK, &
                     2*(2**s)-1,RWORK,INFO )
do i=0,2**s-1
    eval(i)=W(i)
    lneval(i)=-log(W(i))
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.9.1
subroutine RESENTVR(vin,t)
integer,parameter::s=3,s1=1
real*8::vin(0:2**s-1),rdm(0:1,0:1),det,t,con12,con13
integer::site(0:s1-1)
site(0)=0
call PTRVR(s,s1,site,vin,rdm)
call CONVR(s,0,1,vin,con12)
call CONVR(s,0,2,vin,con13)
call DETSR(rdm,2,det)
t=4*det-con12**2-con13**2
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.9.2
subroutine RESENTVC(vin,t)
integer,parameter::s=3,s1=1
complex*16::vin(0:2**s-1),rdm(0:1,0:1)
real*8::t,con12,con13,det
integer::site(0:s1-1)
site(0)=0
call PTRVC(s,s1,site,vin,rdm)
call CONVC(s,0,1,vin,con12)
call CONVC(s,0,2,vin,con13)
call DETHC(rdm,2,det)
t=4*real(det)-con12**2-con13**2
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.9.3
subroutine RESENTDMR(rho,t)
integer,parameter::s=3,s1=1
real*8::rho(0:2**s-1,0:2**s-1),rdm(0:1,0:1),det,t,con12,con13
integer::site(0:s1-1)
site(0)=0
call PTRDMR(s,s1,site,rho,rdm)
call CONDMR(s,0,1,rho,con12)
call CONDMR(s,0,2,rho,con13)
call DETSR(rdm,2,det)
t=4*det-con12**2-con13**2
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!5.9.4
subroutine RESENTDMC(rho,t)
integer,parameter::s=3,s1=1
complex*16::rho(0:2**s-1,0:2**s-1),rdm(0:1,0:1)
real*8::t,con12,con13,det
integer::site(0:s1-1)
site(0)=0
call PTRDMC(s,s1,site,rho,rdm)
call CONDMC(s,0,1,rho,con12)
call CONDMC(s,0,2,rho,con13)
call DETHC(rdm,2,det)
t=4*real(det)-con12**2-con13**2
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHAPTER 6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.1.1
subroutine XHAM(s,B,HB)
implicit none
integer::s,t(0:s-1),y(0:s-1),i,j,d,w
real*8::B(0:s-1),h,e
real*8::HB(0:2**s-1,0:2**s-1)
HB=0.0d0
do i=0,2**s-1,1
    do j=i,2**s-1,1
        call DTOB(j,y,s)
        e=0
        t=y
        do d=0,s-1
            t(d)=1-t(d)
            call BTOD(t,w,s)
            if(w==i) then
                h=B(d)  
            else   
                h=0
            end if
            e=e+h
            t=y
        end do
        HB(i,j)=e 
        HB(j,i)=HB(i,j)  
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.1.2
subroutine YHAM(s,B,HB)
implicit none
integer::s,t(0:s-1),y(0:s-1),i,j,d,w
real*8::B(0:s-1)
complex*16::HB(0:2**s-1,0:2**s-1),phase,h,e
HB=0.0d0
do i=0,2**s-1,1
    do j=i,2**s-1,1
        call DTOB(j,y,s)
        e=0
        t=y
        do d=0,s-1
            if(t(d)==0) then
                t(d)=1
                phase=dcmplx(0,1)
            else   
                t(d)=0
                phase=dcmplx(0,-1)
            end if
            call BTOD(t,w,s)
            if(w==i) then
                h=phase*B(d)
            else   
                h=0
            end if
            e=e+h
            t=y
        end do
        HB(i,j)=e 
        HB(j,i)=conjg(HB(i,j))   
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.1.3
subroutine ZHAM(s,B,HB)
implicit none
integer::s,i,j,v,x(0:s-1)
real*8::B(0:s-1)
real*8::HB(0:2**s-1,0:2**s-1),m1
HB=0.0d0
do i=0,2**s-1,1
    do j=i,2**s-1,1
        call DTOB(i,x,s)
        if (i==j) then    
            do v=0,s-1           
                m1=(1-2*x(v))*B(v)
                HB(i,j)=HB(i,j)+m1
            end do
        else
            HB(i,j)=0
        end if
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.2.1
subroutine  ISINGXHAM(s,r,BC,Jr,HI)
implicit none
integer::s,r,t(0:s-1),y(0:s-1),i,j,d,w
real*8::HI(0:2**s-1,0:2**s-1),e,e1
real*8::h,Jr(0:s-1)
character*1::BC
if (r>s) then
    print*, 'r should be less than s'
    stop
end if
HI=0.0d0
do i=0,2**s-1,1
    do j=i,2**s-1,1
        call DTOB(j,y,s)
        e=0
        e1=0
        t=y
        do d=0,s-1
            if(d.le.s-r-1) then  
                t(d)=1-t(d)
                t(d+r)=1-t(d+r)
                call BTOD(t,w,s)
                if(w==i) then
                    h=Jr(d) 
                else   
                    h=0
                end if
                e=e+h
                t=y
            else
                if (BC=='y' ) then
                    t(d)=1-t(d)
                    t(d+r-s)=1-t(d+r-s)
                    call BTOD(t,w,s)
                    if(w==i) then
                        h=Jr(d)  
                    else   
                        h=0
                    end if
                    e1=e1+h
                    t=y  
                end if
            end if
        end do
        HI(i,j)=e+e1
        HI(j,i)=HI(i,j)
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.2.2
subroutine  ISINGYHAM(s,r,BC,Jr,HI)
implicit none
integer*4::s,r,t(0:s-1),y(0:s-1),d,w,i,j
real*8::HI(0:2**s-1,0:2**s-1),e,e1,Jr(0:s-1)
complex*16::i2,p1,p2,p3,p4,h
character*1::BC
if (r>s) then
    print*, 'r should be less than s'
    stop
end if
HI=0.0d0
i2=dcmplx(0,1)
do i=0,2**s-1,1
    do j=i,2**s-1,1
        call DTOB(j,y,s)
        e=0
        e1=0
        t=y
        do d=0,s-1
            if(d.le.s-r-1) then  
                p1=i2*((-1)**(t(d)))
                p2=i2*((-1)**(t(d+r)))
                t(d)=1-t(d)
                t(d+r)=1-t(d+r)
                call BTOD(t,w,s)
                if(w==i) then
                    h=p1*p2*Jr(d)  
                else   
                    h=0
                end if
                e=e+real(h)
                t=y
            else
                if (BC=='y') then
                    p3=i2*((-1)**(t(d)))
                    p4=i2*((-1)**(t(d+r-s)))
                    t(d)=1-t(d)
                    t(d+r-s)=1-t(d+r-s)
                    call BTOD(t,w,s)
                    if(w==i) then
                        h=p3*p4*Jr(d)
                    else   
                        h=0
                    end if
                    e1=e1+real(h)
                    t=y 
                end if
            end if
        end do
        HI(i,j)=e+e1
        HI(j,i)=HI(i,j)
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.2.3
subroutine  ISINGZHAM(s,r,BC,Jr,HI)
implicit none
integer*4::s,r,i,j,v,x(0:s-1)
real*8::HI(0:2**s-1,0:2**s-1),m1,Jr(0:s-1)
character*1::BC
if (r>s) then
    print*, 'r should be less than s'
    stop
end if
HI=0.0d0
do i=0,2**s-1,1
    do j=i,2**s-1,1
        HI(i,j)=0.0d0
        call DTOB(j,x,s)
        if (i==j) then    
            do v=0,s-1
                if(v.le.s-r-1) then         
                    m1=(1-2*x(v))*(1-2*x(v+r))*Jr(v)
                    HI(i,j)=HI(i,j)+m1
                else
                    if (BC=='y') then
                        m1=(1-2*x(v))*(1-2*x(v+r-s))*Jr(v)
                        HI(i,j)=HI(i,j)+m1
                    end if 
                end if 
            end do
        end if
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.3
subroutine HEIHAM(s,r,BC,Jr,H)
implicit none 
integer::s,r,i,j,y(0:s-1),z(0:s-1),o,w,temp
real*8,dimension(0:2**s-1,0:2**s-1)::H
real*8::c,h1,h2,h3,h4,h5
character*1::BC
real*8::Jr(0:s-1),summ
summ=0.0d0
if (BC=='y') then
    do i=0,s-1
        summ=summ+Jr(i)
    end do
end if
if (BC=='n') then
    do i=0,s-r-1
        summ=summ+Jr(i)
    end do
end if
do i=0,2**s-1,1
    do j=i,2**s-1,1
        call DTOB(j,y,s)
        h1=0
        h2=0
        h3=0
        h4=0
        do o=0,s-1,1
            z=y
            if(o.le.s-r-1) then 
                temp=z(o)
                z(o)=z(o+r)    
                z(o+r)=temp
                call BTOD(z,w,s)
                if(w==i) then
                    c=Jr(o)
                else
                    c=0
                end if
                h1=h1+c
            else
                if (BC=='y') then
                    temp=z(o)
                    z(o)=z(o+r-s)
                    z(o+r-s)=temp
                    call BTOD(z,w,s)
                    if(w==i) then
                        c=Jr(o)
                    else
                        c=0
                    end if
                    h2=h2+c
                end if
            end if
        end do 
        if(i==j)then
            h3=summ
        else 
            h3=0
        end if
        h4=h1+h2
        H(i,j)=4.0d0*(h4/2.0d0-h3/4.0d0)
        H(j,i)=H(i,j)
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.4.1
subroutine DMHAMX(s,r,BC,Dr,HDM)
implicit none
integer::s,r
integer::x,i,j,t(0:s-1),y(0:s-1),w1
real*8::Dr(0:s-1),m
character*1::BC
complex*16::i2,phase,h,HDM(0:2**s-1,0:2**s-1),e,e1,h1,h2,e2,e3
i2=dcmplx(0,1)
do i=0,2**s-1
    do j=i,2**s-1
        call DTOB(j,y,s)
        e=0
        e1=0
        t=y
        e2=0
        e3=0
        do x=0,s-1
            if(x.le.s-r-1) then  
                phase=i2*((-1)**(t(x)))
                t(x)=1-t(x)
                m=1-2*t(x+r)
                call BTOD(t,w1,s)
                if(w1==i) then
                    h=m*phase*Dr(x)
                else   
                    h=0
                end if
                e=e+h
                t=y
            else
                if (BC=='y') then
                    phase=i2*((-1.0)**(t(x)))
                    t(x)=1-t(x)
                    m=1-2*t(x+r-s)
                    call BTOD(t,w1,s)
                    if(w1==i) then
                        h=m*phase*Dr(x)
                    else   
                        h=0
                    end if
                    e1=e1+h
                    t=y
                end if
            end if
        end do
        h1=e+e1
        do x=0,s-1
            if(x.le.s-r-1) then
                phase=i2*((-1.0)**(t(x+r)))
                m=1-2*t(x)  
                t(x+r)=1-t(x+r)
                call BTOD(t,w1,s)
                if(w1==i) then
                    h=m*phase*Dr(x) 
                else   
                    h=0
                end if
                e2=e2+h
                t=y
            else
                if (BC=='y') then
                    phase=i2*((-1)**(t(x+r-s)))
                    t(x+r-s)=1-t(x+r-s)
                    m=1-2*t(x)
                    call BTOD(t,w1,s)
                    if(w1==i) then
                        h=m*phase*Dr(x)
                    else   
                        h=0
                    end if
                    e3=e3+h
                    t=y
                end if
            end if
        end do
        h2=e2+e3
        HDM(i,j)=h1-h2
        HDM(j,i)=conjg(HDM(i,j))
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.4.2
subroutine DMHAMY(s,r,BC,Dr,HDM)
implicit none
integer::s,r
integer::x,i,j,t(0:s-1),y(0:s-1),w1
real*8::Dr(0:s-1),m
character*1::BC
real*8::i2,h,HDM(0:2**s-1,0:2**s-1),e,e1,h1,h2,e2,e3
do i=0,2**s-1
    do j=i,2**s-1
        call DTOB(j,y,s)
        e=0
        e1=0
        t=y
        e2=0
        e3=0
        do x=0,s-1
            if(x.le.s-r-1) then  
                t(x)=1-t(x)
                m=1-2*t(x+r)
                call BTOD(t,w1,s)
                if(w1==i) then
                    h=m*Dr(x)
                else   
                    h=0
                end if
                e=e+h
                t=y
            else
                if (BC=='y') then
                    t(x)=1-t(x)
                    m=1-2*t(x+r-s)
                    call BTOD(t,w1,s)
                    if(w1==i) then
                        h=m*Dr(x)
                    else   
                        h=0
                    end if
                    e1=e1+h
                    t=y
                end if
            end if
        end do
        h1=e+e1
        do x=0,s-1
            if(x.le.s-r-1) then
                m=1-2*t(x)  
                t(x+r)=1-t(x+r)
                call BTOD(t,w1,s)
                if(w1==i) then
                    h=m*Dr(x)
                else   
                    h=0
                end if
                e2=e2+h
                t=y
            else
                if (BC=='y') then
                    t(x+r-s)=1-t(x+r-s)
                    m=1-2*t(x)
                    call BTOD(t,w1,s)
                    if(w1==i) then
                        h=m*Dr(x)
                    else   
                        h=0
                    end if
                    e3=e3+h
                    t=y
                end if
            end if
        end do
        h2=e2+e3
        HDM(i,j)=(h2-h1)
        HDM(j,i)=HDM(i,j)
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!6.4.3
subroutine DMHAMZ(s,r,BC,Dr,HDM)
implicit none
integer::s,r
integer::x,i,j,t(0:s-1),y(0:s-1),w1
real*8::Dr(0:s-1),m
character*1::BC
complex*16::i2,phase,h,HDM(0:2**s-1,0:2**s-1),e,e1,h1,h2,e2,e3
i2=dcmplx(0,1)
do i=0,2**s-1
    do j=i,2**s-1
        call DTOB(j,y,s)
        e=0
        e1=0
        t=y
        e2=0
        e3=0
        do x=0,s-1
            if(x.le.s-r-1) then  
                phase=i2*((-1)**(t(x+r)))
                t(x)=1-t(x)
                t(x+r)=1-t(x+r)
                call BTOD(t,w1,s)
                if(w1==i) then
                    h=phase*Dr(x)
                else   
                    h=0
                end if
                e=e+h
                t=y
            else
                if (BC=='y') then
                    phase=i2*((-1.0)**(t(x+r-s)))
                    t(x)=1-t(x)
                    t(x+r-s)=1-t(x+r-s)
                    call BTOD(t,w1,s)
                    if(w1==i) then
                        h=phase*Dr(x)
                    else   
                        h=0
                    end if
                    e1=e1+h
                    t=y
                end if
            end if
        end do
        h1=e+e1
        do x=0,s-1
            if(x.le.s-r-1) then
                phase=i2*((-1.0)**(t(x)))
                t(x)=1-t(x)  
                t(x+r)=1-t(x+r)
                call BTOD(t,w1,s)
                if(w1==i) then
                    h=phase*Dr(x)
                else   
                    h=0
                end if
                e2=e2+h
                t=y
            else
                if (BC=='y') then
                    phase=i2*((-1)**(t(x)))
                    t(x+r-s)=1-t(x+r-s)
                    t(x)=1-t(x)
                    call BTOD(t,w1,s)
                    if(w1==i) then
                        h=phase*Dr(x)
                    else   
                        h=0
                    end if
                    e3=e3+h
                    t=y
                end if
            end if
        end do
        h2=e2+e3
        HDM(i,j)=(h1-h2)
        HDM(j,i)=conjg(HDM(i,j))
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHAPTER 7!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.1.1
subroutine RANDDISTGAUSSR(d1,d2,mu,sigma,state)
implicit none
integer::d1,d2,i,j
real*8::state(0:d1-1,0:d2-1),rn1,rn2,pi,mu,sigma
pi=4.0d0*atan(1.0d0)
do i=0,d1-1
    do j=0,d2-1
        call random_number(rn1)
        call random_number(rn2)
        state(i,j)=mu+sigma*sqrt(-2*log(rn1))*cos(2*pi*rn2)
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.1.2
subroutine RANDDISTGAUSSC(d1,d2,mu,sigma,state)
implicit none
integer::d1,d2,i,j
complex*16::state(0:d1-1,0:d2-1)
real*8::rn1,rn2,pi,n1,n2,rn3,rn4,mu,sigma
pi=4.0d0*atan(1.0d0)
do i=0,d1-1
    do j=0,d2-1
        call random_number(rn1)
        call random_number(rn2)
        call random_number(rn3)
        call random_number(rn4)
        n1=mu+sigma*sqrt(-2*log(rn1))*cos(2*pi*rn2)
        n2=mu+sigma*sqrt(-2*log(rn3))*cos(2*pi*rn4)
        state(i,j)=dcmplx(n1,n2)
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.2.1
subroutine RANDDISTUNIR(a,b,d1,d2,state)
implicit none
integer::d1,d2
real*8::state(0:d1-1,0:d2-1),a,b
call random_number(state)
state=a+(b-a)*state
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.2.2
subroutine RANDDISTUNIC(a,b,d1,d2,state)
implicit none
integer::d1,d2,i,j
real*8::state1(0:d1-1,0:d2-1),state2(0:d1-1,0:d2-1),a,b
complex*16::state(0:d1-1,0:d2-1)
call random_number(state1)
call random_number(state2)
state1=a+(b-a)*state1
state2=a+(b-a)*state2
do i=0,d1-1
    do j=0,d2-1
        state(i,j)=dcmplx(state1(i,j),state2(i,j))
    end do
end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.3
subroutine RANDSR(d,dist,a,b,mu,sigma,state)
implicit none
integer::i,d
real*8::state(0:d-1,0:d-1),a,b,mu,sigma
character*3::dist
if (dist=='GAU') then
    call RANDDISTGAUSSR(d,d,mu,sigma,state)
end if
if (dist=='UNI') then
    call RANDDISTUNIR(a,b,d,d,state)
end if
state=(state+transpose(state))/2.d0
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.4
subroutine RANDHC(d,dist,a,b,mu,sigma,state)
implicit none
integer::i,d
complex*16::state(0:d-1,0:d-1)
real*8::a,b,mu,sigma
character*3::dist
if (dist=='GAU') then
    call RANDDISTGAUSSC(d,d,mu,sigma,state)
end if
if (dist=='UNI') then
    call RANDDISTUNIC(a,b,d,d,state)
end if
state=(state+transpose(conjg(state)))/2.d0
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.5.1
subroutine RANDORTHOGONAL(d1,d2,dist,a,b,mu,sigma,state)
implicit none
integer::i,d1,d2
real*8::state(0:d1-1,0:d2-1),a,b,mu,sigma,R(0:d1-1,0:d2-1)
character*3::dist
if (dist=='GAU') then
    call RANDDISTGAUSSR(d1,d2,mu,sigma,state)
end if
if (dist=='UNI') then
    call RANDDISTUNIR(a,b,d1,d2,state)
end if
call QRMR(d1,d2,state,R)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.5.2
subroutine RANDUNITARY(d1,d2,dist,a,b,mu,sigma,state)
implicit none
integer::i,d1,d2
complex*16::state(0:d1-1,0:d2-1),R(0:d1-1,0:d2-1)
real*8::a,b,mu,sigma
character*3::dist
if (dist=='GAU') then
    call RANDDISTGAUSSC(d1,d2,mu,sigma,state)
end if
if (dist=='UNI') then
    call RANDDISTUNIC(a,b,d1,d2,state)
end if
call QRMC(d1,d2,state,R)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.6
subroutine GINIBRER(d,state)
implicit none
integer::d
real*8::state(0:d-1,0:d-1),mu,sigma
mu=0.0d0
sigma=1.0d0
call RANDDISTGAUSSR(d,d,mu,sigma,state)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.7
subroutine GINIBREC(d,state)
implicit none
integer::d
complex*16::state(0:d-1,0:d-1)
real*8::mu,sigma
mu=0.0d0
sigma=1.0d0/sqrt(2.0d0)
call RANDDISTGAUSSC(d,d,mu,sigma,state)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.8.1
subroutine WISHARTSR(d,state)
implicit none
integer::d
real*8::state(0:d-1,0:d-1)
call GINIBRER(d,state)
state=matmul(state,transpose(state))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.8.2
subroutine WISHARTHC(d,state)
implicit none
integer::d
complex*16::state(0:d-1,0:d-1)
call GINIBREC(d,state)
state=matmul(state,transpose(conjg(state)))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.9
subroutine RANDPVEC(d,state)
implicit none
integer::d,i
real*8::state(0:d-1),norm
character*3::dist
call random_number(state)
norm=0.0d0
do i=0,d-1
    norm=norm+state(i)
end do
state=state/norm
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.10.1
subroutine RANDVR(d,dist,a,b,mu,sigma,state)
implicit none
integer::i,d
real*8::state1(0:d-1,0:0),a,b,state(0:d-1),norm,mu,sigma
character*3::dist
if (dist=='GAU') then
    call RANDDISTGAUSSR(d,1,mu,sigma,state1)
end if
if (dist=='UNI') then
    call RANDDISTUNIR(a,b,d,1,state1)
end if
norm=0.0d0
do i=0,d-1,1
    state(i)=state1(i,0)
end do
call NORMALIZEVR(state,d)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.10.2
subroutine RANDVC(d,dist,a,b,mu,sigma,state)
implicit none
integer::i,d
complex*16::state1(0:d-1,0:0),state(0:d-1)
real*8::a,b,mu,sigma
character*3::dist
if (dist=='GAU') then
    call RANDDISTGAUSSC(d,1,mu,sigma,state1)
end if
if (dist=='UNI') then
    call RANDDISTUNIC(a,b,d,1,state1)
end if
do i=0,d-1,1
    state(i)=state1(i,0)
end do
call NORMALIZEVC(state,d)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.11.1
subroutine RANDDMR(d,state)
integer::d
real*8::state(0:d-1,0:d-1),mu,sigma,trace
call WISHARTSR(d,state)
trace=0.0d0
do i=0,d-1,1
    trace=trace+state(i,i)
end do
state=state/trace
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!7.11.2
subroutine RANDDMC(d,state)
integer::d
complex*16::state(0:d-1,0:d-1),trace
call WISHARTHC(d,state)
trace=0.0d0
do i=0,d-1,1
trace=trace+state(i,i)
end do
state=state/trace
end subroutine


