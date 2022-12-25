! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes for chapter 3.

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



