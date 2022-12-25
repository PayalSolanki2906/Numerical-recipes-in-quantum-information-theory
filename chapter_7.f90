! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes for chapter 7.

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


