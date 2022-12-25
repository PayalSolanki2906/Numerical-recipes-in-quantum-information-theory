! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes for chapter 4.

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


