! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes for chapter 5.

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


