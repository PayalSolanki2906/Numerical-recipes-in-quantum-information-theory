! Source codes for Numerical Recipes in Quantum Information Theory and Quantum Computing: 
!An Adventure in FORTRAN 90 -by M. S. Ramkarthik and Payal D. Solanki

!This file contains all the source codes for chapter 6.

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


