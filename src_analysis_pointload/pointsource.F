      subroutine pointf(xx,rf,strike,dip,rake,ftmp)
      real*8 r1,r2,r3
      real*8 xx(4,3),rm(3,3),ftmp(30)
      real*8 DIPt,RAKt,STRt,pi,strike,dip,rake,rf(3)
c
      pi=4.*atan(1.)
c
      r1=rf(1)
      r2=rf(2)
      r3=rf(3)
c
c strike
      STRt=pi/180.*  strike
c dip
      DIPt=pi/180.*  dip
c rake
      RAKt=pi/180.*  rake
c
      call coemoment(DIPt,RAKt,STRt,rm)
c
      call point(xx,r1,r2,r3,rm,ftmp)
c
      end
c_______________________________________________________________________
      subroutine point(xx,r1,r2,r3,rm,ftmp)
      real*8 dfaix(10,3),rm(3,3),ftmp(30)
      real*8 r1,r2,r3
      real*8 detjmat,xx(4,3)
      real*8 x2,x3,x4,y2,y3,y4,z2,z3,z4
      real*8 ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33
      real*8 invj(3,3),dfair(10,3)
c
      x2=xx(2,1)-xx(1,1)
      x3=xx(3,1)-xx(1,1)
      x4=xx(4,1)-xx(1,1)
      y2=xx(2,2)-xx(1,2)
      y3=xx(3,2)-xx(1,2)
      y4=xx(4,2)-xx(1,2)
      z2=xx(2,3)-xx(1,3)
      z3=xx(3,3)-xx(1,3)
      z4=xx(4,3)-xx(1,3)
c
      ij11=-(y4*z3) + y3*z4
      ij12=y4*z2-y2*z4
      ij13=-(y3*z2) + y2*z3
      ij21=x4*z3-x3*z4
      ij22=-(x4*z2) + x2*z4
      ij23=x3*z2-x2*z3
      ij31=-(x4*y3) + x3*y4
      ij32=x4*y2-x2*y4
      ij33=-(x3*y2) + x2*y3
c
      detjmat=-x4*y3*z2 + x3*y4*z2 +
     - x4*y2*z3-x2*y4*z3-x3*y2*z4 + x2*y3*z4
c
      invj(1,1)=ij11
      invj(1,2)=ij12
      invj(1,3)=ij13
      invj(2,1)=ij21
      invj(2,2)=ij22
      invj(2,3)=ij23
      invj(3,1)=ij31
      invj(3,2)=ij32
      invj(3,3)=ij33
      invj=invj/detjmat
c
      call compdfair(r1,r2,r3,dfair)
c
      do i=1,10
      do j=1,3
      dfaix(i,j)=
     - invj(j,1)*dfair(i,1)+invj(j,2)*dfair(i,2)+invj(j,3)*dfair(i,3)
      enddo
      enddo
c
      do i=1,10
      do k=1,3
      tmp=0.
      do j=1,3
      tmp=tmp+dfaix(i,j)*rm(k,j)
      enddo
      ftmp(3*(i-1)+k)=tmp
      enddo
      enddo
c
      end
c_______________________________________________________________________
      subroutine compdfair(r1,r2,r3,dfair)
      real*8 r1,r2,r3,dfair(10,3)
c
      dfair(1,1)=1 - 4*(1 - r1 - r2 - r3)
      dfair(2,1)=-1 + 4*r1
      dfair(3,1)=0
      dfair(4,1)=0
      dfair(5,1)=-4*r1 + 4*(1 - r1 - r2 - r3)
      dfair(6,1)=4*r2
      dfair(7,1)=-4*r2
      dfair(8,1)=-4*r3
      dfair(9,1)=4*r3
      dfair(10,1)=0
      dfair(1,2)=1 - 4*(1 - r1 - r2 - r3)
      dfair(2,2)=0
      dfair(3,2)=-1 + 4*r2
      dfair(4,2)=0
      dfair(5,2)=-4*r1
      dfair(6,2)=4*r1
      dfair(7,2)=-4*r2 + 4*(1 - r1 - r2 - r3)
      dfair(8,2)=-4*r3
      dfair(9,2)=0
      dfair(10,2)=4*r3
      dfair(1,3)=1 - 4*(1 - r1 - r2 - r3)
      dfair(2,3)=0
      dfair(3,3)=0
      dfair(4,3)=-1 + 4*r3
      dfair(5,3)=-4*r1
      dfair(6,3)=0
      dfair(7,3)=-4*r2
      dfair(8,3)=4*(1 - r1 - r2 - r3) - 4*r3
      dfair(9,3)=4*r1
      dfair(10,3)=4*r2
c
      end
c_______________________________________________________________________
      subroutine coemoment(DIP,RAK,STR,coem)
      real*8 coem(3,3),DIP,RAK,STR
c
      coem(1,1)=2*Cos(STR)*Sin(DIP)*
     -  (-(Cos(DIP)*Cos(STR)*Sin(RAK)) + Cos(RAK)*Sin(STR))
      coem(2,2)=-2*Sin(DIP)*Sin(STR)*
     -  (Cos(RAK)*Cos(STR) + Cos(DIP)*Sin(RAK)*Sin(STR))
      coem(3,3)=2*Cos(DIP)*Sin(DIP)*Sin(RAK)
      coem(1,2)=Sin(DIP)*(Cos(RAK)*Cos(2*STR) +
     -    Cos(DIP)*Sin(RAK)*Sin(2*STR))
      coem(1,3)=-(Cos(DIP)**2*Cos(STR)*Sin(RAK)) +
     -  Cos(STR)*Sin(DIP)**2*Sin(RAK) + Cos(DIP)*Cos(RAK)*Sin(STR)
      coem(2,3)=Cos(DIP)*Cos(RAK)*Cos(STR) +
     -   Cos(DIP)**2*Sin(RAK)*Sin(STR) -
     -   Sin(DIP)**2*Sin(RAK)*Sin(STR)
      coem(2,1)=coem(1,2)
      coem(3,2)=coem(2,3)
      coem(3,1)=coem(1,3)
      end
