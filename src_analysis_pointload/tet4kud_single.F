      subroutine tet4kus(young,gxy,gyz,gxz,rnyu,xx,u,bdbu,detj)
      implicit none
      REAL_4 young,rnyu,cor,u(12),bdbu(12),xl(4,3),xx(4,3)
      REAL_4 gxy,gyz,gxz
      REAL_4 bu(6),dbu(6),d1,d2,d3,bmc(3),detjt,tmp1,tmp2
      REAL_4 detj,ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33
      integer i,ii
c
      tmp1=1-rnyu
      tmp2=1+rnyu
      d1=(tmp1*young)/((1-2*rnyu)*tmp2)
      d2=d1*rnyu/tmp1
      d3=young/(2.*tmp2)
c
      do i=2,4
      do ii=1,3
      xl(i,ii)=xx(i,ii)-xx(1,ii)
      enddo
      enddo
c
      detj=-(xl(2,3)*xl(3,2)*xl(4,1))+xl(2,2)*xl(3,3)*xl(4,1)+
     -  xl(2,3)*xl(3,1)*xl(4,2)-xl(2,1)*xl(3,3)*xl(4,2)-
     -  xl(2,2)*xl(3,1)*xl(4,3)+xl(2,1)*xl(3,2)*xl(4,3)
c
      ij11=-(xl(3,3)*xl(4,2))+xl(3,2)*xl(4,3)
      ij12=xl(2,3)*xl(4,2)-xl(2,2)*xl(4,3)
      ij13=-(xl(2,3)*xl(3,2))+xl(2,2)*xl(3,3)
      ij21=xl(3,3)*xl(4,1)-xl(3,1)*xl(4,3)
      ij22=-(xl(2,3)*xl(4,1))+xl(2,1)*xl(4,3)
      ij23=xl(2,3)*xl(3,1)-xl(2,1)*xl(3,3)
      ij31=-(xl(3,2)*xl(4,1))+xl(3,1)*xl(4,2)
      ij32=xl(2,2)*xl(4,1)-xl(2,1)*xl(4,2)
      ij33=-(xl(2,2)*xl(3,1))+xl(2,1)*xl(3,2)
c
      bmc(1)=-(ij11+ij12+ij13)
      bmc(2)=-(ij21+ij22+ij23)
      bmc(3)=-(ij31+ij32+ij33)
c
      bu(1)=bmc(1)*u(1) + ij11*u(4) + ij12*u(7) + ij13*u(10)
      bu(2)=bmc(2)*u(2) + ij21*u(5) + ij22*u(8) + ij23*u(11)
      bu(3)=bmc(3)*u(3) + ij31*u(6) + ij32*u(9) + ij33*u(12)
      bu(4)=bmc(2)*u(1) + bmc(1)*u(2) + ij21*u(4) + ij11*u(5) +
     -  ij22*u(7) + ij12*u(8) + ij23*u(10) + ij13*u(11)
      bu(5)= bmc(3)*u(2) + bmc(2)*u(3) + ij31*u(5) + ij21*u(6) +
     -  ij32*u(8) + ij22*u(9) + ij33*u(11) + ij23*u(12)
      bu(6)=bmc(3)*u(1) + bmc(1)*u(3) + ij31*u(4) + ij11*u(6) +
     -  ij32*u(7) + ij12*u(9) + ij33*u(10) + ij13*u(12)
c
      dbu(1)=d1*bu(1) + d2*(bu(2) + bu(3))
      dbu(2)=d1*bu(2) + d2*(bu(1) + bu(3))
      dbu(3)=d2*(bu(1) + bu(2)) + d1*bu(3)
      dbu(4)=gxy*bu(4)
      dbu(5)=gyz*bu(5)
      dbu(6)=gxz*bu(6)
c
      bdbu(1)=bmc(1)*dbu(1) + bmc(2)*dbu(4) + bmc(3)*dbu(6)
      bdbu(2)=bmc(2)*dbu(2) + bmc(1)*dbu(4) + bmc(3)*dbu(5)
      bdbu(3)=bmc(3)*dbu(3) + bmc(2)*dbu(5) + bmc(1)*dbu(6)
      bdbu(4)=ij11*dbu(1) + ij21*dbu(4) + ij31*dbu(6)
      bdbu(5)=ij21*dbu(2) + ij11*dbu(4) + ij31*dbu(5)
      bdbu(6)=ij31*dbu(3) + ij21*dbu(5) + ij11*dbu(6)
      bdbu(7)=ij12*dbu(1) + ij22*dbu(4) + ij32*dbu(6)
      bdbu(8)=ij22*dbu(2) + ij12*dbu(4) + ij32*dbu(5)
      bdbu(9)=ij32*dbu(3) + ij22*dbu(5) + ij12*dbu(6)
      bdbu(10)=ij13*dbu(1) + ij23*dbu(4) + ij33*dbu(6)
      bdbu(11)=ij23*dbu(2) + ij13*dbu(4) + ij33*dbu(5)
      bdbu(12)=ij33*dbu(3) + ij23*dbu(5) + ij13*dbu(6)
c
      detjt=1./detj*0.16666666666666667
c
      do i=1,12
      bdbu(i)=bdbu(i)*detjt
      enddo
c
      end
