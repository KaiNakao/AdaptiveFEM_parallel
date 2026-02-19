c 20090806 by Tsuyoshi Ichimura
c mix. disp:heviside, stress:constant
c stifness matrix: node i, componet j -> 3*(i-1)+j
      subroutine nontetra4(young,gxy,gyz,gxz,rnyu,xx,ke)
      real*8 young,gxy,gyz,gxz,rnyu,xx(4,3),ke(12,12)
c
      real*8 detj,ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33
      real*8 d1,d2,d3,d4,d5,d6,d7,d8,d9
c
c
c      d1=((1-nyu)*young)/((1-2*nyu)*(1+nyu))
c      d2=(nyu*young)/((1-2*nyu)*(1 + nyu))
c      d3=((1-nyu)*young)/((1-2*nyu)*(1+nyu))
c      d4=(nyu*young)/((1-2*nyu)*(1+nyu))
c      d5=((1-nyu)*young)/((1-2*nyu)*(1+nyu))
c      d6=(nyu*young)/((1-2*nyu)*(1+nyu))
c      d7=young/(2.*(1+nyu))
c      d8=young/(2.*(1+nyu))
c      d9=young/(2.*(1+nyu))
      d7=gxy
      d8=gyz
      d9=gxz
c
      d1=((1-rnyu)*young)/((1-2*rnyu)*(1+rnyu))
	d2=(rnyu*young)/((1-2*rnyu)*(1+rnyu))
	d3=young/(2.*(1+rnyu))
c
      do i=2,4
      do ii=1,3
      xx(i,ii)=xx(i,ii)-xx(1,ii)
      enddo
      enddo
      do ii=1,3
      xx(1,ii)=0.
      enddo
c
      detj=-(xx(2,3)*xx(3,2)*xx(4,1))+xx(2,2)*xx(3,3)*xx(4,1)+
     -  xx(2,3)*xx(3,1)*xx(4,2)-xx(2,1)*xx(3,3)*xx(4,2)-
     -  xx(2,2)*xx(3,1)*xx(4,3)+xx(2,1)*xx(3,2)*xx(4,3)
c
      ij11=-(xx(3,3)*xx(4,2))+xx(3,2)*xx(4,3)
	ij12=xx(2,3)*xx(4,2)-xx(2,2)*xx(4,3)
	ij13=-(xx(2,3)*xx(3,2))+xx(2,2)*xx(3,3)
	ij21=xx(3,3)*xx(4,1)-xx(3,1)*xx(4,3)
	ij22=-(xx(2,3)*xx(4,1))+xx(2,1)*xx(4,3)
	ij23=xx(2,3)*xx(3,1)-xx(2,1)*xx(3,3)
	ij31=-(xx(3,2)*xx(4,1))+xx(3,1)*xx(4,2)
	ij32=xx(2,2)*xx(4,1)-xx(2,1)*xx(4,2)
	ij33=-(xx(2,2)*xx(3,1))+xx(2,1)*xx(3,2)
c
      ke(1,1)=d1*(ij11+ij12+ij13)**2+d7*(ij21+ij22+ij23)**2+
     -  d9*(ij31+ij32+ij33)**2
      ke(2,1)=(d2+d7)*(ij11+ij12+ij13)*(ij21+ij22+ij23)
      ke(3,1)=(d2+d9)*(ij11+ij12+ij13)*(ij31+ij32+ij33)
      ke(4,1)=-(d1*ij11*(ij11+ij12+ij13)) -
     -  d7*ij21*(ij21+ij22+ij23) - d9*ij31*(ij31+ij32+ij33)
      ke(5,1)=-(d2*(ij11+ij12+ij13)*ij21) - d7*ij11*(ij21+ij22+ij23)
      ke(6,1)=-(d2*(ij11+ij12+ij13)*ij31) - d9*ij11*(ij31+ij32+ij33)
      ke(7,1)=-(d1*ij12*(ij11+ij12+ij13)) -
     -  d7*ij22*(ij21+ij22+ij23) - d9*ij32*(ij31+ij32+ij33)
      ke(8,1)=-(d2*(ij11+ij12+ij13)*ij22) - d7*ij12*(ij21+ij22+ij23)
      ke(9,1)=-(d2*(ij11+ij12+ij13)*ij32) - d9*ij12*(ij31+ij32+ij33)
      ke(10,1)=-(d1*ij13*(ij11+ij12+ij13)) -
     -  d7*ij23*(ij21+ij22+ij23) - d9*ij33*(ij31+ij32+ij33)
      ke(11,1)=-(d2*(ij11+ij12+ij13)*ij23) - d7*ij13*(ij21+ij22+ij23)
      ke(12,1)=-(d2*(ij11+ij12+ij13)*ij33) - d9*ij13*(ij31+ij32+ij33)
      ke(2,2)=d7*(ij11+ij12+ij13)**2+d1*(ij21+ij22+ij23)**2+
     -  d8*(ij31+ij32+ij33)**2
      ke(3,2)=(d2+d8)*(ij21+ij22+ij23)*(ij31+ij32+ij33)
      ke(4,2)=-(d7*(ij11+ij12+ij13)*ij21) - d2*ij11*(ij21+ij22+ij23)
      ke(5,2)=-(d7*ij11*(ij11+ij12+ij13)) -
     -  d1*ij21*(ij21+ij22+ij23) - d8*ij31*(ij31+ij32+ij33)
      ke(6,2)=-(d2*(ij21+ij22+ij23)*ij31) - d8*ij21*(ij31+ij32+ij33)
      ke(7,2)=-(d7*(ij11+ij12+ij13)*ij22) - d2*ij12*(ij21+ij22+ij23)
      ke(8,2)=-(d7*ij12*(ij11+ij12+ij13)) -
     -  d1*ij22*(ij21+ij22+ij23) - d8*ij32*(ij31+ij32+ij33)
      ke(9,2)=-(d2*(ij21+ij22+ij23)*ij32) - d8*ij22*(ij31+ij32+ij33)
      ke(10,2)=-(d7*(ij11+ij12+ij13)*ij23) - d2*ij13*(ij21+ij22+ij23)
      ke(11,2)=-(d7*ij13*(ij11+ij12+ij13)) -
     -  d1*ij23*(ij21+ij22+ij23) - d8*ij33*(ij31+ij32+ij33)
      ke(12,2)=-(d2*(ij21+ij22+ij23)*ij33) - d8*ij23*(ij31+ij32+ij33)
      ke(3,3)=d9*(ij11+ij12+ij13)**2+d8*(ij21+ij22+ij23)**2+
     -  d1*(ij31+ij32+ij33)**2
      ke(4,3)=-(d9*(ij11+ij12+ij13)*ij31) - d2*ij11*(ij31+ij32+ij33)
      ke(5,3)=-(d8*(ij21+ij22+ij23)*ij31) - d2*ij21*(ij31+ij32+ij33)
      ke(6,3)=-(d9*ij11*(ij11+ij12+ij13)) -
     -  d8*ij21*(ij21+ij22+ij23) - d1*ij31*(ij31+ij32+ij33)
      ke(7,3)=-(d9*(ij11+ij12+ij13)*ij32) - d2*ij12*(ij31+ij32+ij33)
      ke(8,3)=-(d8*(ij21+ij22+ij23)*ij32) - d2*ij22*(ij31+ij32+ij33)
      ke(9,3)=-(d9*ij12*(ij11+ij12+ij13)) -
     -  d8*ij22*(ij21+ij22+ij23) - d1*ij32*(ij31+ij32+ij33)
      ke(10,3)=-(d9*(ij11+ij12+ij13)*ij33) - d2*ij13*(ij31+ij32+ij33)
      ke(11,3)=-(d8*(ij21+ij22+ij23)*ij33) - d2*ij23*(ij31+ij32+ij33)
      ke(12,3)=-(d9*ij13*(ij11+ij12+ij13)) -
     -  d8*ij23*(ij21+ij22+ij23) - d1*ij33*(ij31+ij32+ij33)
      ke(4,4)=d1*ij11**2+d7*ij21**2+d9*ij31**2
      ke(5,4)=(d2+d7)*ij11*ij21
      ke(6,4)=(d2+d9)*ij11*ij31
      ke(7,4)=d1*ij11*ij12+d7*ij21*ij22+d9*ij31*ij32
      ke(8,4)=d7*ij12*ij21+d2*ij11*ij22
      ke(9,4)=d9*ij12*ij31+d2*ij11*ij32
      ke(10,4)=d1*ij11*ij13+d7*ij21*ij23+d9*ij31*ij33
      ke(11,4)=d7*ij13*ij21+d2*ij11*ij23
      ke(12,4)=d9*ij13*ij31+d2*ij11*ij33
      ke(5,5)=d7*ij11**2+d1*ij21**2+d8*ij31**2
      ke(6,5)=(d2+d8)*ij21*ij31
      ke(7,5)=d2*ij12*ij21+d7*ij11*ij22
      ke(8,5)=d7*ij11*ij12+d1*ij21*ij22+d8*ij31*ij32
      ke(9,5)=d8*ij22*ij31+d2*ij21*ij32
      ke(10,5)=d2*ij13*ij21+d7*ij11*ij23
      ke(11,5)=d7*ij11*ij13+d1*ij21*ij23+d8*ij31*ij33
      ke(12,5)=d8*ij23*ij31+d2*ij21*ij33
      ke(6,6)=d9*ij11**2+d8*ij21**2+d1*ij31**2
      ke(7,6)=d2*ij12*ij31+d9*ij11*ij32
      ke(8,6)=d2*ij22*ij31+d8*ij21*ij32
      ke(9,6)=d9*ij11*ij12+d8*ij21*ij22+d1*ij31*ij32
      ke(10,6)=d2*ij13*ij31+d9*ij11*ij33
      ke(11,6)=d2*ij23*ij31+d8*ij21*ij33
      ke(12,6)=d9*ij11*ij13+d8*ij21*ij23+d1*ij31*ij33
      ke(7,7)=d1*ij12**2+d7*ij22**2+d9*ij32**2
      ke(8,7)=(d2+d7)*ij12*ij22
      ke(9,7)=(d2+d9)*ij12*ij32
      ke(10,7)=d1*ij12*ij13+d7*ij22*ij23+d9*ij32*ij33
      ke(11,7)=d7*ij13*ij22+d2*ij12*ij23
      ke(12,7)=d9*ij13*ij32+d2*ij12*ij33
      ke(8,8)=d7*ij12**2+d1*ij22**2+d8*ij32**2
      ke(9,8)=(d2+d8)*ij22*ij32
      ke(10,8)=d2*ij13*ij22+d7*ij12*ij23
      ke(11,8)=d7*ij12*ij13+d1*ij22*ij23+d8*ij32*ij33
      ke(12,8)=d8*ij23*ij32+d2*ij22*ij33
      ke(9,9)=d9*ij12**2+d8*ij22**2+d1*ij32**2
      ke(10,9)=d2*ij13*ij32+d9*ij12*ij33
      ke(11,9)=d2*ij23*ij32+d8*ij22*ij33
      ke(12,9)=d9*ij12*ij13+d8*ij22*ij23+d1*ij32*ij33
      ke(10,10)=d1*ij13**2+d7*ij23**2+d9*ij33**2
      ke(11,10)=(d2+d7)*ij13*ij23
      ke(12,10)=(d2+d9)*ij13*ij33
      ke(11,11)=d7*ij13**2+d1*ij23**2+d8*ij33**2
      ke(12,11)=(d2+d8)*ij23*ij33
      ke(12,12)=d9*ij13**2+d8*ij23**2+d1*ij33**2
c
      do i=1,12
	do j=1,i
      ke(i,j)=ke(i,j)/6./detj
      enddo
      enddo
c
      do i=2,12
	do j=1,i-1
      ke(j,i)=ke(i,j)
      enddo
      enddo
c
      end
c
      subroutine nontetra4_KG(rKc,rGc,xx,ke)
      implicit none
      real*8 rKc,rGc
      real*8 young,gxy,gyz,gxz,rnyu,xx(4,3),ke(12,12)
c
      real*8 detj,ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33
      real*8 d1,d2,d7,d8,d9
      integer i,j,ii
      d1=1/3.*(4*rGc+3*rKc)
      d2=1/3.*(-2*rGc+3*rKc)
      d7=rGc
      d8=rGc
      d9=rGc
c
      do i=2,4
      do ii=1,3
      xx(i,ii)=xx(i,ii)-xx(1,ii)
      enddo
      enddo
      do ii=1,3
      xx(1,ii)=0.
      enddo
c
      detj=-(xx(2,3)*xx(3,2)*xx(4,1))+xx(2,2)*xx(3,3)*xx(4,1)+
     -  xx(2,3)*xx(3,1)*xx(4,2)-xx(2,1)*xx(3,3)*xx(4,2)-
     -  xx(2,2)*xx(3,1)*xx(4,3)+xx(2,1)*xx(3,2)*xx(4,3)
c
      ij11=-(xx(3,3)*xx(4,2))+xx(3,2)*xx(4,3)
	ij12=xx(2,3)*xx(4,2)-xx(2,2)*xx(4,3)
	ij13=-(xx(2,3)*xx(3,2))+xx(2,2)*xx(3,3)
	ij21=xx(3,3)*xx(4,1)-xx(3,1)*xx(4,3)
	ij22=-(xx(2,3)*xx(4,1))+xx(2,1)*xx(4,3)
	ij23=xx(2,3)*xx(3,1)-xx(2,1)*xx(3,3)
	ij31=-(xx(3,2)*xx(4,1))+xx(3,1)*xx(4,2)
	ij32=xx(2,2)*xx(4,1)-xx(2,1)*xx(4,2)
	ij33=-(xx(2,2)*xx(3,1))+xx(2,1)*xx(3,2)
c
      ke(1,1)=d1*(ij11+ij12+ij13)**2+d7*(ij21+ij22+ij23)**2+
     -  d9*(ij31+ij32+ij33)**2
      ke(2,1)=(d2+d7)*(ij11+ij12+ij13)*(ij21+ij22+ij23)
      ke(3,1)=(d2+d9)*(ij11+ij12+ij13)*(ij31+ij32+ij33)
      ke(4,1)=-(d1*ij11*(ij11+ij12+ij13)) -
     -  d7*ij21*(ij21+ij22+ij23) - d9*ij31*(ij31+ij32+ij33)
      ke(5,1)=-(d2*(ij11+ij12+ij13)*ij21) - d7*ij11*(ij21+ij22+ij23)
      ke(6,1)=-(d2*(ij11+ij12+ij13)*ij31) - d9*ij11*(ij31+ij32+ij33)
      ke(7,1)=-(d1*ij12*(ij11+ij12+ij13)) -
     -  d7*ij22*(ij21+ij22+ij23) - d9*ij32*(ij31+ij32+ij33)
      ke(8,1)=-(d2*(ij11+ij12+ij13)*ij22) - d7*ij12*(ij21+ij22+ij23)
      ke(9,1)=-(d2*(ij11+ij12+ij13)*ij32) - d9*ij12*(ij31+ij32+ij33)
      ke(10,1)=-(d1*ij13*(ij11+ij12+ij13)) -
     -  d7*ij23*(ij21+ij22+ij23) - d9*ij33*(ij31+ij32+ij33)
      ke(11,1)=-(d2*(ij11+ij12+ij13)*ij23) - d7*ij13*(ij21+ij22+ij23)
      ke(12,1)=-(d2*(ij11+ij12+ij13)*ij33) - d9*ij13*(ij31+ij32+ij33)
      ke(2,2)=d7*(ij11+ij12+ij13)**2+d1*(ij21+ij22+ij23)**2+
     -  d8*(ij31+ij32+ij33)**2
      ke(3,2)=(d2+d8)*(ij21+ij22+ij23)*(ij31+ij32+ij33)
      ke(4,2)=-(d7*(ij11+ij12+ij13)*ij21) - d2*ij11*(ij21+ij22+ij23)
      ke(5,2)=-(d7*ij11*(ij11+ij12+ij13)) -
     -  d1*ij21*(ij21+ij22+ij23) - d8*ij31*(ij31+ij32+ij33)
      ke(6,2)=-(d2*(ij21+ij22+ij23)*ij31) - d8*ij21*(ij31+ij32+ij33)
      ke(7,2)=-(d7*(ij11+ij12+ij13)*ij22) - d2*ij12*(ij21+ij22+ij23)
      ke(8,2)=-(d7*ij12*(ij11+ij12+ij13)) -
     -  d1*ij22*(ij21+ij22+ij23) - d8*ij32*(ij31+ij32+ij33)
      ke(9,2)=-(d2*(ij21+ij22+ij23)*ij32) - d8*ij22*(ij31+ij32+ij33)
      ke(10,2)=-(d7*(ij11+ij12+ij13)*ij23) - d2*ij13*(ij21+ij22+ij23)
      ke(11,2)=-(d7*ij13*(ij11+ij12+ij13)) -
     -  d1*ij23*(ij21+ij22+ij23) - d8*ij33*(ij31+ij32+ij33)
      ke(12,2)=-(d2*(ij21+ij22+ij23)*ij33) - d8*ij23*(ij31+ij32+ij33)
      ke(3,3)=d9*(ij11+ij12+ij13)**2+d8*(ij21+ij22+ij23)**2+
     -  d1*(ij31+ij32+ij33)**2
      ke(4,3)=-(d9*(ij11+ij12+ij13)*ij31) - d2*ij11*(ij31+ij32+ij33)
      ke(5,3)=-(d8*(ij21+ij22+ij23)*ij31) - d2*ij21*(ij31+ij32+ij33)
      ke(6,3)=-(d9*ij11*(ij11+ij12+ij13)) -
     -  d8*ij21*(ij21+ij22+ij23) - d1*ij31*(ij31+ij32+ij33)
      ke(7,3)=-(d9*(ij11+ij12+ij13)*ij32) - d2*ij12*(ij31+ij32+ij33)
      ke(8,3)=-(d8*(ij21+ij22+ij23)*ij32) - d2*ij22*(ij31+ij32+ij33)
      ke(9,3)=-(d9*ij12*(ij11+ij12+ij13)) -
     -  d8*ij22*(ij21+ij22+ij23) - d1*ij32*(ij31+ij32+ij33)
      ke(10,3)=-(d9*(ij11+ij12+ij13)*ij33) - d2*ij13*(ij31+ij32+ij33)
      ke(11,3)=-(d8*(ij21+ij22+ij23)*ij33) - d2*ij23*(ij31+ij32+ij33)
      ke(12,3)=-(d9*ij13*(ij11+ij12+ij13)) -
     -  d8*ij23*(ij21+ij22+ij23) - d1*ij33*(ij31+ij32+ij33)
      ke(4,4)=d1*ij11**2+d7*ij21**2+d9*ij31**2
      ke(5,4)=(d2+d7)*ij11*ij21
      ke(6,4)=(d2+d9)*ij11*ij31
      ke(7,4)=d1*ij11*ij12+d7*ij21*ij22+d9*ij31*ij32
      ke(8,4)=d7*ij12*ij21+d2*ij11*ij22
      ke(9,4)=d9*ij12*ij31+d2*ij11*ij32
      ke(10,4)=d1*ij11*ij13+d7*ij21*ij23+d9*ij31*ij33
      ke(11,4)=d7*ij13*ij21+d2*ij11*ij23
      ke(12,4)=d9*ij13*ij31+d2*ij11*ij33
      ke(5,5)=d7*ij11**2+d1*ij21**2+d8*ij31**2
      ke(6,5)=(d2+d8)*ij21*ij31
      ke(7,5)=d2*ij12*ij21+d7*ij11*ij22
      ke(8,5)=d7*ij11*ij12+d1*ij21*ij22+d8*ij31*ij32
      ke(9,5)=d8*ij22*ij31+d2*ij21*ij32
      ke(10,5)=d2*ij13*ij21+d7*ij11*ij23
      ke(11,5)=d7*ij11*ij13+d1*ij21*ij23+d8*ij31*ij33
      ke(12,5)=d8*ij23*ij31+d2*ij21*ij33
      ke(6,6)=d9*ij11**2+d8*ij21**2+d1*ij31**2
      ke(7,6)=d2*ij12*ij31+d9*ij11*ij32
      ke(8,6)=d2*ij22*ij31+d8*ij21*ij32
      ke(9,6)=d9*ij11*ij12+d8*ij21*ij22+d1*ij31*ij32
      ke(10,6)=d2*ij13*ij31+d9*ij11*ij33
      ke(11,6)=d2*ij23*ij31+d8*ij21*ij33
      ke(12,6)=d9*ij11*ij13+d8*ij21*ij23+d1*ij31*ij33
      ke(7,7)=d1*ij12**2+d7*ij22**2+d9*ij32**2
      ke(8,7)=(d2+d7)*ij12*ij22
      ke(9,7)=(d2+d9)*ij12*ij32
      ke(10,7)=d1*ij12*ij13+d7*ij22*ij23+d9*ij32*ij33
      ke(11,7)=d7*ij13*ij22+d2*ij12*ij23
      ke(12,7)=d9*ij13*ij32+d2*ij12*ij33
      ke(8,8)=d7*ij12**2+d1*ij22**2+d8*ij32**2
      ke(9,8)=(d2+d8)*ij22*ij32
      ke(10,8)=d2*ij13*ij22+d7*ij12*ij23
      ke(11,8)=d7*ij12*ij13+d1*ij22*ij23+d8*ij32*ij33
      ke(12,8)=d8*ij23*ij32+d2*ij22*ij33
      ke(9,9)=d9*ij12**2+d8*ij22**2+d1*ij32**2
      ke(10,9)=d2*ij13*ij32+d9*ij12*ij33
      ke(11,9)=d2*ij23*ij32+d8*ij22*ij33
      ke(12,9)=d9*ij12*ij13+d8*ij22*ij23+d1*ij32*ij33
      ke(10,10)=d1*ij13**2+d7*ij23**2+d9*ij33**2
      ke(11,10)=(d2+d7)*ij13*ij23
      ke(12,10)=(d2+d9)*ij13*ij33
      ke(11,11)=d7*ij13**2+d1*ij23**2+d8*ij33**2
      ke(12,11)=(d2+d8)*ij23*ij33
      ke(12,12)=d9*ij13**2+d8*ij23**2+d1*ij33**2
c
      do i=1,12
	do j=1,i
      ke(i,j)=ke(i,j)/6./detj
      enddo
      enddo
c
      do i=2,12
	do j=1,i-1
      ke(j,i)=ke(i,j)
      enddo
      enddo
c
      end
c_______________________________________________________________________
      subroutine tetra4me(rho,xx,me)
      real*8 rho,xx(4,3),me(4,4)
c
      real*8 detj
c
      do i=2,4
      do ii=1,3
      xx(i,ii)=xx(i,ii)-xx(1,ii)
      enddo
      enddo
      do ii=1,3
      xx(1,ii)=0.
      enddo
c
      detj=-(xx(2,3)*xx(3,2)*xx(4,1))+xx(2,2)*xx(3,3)*xx(4,1)+
     -  xx(2,3)*xx(3,1)*xx(4,2)-xx(2,1)*xx(3,3)*xx(4,2)-
     -  xx(2,2)*xx(3,1)*xx(4,3)+xx(2,1)*xx(3,2)*xx(4,3)
c
      me=0.
      do i=1,4
      me(i,i)=detj*rho/6. /4.
      enddo
      end
