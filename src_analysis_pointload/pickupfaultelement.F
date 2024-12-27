c_______________________________________________________________________
      subroutine pickupelem(n,ne,coor,cny,cort,iet,r)
      implicit none
      integer n,ne,iet
      real*8 coor(3,n)
      integer cny(10,ne)
      real*8 cort(3),vec(3,3),vol,r(3)
      real*8 vol1,vol2,vol3,vol4,volt
      integer ie,ii,i1,i2,i3,i4
      real*8 damy

      do ie=1,ne
      i1=cny(1,ie)
      i2=cny(2,ie)
      i3=cny(3,ie)
      i4=cny(4,ie)
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-coor(ii,i1)
      vec(2,ii)=coor(ii,i3)-coor(ii,i1)
      vec(3,ii)=coor(ii,i4)-coor(ii,i1)
      enddo
      call tet4volume(vec,vol)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-coor(ii,i1)
      vec(2,ii)=coor(ii,i3)-coor(ii,i1)
      vec(3,ii)=cort(ii)-coor(ii,i1)
      enddo
      call tet4volume(vec,vol1)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-coor(ii,i1)
      vec(2,ii)=cort(ii)-coor(ii,i1)
      vec(3,ii)=coor(ii,i4)-coor(ii,i1)
      enddo
      call tet4volume(vec,vol2)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-cort(ii)
      vec(2,ii)=coor(ii,i3)-cort(ii)
      vec(3,ii)=coor(ii,i4)-cort(ii)
      enddo
      call tet4volume(vec,vol3)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i3)-coor(ii,i1)
      vec(2,ii)=coor(ii,i4)-coor(ii,i1)
      vec(3,ii)=cort(ii)-coor(ii,i1)
      enddo
      call tet4volume(vec,vol4)
c
      volt=vol1+vol2+vol3+vol4
c
      if(volt/vol.le.1.0000005)then
      iet=ie
      r(1)=vol4/vol
      r(2)=vol2/vol
      r(3)=vol1/vol
      goto 201
      endif
c
      enddo
      iet=0
!      write(6,*) 'failed'
!      stop

  201 damy=0.

      end
c_______________________________________________________________________
      subroutine tet4volume(vec,vol)
      real*8 vec(3,3)  ! vec(i,j): i-th vector, j-th component
      real*8 vol,outer(3)
c
c outer: vec(1,*) vec(2,*)
      outer(1)=-(vec(1,3)*vec(2,2)) + vec(1,2)*vec(2,3)
      outer(2)=vec(1,3)*vec(2,1) - vec(1,1)*vec(2,3)
      outer(3)=-(vec(1,2)*vec(2,1)) + vec(1,1)*vec(2,2)
c inner: outer vec(3,*)
      vol=1/6.*
     -abs(outer(1)*vec(3,1)+outer(2)*vec(3,2)+outer(3)*vec(3,3))
c
      end
