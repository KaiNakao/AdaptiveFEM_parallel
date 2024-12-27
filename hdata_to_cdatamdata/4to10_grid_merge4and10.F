c_______________________________________________________________________
      subroutine convert_4to10_grid_merge4and10
     - (n4,ne,coorsize,coor,conn,nl,dsg,
     - nxg,nyg,nzg,co,n10,flag)
      implicit none
      integer*8 n4,ne,nl,nxg,nyg,nzg,n10,flag,coorsize
      integer*8 conn(ne*11)
      real*8 coor(coorsize*3)
      real*8 dsg,co(3)
      integer*8 i,j,ii,k,ix,iy,iz,k1,j1,i1,nm1,nm2,i2
      integer*8 cnyele(5:10,2)
      real*8 coorlocal(10,3),cri,damy
      integer*8 maxcount
      integer*4 listnum(nxg,nyg,nzg)
      integer*8 listn(nl,nxg,nyg,nzg)
      integer*8 istart, iend, jstart, jend, kstart, kend
      real*8 xcoor,ycoor,zcoor
      integer*8 ntmp

      flag=1
      n10=n4
      cri=0.001

      cnyele(5,1)=1
      cnyele(5,2)=2
      cnyele(6,1)=2
      cnyele(6,2)=3
      cnyele(7,1)=1
      cnyele(7,2)=3
      cnyele(8,1)=1
      cnyele(8,2)=4
      cnyele(9,1)=2
      cnyele(9,2)=4
      cnyele(10,1)=4
      cnyele(10,2)=3

      listnum=0
      
! add tet4 nodes to grid
      ntmp=0
      do k=1,n4
      xcoor=coor(3*(k-1)+1)
      ycoor=coor(3*(k-1)+2)
      zcoor=coor(3*(k-1)+3)
      
      ix=int((xcoor-co(1))/dsg)+1
      if(ix + 1 - ((xcoor-co(1))/dsg + 1) < 10e-4) then
        istart = ix
        iend   = ix + 1
      elseif ((xcoor-co(1))/dsg + 1 - ix < 10e-4) then
        istart = ix - 1
        iend   = ix
      else
        istart = ix
        iend   = ix 
      endif

      iy=int((ycoor-co(2))/dsg)+1
      if(iy + 1 - ((ycoor-co(2))/dsg + 1) < 10e-4) then
        jstart = iy
        jend   = iy + 1
      elseif ((ycoor-co(2))/dsg + 1 - iy < 10e-4) then
        jstart = iy - 1
        jend   = iy 
      else
        jstart = iy
        jend   = iy 
      endif

      iz=int((zcoor-co(3))/dsg)+1   
      if(iz + 1 - ((zcoor-co(3))/dsg + 1) < 10e-4) then
        kstart = iz
        kend   = iz + 1
      elseif ((zcoor-co(3))/dsg + 1 - iz < 10e-4) then
        kstart = iz - 1
        kend   = iz 
      else
        kstart = iz
        kend   = iz 
      endif

      do k1=kstart,kend !iz-1,iz+1
      do j1=jstart,jend !iy-1,iy+1
      do i1=istart,iend !ix-1,ix+1
      nm1=listnum(i1,j1,k1)

      if(nm1.ne.0)then
      do i2=1,nm1
      nm2=listn(i2,i1,j1,k1)
      if(
     - (abs(xcoor-coor(3*(nm2-1)+1)).lt.cri).and.
     - (abs(ycoor-coor(3*(nm2-1)+2)).lt.cri).and.
     - (abs(zcoor-coor(3*(nm2-1)+3)).lt.cri))then
      write(*,*) 'nodes with same position in n4-nodeset'
      stop
!      goto 104
      endif
      enddo   
      endif
      enddo
      enddo
      enddo
      ntmp=ntmp+1
      if(listnum(ix,iy,iz)+1.gt.nl)then
      flag=0      
      goto 102
      endif
      listnum(ix,iy,iz)=listnum(ix,iy,iz)+1
      listn(listnum(ix,iy,iz),ix,iy,iz)=ntmp
!  104 damy=0.
      enddo
! end add tet4 nodes to grid
      do i=1,ne
      do j=1,4
      do ii=1,3
      coorlocal(j,ii)=coor(3*(conn(11*(i-1)+j)-1)+ii)
      enddo
      enddo
      do ii=1,3
      do k=5,10
      coorlocal(k,ii)=
     - 0.5*(coorlocal(cnyele(k,1),ii)+coorlocal(cnyele(k,2),ii))
      enddo
      enddo

      do k=5,10
      ix=int((coorlocal(k,1)-co(1))/dsg)+1
      if(ix + 1 - ((coorlocal(k,1)-co(1))/dsg + 1) < 10e-4) then
        istart = ix
        iend   = ix + 1
      elseif ((coorlocal(k,1)-co(1))/dsg + 1 - ix < 10e-4) then
        istart = ix - 1
        iend   = ix
      else
        istart = ix
        iend   = ix 
      endif

      iy=int((coorlocal(k,2)-co(2))/dsg)+1
      if(iy + 1 - ((coorlocal(k,2)-co(2))/dsg + 1) < 10e-4) then
        jstart = iy
        jend   = iy + 1
      elseif ((coorlocal(k,2)-co(2))/dsg + 1 - iy < 10e-4) then
        jstart = iy - 1
        jend   = iy 
      else
        jstart = iy
        jend   = iy 
      endif

      iz=int((coorlocal(k,3)-co(3))/dsg)+1   
      if(iz + 1 - ((coorlocal(k,3)-co(3))/dsg + 1) < 10e-4) then
        kstart = iz
        kend   = iz + 1
      elseif ((coorlocal(k,3)-co(3))/dsg + 1 - iz < 10e-4) then
        kstart = iz - 1
        kend   = iz 
      else
        kstart = iz
        kend   = iz 
      endif

      do k1=kstart,kend !iz-1,iz+1
      do j1=jstart,jend !iy-1,iy+1
      do i1=istart,iend !ix-1,ix+1
      nm1=listnum(i1,j1,k1)

      if(nm1.ne.0)then
      do i2=1,nm1
      nm2=listn(i2,i1,j1,k1)
      if(
     - (abs(coorlocal(k,1)-coor(3*(nm2-1)+1)).lt.cri).and.
     - (abs(coorlocal(k,2)-coor(3*(nm2-1)+2)).lt.cri).and.
     - (abs(coorlocal(k,3)-coor(3*(nm2-1)+3)).lt.cri))then
      conn(11*(i-1)+k)=nm2
      goto 101
      endif
      enddo   
      endif
      enddo
      enddo
      enddo
      call check_increment(n10,coorsize)
      conn(11*(i-1)+k)=n10
      do ii=1,3
      coor(3*(n10-1)+ii)=coorlocal(k,ii)
      enddo

      if(listnum(ix,iy,iz)+1.gt.nl)then
      flag=0      
      goto 102
      endif
      listnum(ix,iy,iz)=listnum(ix,iy,iz)+1
      listn(listnum(ix,iy,iz),ix,iy,iz)=n10
  101 damy=0.
      enddo
      enddo

      maxcount=0
      do iz=1,nzg
      do iy=1,nyg
      do ix=1,nxg
      if(maxcount.lt.listnum(ix,iy,iz))then
      maxcount=listnum(ix,iy,iz)
      endif
      enddo
      enddo
      enddo
#ifdef OUTPUT_LOG
      write(*,*) 'nl,maxcount',nl,maxcount
#endif
  102 damy=0.

      end
c_______________________________________________________________________
