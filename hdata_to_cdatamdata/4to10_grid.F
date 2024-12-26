!      implicit none
!      real*8 cri
!      integer n,ne,coornewsize
!
!      open(11,file='./data/tet4setting.dat',status='old')
!      read(11,*) 
!      read(11,*) n
!      read(11,*) 
!      read(11,*) ne
!      close(11)
!
!!--- coornewsize: number of nodes in tet10 mesh
!!----- at pessimistic case, coornewsize=ne*10 (none of the nodes are shared among elements)
!!      coornewsize=ne*10
!!----- normally, we can assume that at least 3 elements share one node in average
!      coornewsize=(ne*10)/3
!!---
!      cri=0.001
!
!      call main(n,ne,cri,coornewsize)
!
!      end
c_______________________________________________________________________
!      subroutine main(n,ne,cri,coornewsize)
!      use omp_lib
!      implicit none
!      integer n,ne,ii,nl,i,nxg,nyg,nzg,nnew
!      integer conn4(5*ne)
!      real*8 coor4(3*n)
!      integer conn10(11*ne)
!      integer fl,coornewsize
!      real*8 coor10(3*coornewsize)
!      real*8 cri,cmm(2,3),dsg,co(3)
!      real*8 t0,t1,t2,t3,t4,t5,tt0
!
!      t0=omp_get_wtime()
!      write(*,*) 'starting 4to10_grid_pickupBCnode_hdf.exe'
!      call read_tet4_geometry_hdf(n,ne,coor4,conn4)
!      t1=omp_get_wtime()
!      write(*,*) 'read_tet4_geometry_hdf() took',t1-t0
!
!      do ii=1,3
!      cmm(1,ii)=coor4(ii)
!      cmm(2,ii)=coor4(ii)
!      enddo
!      do i=1,n
!      do ii=1,3
!      if(cmm(1,ii).ge.coor4(3*(i-1)+ii))then
!      cmm(1,ii)=coor4(3*(i-1)+ii)
!      endif
!      if(cmm(2,ii).le.coor4(3*(i-1)+ii))then
!      cmm(2,ii)=coor4(3*(i-1)+ii)
!      endif
!      enddo
!      enddo
!
!      open(21,file='./data/modeling_setting.dat',status='old')
!      read(21,*)	
!      read(21,*)	
!      read(21,*)	
!      read(21,*) dsg	
!      close(21)
!!      dsg=dsg*2
!
!      nxg=int((cmm(2,1)-cmm(1,1))/dsg)+1 +2 +2
!      nyg=int((cmm(2,2)-cmm(1,2))/dsg)+1 +2 +2
!      nzg=int((cmm(2,3)-cmm(1,3))/dsg)+1 +2 +2
!      write(*,*) 'grids',nxg,nyg,nzg
!      do ii=1,3
!      co(ii)=cmm(1,ii)-2*dsg
!      enddo
!      nl=20 !!! DEFAULT
!      t2=omp_get_wtime()
!      write(*,*) 'prep() took',t2-t1
!
! 501  fl=1
!      nl=nl+10 !!! DEFAULT
!      write(6,*) 'max num of grid list', nl
!      tt0=omp_get_wtime()
!      call make4to10_grid
!     - (n,nnew,ne,conn4,conn10,coor4,coor10,cri,nxg,nyg,nzg,nl,
!     - dsg,co,fl,coornewsize)
!      write(*,*) 'make4to10_grid() loop took',omp_get_wtime()-tt0
!      if(fl.eq.0)then
!      goto 501
!      endif
!      t3=omp_get_wtime()
!      write(*,*) 'make4to10_grid() took',t3-t2
!
!      open(11,file='./data/setting.dat',status='unknown')
!      write(11,*) 'num of node'
!      write(11,*) nnew
!      write(11,*) 'num of tet elem'
!      write(11,*) ne
!      close(11)
!
!      call pickupBCnode_lib(nnew,coornewsize,coor10)
!      t4=omp_get_wtime()
!      write(*,*) 'pickupBCnode_lib() took',t4-t3
!      
!      call write_tet10_geometry(nnew,ne,coornewsize,coor10,conn10)
!      t5=omp_get_wtime()
!      write(*,*) 'write_tet10_geometry() took',t5-t4
!      write(*,*) '4to10_grid_pickupBCnode_hdf.exe took',t5-t0
!
!      end
c_______________________________________________________________________
      subroutine write_tet10_geometry(n,ne,coorsize,coor10,conn10)
      implicit none
      integer*8 n,ne,coorsize,isize
      integer*8 conn10(11*ne)
      real*8 coor10(3*coorsize)
      character filename*50
      character dataname*50
      integer*8 fid,coorbufsize,connbufsize

      write(filename,'(a20)')
     - './data/tet10model.h5'
      isize=20
      call set_null_char(filename,isize)
      fid=0
      call hdf_create_file(filename,fid)

      write(dataname,'(a5)') '/coor'
      isize=5
      call set_null_char(dataname,isize)
      coorbufsize=n*3
      call hdf_write_double_array(fid,dataname,coor10,coorbufsize)

      write(dataname,'(a5)') '/conn'
      isize=5
      call set_null_char(dataname,isize)
      connbufsize=ne*11
      call hdf_write_long_array(fid,dataname,conn10,connbufsize)

      call hdf_close_file(fid)

      end
c_______________________________________________________________________
!      subroutine make4to10_grid
!     - (n,nnew,ne,conn4,conn10,coor4,coor10,cri,nxg,nyg,nzg,nl,
!     - dsg,co,fl,coornewsize)
!      implicit none
!      integer i,j,nnew,ii,k,ix,iy,iz,k1,j1,i1,nm1,nm2,i2
!      integer n,ne,nxg,nyg,nzg,nl
!      integer conn4(5*ne),conn10(11*ne),cnyele(5:10,2),coornewsize
!      real*8 coor4(3*n),coor10(3*coornewsize)
!      real*8 coorlocal(10,3),cri,dsg,co(3),damy
!      integer fl,maxcount
!      integer*4 listnum(nxg,nyg,nzg)
!      integer listn(nl,nxg,nyg,nzg)
!      integer istart, iend, jstart, jend, kstart, kend
!
!      do i=1,ne
!      do j=1,4
!      conn10(11*(i-1)+j)=conn4(5*(i-1)+j)
!      enddo
!      conn10(11*(i-1)+11)=conn4(5*(i-1)+5)
!      enddo
!      do i=1,3*n
!      coor10(i)=coor4(i)
!      enddo
!      nnew=n
!
!      cnyele(5,1)=1
!      cnyele(5,2)=2
!      cnyele(6,1)=2
!      cnyele(6,2)=3
!      cnyele(7,1)=1
!      cnyele(7,2)=3
!      cnyele(8,1)=1
!      cnyele(8,2)=4
!      cnyele(9,1)=2
!      cnyele(9,2)=4
!      cnyele(10,1)=4
!      cnyele(10,2)=3
!
!      listnum=0
!
!! add tet4 nodes to grid
!!--- commented out because of error during analysis
!!      do i=1,n
!!      ix=int((coor(i,1)-co(1))/dsg)+1
!!      iy=int((coor(i,2)-co(2))/dsg)+1
!!      iz=int((coor(i,3)-co(3))/dsg)+1	
!!      if(listnum(ix,iy,iz)+1.gt.nl)then
!!      fl=0      
!!      goto 102
!!      endif
!!      listnum(ix,iy,iz)=listnum(ix,iy,iz)+1
!!      listn(ix,iy,iz,listnum(ix,iy,iz))=i
!!      enddo
!!---
!	
!      do i=1,ne
!      do j=1,4
!      do ii=1,3
!      coorlocal(j,ii)=coor10(3*(conn10(11*(i-1)+j)-1)+ii)
!      enddo
!      enddo
!      do ii=1,3
!      do k=5,10
!      coorlocal(k,ii)=
!     - 0.5*(coorlocal(cnyele(k,1),ii)+coorlocal(cnyele(k,2),ii))
!      enddo
!      enddo
!
!      do k=5,10
!      ix=int((coorlocal(k,1)-co(1))/dsg)+1
!      if(ix + 1 - ((coorlocal(k,1)-co(1))/dsg + 1) < 10e-4) then
!        istart = ix
!        iend   = ix + 1
!      elseif ((coorlocal(k,1)-co(1))/dsg + 1 - ix < 10e-4) then
!        istart = ix - 1
!        iend   = ix
!      else
!        istart = ix
!        iend   = ix 
!      endif
!
!      iy=int((coorlocal(k,2)-co(2))/dsg)+1
!      if(iy + 1 - ((coorlocal(k,1)-co(1))/dsg + 1) < 10e-4) then
!        jstart = iy
!        jend   = iy + 1
!      elseif ((coorlocal(k,1)-co(1))/dsg + 1 - iy < 10e-4) then
!        jstart = iy - 1
!        jend   = iy 
!      else
!        jstart = iy
!        jend   = iy 
!      endif
!
!      iz=int((coorlocal(k,3)-co(3))/dsg)+1   
!      if(iz + 1 - ((coorlocal(k,1)-co(1))/dsg + 1) < 10e-4) then
!        kstart = iz
!        kend   = iz + 1
!      elseif ((coorlocal(k,1)-co(1))/dsg + 1 - iz < 10e-4) then
!        kstart = iz - 1
!        kend   = iz 
!      else
!        kstart = iz
!        kend   = iz 
!      endif
!
!      do k1=kstart,kend !iz-1,iz+1
!      do j1=jstart,jend !iy-1,iy+1
!      do i1=istart,iend !ix-1,ix+1
!      nm1=listnum(i1,j1,k1)
!
!      if(nm1.ne.0)then
!      do i2=1,nm1
!      nm2=listn(i2,i1,j1,k1)
!      if(
!     - (abs(coorlocal(k,1)-coor10(3*(nm2-1)+1)).lt.cri).and.
!     - (abs(coorlocal(k,2)-coor10(3*(nm2-1)+2)).lt.cri).and.
!     - (abs(coorlocal(k,3)-coor10(3*(nm2-1)+3)).lt.cri))then
!      conn10(11*(i-1)+k)=nm2
!      goto 101
!      endif
!      enddo   
!      endif
!      enddo
!      enddo
!      enddo
!      call check_increment(nnew,coornewsize)
!      conn10(11*(i-1)+k)=nnew
!      do ii=1,3
!      coor10(3*(nnew-1)+ii)=coorlocal(k,ii)
!      enddo
!
!      if(listnum(ix,iy,iz)+1.gt.nl)then
!      fl=0      
!      goto 102
!      endif
!      listnum(ix,iy,iz)=listnum(ix,iy,iz)+1
!      listn(listnum(ix,iy,iz),ix,iy,iz)=nnew
!  101 damy=0.
!      enddo
!      enddo
!  102 damy=0.
!
!      maxcount=0
!      do iz=1,nzg
!      do iy=1,nyg
!      do ix=1,nxg
!      if(maxcount.lt.listnum(ix,iy,iz))then
!      maxcount=listnum(ix,iy,iz)
!      endif
!      enddo
!      enddo
!      enddo
!      write(*,*) 'nl,maxcount',nl,maxcount
!
!      end
c_______________________________________________________________________
      subroutine read_tet4_geometry_hdf_open(n,ne,coorsize,coor4,conn4)
      implicit none
      integer*8 n,ne,coorsize,isize
      real*8 coor4(3*coorsize)
      integer*8 conn4(5*ne)
      character filename*50
      character dataname*50
      integer*8 bufsize,bufcount
      integer*8 fid      

      write(filename,'(a17)') 'data/tet4model.h5'
      isize=17
      call set_null_char(filename,isize)
      fid=0
      call hdf_open_file(filename,fid)

      write(dataname,'(a5)') '/coor'
      isize=5
      call set_null_char(dataname,isize)
      bufsize=n*3
      call hdf_read_double_array(fid,dataname,coor4,bufsize,bufcount)

      write(dataname,'(a5)') '/conn'
      isize=5
      call set_null_char(dataname,isize)
      bufsize=ne*5
      call hdf_read_long_array(fid,dataname,conn4,bufsize,bufcount)

      call hdf_close_file(fid)

      end
c_______________________________________________________________________
!      subroutine make4to10(n,nnew,ne,cny,cnynew,coor,coornew,cri)
!      integer cny(ne,4),cnynew(ne,10),cnyele(5:10,2)
!      real*8 coor(n,3),coornew(ne*10,3)
!      real*8 coorlocal(10,3)
!      real*8 cri
!
!c initial setting
!      do i=1,ne
!      do j=1,10
!      cnynew(i,j)=0
!      enddo
!      enddo
!      do i=1,ne
!      do j=1,4
!      cnynew(i,j)=cny(i,j)
!      enddo
!      enddo
!      do i=1,3*n
!      do j=1,3
!      coornew(i,j)=0
!      enddo
!      enddo
!      do i=1,n
!      do j=1,3
!      coornew(i,j)=coor(i,j)
!      enddo
!      enddo
!      nnew=n
!c
!      cnyele(5,1)=1
!      cnyele(5,2)=2
!      cnyele(6,1)=2
!      cnyele(6,2)=3
!      cnyele(7,1)=1
!      cnyele(7,2)=3
!      cnyele(8,1)=1
!      cnyele(8,2)=4
!      cnyele(9,1)=2
!      cnyele(9,2)=4
!      cnyele(10,1)=4
!      cnyele(10,2)=3
!c
!      do i=1,ne
!c
!      do j=1,4
!      do ii=1,3
!      coorlocal(j,ii)=coornew(cnynew(i,j),ii)
!      enddo
!      enddo
!      do ii=1,3
!      do k=5,10
!      coorlocal(k,ii)=
!     -0.5*(coorlocal(cnyele(k,1),ii)+coorlocal(cnyele(k,2),ii))
!      enddo
!      enddo
!c
!      do k=5,10
!      do j=1,nnew
!      if(
!     - (abs(coorlocal(k,1)-coornew(j,1)).lt.cri).and.
!     - (abs(coorlocal(k,2)-coornew(j,2)).lt.cri).and.
!     - (abs(coorlocal(k,3)-coornew(j,3)).lt.cri))then
!      cnynew(i,k)=j
!      goto 101
!      endif
!      enddo
!      nnew=nnew+1
!      cnynew(i,k)=nnew
!      do ii=1,3
!      coornew(nnew,ii)=coorlocal(k,ii)
!      enddo
!
!  101 write(6,*) nnew
!      enddo
!c
!      enddo
!c
!      end
c_______________________________________________________________________
      subroutine check_increment(count,size)
      implicit none
      integer*8 count,size

      if(count.gt.(size-1))then
      write(*,*)'accessing out of array',size
      stop
      endif
      count=count+1

      end
c_______________________________________________________________________
      subroutine read_settings(n,ne,kd,dsg)
      implicit none
      integer*8 n,ne,kd
      real*8 dsg

      open(11,file='./data/tet4setting.dat',status='old')
      read(11,*) 
      read(11,*) n
      read(11,*) 
      read(11,*) ne
      read(11,*) 
      read(11,*) kd
      close(11)
      
      open(21,file='./data/modeling_setting.dat',status='old')
      read(21,*)	
      read(21,*)	
      read(21,*)	
      read(21,*) dsg	
      close(21)
!      dsg=dsg*2

      end
c_______________________________________________________________________
      subroutine set_gridsize(n4,coorsize,coor4,dsg,nxg,nyg,nzg,co)
      implicit none
      integer*8 n4,nxg,nyg,nzg,i,ii,coorsize
      real*8 co(3),cmm(2,3),dsg
      real*8 coor4(coorsize*3)

      do ii=1,3
      cmm(1,ii)=coor4(ii)
      cmm(2,ii)=coor4(ii)
      enddo
      do i=1,n4
      do ii=1,3
      if(cmm(1,ii).ge.coor4(3*(i-1)+ii))then
      cmm(1,ii)=coor4(3*(i-1)+ii)
      endif
      if(cmm(2,ii).le.coor4(3*(i-1)+ii))then
      cmm(2,ii)=coor4(3*(i-1)+ii)
      endif
      enddo
      enddo

      nxg=int((cmm(2,1)-cmm(1,1))/dsg)+1 +2 +2
      nyg=int((cmm(2,2)-cmm(1,2))/dsg)+1 +2 +2
      nzg=int((cmm(2,3)-cmm(1,3))/dsg)+1 +2 +2
#ifdef OUTPUT_LOG
      write(*,'(a7,i8,i8,i8)') 'grids: ',nxg,nyg,nzg
#endif
      do ii=1,3
      co(ii)=cmm(1,ii)-2*dsg
      enddo

      end
c_______________________________________________________________________
      subroutine write_settings(n10,ne,kd)
      implicit none
      integer*8 n10,ne,kd

      open(11,file='./data/setting.dat',status='unknown')
      write(11,*) 'num of node'
      write(11,*) n10
      write(11,*) 'num of tet elem'
      write(11,*) ne
      write(11,*) 'num of material properties'
      write(11,*) kd
      close(11)

      end
c_______________________________________________________________________
      subroutine convert_4to10_grid(n4,ne,coorsize,coor,conn,nl,dsg,
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
      integer*4 listn(nl,nxg,nyg,nzg)
      integer*8 istart, iend, jstart, jend, kstart, kend

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
      write(*,*) 'nl,maxcnt,x,y,z',nl,maxcount,nxg,nyg,nzg
#endif
  102 damy=0.

      end
c_______________________________________________________________________
