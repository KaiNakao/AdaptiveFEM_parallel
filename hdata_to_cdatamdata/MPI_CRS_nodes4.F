
#define CRS_NC 800
!#define CRS_NC 50

c_______________________________________________________________________
      subroutine CRS_nodes4_lib(n,ne,fid)
      implicit none
      integer*4 n,ne,i,ie,j1,inum,ntm,j,k,i1
      integer*8 fid
      real*8 damy
 !     parameter(nc=200) 
      integer*4 nc
      real*8 coor(n,3)
      integer*4 num(ne),cny(ne,4)
!      integer mat(n,nc),matnum(n)
      integer*4 mat(n,CRS_NC),matnum(n)
      integer*4 bufsize
      
      nc=CRS_NC
c
      call read_tet4_geometry_hdf(n,ne,coor,cny,num,fid)
      matnum=0
      mat=0

      do ie=1,ne     
      do i=1,4
      i1=cny(ie,i)   
      do j=1,4
      j1=cny(ie,j)
      if(matnum(i1).eq.0)then
      mat(i1,1)=j1
      matnum(i1)=1
      else
      do k=1,matnum(i1)
      if(mat(i1,k).eq.j1)then
      goto 110
      endif
      enddo
      if(matnum(i1).eq.nc-1)then
      write(*,*)'error: enlarge nc'
      stop
      endif
      matnum(i1)=matnum(i1)+1
      mat(i1,matnum(i1))=j1
  110 damy=0
      endif      
      enddo 
      enddo 
      enddo

      inum=0
      do i=1,n
      inum=inum+matnum(i)
      enddo

      ntm=0
      do i=1,n
      ntm=max(ntm,matnum(i))
      enddo
      bufsize=2+n
      do i=1,n
      bufsize=bufsize+matnum(i)
      enddo
      call write_hdf_crs_nodes(ntm,inum,matnum,mat,bufsize,n,nc,fid)

      end
c_______________________________________________________________________
      subroutine write_hdf_crs_nodes
     -(ntm,inum,matnum,mat,bufsize,n,nc,fid)
      implicit none
      integer*4 ntm,inum,bufsize,i,j,n,nc
      integer*4 buf(bufsize)
      character dataname*50
      integer*4 mat(n,nc),matnum(n)
      integer*8 fid,bufcount
      
      buf(1)=ntm
      buf(2)=inum
      do i=1,n
      buf(2+i)=matnum(i)
      enddo
      bufcount=2+n
      do i=1,n
      do j=1,matnum(i)
      bufcount=bufcount+1
      buf(bufcount)=mat(i,j)
      enddo
      enddo
      
      write(dataname,'(a10)') '/crs_nodes'
      call set_null_char(dataname,10)
      call hdf_write_int_array
     -(fid,dataname,buf,bufcount)
     
      end
c_______________________________________________________________________

