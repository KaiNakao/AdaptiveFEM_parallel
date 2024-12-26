      subroutine read_tet10_geometry_hdf(n,ne,coor,cny,num,fid)
      implicit none
      integer*4 n,ne,i,j
      real*8 coor(n,3)
      integer*4 cny(ne,10),num(ne)
      character dataname*50
      integer*4 connbuf(ne*11)
      real*8 coorbuf(n*3)
      integer*8 bufsize,bufcount,fid,tmp
 
      write(dataname,'(a5)') '/coor'
      tmp=5
      call set_null_char(dataname,tmp)
      bufsize=n*3
      call hdf_read_double_array(fid,dataname,coorbuf,bufsize,bufcount)
      
      write(dataname,'(a5)') '/conn'
      call set_null_char(dataname,tmp)
      bufsize=ne*11
      call hdf_read_int_array(fid,dataname,connbuf,bufsize,bufcount)

      do i=1,n
      do j=1,3
      coor(i,j)=coorbuf(3*(i-1)+j)
      enddo
      enddo
            
      do i=1,ne
      do j=1,10
      cny(i,j)=connbuf(11*(i-1)+j)
      enddo
      num(i)=connbuf(11*(i-1)+11)
      enddo
      
      end
