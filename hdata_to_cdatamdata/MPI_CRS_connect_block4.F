c_______________________________________________________________________
      subroutine CRS_connect_block4_lib(n,ne,ib,fid)
      implicit none
      integer*4 n,ne,ib
      integer*4 ntm,nab
      integer*8 fid,one,zero,intbufsize,realbufsize
      character dataname*50

      one=1
      zero=0
      
      write(dataname,'(a10)') '/crs_nodes'
      call set_null_char(dataname,10)
      call hdf_read_int_array_part(fid,dataname,zero,one,ntm)
c
      write(dataname,'(a11)') '/ABCsetting'
      call set_null_char(dataname,11)
      call hdf_read_int_array_part(fid,dataname,zero,one,nab)
c
      intbufsize=n*100 !n*100
      realbufsize=n*40 !n*30
      
      call main_CRS_connect_block4
     - (n,ne,ntm,ib,nab,fid,intbufsize,realbufsize)
      
      end
c_______________________________________________________________________
      subroutine main_CRS_connect_block4
     - (n,ne,ntm,ib,nab,fid,intbufsize,realbufsize)
      implicit none
      integer*4 n,ne,ntm,ib,nab
      integer*4 i,j,j2,j1,ie,i1,i2,i3,i4,i5,itmp,i_val
      real*8 damy
      real*8 coor(n,3)
      integer*4 num(ne),cny(ne,4)
!      real*8 abc(nab,3,6,6)
      integer*4 nabc(nab,6)
      integer*4 ibc(ib)
      integer*4 numcrs(n,ntm),liscrs(n)
      integer*4 liscrssum(n)
      character dataname*50
      integer*8 intbufsize,realbufsize,bufcount,fid
      integer*4 intbuf(intbufsize)
!      real*8 realbuf(realbufsize)
      

      call read_tet4_geometry_hdf(n,ne,coor,cny,num,fid)
      call MPI_ABC4_read_BC(ib,ibc,fid)

      write(dataname,'(a8)') '/ABCnode'
      call set_null_char(dataname,8)
      call hdf_read_int_array(fid,dataname,
     - intbuf,intbufsize,bufcount)
      do i=1,nab
      do j=1,3
      nabc(i,j)=intbuf(3*(i-1)+j)
      enddo
      enddo
      
!      write(dataname,'(a7)') '/ABCcoe'
!      call set_null_char(dataname,7)
!      call hdf_read_double_array(fid,dataname,
!     - realbuf,realbufsize,bufcount)
!      bufcount=0
!      do i=1,nab
!      do ii=1,3
!      do j2=1,3
!      do j1=1,3
!      bufcount=bufcount+1
!      abc(i,ii,j1,j2)=realbuf(bufcount)
!      enddo
!      enddo
!      enddo
!      enddo

      numcrs=0
      write(dataname,'(a10)') '/crs_nodes'
      call set_null_char(dataname,10)
      call hdf_read_int_array(fid,dataname,
     - intbuf,intbufsize,bufcount)
      do i=1,n
      liscrs(i)=intbuf(i+2)
      enddo
      bufcount=2+n
      do i=1,n
      do j=1,liscrs(i)
      bufcount=bufcount+1
      numcrs(i,j)=intbuf(bufcount)
      enddo
      enddo
c
      liscrssum=0
      do i=1,n-1
      liscrssum(i+1)=liscrssum(i)+liscrs(i)
      enddo
      bufcount=1
      intbuf(1)=4*4*ne
      
      do ie=1,ne
      do i1=1,4
      i2=cny(ie,i1)
      do i3=1,4
      i4=cny(ie,i3)
      do i5=1,liscrs(i2)
      if(numcrs(i2,i5).eq.i4)then
      itmp=i5
      goto 252
      endif      
      enddo    
  252 damy=0
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=liscrssum(i2)+itmp
      enddo
      enddo
      enddo
      write(dataname,'(a12)') '/crs_connect'
      call set_null_char(dataname,12)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      bufcount=1
      intbuf(1)=3*3*nab
      
      do i=1,nab
      do j2=1,3
      i2=nabc(i,j2)
      do j1=1,3
      i4=nabc(i,j1)
      do i5=1,liscrs(i2)
      if(numcrs(i2,i5).eq.i4)then
      itmp=i5
      goto 251
      endif      
      enddo    
  251 damy=0     
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=liscrssum(i2)+itmp
	enddo
      enddo
      enddo
      
      write(dataname,'(a16)') '/crs_ABC_connect'
      call set_null_char(dataname,16)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      intbuf(1)=n+1
      i_val=0
      intbuf(2)=1
      do i=1,n
      i_val=i_val+liscrs(i)
      intbuf(2+i)=i_val+1
      enddo
      
      bufcount=2+n-1
      call check_increment(bufcount,intbufsize)
      write(dataname,'(a8)') '/crs_ptr'
      call set_null_char(dataname,8)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)
   
      intbuf(1)=i_val
      bufcount=1
      do i=1,n
      do j=1,liscrs(i)
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=numcrs(i,j)
      enddo
      enddo
      write(dataname,'(a8)') '/crs_ind'
      call set_null_char(dataname,8)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      bufcount=0
      do i=1,n
      do j=1,liscrs(i)
      if(i.eq.numcrs(i,j))then
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=liscrssum(i)+j
      endif 
      enddo
      enddo     
      write(dataname,'(a15)') '/crs_block_list'
      call set_null_char(dataname,15)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      end	   
c_______________________________________________________________________
