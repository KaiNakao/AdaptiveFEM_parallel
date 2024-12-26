#define PAIRSIZE 80
c_______________________________________________________________________
      subroutine convert_4to10_confirming_mesh
     - (n4,ne,coorsize,coor,conn,n10)
      implicit none
      integer*8 n4,ne,n10,coorsize
      integer*8 conn(ne*11)
      real*8 coor(coorsize*3)
      integer*8 i,j,k,ii,i1,i2
      integer*8 cnyele(5:10,2)
      real*8 damy
      integer*8 pairs(2,PAIRSIZE,n4),pairs_size(n4)

      n10=n4
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
      pairs_size=0

      do i=1,ne
      do j=5,10
      i1=conn(11*(i-1)+cnyele(j,1))
      i2=conn(11*(i-1)+cnyele(j,2))
      if(i1.lt.i2)then
      ii=i1
      i1=i2
      i2=ii
      endif
      do k=1,pairs_size(i1)
      if(pairs(1,k,i1).eq.i2)then
        conn(11*(i-1)+j)=pairs(2,k,i1)
        goto 123
      endif
      enddo
      if(pairs_size(i1).eq.PAIRSIZE)then
        write(*,*) 'PAIRSIZE too small'
        stop
      endif
      pairs_size(i1)=pairs_size(i1)+1
      pairs(1,pairs_size(i1),i1)=i2
      call check_increment(n10,coorsize)
      coor(3*n10-2)=0.5*(coor(3*i1-2)+coor(3*i2-2))
      coor(3*n10-1)=0.5*(coor(3*i1-1)+coor(3*i2-1))
      coor(3*n10  )=0.5*(coor(3*i1  )+coor(3*i2  ))
      pairs(2,pairs_size(i1),i1)=n10
      conn(11*(i-1)+j)=n10
123   damy=0.0 
      enddo
      enddo
#ifdef OUTPUT_LOG
      j=0
      do i=1,n4
      j=max(j,pairs_size(n4))
      enddo
      write(*,*) 'max pairs_size',j,PAIRSIZE
#endif

      end
c_______________________________________________________________________
