#define NSIZE 50
c_______________________________________________________________________
      subroutine ABC4_lib(n,ne,kd,ib,fid)
      implicit none
      integer*4 n,ne,kd,ib,i,ie,in,ine
      integer*8 fid,isize
      real*8 coor(n,3),rmat(kd,10)
      integer*4 cny(ne,4),num(ne),ibc(ib),nf(n),ino(6)
      real*8 co1,co2,cm(3,3),xmin,ymin
      character dataname*50
      integer*4 buf10(n*NSIZE)
      real*8 buf11(n*NSIZE)
      real*8 buf12(n*NSIZE)
      real*8 buf13(n*NSIZE)
      integer*4 buf14(n*NSIZE)
      integer*8 bufs(10:14),bufc(10:14),one
      real*8 domainsetting(5),cri

      open(42,file='./data/modeldomain.dat',status='unknown')
      read(42,*)
      read(42,*) domainsetting(1)
      read(42,*)
      read(42,*) domainsetting(2)
      read(42,*)
      read(42,*) domainsetting(3)
      read(42,*)
      read(42,*) domainsetting(4)
      read(42,*)
      read(42,*) domainsetting(5)
      close(42)
      cri=0.001

      bufs(10)=n*NSIZE
      bufs(11)=n*NSIZE
      bufs(12)=n*NSIZE
      bufs(13)=n*NSIZE
      bufs(14)=n*NSIZE

      bufc=0
      one=1
      
      call read_tet4_geometry_hdf(n,ne,coor,cny,num,fid)
      call MPI_ABC4_read_material(kd,rmat)
      call MPI_ABC4_read_BC(ib,ibc,fid)
      call MPI_ABC4_def_weight(cm)
c 
      xmin=1.0e+10
      ymin=1.0e+10
      do i=1,n
      xmin=dmin1(xmin,coor(i,1))
      ymin=dmin1(ymin,coor(i,2))
      enddo
c
      ine=0
c_______________________________________________________________________
      
      nf=0
      do i=1,n
      if(abs(domainsetting(1)-coor(i,1)).le.cri)then        
      nf(i)=1
      endif
      enddo
      
      do ie=1,ne
      co2=rmat(num(ie),3)*rmat(num(ie),2)
      co1=rmat(num(ie),3)*rmat(num(ie),1)-co2
c 1 5 2 8 9 4: 1 2 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 2 6 3 9 10 4: 2 3 4
      ino(1)=cny(ie,2)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 7 3 8 10 4: 3 1 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 5 2 7 6 3: 2 1 3
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,3)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
      enddo
c_______________________________________________________________________
      nf=0
      do i=1,n
      if(abs(domainsetting(2)-coor(i,1)).le.cri)then        
      nf(i)=1
      endif
      enddo
      
      do ie=1,ne
      co2=rmat(num(ie),3)*rmat(num(ie),2)
      co1=rmat(num(ie),3)*rmat(num(ie),1)-co2
c 1 5 2 8 9 4: 1 2 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 2 6 3 9 10 4: 2 3 4
      ino(1)=cny(ie,2)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 7 3 8 10 4: 3 1 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 5 2 7 6 3: 2 1 3
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,3)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
      enddo
c_______________________________________________________________________
      nf=0
      do i=1,n
      if(abs(domainsetting(3)-coor(i,2)).le.cri)then        
      nf(i)=1
      endif
      enddo
      
      do ie=1,ne
      co2=rmat(num(ie),3)*rmat(num(ie),2)
      co1=rmat(num(ie),3)*rmat(num(ie),1)-co2
c 1 5 2 8 9 4: 1 2 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 2 6 3 9 10 4: 2 3 4
      ino(1)=cny(ie,2)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 7 3 8 10 4: 3 1 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 5 2 7 6 3: 2 1 3
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,3)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
      enddo
c_______________________________________________________________________
      nf=0
      do i=1,n
      if(abs(domainsetting(4)-coor(i,2)).le.cri)then        
      nf(i)=1
      endif
      enddo
      
      do ie=1,ne
      co2=rmat(num(ie),3)*rmat(num(ie),2)
      co1=rmat(num(ie),3)*rmat(num(ie),1)-co2
c 1 5 2 8 9 4: 1 2 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 2 6 3 9 10 4: 2 3 4
      ino(1)=cny(ie,2)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 7 3 8 10 4: 3 1 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 5 2 7 6 3: 2 1 3
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,3)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
      enddo
c_______________________________________________________________________
      nf=0
      do i=1,n
      if(abs(domainsetting(5)-coor(i,3)).le.cri)then        
      nf(i)=1
      endif
      enddo
      
      do ie=1,ne
      co2=rmat(num(ie),3)*rmat(num(ie),2)
      co1=rmat(num(ie),3)*rmat(num(ie),1)-co2
c 1 5 2 8 9 4: 1 2 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 2 6 3 9 10 4: 2 3 4
      ino(1)=cny(ie,2)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 7 3 8 10 4: 3 1 4
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,3)
      ino(3)=cny(ie,4)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
c 1 5 2 7 6 3: 2 1 3
      ino(1)=cny(ie,1)
      ino(2)=cny(ie,2)
      ino(3)=cny(ie,3)
      in=nf(ino(1))+nf(ino(2))+nf(ino(3))
      if(in.eq.3)then
      call MPI_ABC4_add_abc(n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      endif
      enddo
c_______________________________________________________________________

c
      write(dataname,'(a8)') '/ABCnode'
      isize=8
      call set_null_char(dataname,isize)
      call hdf_write_int_array(fid,dataname,buf10,bufc(10))
     
      write(dataname,'(a7)') '/ABCcoe'
      isize=7
      call set_null_char(dataname,isize)
      call hdf_write_double_array(fid,dataname,buf11,bufc(11))
     
      write(dataname,'(a10)') '/Stresscoe'
      isize=10
      call set_null_char(dataname,isize)
      call hdf_write_double_array(fid,dataname,buf12,bufc(12))
     
      write(dataname,'(a7)') '/normal'
      isize=7
      call set_null_char(dataname,isize)
      call hdf_write_double_array(fid,dataname,buf13,bufc(13))
     
      write(dataname,'(a8)') '/abcelem'
      isize=8
      call set_null_char(dataname,isize)
      call hdf_write_int_array(fid,dataname,buf14,bufc(14))

      write(dataname,'(a11)') '/ABCsetting'
      isize=11
      call set_null_char(dataname,isize)
      call hdf_write_int_array(fid,dataname,ine,one,one)

      end
c_______________________________________________________________________
      subroutine MPI_ABC4_add_abc
     - (n,coor,ino,co1,co2,cm,ine,xmin,ymin,ie2,
     - bufs,bufc,buf10,buf11,buf12,buf13,buf14)
      implicit none
      integer*4 n,ie2,j,i,ii,ine
      real*8 xx(3,3),area,norm(3),coor(n,3),co1,co2,cm(3,3)
      real*8 xmin,ymin
      integer*4 ino(6)
      integer*4 buf10(n*NSIZE)
      real*8 buf11(n*NSIZE)
      real*8 buf12(n*NSIZE)
      real*8 buf13(n*NSIZE)
      integer*4 buf14(n*NSIZE)
      integer*8 bufs(10:14),bufc(10:14)
c
      call check_increment(bufc(14),bufs(14))
      buf14(bufc(14))=ie2
      do j=1,3
      do ii=1,3
      xx(j,ii)=coor(ino(j),ii)
      enddo
      enddo
      call MPI_ABC4_calc_area_norm(xx,area,norm)
      if(abs(abs(norm(1)+norm(2)+norm(3))-1).le.0.1)then
      ine=ine+1
      do i=1,3
      call check_increment(bufc(10),bufs(10))
      buf10(bufc(10))=ino(i)
      enddo
      do ii=1,3
      do j=1,3
      do i=1,3
      call check_increment(bufc(11),bufs(11))
      buf11(bufc(11))=cm(i,j)*(co2+norm(ii)*co1)*abs(area)*2
      enddo
      enddo
      enddo
      do j=1,3
      do i=1,3
      call check_increment(bufc(12),bufs(12))
      buf12(bufc(12))=cm(i,j)*abs(area)*2
      enddo
      enddo
c
      if(norm(3).gt.0.9)then
      norm(3)=-norm(3)
      endif
      if((norm(1).gt.0.9).and.(abs(xmin-coor(ino(1),1))).le.0.01)then
      norm(1)=-norm(1)
      endif
      if((norm(2).gt.0.9).and.(abs(ymin-coor(ino(1),2))).le.0.01)then
      norm(2)=-norm(2)
      endif
c
      do i=1,3
      call check_increment(bufc(13),bufs(13))
      buf13(bufc(13))=norm(i)
      enddo
      endif
      
      end
c_______________________________________________________________________
      subroutine MPI_ABC4_calc_area_norm(xx,area,norm)
      implicit none
      integer*4 i,ii
      real*8 xx(3,3)
      real*8 inner,amag,bmag,cos,sin,area,norm(3)
c
      do i=2,3
      do ii=1,3
      xx(i,ii)=xx(i,ii)-xx(1,ii)
      enddo
      enddo
c
      inner=0.
      do ii=1,3
      inner=inner+xx(2,ii)*xx(3,ii)
      enddo
      amag=0.
      do ii=1,3
      amag=amag+xx(2,ii)*xx(2,ii)
      enddo
      amag=sqrt(amag)
      bmag=0.
      do ii=1,3
      bmag=bmag+xx(3,ii)*xx(3,ii)
      enddo
      bmag=sqrt(bmag)
c
      cos=inner/amag/bmag
      sin=sqrt(1-cos**2)
c
      area=amag*bmag*sin/2.
c
      norm(1)=-(xx(2,3)*xx(3,2)) + xx(2,2)*xx(3,3)
      norm(2)=xx(2,3)*xx(3,1) - xx(2,1)*xx(3,3)
      norm(3)= -(xx(2,2)*xx(3,1)) + xx(2,1)*xx(3,2)
c
      amag=sqrt(norm(1)*norm(1)+norm(2)*norm(2)+norm(3)*norm(3))
      do i=1,3
      norm(i)=abs(norm(i)/amag)
      enddo
c
      end
c_______________________________________________________________________
      subroutine MPI_ABC4_read_BC(ib,ibc,fid)
      implicit none
      integer*4 ibc(ib),i,ib,intbuf(ib+1)
      integer*8 fid,bufsize,bufcount,isize
      character dataname*50

      write(dataname,'(a3)') '/BC'
      isize=3
      call set_null_char(dataname,isize)
      bufsize=ib+1
      call hdf_read_int_array(fid,dataname,intbuf,bufsize,bufcount)
      do i=1,ib
      ibc(i)=intbuf(i+1)
      enddo      
      
      end
c_______________________________________________________________________
      subroutine MPI_ABC4_read_material(kd,rmat)
      implicit none
      integer*4 kd,i,j
      real*8 rmat(kd,10)
      
      print *, "kd: ", kd
      open(42,file='./data/material.dat',status='unknown')
      do j=1,kd
      read(42,*)
      do i=1,10
      read(42,*) rmat(j,i)
      enddo
      enddo
      close(42)
      
      end
c_______________________________________________________________________
      subroutine MPI_ABC4_def_weight(sm)
      implicit none
      real*8 sm(3,3)
      
c 1,2,3,5,6,7
      sm(1,1)=1/6.0
      sm(1,2)=0.0
      sm(1,3)=0.0
      sm(2,1)=0.0
      sm(2,2)=1/6.
      sm(2,3)=0.0
      sm(3,1)=0.0
      sm(3,2)=0.0
      sm(3,3)=1/6.0
c 150115 tmp           
      sm(1,1)=1/12.
      sm(1,2)=1/24.
      sm(1,3)=1/24.
      sm(2,1)=1/24.
      sm(2,2)=1/12.
      sm(2,3)=1/24.
      sm(3,1)=1/24.
      sm(3,2)=1/24.
      sm(3,3)=1/12.
c 150115 tmp       
      end
c_______________________________________________________________________
