      subroutine check_macro_parameters
      implicit none
      logical :: something_wrong
      something_wrong = .false.
      write(6, *) "npc:", npc
      write(6, *) "nkl", nkl
      write(6, *) "nsta", nsta
      if(npc == 1)then
        if(nsta /= 0)then
          write(6, *) "nsta must be 0, but is",nsta
          something_wrong = .true.
        endif
        if(nkl /= 0)then
          write(6, *) "nkl must be 0, but is",nkl
          something_wrong = .true.
        endif
      elseif(npc == 21)then
        if(nsta /= 5)then
          write(6, *) "nsta must be 5, but is",nsta
          something_wrong = .true.
        endif
        if(nkl /= 5)then
          write(6, *) "nkl must be 5, but is",nkl
          something_wrong = .true.
        endif
      elseif(npc == 288)then
        if(nsta /= 5)then
          write(6, *) "nsta must be 5, but is",nsta
          something_wrong = .true.
        endif
        if(nkl /= 10)then
          write(6, *) "nkl must be 10, but is",nkl
          something_wrong = .true.
        endif
      else
        something_wrong = .true.
      endif
      if(something_wrong)then
        stop "Something is wrong with the macro parameters"
      endif
      end subroutine check_macro_parameters
c_______________________________________________________________________
      subroutine compute_nsize(n,n_size)
      implicit none
      INCLUDE 'mpif.h'
      integer*4 ierr,n,n_size,i4_1
      i4_1=1
!      call MPI_ALLREDUCE(n,n_size,i4_1,MPI_INTEGER,
!     - MPI_MAX,MPI_COMM_WORLD,ierr)
      n_size=n
      end
c_______________________________________________________________________

      subroutine write_tet_vtk_s
     - (n_size,ne_size,n,ne,tettype,cny,coor,u,im,outtype)
      implicit none
      integer im,i,tettype,outtype
      integer n,ne,n_size,ne_size
      REAL_4 u(3*n_size)
      REAL_4 coor(3,n_size)
      integer cny(tettype,ne_size)
      character filename*50

      if(tettype.eq.4)then
      write(filename,'(a7,i6.6,a1,i6.6,a4)'),
     - './vtk4_',outtype,'.',im,'.vtk'
      else
      write(filename,'(a8,i6.6,a1,i6.6,a4)'),
     - './vtk10_',outtype,'.',im,'.vtk'
      endif
      open(109,file=filename,status='unknown')
      write(109,'(a26)') '# vtk DataFile Version 2.0'
      write(109,'(a4)')  'file'
      write(109,'(a5)')  'ASCII'
      write(109,'(a25)') 'DATASET UNSTRUCTURED_GRID'

      write(109,'(a7,i9.9,a6)') 'POINTS ',n,' float'
      do i=1,n
      write(109,*) coor(1,i),coor(2,i),coor(3,i)
      enddo

      write(109,'(a6,i9.9,a1,i9.9)') 'CELLS ',ne,' ',ne*(tettype+1)
      if(tettype.eq.4)then
      do i=1,ne
      write(109,*) '4 ',cny(1,i)-1,cny(2,i)-1,
     - cny(3,i)-1,cny(4,i)-1
      enddo
      else
      do i=1,ne
      write(109,'(a,10i)') '10 ',cny(1,i)-1,cny(2,i)-1,
     - cny(3,i)-1,cny(4,i)-1,
     - cny(5,i)-1,cny(6,i)-1,
     - cny(7,i)-1,cny(8,i)-1,
     - cny(9,i)-1,cny(10,i)-1
      enddo
      endif

      write(109,'(a11,i9.9)') 'CELL_TYPES ',ne
      if(tettype.eq.4)then
      do i=1,ne
      write(109,'(a2)') '10'
      enddo
      else
      do i=1,ne
      write(109,'(a2)') '24'
      enddo
      endif

      write(109,'(a11,i9.9)') 'POINT_DATA ',n
      write(109,'(a26)') 'VECTORS displacement float'
      do i=1,n
      write(109,*) u(3*(i-1)+1),u(3*(i-1)+2),u(3*(i-1)+3)
      enddo

      close(109)
      end
c_______________________________________________________________________
      subroutine write_tet_tri_vtk_s
     - (n_size,ne,nei,cny,cnyi,coor,u,im,outtype)
      implicit none
      integer im,i,tettype,outtype
      integer nei,ne,n_size
      REAL_4 u(3*n_size)
      REAL_4 coor(3,n_size)
      integer cny(10,ne)
      integer cnyi(12,nei)
      character filename*50

      write(filename,'(a8,i6.6,a1,i6.6,a4)'),
     - './vtk12_',outtype,'.',im,'.vtk'
      open(109,file=filename,status='unknown')
      write(109,'(a26)') '# vtk DataFile Version 2.0'
      write(109,'(a4)')  'file'
      write(109,'(a5)')  'ASCII'
      write(109,'(a25)') 'DATASET UNSTRUCTURED_GRID'

      write(109,'(a7,i9.9,a6)') 'POINTS ',n_size,' float'
      do i=1,n_size
      write(109,*) coor(1,i),coor(2,i),coor(3,i)
      enddo

#ifdef OUTPUTALL
      write(109,'(a6,i9.9,a1,i9.9)') 'CELLS ',ne+nei,' ',ne*11+nei*7
      do i=1,ne
      write(109,*) '10 ',cny(1,i)-1,cny(2,i)-1,
     - cny(3,i)-1,cny(4,i)-1,
     - cny(5,i)-1,cny(6,i)-1,
     - cny(7,i)-1,cny(8,i)-1,
     - cny(9,i)-1,cny(10,i)-1
      enddo
#else
      write(109,'(a6,i9.9,a1,i9.9)') 'CELLS ',nei,' ',nei*7
#endif
      do i=1,nei
      write(109,*) '6 ',cnyi(7,i)-1,cnyi(8,i)-1,
     - cnyi(9,i)-1,cnyi(10,i)-1,
     - cnyi(11,i)-1,cnyi(12,i)-1
      enddo

#ifdef OUTPUTALL
      write(109,'(a11,i9.9)') 'CELL_TYPES ',ne+nei
      do i=1,ne
      write(109,'(a2)') '24'
      enddo
#else
      write(109,'(a11,i9.9)') 'CELL_TYPES ',nei
#endif
      do i=1,nei
      write(109,'(a2)') '22'
      enddo

      write(109,'(a11,i9.9)') 'POINT_DATA ',n_size
      write(109,'(a26)') 'VECTORS displacement float'
      do i=1,n_size
      write(109,*) u(3*(i-1)+1),u(3*(i-1)+2),u(3*(i-1)+3)
      enddo

      close(109)
      end
c_______________________________________________________________________
      subroutine write_tet4_tri_vtk_s
     - (n_size,ne,nei,cny,cnyi,coor,u,im,outtype)
      implicit none
      integer im,i,tettype,outtype
      integer nei,ne,n_size
      REAL_4 u(3*n_size)
      REAL_4 coor(3,n_size)
      integer cny(4,ne)
      integer cnyi(6,nei)
      character filename*50

      write(filename,'(a7,i6.6,a1,i6.6,a4)'),
     - './vtk6_',outtype,'.',im,'.vtk'
      open(109,file=filename,status='unknown')
      write(109,'(a26)') '# vtk DataFile Version 2.0'
      write(109,'(a4)')  'file'
      write(109,'(a5)')  'ASCII'
      write(109,'(a25)') 'DATASET UNSTRUCTURED_GRID'

      write(109,'(a7,i9.9,a6)') 'POINTS ',n_size,' float'
      do i=1,n_size
      write(109,*) coor(1,i),coor(2,i),coor(3,i)
      enddo

#ifdef OUTPUTALL
      write(109,'(a6,i9.9,a1,i9.9)') 'CELLS ',ne+nei,' ',ne*5+nei*4
      do i=1,ne
      write(109,*)'4 ',cny(1,i)-1,cny(2,i)-1,cny(3,i)-1,cny(4,i)-1
      enddo
#else
      write(109,'(a6,i9.9,a1,i9.9)') 'CELLS ',nei,' ',nei*4
#endif
      do i=1,nei
      write(109,*)'3 ',cnyi(4,i)-1,cnyi(5,i)-1,cnyi(6,i)-1
      enddo

#ifdef OUTPUTALL
      write(109,'(a11,i9.9)') 'CELL_TYPES ',ne+nei
      do i=1,ne
      write(109,'(a2)') '10'
      enddo
#else
      write(109,'(a11,i9.9)') 'CELL_TYPES ',nei
#endif
      do i=1,nei
      write(109,'(a2)') '5'
      enddo

      write(109,'(a11,i9.9)') 'POINT_DATA ',n_size
      write(109,'(a26)') 'VECTORS displacement float'
      do i=1,n_size
      write(109,*) u(3*(i-1)+1),u(3*(i-1)+2),u(3*(i-1)+3)
      enddo

      close(109)
      end
c_______________________________________________________________________
      subroutine setmmat_s(me)
      REAL_4 me(10,10)

      me(1,1)=0.002380952380952381
      me(2,1)=0.0003968253968253968
      me(2,2)=0.002380952380952381
      me(3,1)=0.0003968253968253968
      me(3,2)=0.0003968253968253968
      me(3,3)=0.002380952380952381
      me(4,1)=0.0003968253968253968
      me(4,2)=0.0003968253968253968
      me(4,3)=0.0003968253968253968
      me(4,4)=0.002380952380952381
      me(5,1)=-0.0015873015873015873
      me(5,2)=-0.0015873015873015873
      me(5,3)=-0.002380952380952381
      me(5,4)=-0.002380952380952381
      me(5,5)=0.012698412698412698
      me(6,1)=-0.002380952380952381
      me(6,2)=-0.0015873015873015873
      me(6,3)=-0.0015873015873015873
      me(6,4)=-0.002380952380952381
      me(6,5)=0.006349206349206349
      me(6,6)=0.012698412698412698
      me(7,1)=-0.0015873015873015873
      me(7,2)=-0.002380952380952381
      me(7,3)=-0.0015873015873015873
      me(7,4)=-0.002380952380952381
      me(7,5)=0.006349206349206349
      me(7,6)=0.006349206349206349
      me(7,7)=0.012698412698412698
      me(8,1)=-0.0015873015873015873
      me(8,2)=-0.002380952380952381
      me(8,3)=-0.002380952380952381
      me(8,4)=-0.0015873015873015873
      me(8,5)=0.006349206349206349
      me(8,6)=0.0031746031746031746
      me(8,7)=0.006349206349206349
      me(8,8)=0.012698412698412698
      me(9,1)=-0.002380952380952381
      me(9,2)=-0.0015873015873015873
      me(9,3)=-0.002380952380952381
      me(9,4)=-0.0015873015873015873
      me(9,5)=0.006349206349206349
      me(9,6)=0.006349206349206349
      me(9,7)=0.0031746031746031746
      me(9,8)=0.006349206349206349
      me(9,9)=0.012698412698412698
      me(10,1)=-0.002380952380952381
      me(10,2)=-0.002380952380952381
      me(10,3)=-0.0015873015873015873
      me(10,4)=-0.0015873015873015873
      me(10,5)=0.0031746031746031746
      me(10,6)=0.006349206349206349
      me(10,7)=0.006349206349206349
      me(10,8)=0.006349206349206349
      me(10,9)=0.006349206349206349
      me(10,10)=0.012698412698412698
      do j=1,10
      do i=1,j-1
      me(i,j)=me(j,i)
      enddo
      enddo

      end
c_______________________________________________________________________
      subroutine setmmat_d(me)
      REAL*8 me(10,10)

      me(1,1)=0.002380952380952381
      me(2,1)=0.0003968253968253968
      me(2,2)=0.002380952380952381
      me(3,1)=0.0003968253968253968
      me(3,2)=0.0003968253968253968
      me(3,3)=0.002380952380952381
      me(4,1)=0.0003968253968253968
      me(4,2)=0.0003968253968253968
      me(4,3)=0.0003968253968253968
      me(4,4)=0.002380952380952381
      me(5,1)=-0.0015873015873015873
      me(5,2)=-0.0015873015873015873
      me(5,3)=-0.002380952380952381
      me(5,4)=-0.002380952380952381
      me(5,5)=0.012698412698412698
      me(6,1)=-0.002380952380952381
      me(6,2)=-0.0015873015873015873
      me(6,3)=-0.0015873015873015873
      me(6,4)=-0.002380952380952381
      me(6,5)=0.006349206349206349
      me(6,6)=0.012698412698412698
      me(7,1)=-0.0015873015873015873
      me(7,2)=-0.002380952380952381
      me(7,3)=-0.0015873015873015873
      me(7,4)=-0.002380952380952381
      me(7,5)=0.006349206349206349
      me(7,6)=0.006349206349206349
      me(7,7)=0.012698412698412698
      me(8,1)=-0.0015873015873015873
      me(8,2)=-0.002380952380952381
      me(8,3)=-0.002380952380952381
      me(8,4)=-0.0015873015873015873
      me(8,5)=0.006349206349206349
      me(8,6)=0.0031746031746031746
      me(8,7)=0.006349206349206349
      me(8,8)=0.012698412698412698
      me(9,1)=-0.002380952380952381
      me(9,2)=-0.0015873015873015873
      me(9,3)=-0.002380952380952381
      me(9,4)=-0.0015873015873015873
      me(9,5)=0.006349206349206349
      me(9,6)=0.006349206349206349
      me(9,7)=0.0031746031746031746
      me(9,8)=0.006349206349206349
      me(9,9)=0.012698412698412698
      me(10,1)=-0.002380952380952381
      me(10,2)=-0.002380952380952381
      me(10,3)=-0.0015873015873015873
      me(10,4)=-0.0015873015873015873
      me(10,5)=0.0031746031746031746
      me(10,6)=0.006349206349206349
      me(10,7)=0.006349206349206349
      me(10,8)=0.006349206349206349
      me(10,9)=0.006349206349206349
      me(10,10)=0.012698412698412698
      do j=1,10
      do i=1,j-1
      me(i,j)=me(j,i)
      enddo
      enddo

      end

c_______________________________________________________________________
      subroutine map4to10_init
     & (im,ds,
     & n,ne,ni,nei,coor,cny,cnyi,
     & n4,ne4,ni4,nei4,coor4,cny4,cnyi4,
     & map_c_to_f)
      implicit none
      integer n,ne,n4,ne4,map_c_to_f(2,n+ni),im,nei,nei4,ni,ni4
      integer cny(10,ne),cny4(4,ne4),cnyi(12,nei),cnyi4(6,nei4)
      real*8 coor(3,n+ni),coor4(3,n4+ni4),ds
      integer i,j,ie,ii,flag(n+ni),flagarr(n4+ni4),jj,ie4,cnyele(2,12)
      real*8 x,y,z,tmp
      map_c_to_f=0

      if(im.eq.0)then
      write(*,*) 'map4to10_init'
      endif
      flag=0

      do ie=1,ne
      flag(cny(1,ie))=1
      flag(cny(2,ie))=1
      flag(cny(3,ie))=1
      flag(cny(4,ie))=1
      enddo
      j=0
      do i=1,n
      if(flag(i).eq.1)then
      j=j+1
      flagarr(j)=i
      endif
      enddo
      if(j.ne.n4)then
      write(*,*) 'something wrong?',j,n4
      stop
      endif

      do ie=1,ne
      do ii=1,4
      i=cny4(ii,ie)
      j=cny(ii,ie)
      tmp=(coor4(1,i)-coor(1,j))**2+
     &    (coor4(2,i)-coor(2,j))**2+
     &    (coor4(3,i)-coor(3,j))**2
      if(tmp.lt.(0.01*ds)**2)then
        map_c_to_f(1,j)=i
        map_c_to_f(2,j)=i
      else
        write(*,*) 'something wrong??? 1'
        stop
      endif
      enddo
      enddo

      cnyele=0
      cnyele(1,5)=1
      cnyele(2,5)=2
      cnyele(1,6)=2
      cnyele(2,6)=3
      cnyele(1,7)=1
      cnyele(2,7)=3
      cnyele(1,8)=1
      cnyele(2,8)=4
      cnyele(1,9)=2
      cnyele(2,9)=4
      cnyele(1,10)=4
      cnyele(2,10)=3
      do ie=1,ne
      do j=5,10
        map_c_to_f(1,cny(j,ie))=map_c_to_f(1,cny(cnyele(1,j),ie))
        map_c_to_f(2,cny(j,ie))=map_c_to_f(1,cny(cnyele(2,j),ie))
      enddo
      enddo

      if(nei.ne.nei4)then
      write(*,*) 'something wrong in nei,nei4',nei,nei4
      stop
      endif

      do ie=1,nei

      do ii=1,3
      i=cnyi4(ii,ie)
      j=cnyi(ii,ie)
      tmp=(coor4(1,i)-coor(1,j))**2+
     &    (coor4(2,i)-coor(2,j))**2+
     &    (coor4(3,i)-coor(3,j))**2
      if(tmp.lt.(0.01*ds)**2)then
c        map_c_to_f(1,j)=i
c        map_c_to_f(2,j)=i
      else
        write(*,*) 'something wrong??? 2'
        stop
      endif
      enddo

      do ii=4,6
      i=cnyi4(ii,ie)
      j=cnyi(ii+3,ie)
      tmp=(coor4(1,i)-coor(1,j))**2+
     &    (coor4(2,i)-coor(2,j))**2+
     &    (coor4(3,i)-coor(3,j))**2
      if(tmp.lt.(0.01*ds)**2)then
        map_c_to_f(1,j)=i
        map_c_to_f(2,j)=i
      else
        write(*,*) 'something wrong??? 3',ie
        write(*,*) coor4(1,i),coor(1,j)
        write(*,*) coor4(2,i),coor(2,j)
        write(*,*) coor4(3,i),coor(3,j)
        stop
      endif
      enddo

      enddo

c   3
c   6   5
c   1  4  2
      cnyele=0
      cnyele(1,4)=1
      cnyele(2,4)=2
      cnyele(1,5)=2
      cnyele(2,5)=3
      cnyele(1,6)=3
      cnyele(2,6)=1
      cnyele(1,10)=7
      cnyele(2,10)=8
      cnyele(1,11)=8
      cnyele(2,11)=9
      cnyele(1,12)=9
      cnyele(2,12)=7
      do ie=1,nei
      do j=4,6
c        map_c_to_f(1,cnyi(j,ie))=map_c_to_f(1,cnyi(cnyele(1,j),ie))
c        map_c_to_f(2,cnyi(j,ie))=map_c_to_f(1,cnyi(cnyele(2,j),ie))
      enddo
      do j=10,12
        map_c_to_f(1,cnyi(j,ie))=map_c_to_f(1,cnyi(cnyele(1,j),ie))
        map_c_to_f(2,cnyi(j,ie))=map_c_to_f(1,cnyi(cnyele(2,j),ie))
      enddo
      enddo

      end
c_______________________________________________________________________
      subroutine force_to_nodal_force_s(n,ne,cny,rf,r)
      implicit none
      integer n,ne,cny(10,ne)
      REAL_4 rf(npc,3,n),r(npc,3,n)
      REAL_4 me(10,10)
      integer ie,i,j,k,l
      REAL_4 tmp,arayi(10)

      call setmmat_s(me)
      r=0.0
      do ie=1,ne
        do i=1,3
          do j=1,npc
            do k=1,10
              arayi(k)=rf(j,i,cny(k,ie))
            enddo
            do k=1,10
              tmp=0.0
              do l=1,10
                tmp=tmp+me(k,l)*arayi(l)
              enddo
              r(j,i,cny(k,ie))=tmp
            enddo
          enddo
        enddo
      enddo
      end
c_______________________________________________________________________
      subroutine makmatlist(kd,rmat,younglst)
      real*8 rmat(10,kd)
      real*8 younglst(2,kd)
      real*8 c1tmp,c2tmp,rhtmp,ygtmp,rntmp
      do it=1,kd
      c1tmp=rmat(1,it)
      c2tmp=rmat(2,it)
      rhtmp=rmat(3,it)
      call mak_E_nyu(c1tmp,c2tmp,rhtmp,ygtmp,rntmp)
      younglst(1,it)=ygtmp
      younglst(2,it)=rntmp
      enddo
      end
c_______________________________________________________________________
      subroutine mak_E_nyu(c1,c2,rho,E,nyu)
      real*8 c1,c2,rho
      real*8 E,nyu
      E=((3*c1**2*c2**2-4*c2**4)*rho)/(c1**2-c2**2)
      nyu=(c1**2-2*c2**2)/(2*c1**2-2*c2**2)
c      E=c1**2*rho-4/3.*c2**2*rho
c	nyu=c2**2*rho
      end
c_______________________________________________________________________
      subroutine set_fault_input
     & (im,kd,mpibuf,younglst,
     & n_size,ne_size,n,ne,coor,cny,num,
     & nad,nadto,mpiad,mpiadlist,mpl,
     & rv
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      implicit none
      INCLUDE 'mpif.h'
      integer im,kd
      integer n_size,ne_size,n,ne
      integer cny(10,ne_size),num(ne_size)
      real*8 coor(3,n_size),younglst(2,kd)
      integer nad,nadto
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real*8 mpibuf(2*3*npc*(nadto+1))
      integer ierr,nfp,nxsplit,nysplit
      real*8 dxsplit,dysplit,originsplit(2)
      real*8 rv(3*n_size*npc)

      if(im.eq.0)then
      open(42,file='./data/faultinput_grid.dat',status='unknown')
      read(42,*) !nfp
      read(42,*) nfp
      read(42,*) !nx,ny
      read(42,*) nxsplit,nysplit
      read(42,*) !dx,dy
      read(42,*) dxsplit,dysplit
      read(42,*) !origin(x),origin(y)
      read(42,*) originsplit(1),originsplit(2)
      close(42)

      endif
      call MPI_BCAST(nfp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nxsplit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nysplit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dxsplit,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dysplit,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(originsplit,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      call set_fault_input_sub
     & (im,kd,mpibuf,younglst,
     & n_size,ne_size,n,ne,coor,cny,num,
     & nad,nadto,mpiad,mpiadlist,mpl,
     & rv,
     & nfp,nxsplit,nysplit,dxsplit,dysplit,originsplit
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      end
c_______________________________________________________________________
      subroutine set_fault_input_sub
     & (im,kd,mpibuf,younglst,
     & n_size,ne_size,n,ne,coor,cny,num,
     & nad,nadto,mpiad,mpiadlist,mpl,
     & rv,
     & nfp,nx,ny,dx,dy,origin
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      implicit none
      INCLUDE 'mpif.h'

      integer im,kd
      integer n_size,ne_size,n,ne
      integer cny(10,ne_size),num(ne_size)
      real*8 coor(3,n_size),younglst(2,kd)
      integer nad,nadto
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real*8 mpibuf(2*3*npc*(nadto+1))
      REAL_4 mpibufs(2*3*npc*(nadto+1))
      integer ierr,nfp,nx,ny
      real*8 dx,dy,origin(2)
      integer flag1(n_size),flag2(n_size)
      REAL_4 utmp(3*npc,n_size)
      real*8 slip(3,nx+1,ny+1),ratiox,ratioy,x,y,z
      real*8 maxx,maxy,dsx,dsy
      real*8 rv(3*n_size*npc)
      real*8 up(3,n_size),rnyu,young,gxy,gyz,gxz
      integer ie,in,j,i,i1,j1,i2,i3,ix,iy
      real*8 xx(4,3),arayi(30),arayo(30),detjmat

      flag1=0
      flag2=0
      do ie=1,ne
      in=num(ie)
      if(in.le.nfp)then
        do j=1,10
        flag1(cny(j,ie))=1
        enddo
      endif
      if(in.gt.nfp)then
        do j=1,10
        flag2(cny(j,ie))=1
        enddo
      endif
      enddo

      utmp=0.0
      do i=1,n
      utmp(1,i)=flag1(i)
      enddo
      call sync_parallel_block_s
     & (n_size,n,utmp,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibufs
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      do i=1,n
      flag1(i)=utmp(1,i)
      enddo

      utmp=0.0
      do i=1,n
      utmp(1,i)=flag2(i)
      enddo
      call sync_parallel_block_s
     & (n_size,n,utmp,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibufs
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      do i=1,n
      flag2(i)=utmp(1,i)
      enddo

      if(im.eq.0)then

      open(42,file='./data/dspli_10001_cart.dat',status='unknown')
      read(42,*) !nfp
      read(42,*) !nfp
      read(42,*) !nx,ny
      read(42,*) !nx,ny
      read(42,*) !dx,dy
      read(42,*) !dx,yd
      read(42,*) !origin(x),origin(y)
      read(42,*) !origin(1),origin(2)
      do iy=1,ny+1
      do ix=1,nx+1
      read(42,*) x,y,z,slip(1,ix,iy),slip(2,ix,iy),slip(3,ix,iy)
      enddo
      enddo
      close(42)
      endif
      call MPI_BCAST
     - (slip,3*(nx+1)*(ny+1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      maxx = origin(1) + dx*nx
      maxy = origin(2) + dy*ny

      dsx=nx*dx
      dsy=ny*dy

      up=0.0

      do i=1,n
      if((flag1(i).gt.0.5).and.(flag2(i).gt.0.5) .and.
     &   coor(1,i) >= origin(1) .and. coor(1,i) <= maxx .and.
     &   coor(2,i) >= origin(2) .and. coor(2,i) <= maxy) then

      ! found node on fault plane
      x=coor(1,i)-origin(1)
      y=coor(2,i)-origin(2)
!      itx=x/dsx
!      ity=y/dsy
!      x=x-itx*dsx
!      y=y-ity*dsy
      if( (x.lt.0.0).or.(x.gt.dsx))then
        write(*,*) 'something wrong in fault setting? x'
        stop
      endif
      if( (y.lt.0.0).or.(y.gt.dsy))then
        write(*,*) 'something wrong in fault setting? y'
        stop
      endif

      ix=x/dx+1.000001
      iy=y/dy+1.000001
      ratiox=(x-(ix-1)*dx)/dx
      ratioy=(y-(iy-1)*dy)/dy

      if(ix.eq.nx+1)then
        ix=nx
        ratiox=1.0
      endif
      if(iy.eq.ny+1)then
        iy=ny
        ratioy=1.0
      endif

      if((ratiox.lt.0.0).or.(ratiox.gt.1.0))then
        write(*,*) 'something wrong in fault setting? ratiox',ratiox,x
!        stop
        goto 528
      endif
      if((ratioy.lt.0.0).or.(ratioy.gt.1.0))then
        write(*,*) 'something wrong in fault setting? ratioy',ratioy,y
!        stop
        goto 528
      endif
      if((ix.gt.nx).or.(ix.lt.1))then
        write(*,*) 'something wrong in fault setting? ix',ix
        stop
      endif
      if((iy.gt.ny).or.(iy.lt.1))then
        write(*,*) 'something wrong in fault setting? iy',iy
        stop
      endif
      do j=1,3
      up(j,i)=
     - slip(j,ix  ,iy  )*(1.0-ratiox)*(1.0-ratioy)+
     - slip(j,ix+1,iy  )*(ratiox    )*(1.0-ratioy)+
     - slip(j,ix  ,iy+1)*(1.0-ratiox)*(ratioy    )+
     - slip(j,ix+1,iy+1)*(ratiox    )*(ratioy    )
      enddo
      endif
528   x=0
      enddo

      up=up*0.5d0
      rv=0.0
      do ie=1,ne
        in=num(ie)
        young=younglst(1,in)
        rnyu=younglst(2,in)
        gxy=young/(2.*(1+rnyu))
        gyz=young/(2.*(1+rnyu))
        gxz=young/(2.*(1+rnyu))
        do i1=1,4
        do j1=1,3
        xx(i1,j1)=coor(j1,cny(i1,ie))
        enddo
        enddo
        if(in.le.nfp)then
          do i1=1,10
          do j1=1,3
          arayi(3*(i1-1)+j1)=up(j1,cny(i1,ie))
          enddo
          enddo
        else
          do i1=1,10
          do j1=1,3
          arayi(3*(i1-1)+j1)=-up(j1,cny(i1,ie))
          enddo
          enddo
        endif
        arayo=0.
        call nontet10kud(young,gxy,gyz,gxz,rnyu,xx,arayi,arayo,detjmat)
        do i1=1,10
        do j1=1,3
        rv(npc*3*(cny(i1,ie)-1)+1+(j1-1)*npc)=
     &  rv(npc*3*(cny(i1,ie)-1)+1+(j1-1)*npc)+arayo(3*(i1-1)+j1)
        enddo
        enddo
      enddo
c
      end
c_______________________________________________________________________
      subroutine setduplixyz(n,ni,im,mpibuf,
     & nad,nadto,mpiad,mpiadlist,mpl,duplixyz)
      implicit none
      INCLUDE 'mpif.h'
      integer n,ni,im,nad,nadto
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      !real*8 mpibuf(2*3*npc*(nadto+1)*2)
      real*8 mpibuf(2*3*npc*(nadto+1))
      real*8 duplixyz(3*(n+ni)),up(3*(n+ni)*npc)
      integer i,j,ierr

      up=0.0
      do i=1,n+ni
      do j=1,3
      up(3*npc*(i-1)+j)=1
      enddo
      enddo
      call sync_parallel_block_d
     & (n+ni,n+ni,
     & up,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      do i=1,n+ni
      do j=1,3
      duplixyz(3*(i-1)+j)=1.0/up(3*npc*(i-1)+j)
      enddo
      enddo
      end
c_______________________________________________________________________
      subroutine sort_bubble(num,a)
      implicit none
      integer num,a(num)
      integer b(num)
      integer i,j,w
      b=a
      do i=1,num-1
      do j=num,i+1,-1
        if (b(j) >= b(j-1)) cycle
        w     =b(j-1)
        b(j-1)=b(j)
        b(j)  =w
      enddo
      enddo
      a=b
      end

c_______________________________________________________________________
      subroutine setcnysep(im,n,ne,n4,num,cny,cny4,
     & nsep,nsep4,numsep,
     & nsep_n_ptr,nsep_n_ptr4,nsep_ne_ptr,nsepnmap,nsepnmap4,
     & cnysep_read,cnysep_write,cnysep_read4,cnysep_write4)
      implicit none
      integer n,ne,n4,im,num(ne)
      integer cny(10,ne),cny4(4,ne)
      integer cnysep_read(10,ne),cnysep_write(10,ne)
      integer cnysep_read4(4,ne),cnysep_write4(4,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep_n_ptr4(nsta+2)
      integer nsep,nsepnmap(2*n)
      integer nsep4,nsepnmap4(2*n4)

      integer flag4(n4),flag(n),numtmp(ne)
      integer ikd,i,j,ie,nesep

c make separate temporary buffers  (move to main function?)
      nsep=0
      nsep4=0
      nesep=0
      do ie=1,ne
        if(num(ie).gt.nsta)then
          numtmp(ie)=nsta+1
        else
          numtmp(ie)=num(ie)
        endif
      enddo
c
      do ikd=1,nsta+1

      nsep_n_ptr(ikd)=nsep
      nsep_n_ptr4(ikd)=nsep4
      nsep_ne_ptr(ikd)=nesep
      flag=0
      flag4=0
      do ie=1,ne
        if(numtmp(ie).eq.ikd)then
          do j=1,10
            if(flag(cny(j,ie)).eq.0)then
              nsep=nsep+1
              flag(cny(j,ie))=nsep
              nsepnmap(nsep)=cny(j,ie)
            endif
          enddo
          do j=1,4
            if(flag4(cny4(j,ie)).eq.0)then
              nsep4=nsep4+1
              flag4(cny4(j,ie))=nsep4
              nsepnmap4(nsep4)=cny4(j,ie)
            endif
          enddo
        endif
      enddo
      do ie=1,ne
        if(numtmp(ie).eq.ikd)then
          nesep=nesep+1
          numsep(nesep)=num(ie)
          do j=1,10
            cnysep_read(j,nesep)=cny(j,ie)
            cnysep_write(j,nesep)=flag(cny(j,ie))
          enddo
           do j=1,4
            cnysep_read4(j,nesep)=cny4(j,ie)
            cnysep_write4(j,nesep)=flag4(cny4(j,ie))
          enddo
        endif
      enddo
      enddo
c
      nsep_n_ptr(nsta+2)=nsep
      nsep_n_ptr4(nsta+2)=nsep4
      nsep_ne_ptr(nsta+2)=nesep
c
      if(nesep.ne.ne)then
      write(*,*) 'something wrong with nesep',ne,nesep
      endif
      end
c_______________________________________________________________________

      subroutine remap_material_id(im,kd,ne,num)
        use mpi
        implicit none
        integer, intent(in) :: im, kd, ne
        integer, intent(inout) :: num(ne)
        integer :: material_mapping(kd)
        logical :: remap_mid
        integer :: ii,ie,ierr
        if(im == 0)then
          remap_mid = .false.
          open(42, file = 'data/material_mapping.dat',
     &             status = "old", err = 999)
          read(42,*)
          do ii = 1, kd
            read(42,*) material_mapping(ii)
          enddo
          remap_mid = .true.
          write(6, *) "remap material id"
          write(6, '(100i4)') material_mapping(:)
999       close(42)
        endif
        call MPI_BCAST(remap_mid,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

        if(.not. remap_mid)return
        call MPI_BCAST(material_mapping,kd,MPI_INTEGER,
     &                 0,MPI_COMM_WORLD,ierr)
        do ie = 1, ne
          num(ie) = material_mapping(num(ie))
        enddo
      end subroutine remap_material_id
c_______________________________________________________________________
      ! subroutines for debugging
      subroutine set_unit_uval(n,uval)
      implicit none
      integer, intent(in) :: n
      real(8), intent(out) :: uval(6,n)
      integer :: i
      real(8), parameter :: uvaldiag = 1d-14
      do i = 1,n
        uval(1,i) = uvaldiag
        uval(2,i) = uvaldiag
        uval(3,i) = uvaldiag
        uval(4,i) = 0d0
        uval(5,i) = 0d0
        uval(6,i) = 0d0
      enddo
      end subroutine set_unit_uval

      subroutine set_dummy_force(n,ni,coor,p)
      implicit none
      integer, intent(in) :: n,ni
      real(8), intent(in) :: coor(3,n)
      real(8), intent(out) :: p(npc,3,n+ni)
      integer :: i,j,ipc
      real(8) :: x,y,z
      integer :: npc_,icase
      icase = 2
      npc_ = min(npc,286)
      p = 0d0
      do i = 1, n
      do j = 1, 3
      do ipc = 1, npc_
        x = coor(1,i)
        y = coor(2,i)
        z = coor(3,i)
        select case(icase)
        case(1)
        p(ipc,j,i) = sin(x*0.11 + y*0.12 + z*0.13+j+ipc*3)
        case(2)
        p(1,j,i) = sin(x*0.11 + y*0.12 + z*0.13+j)
        end select
      enddo
      enddo
      enddo
      end subroutine set_dummy_force

      subroutine print_norm_s(npc_,n,ipc_sta,ipc_end,x,name)
      implicit none
      integer, intent(in) :: npc_, n, ipc_sta, ipc_end
      REAL_4, intent(in) :: x(npc_,3,n)
      character(len=*), intent(in) :: name
      real(8) x2
      integer :: ipc,j1,i
      x2 = 0d0
      do i = 1, n
      do j1 = 1, 3
      do ipc = ipc_sta, ipc_end
        x2 = x2 + x(ipc,j1,i)**2
      enddo
      enddo
      enddo
      write(6,*) name, x2
      end subroutine print_norm_s

      subroutine print_norm_d(npc_,n,ipc_sta,ipc_end,x,name)
      implicit none
      integer, intent(in) :: npc_, n, ipc_sta, ipc_end
      real(8), intent(in) :: x(npc_,3,n)
      character(len=*), intent(in) :: name
      real(8) x2
      integer :: ipc,j1,i
      x2 = 0d0
      do i = 1, n
      do j1 = 1, 3
      do ipc = ipc_sta, ipc_end
        x2 = x2 + x(ipc,j1,i)**2
      enddo
      enddo
      enddo
      write(6,*) name, x2
      end subroutine print_norm_d

      subroutine print_problem_size(im,n,ne,nadto,ni,n4,nadto4,ni4,
     &  niw,nir,niw4,nir4)
        implicit none
        integer, intent(in) :: im,n,ne,nadto,ni,n4,nadto4,ni4
        integer, intent(in) :: niw,nir,niw4,nir4

        if(im == 0)then
          write(6,*) "problem size"
          write(6,'(a7,3a12)') "","min", "max", "average"
        endif
        call print_problem_size_sub(im,"n",n)
        call print_problem_size_sub(im,"ne",ne)
        call print_problem_size_sub(im,"nadto",nadto)
        call print_problem_size_sub(im,"ni",ni)
        call print_problem_size_sub(im,"niw",niw)
        call print_problem_size_sub(im,"nir",nir)
        call print_problem_size_sub(im,"n4",n4)
        call print_problem_size_sub(im,"nadto4",nadto4)
        call print_problem_size_sub(im,"ni4",ni4)
        call print_problem_size_sub(im,"niw4",niw4)
        call print_problem_size_sub(im,"nir4",nir4)
      end subroutine print_problem_size

      subroutine print_problem_size_sub(im,name,m)
        use mpi
        implicit none
        integer, intent(in) :: im, m
        character(len=*), intent(in) :: name
        integer :: n_mpi, m_max, m_min, ierr
        integer(8) :: m_avg
        m_avg = m
        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ierr)
        call MPI_Allreduce(m,m_min,1,MPI_INTEGER,MPI_MIN,
     &                     MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(m,m_max,1,MPI_INTEGER,MPI_MAX,
     &                     MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,m_avg,1,MPI_INTEGER8,MPI_SUM,
     &                     MPI_COMM_WORLD,ierr)
        if(im == 0)then
          write(6,'(a7,3i)') name, m_min, m_max, nint(1d0*m_avg/n_mpi)
        endif
      end subroutine print_problem_size_sub

      subroutine reorder_nodes
     & (im,n,ne,tettype,nadto,mpiadlist,dupli,coor,cny,old2new)
      implicit none
      integer, intent(in) :: im,n,ne,tettype,nadto
      integer, intent(inout) :: mpiadlist(nadto),cny(tettype,ne)
      real(8), intent(inout) :: dupli(n), coor(3,n)
      integer, intent(out) :: old2new(n)
      integer, parameter :: flag_max = 2
      integer :: iflag,i,j,inew,ie
      integer :: node_flag(n),new2old(n)
      integer :: node_ptr(flag_max)
      real(8) :: dupli_tmp(n), coor_tmp(3,n)

      if(im == 0)then
        write(6,'(a,i0,a)') "reorder_nodes (tettype: ", tettype,")"
      endif
      node_flag(:) = 2
      do j = 1, nadto
        i = mpiadlist(j)
        node_flag(i) = 1
      enddo

      ! node_flag が小さい順に並べなおす
      node_ptr(:) = 0
      do i = 1, n
        iflag = node_flag(i)
        if(iflag == flag_max) cycle
        node_ptr(iflag+1) = node_ptr(iflag+1) + 1
      enddo
      do iflag = 2, flag_max
        node_ptr(iflag) = node_ptr(iflag) + node_ptr(iflag-1)
      enddo

      new2old(:) = -1
      do i = 1, n
        iflag = node_flag(i)
        node_ptr(iflag) = node_ptr(iflag) + 1
        inew = node_ptr(iflag)
        new2old(inew) = i
        old2new(i) = inew
      enddo

      coor_tmp = coor
      dupli_tmp = dupli
      do i = 1, n
        inew = old2new(i)
        coor(:,inew) = coor_tmp(:,i)
        dupli(inew) = dupli_tmp(i)
      enddo
      do j = 1, nadto
        i = mpiadlist(j)
        inew = old2new(i)
        mpiadlist(j) = inew
      enddo
      do ie = 1, ne
        do j = 1, tettype
          i = cny(j,ie)
          inew = old2new(i)
          cny(j,ie) = inew
        enddo
      enddo

      end subroutine reorder_nodes

      subroutine reorder_outflag(im,n,old2new,n2dvisl,u2dflagl)
        implicit none
        integer, intent(in) :: im,n,n2dvisl
        integer, intent(in) :: old2new(n)
        integer, intent(inout) :: u2dflagl(n2dvisl)
        integer :: i,j,inew
        do j = 1, n2dvisl
          i = u2dflagl(j)
          inew = old2new(i)
          u2dflagl(j) = inew
        enddo
      end subroutine reorder_outflag
