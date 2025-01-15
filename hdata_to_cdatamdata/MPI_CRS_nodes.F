
!~      INCLUDE 'mpif.h'
!~      integer*4 myrank,ncpu,ierr
!~	integer im,mp
!~	character filename*50
!~
!~      call MPI_INIT(ierr)
!~      call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
!~      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
!~	
!~      open(10,file='./data/para_setting.dat',status='unknown')
!~      read(10,*) !'num of MPI thread'
!~      read(10,*) mp
!~c      read(10,*) !'num of OpenMP thread'
!~c      read(10,*) np
!~      close(10)
!~
!~      if(0.NE.MOD(mp,ncpu))then
!~      write(*,*) 'mp is not divisible by ncpu'
!~	call exit()
!~	endif
!~      do im=myrank*(mp/ncpu)+1,(myrank+1)*(mp/ncpu)
!~      write(*,*) myrank,im
!~!      do im=1,mp 
!~
!~      write(filename,'(a8,i6.6,a12)')
!~     - './mdata/',im-1,'.setting.dat'
!~      open(10,file=filename,status='unknown')
!~      read(10,*) !'num of node' 
!~      read(10,*)  n
!~      read(10,*) !'num of elem' 
!~      read(10,*)  ne
!~      read(10,*) !'num of materials' 
!~      read(10,*)  kd
!~      close(10)
!~c
!~      call main(n,ne,kd,im)
!~      
!~      enddo
!~c
!~      call MPI_FINALIZE(ierr)
!~      end
!~c_______________________________________________________________________
 !     subroutine main(n,ne,kd,im)
!#define NC_CRSNODESLIB 200
#define NC_CRSNODESLIB 800
!      parameter(nc=200)
      subroutine CRS_nodes_lib(n,ne,fid)
      implicit none
      integer*4 ne,n,ie,i,j,k,i1,j1,inum,ntm,nc
      integer*8 fid
      real*8 damy
      real*8 coor(n,3)
      integer*4 num(ne),cny(ne,10)
      integer*4 mat(n,NC_CRSNODESLIB),matnum(n)
      integer*4 bufsize

      nc=NC_CRSNODESLIB
c
      call read_tet10_geometry_hdf(n,ne,coor,cny,num,fid)
      matnum=0
      mat=0

      do ie=1,ne
      do i=1,10
      i1=cny(ie,i)   
      do j=1,10
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

!~      do i=1,n
!~      ntm=max(ntm,matnum(i))
!~      enddo
!~c
!~      write(filename,'(a8,i6.6,a14)')
!~     - './mdata/',im,'.crs_nodes.dat'
!~      open(42,file=filename,status='unknown')
!~c      open(42,file='./data/crs_nodes.dat',status='unknown')
!~      write(42,*) 'num of temp variable'
!~      write(42,*) ntm
!~      write(42,*) 'total num of components'
!~      write(42,*) inum
!~      write(42,*) 'num of components for each node'     
!~      do i=1,n
!~      write(42,*) matnum(i)
!~      enddo
!~      write(42,*) 'components num for each node'     
!~      do i=1,n
!~      do j=1,matnum(i)
!~      write(42,*) mat(i,j)
!~      enddo
!~      enddo
!~      close(42) 

      end
!~c_______________________________________________________________________
!~      subroutine read_geometry(n,ne,coor,cny,num,im)
!~      character filename*50
!~	real*8 coor(n,3)
!~	integer cny(ne,10),num(ne)
!~
!~      write(filename,'(a8,i6.6,a20)')
!~     - './mdata/',im-1,'.geometry.flavia.msh'
!~      open(42,file=filename,status='unknown')
!~      read(42,*) 
!~      read(42,*) 
!~      read(42,*) 
!~      do i=1,n
!~      read(42,'(i10,3f16.6)') 
!~     -         itmp,coor(itmp,1),coor(itmp,2),coor(itmp,3)
!~      enddo
!~      read(42,*) 
!~      read(42,*) 
!~      read(42,*) 
!~      read(42,*) 
!~      do i=1,ne
!~      read(42,200) itmp,cny(i,1),cny(i,2),cny(i,3),cny(i,4),
!~     -                  cny(i,5),cny(i,6),cny(i,7),cny(i,8),
!~     -                  cny(i,9),cny(i,10),
!~     -                  num(i)  
!~      enddo
!~      close(42)
!~  200 format(12i10)
!~
!~      end
