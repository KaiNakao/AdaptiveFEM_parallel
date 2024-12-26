# define NBLOCKS 8

#ifdef TESTSETFAULTPOINT

!#define DIR_MDATA_CDATA "/tmp/u00138/"
#define DIR_MDATA_CDATA "./"
c_______________________________________________________________________
      implicit none
      INCLUDE 'mpif.h'
      integer im,ncpu,np,kd,nfa,n2dvisl,n2dvisa
      integer n_size,ne_size,n,ne,ib,nab,nad,nadto
      integer n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,imtmp

      ! input/output to hdf_lib libarary is in 8 byte integers
      integer*8 bufls,bufds,fidc,fidm
      integer ierr

      fidm=0
      fidc=1

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, im, ierr)

      if(im.eq.0)then
      write(6, *) "read mdata, cdata from ", DIR_MDATA_CDATA
      np=1
      open(10,file='./data/setting.dat',status='old')
      read(10,*) !'num of node'
      read(10,*) ! n
      read(10,*) !'num of elem'
      read(10,*) ! ne
      read(10,*) !'num of elem'
      read(10,*)  kd
      close(10)
!      call check_macro_parameters
      endif
      call MPI_BCAST(np,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      do imtmp=0,9215 
!      do imtmp=0,9000
      write(*,*) 'processing rank',imtmp
      
      call read_local_setting_hdf(imtmp,kd,
     & n,ne,ib,nab,nad,nadto,fidm,0,DIR_MDATA_CDATA)
      call compute_nsize(n,n_size)
      call compute_nsize(ne,ne_size)
      call read_local_setting_hdf(imtmp,kd,
     & n4,ne4,ib4,nab4,nad4,nadto4,fidc,1,DIR_MDATA_CDATA)
      call compute_nsize(n4,n4_size)
      call compute_nsize(ne4,ne4_size)
      call read_outputflag_settings
     & (0,imtmp,n2dvisl,n2dvisa,DIR_MDATA_CDATA)
      bufds=n*30
      bufls=n*30
c
      call test_compute_ni(imtmp,ncpu,np,kd,nfa,n2dvisl,
     & n_size,ne_size,n,ne,ib,nab,nad,nadto,fidm,
     & n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,fidc,
     & bufds,bufls)
     
      call hdf_close_file(fidc)
      call hdf_close_file(fidm)
     
      enddo
      call MPI_FINALIZE(ierr)
      end
c_______________________________________________________________________
      subroutine test_compute_ni(im,ncpu,np,kd,nfa,n2dvisl,
     & n_size,ne_size,n,ne,ib,nab,nad,nadto,fidm,
     & n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,fidc,
     & bufds,bufls)
      implicit none
      integer im,ncpu,np,kd,nfa,n2dvisl,n2dvisa
      integer n_size,ne_size,n,ne,ib,nab,nad,nadto
      integer n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4
      integer*8 bufls,bufds,fidc,fidm

      
      real*8 coor(3,n_size),ds
      integer cny(10,ne_size),num(ne_size)
      real*8 dupli(n_size)
      integer ibi(ib),nabi(6,nab),npl(np+1),abcelem(nab)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real*8 rv(3*n*npc)
      
      
      open(61,file='./data/modeling_setting.dat',status='unknown')
      read(61,*) !'nx, ny'
      read(61,*) !'nx, ny'
      read(61,*) !'ds'
      read(61,*) ds
      close(61)
      
      call read_local_data_hdf
     - (n_size,ne_size,n,ne,ib,nab,nad,nadto,
     - np,bufls,bufds,mpiad,mpiadlist,mpl,dupli,
     - coor,num,cny,ibi,nabi,abcelem,fidm,10,6)
     
      call set_fault_point
     -(n,n_size,ne,ne_size,0,0,im,rv,coor,cny,1,ds)
      end
#endif
c_______________________________________________________________________
      subroutine set_fault_point
     -(n,n_size,ne,ne_size,ni,ni_size,im,rv,coor,cny,iunit,ds)
      implicit none
      INCLUDE 'mpif.h'
      integer n,n_size,ne,ne_size,ni,im,ierr,ni_size,nfa
      real*8 rv(3*(n+ni)*npc),coor(3,n),ds
      integer cny(10,ne)
      integer ki,i,ii,iunit
      character*50 filename

#ifdef TESTSETFAULTPOINT
      if(.true.)then
#else
      if(im.eq.0)then
#endif
#ifdef UNITFAULT
      write(filename,'(a32,i4.4,a4)')
     - './data/faultpara/faultpara_gfmod',iunit,'.dat'
      open(41,file='./data/faultpara.dat',status='old')
#else
      write(filename,*)'./data/faultpara.dat'
      open(41,file='./data/faultpara.dat',status='old')
#endif
      read(41,*)
      read(41,*) nfa
      close(41)
      endif
      call MPI_BCAST(nfa,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call set_fault_point_main
     -(n,n_size,ne,ne_size,ni,ni_size,im,nfa,ierr,
     - rv,coor,cny,iunit,ds)

      end
c_______________________________________________________________________
      subroutine set_fault_point_main
     -(n,n_size,ne,ne_size,ni,ni_size,im,nfa,ierr,
     - rv,coor,cny,iunit,ds)
      implicit none
      INCLUDE 'mpif.h'
      integer n,n_size,ne,ne_size,ni,im,nfa,ierr,ni_size
      real*8 rv(3*(n+ni)*npc),coor(3,n),ds
      integer cny(10,ne)
      real*8 xx(4,3),strike,dip,rake,ftmp(30)
      real*8 corf(3,nfa),faultp(nfa,4),rfb(3,nfa),rf(3)
      integer ietb(nfa),ietbtmp(nfa),ietbim(nfa),iet
      integer ki,i,ii,iunit
      real*8 t1,t2,tmp
      character*50 filename

      rv=0.

      tmp=(8*NBLOCKS**3+8*10+4*3)*nfa/1000.0/1000.0
#ifdef TESTSETFAULTPOINT
      if(.true.)then
#else
      if( im.le.1 )then
!      if(im.eq.1)then
#endif
      write(*,*) 'nfa',nfa
      write(*,*) 'nfa_NBLOCKS',NBLOCKS
      write(*,*) 'memory use(MiB)',tmp
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t1=MPI_WTIME()

#ifdef TESTSETFAULTPOINT
      if(.true.)then
#else
      if(im.eq.0)then
#endif
#ifdef UNITFAULT
      write(filename,'(a32,i4.4,a4)')
     - './data/faultpara/faultpara_gfmod',iunit,'.dat'
      open(41,file='./data/faultpara.dat',status='old')
#else
      write(filename,*)'./data/faultpara.dat'
      open(41,file='./data/faultpara.dat',status='old')
#endif
      read(41,*)
      read(41,*)
      do i=1,nfa
      read(41,*)  ! i,'th point source'
      read(41,*)  corf(1,i)
      read(41,*)  corf(2,i)
      read(41,*)  corf(3,i)
      read(41,*)  faultp(i,1)
      read(41,*)  faultp(i,2)
      read(41,*)  faultp(i,3)
      read(41,*)  faultp(i,4)
      
      corf(1,i)=corf(1,i)+0.0123*ds
      corf(2,i)=corf(2,i)+0.0123*ds
      corf(3,i)=corf(3,i)+0.0123*ds

      print *, "fault coordinate: ", corf(1,i), corf(2,i), corf(3,i)
      
      enddo
      close(41)
      endif
      call MPI_BCAST(corf,nfa*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(faultp,nfa*4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'read and cast fault para took',t2-t1
      t1=t2
      endif
      
      call pickupelem_block(n,ne,nfa,coor,cny,corf,ietb,rfb,im)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'pick up ',t2-t1
      t1=t2
      endif
      
!~      ietbtmp=-1
!~      do ki=1,nfa
!~      if(ietb(ki).ne.0)then
!~      ietbtmp(ki)=im
!~      endif
!~      enddo
!~      call MPI_ALLREDUCE(ietbtmp,ietbim,nfa,MPI_INTEGER,
!~     - MPI_MAX,MPI_COMM_WORLD,ierr)
!~      do ki=1,nfa
!~      if(ietbim(ki).ne.im)then
!~      if(ietb(ki).ne.0)then
!~        write(*,*) "check here in ietb",im,ki,ietb(ki)
!~      endif
!~      ietb(ki)=0
!~      endif
!~      enddo
!~
!~      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!~      t2=MPI_WTIME()
!~      if(im.eq.1)then
!~      write(*,*) 'delete duplicate',t2-t1
!~      t1=t2
!~      endif

      do ki=1,nfa
      if(ietb(ki).ne.0)then
      rf(1)=rfb(1,ki)
      rf(2)=rfb(2,ki)
      rf(3)=rfb(3,ki)
      strike=faultp(ki,1)
      dip=faultp(ki,2)
      rake=faultp(ki,3)
      iet=ietb(ki)
      do i=1,4
      do ii=1,3
      xx(i,ii)=coor(ii,cny(i,iet))
      enddo
      enddo
      
      call check_inside(xx,corf(1,ki),corf(2,ki),corf(3,ki),im,ki,iet)

      call pointf(xx,rf,strike,dip,rake,ftmp)
      ftmp=ftmp*faultp(ki,4)
      do i=1,10
      do ii=1,3
        rv(3*npc*(cny(i,iet)-1)+1+(ii-1)*npc)=
     &  rv(3*npc*(cny(i,iet)-1)+1+(ii-1)*npc)+ftmp(3*(i-1)+ii)
      enddo
      enddo
      endif
      enddo
!      call sync_parallel_block_d
!     &  (n_size+ni_size,n+ni,rv,
!     & nad,nadtoi,mpiadi,mpiadlisti,mpli,ierr,im,mpibuf
!#ifdef MEASURE_TIME_BARRIER
!     - ,mp,timerbufd,timerbufi,ind
!#endif
!     - )

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'sdp into rv',t2-t1
      t1=t2
      endif

      do i=1,nfa
      if(ietb(i).ne.0) ietb(i)=1
      enddo
      call MPI_ALLREDUCE(ietb,ietbtmp,nfa,MPI_INTEGER,
     - MPI_SUM,MPI_COMM_WORLD,ierr)

      if(im.eq.1)then
      write(*,*) 'nfa',nfa,sum(ietbtmp)
      endif

      if(im.eq.1)then
      do i=1,nfa
      if(ietbtmp(i).gt.1) write(*,*)"duplicate fault point",i
      if(ietbtmp(i).eq.0) write(*,*)"not found fault point",i
      enddo      
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'check',t2-t1
      t1=t2
      endif

      end
c_______________________________________________________________________
      subroutine check_inside(xx,x,y,z,im,ki,iet)
      implicit none
      real*8 x,y,z,xx(4,3)
      real*8 cort(3,1),coor(3,4),vec(3,3)
      real*8 volt,vol,vol1,vol2,vol3,vol4
      integer im,ki,iet
      integer ii,j,i1,i2,i3,i4
      
!      write(*,*) 'x,y,z',x,y,z
      
      do ii=1,4
      coor(1,ii)=xx(ii,1)
      coor(2,ii)=xx(ii,2)
      coor(3,ii)=xx(ii,3)
      enddo
      cort(1,1)=x
      cort(2,1)=y
      cort(3,1)=z
      
      
!      write(*,*) 'x,y,z',x,y,z

      i1=1
      i2=2
      i3=3
      i4=4
      j=1

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
      vec(3,ii)=cort(ii,j )-coor(ii,i1)
      enddo
      call tet4volume(vec,vol1)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-coor(ii,i1)
      vec(2,ii)=cort(ii,j )-coor(ii,i1)
      vec(3,ii)=coor(ii,i4)-coor(ii,i1)
      enddo
      call tet4volume(vec,vol2)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-cort(ii,j)
      vec(2,ii)=coor(ii,i3)-cort(ii,j)
      vec(3,ii)=coor(ii,i4)-cort(ii,j)
      enddo
      call tet4volume(vec,vol3)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i3)-coor(ii,i1)
      vec(2,ii)=coor(ii,i4)-coor(ii,i1)
      vec(3,ii)=cort(ii,j )-coor(ii,i1)
      enddo
      call tet4volume(vec,vol4)
c
      volt=vol1+vol2+vol3+vol4
c
!      if(volt/vol.le.1.0000005)then
      if(volt/vol.ge.1.0000005)then
      write(*,*) "error_step2,id",im,ki,iet
      write(*,*) "error_step2,val",volt,vol,volt/vol
      write(*,*) "error_step2,x",x,xx(:,1)
      write(*,*) "error_step2,y",y,xx(:,2)
      write(*,*) "error_step2,z",z,xx(:,3)
      else
      write(*,*) 'found_step2,id',im,ki,iet
      write(*,*) 'found_step2,x',x,xx(:,1)
      write(*,*) 'found_step2,y',y,xx(:,2)
      write(*,*) 'found_step2,z',z,xx(:,3)
      endif

      end
c_______________________________________________________________________
      subroutine pickupelem_block(n,ne,nfa,coor,cny,cort,iet,r,im)
      implicit none
      integer n,ne,iet(nfa),nfa,im
      real*8 coor(3,n)
      integer cny(10,ne)
      real*8 cort(3,nfa),vec(3,3),vol,r(3,nfa)

      real*8 vol1,vol2,vol3,vol4,volt
      integer ie,ii,i1,i2,i3,i4
      real*8 damy
      real*8 eps
      real*8 coormin(3),coormax(3),dblocks(3)
      real*8 coormintmp(3),coormaxtmp(3)
      integer nb_count(NBLOCKS,NBLOCKS,NBLOCKS)
      integer nb(NBLOCKS,NBLOCKS,NBLOCKS,nfa)
      integer ix,iy,iz,itmp,ix1,ix2,iy1,iy2,iz1,iz2,j1,j,i
c 
      do i=1,3
      coormin(i)=coor(i,1)
      coormax(i)=coor(i,1)
      enddo

      do j=1,n
      do i=1,3
      coormin(i)=min(coor(i,j),coormin(i))
      coormax(i)=max(coor(i,j),coormax(i))
      enddo
      enddo

#ifdef TESTSETFAULTPOINT
      if(.true.)then
#else
      if(im.eq.1)then
#endif
      write(*,*) 'coormin',coormin
      write(*,*) 'coormax',coormax
      endif


      eps=1d0

      do i=1,3
      coormax(i)=coormax(i)+eps
      coormin(i)=coormin(i)-eps
      dblocks(i)=(coormax(i)-coormin(i))/NBLOCKS
      enddo

#ifdef TESTSETFAULTPOINT
      if(.true.)then
#else
      if(im.eq.1)then
#endif
      write(*,*) 'dblocks',dblocks
      endif

      nb_count=0
      do i=1,nfa
      ix=int((cort(1,i)-coormin(1))/dblocks(1))+1
      iy=int((cort(2,i)-coormin(2))/dblocks(2))+1
      iz=int((cort(3,i)-coormin(3))/dblocks(3))+1
!      if((im.eq.30001))then
!      write(*,*) 'cort',cort(:,i)
!      write(*,*) 'ix,iy,iz',ix,iy,iz,i
!      write(*,*) 'x',(cort(1,i)-coormin(1))/dblocks(1)
!      write(*,*) 'y',(cort(2,i)-coormin(2))/dblocks(2)
!      write(*,*) 'z',(cort(3,i)-coormin(3))/dblocks(3)
!      endif
      if((ix.ge.1)      .and.(iy.ge.1)      .and.(iz.ge.1).and.
     -   (ix.le.NBLOCKS).and.(iy.le.NBLOCKS).and.(iz.le.NBLOCKS))then
      itmp=nb_count(ix,iy,iz)+1
      nb(ix,iy,iz,itmp)=i
      nb_count(ix,iy,iz)=itmp
      endif
      enddo

#ifdef TESTSETFAULTPOINT
      if(.true.)then
#else
      if(im.eq.1)then
#endif
      write(*,*) 'nb_count'
      do i=1,NBLOCKS
      write(*,*) nb_count(i,:,:)
      enddo
      endif

      iet=0

      do ie=1,ne
      i1=cny(1,ie)
      i2=cny(2,ie)
      i3=cny(3,ie)
      i4=cny(4,ie)
      do ii=1,3
      coormintmp(ii)=
     -   min(coor(ii,i1),coor(ii,i2),coor(ii,i3),coor(ii,i4))
      coormaxtmp(ii)=
     -   max(coor(ii,i1),coor(ii,i2),coor(ii,i3),coor(ii,i4))
      enddo
!      ix1=max(int((coormintmp(1)-coormin(1)-eps)/dblocks(1))+1,1)
!      iy1=max(int((coormintmp(2)-coormin(2)-eps)/dblocks(2))+1,1)
!      iz1=max(int((coormintmp(3)-coormin(3)-eps)/dblocks(3))+1,1)
!      ix2=min(int((coormaxtmp(1)-coormin(1)+eps)/dblocks(1))+1,NBLOCKS)
!      iy2=min(int((coormaxtmp(2)-coormin(2)+eps)/dblocks(2))+1,NBLOCKS)
!      iz2=min(int((coormaxtmp(3)-coormin(3)+eps)/dblocks(3))+1,NBLOCKS)
      ix1=max(int((coormintmp(1)-coormin(1))/dblocks(1))+0,1)
      iy1=max(int((coormintmp(2)-coormin(2))/dblocks(2))+0,1)
      iz1=max(int((coormintmp(3)-coormin(3))/dblocks(3))+0,1)
      ix2=min(int((coormaxtmp(1)-coormin(1))/dblocks(1))+2,NBLOCKS)
      iy2=min(int((coormaxtmp(2)-coormin(2))/dblocks(2))+2,NBLOCKS)
      iz2=min(int((coormaxtmp(3)-coormin(3))/dblocks(3))+2,NBLOCKS)

      do ix=ix1,ix2
      do iy=iy1,iy2
      do iz=iz1,iz2
      do j1=1,nb_count(ix,iy,iz)
      j=nb(ix,iy,iz,j1)
!      if(iet(j).eq.0)then

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
      vec(3,ii)=cort(ii,j )-coor(ii,i1)
      enddo
      call tet4volume(vec,vol1)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-coor(ii,i1)
      vec(2,ii)=cort(ii,j )-coor(ii,i1)
      vec(3,ii)=coor(ii,i4)-coor(ii,i1)
      enddo
      call tet4volume(vec,vol2)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i2)-cort(ii,j)
      vec(2,ii)=coor(ii,i3)-cort(ii,j)
      vec(3,ii)=coor(ii,i4)-cort(ii,j)
      enddo
      call tet4volume(vec,vol3)
c
      do ii=1,3
      vec(1,ii)=coor(ii,i3)-coor(ii,i1)
      vec(2,ii)=coor(ii,i4)-coor(ii,i1)
      vec(3,ii)=cort(ii,j )-coor(ii,i1)
      enddo
      call tet4volume(vec,vol4)
c
      volt=vol1+vol2+vol3+vol4
c
      if(volt/vol.le.1.0000005)then
      
      if(iet(j).ne.0)then
      write(*,*) 'ERROR_FOUND_DUPLICATE_POINT_SOURCE',im,j,iet(j),ie
      endif
      
      iet(j)=ie
c centering
      r(1,j)=0.25!vol4/vol
      r(2,j)=0.25!vol2/vol
      r(3,j)=0.25!vol1/vol
      write(*,*) 'found_step1,id',im,j,ie
      write(*,*) 'found_step1,x',cort(1,j),
     & coor(1,i1),coor(1,i2),coor(1,i3),coor(1,i4)
      write(*,*) 'found_step1,y',cort(2,j),
     & coor(2,i1),coor(2,i2),coor(2,i3),coor(2,i4)
      write(*,*) 'found_step1,z',cort(3,j),
     & coor(3,i1),coor(3,i2),coor(3,i3),coor(3,i4)
      endif
      
!      endif ! if(iet(j).ne.0)then
      
      enddo      
      enddo      
      enddo      
      enddo  
      enddo  

      end

c_______________________________________________________________________
      subroutine output_surface
     -(n2dvisl,n,ni,im,u2dl,up,u2dflagl,coor,iunit)
      implicit none
      include 'mpif.h'
      integer n2dvisl,n,ni,im
      real u2dl(3*n2dvisl*npc)
      real*8 up(3*(n+ni)*npc),coor(3,n)
      integer u2dflagl(n2dvisl)

      real tmpcoor(3*n2dvisl)
      integer ipout,itmp,i,ii,j,n2dvisl2,iunit
      integer fh,mpicount,ierr
      integer ix,iy
      real*8 eps,dsg
      integer mpistat(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) mpifsbuf
      character*50 filename
      real*8 t1,t2
      

      t1=MPI_WTIME()

#if 1
      dsg=1500
      eps=1e0
      ii=0
      do i=1,n2dvisl
      ix=nint(coor(1,u2dflagl(i))/dsg)
      iy=nint(coor(2,u2dflagl(i))/dsg)
      if((abs(coor(1,u2dflagl(i))-ix*dsg).lt.eps).and.
     -   (abs(coor(2,u2dflagl(i))-iy*dsg).lt.eps))then
      ii=ii+1
      endif
      enddo
      n2dvisl2=ii
      ii=0
      do i=1,n2dvisl
      ix=nint(coor(1,u2dflagl(i))/dsg)
      iy=nint(coor(2,u2dflagl(i))/dsg)
      if((abs(coor(1,u2dflagl(i))-ix*dsg).lt.eps).and.
     -   (abs(coor(2,u2dflagl(i))-iy*dsg).lt.eps))then
      ii=ii+1
      tmpcoor(3*(ii-1)+1)=coor(1,u2dflagl(i))
      tmpcoor(3*(ii-1)+2)=coor(2,u2dflagl(i))
      tmpcoor(3*(ii-1)+3)=coor(3,u2dflagl(i))
      do j=1,npc
      u2dl(3*(ii-1)+1+3*(j-1)*n2dvisl2)=
     -up(3*npc*(u2dflagl(i)-1)+0*npc+j)
      u2dl(3*(ii-1)+2+3*(j-1)*n2dvisl2)=
     -up(3*npc*(u2dflagl(i)-1)+1*npc+j)
      u2dl(3*(ii-1)+3+3*(j-1)*n2dvisl2)=
     -up(3*npc*(u2dflagl(i)-1)+2*npc+j)
      enddo
      endif
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'global into outputbuf',t2-t1
      t1=t2
      endif

      itmp=0
#if REFSLIP>0
      itmp=REFSLIP+100
#endif
#ifdef UNITFAULT
      itmp=iunit
#endif
      itmp=0
#ifdef NOUTPUTBLOCK
      ipout=im/NOUTPUTBLOCK
#endif
      if(n2dvisl2>0)then
#ifdef NOUTPUTBLOCK
      write(filename,'(a21,i4.4,a1,i6.6,a6,i6.6,a4)')
     - './2Doutput/proc_block',ipout,'/',im,'.coor.',itmp,'.bin'
#else
      write(filename,'(a11,i6.6,a6,i6.6,a4)')
     - './2Doutput/',im,'.disp.',itmp,'.bin'
#endif
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,
     - MPI_MODE_CREATE+MPI_MODE_WRONLY,
     - MPI_INFO_NULL,fh,ierr)
      mpicount=n2dvisl2*3
      mpifsbuf=mpicount*4
      call MPI_FILE_SET_SIZE(fh,mpifsbuf,ierr)
      call MPI_FILE_WRITE(fh,tmpcoor,mpicount,MPI_REAL4,mpistat,ierr)
      call MPI_FILE_CLOSE(fh,ierr)

#ifdef NOUTPUTBLOCK
      write(filename,'(a21,i4.4,a1,i6.6,a6,i6.6,a4)')
     - './2Doutput/proc_block',ipout,'/',im,'.disp.',itmp,'.bin'
#else
      write(filename,'(a11,i6.6,a6,i6.6,a4)')
     - './2Doutput/',im,'.disp.',itmp,'.bin'
#endif
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,
     - MPI_MODE_CREATE+MPI_MODE_WRONLY,
     - MPI_INFO_NULL,fh,ierr)
      mpicount=n2dvisl2*3*npc
      mpifsbuf=mpicount*4
      call MPI_FILE_SET_SIZE(fh,mpifsbuf,ierr)
      call MPI_FILE_WRITE(fh,u2dl,mpicount,MPI_REAL4,mpistat,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      endif

#else
      do i=1,n2dvisl
      do j=1,npc
      u2dl(3*(i-1)+1+3*(j-1)*n2dvisl)=up(3*npc*(u2dflagl(i)-1)+0*npc+j)
      u2dl(3*(i-1)+2+3*(j-1)*n2dvisl)=up(3*npc*(u2dflagl(i)-1)+1*npc+j)
      u2dl(3*(i-1)+3+3*(j-1)*n2dvisl)=up(3*npc*(u2dflagl(i)-1)+2*npc+j)
      enddo
      enddo

      itmp=0
#if REFSLIP>0
      itmp=REFSLIP+100
#endif
#ifdef fault
      itmp=UNITSLIP
#endif
      itmp=0
#ifdef NOUTPUTBLOCK
      ipout=im/NOUTPUTBLOCK
#endif
      if(n2dvisl>0)then
#ifdef NOUTPUTBLOCK
      write(filename,'(a21,i4.4,a1,i6.6,a6,i6.6,a4)')
     - './2Doutput/proc_block',ipout,'/',im,'.disp.',itmp,'.bin'
#else
      write(filename,'(a11,i6.6,a6,i6.6,a4)')
     - './2Doutput/',im,'.disp.',itmp,'.bin'
#endif
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,
     - MPI_MODE_CREATE+MPI_MODE_WRONLY,
     - MPI_INFO_NULL,fh,ierr)
      mpicount=n2dvisl*3*npc
      mpifsbuf=mpicount*4
      call MPI_FILE_SET_SIZE(fh,mpifsbuf,ierr)
      call MPI_FILE_WRITE(fh,u2dl,mpicount,MPI_REAL4,mpistat,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      endif
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'binary write surface disp',t2-t1
      t1=t2
      endif

      end
c_______________________________________________________________________
      subroutine set_nankai_obs_sub
     - (im,n,n_size,ne,ne_size,ni,ni_size,coor,up,cny,
     -  nad,mpiadi,mpiadlisti,mpli,nadtoi,mpibuf,nadto,
     -  nvis,ivis,iunit)
      implicit none
      include 'mpif.h'
      integer im,n,no,n_size,ne,ne_size,ni,ni_size,ierr
      integer cny(10,ne)
      real*8 up(3*(n+ni)*npc),coor(3,n)
      integer nad,nadto
      integer mpiadi(nad,2),mpiadlisti(2*nadto),mpli(nad+1),nadtoi
      real*8 mpibuf(2*3*npc*(nadto+1)*2)
      integer nvis,ivis(nvis),iunit
      character*50 filename

      if(im.eq.0)then
      open(18,file='./data/rotatedvector.dat',status='unknown')
      no = 0
      do
      read (18, *, end=100) ! ファイル終端ならば999に飛ぶ
      no = no + 1
      end do
100   continue
      close(18)
      write(*,*) "no",no
      endif

      call MPI_BCAST(no,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call set_nankai_obs_main
     - (im,n,n_size,ne,ne_size,ni,ni_size,no,coor,up,cny,
     -  nad,mpiadi,mpiadlisti,mpli,nadtoi,mpibuf,nadto,
     -  nvis,ivis,iunit)

      end
c_______________________________________________________________________
      subroutine set_nankai_obs_main
     - (im,n,n_size,ne,ne_size,ni,ni_size,no,coor,up,cny,
     -  nad,mpiadi,mpiadlisti,mpli,nadtoi,mpibuf,nadto,
     -  nvis,ivis,iunit)
      implicit none
	  include 'mpif.h'
! in
      integer im,n,no,n_size,ne,ne_size,ni,ni_size,iunit
      integer cny(10,ne)
      real*8 up(3*(n+ni)*npc),coor(3,n)
      integer nad,nadto
      integer mpiadi(nad,2),mpiadlisti(2*nadto),mpli(nad+1),nadtoi
      real*8 mpibuf(2*3*npc*(nadto+1)*2)
      integer nvis,ivis(nvis)
!out
      integer nl(no) 
      integer ol(no)
      real wl(6,no)
      integer nnl
      integer ierr,isf,isfn(ne,7)
      real tmpu3dl3(3*no*npc),u3dl3(3*no*npc)
!tmp
      real*8 ds,tmp(3)
      real*8 coorf(2,no)
      integer no_count(no),no_count_all(no)
      integer,allocatable,dimension(:,:) :: xy_count
      integer,allocatable,dimension(:,:,:) :: xy_obs

      integer nx,ny
      integer xmin,xmax,ymin,ymax
      real*8 tmpcoor(2,3),cx,cy,s,t
      integer ix,iy,i,ii,i_vec,j,k,it,ic
      integer c_max,xy_c_tmp,inum

      character*50 filename
      real*8 t1,t2

      integer fh,mpicount
      integer mpistat(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) mpifsbuf

      t1=MPI_WTIME()

      if(im.eq.0)then

      open(21,file='./data/modeling_setting.dat',status='old')
      read(21,*)
      read(21,*)
      read(21,*)
      read(21,*) ds
      close(21)

      open(18,file='./data/rotatedvector.dat',status='unknown')
      do i=1,no
      read(18,*) coorf(1,i),coorf(2,i)
      enddo
      close(18)

      endif

      call MPI_BCAST(coorf(1,1),2*no,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ds,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'load observation data',t2-t1
      t1=t2
      endif

      if(im.eq.0)write(*,*) "no",no

      call estsurfnode
     -(n,n_size,ne,ne_size,ni,ni_size,nvis,cny,ivis,isf,isfn,im,
     - nad,mpiadi,mpiadlisti,mpli,nadtoi,mpibuf,nadto)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) 'est surface node',t2-t1
      t1=t2
      endif


      xmin=int(coorf(1,1)/ds)  
      xmax=int(coorf(1,1)/ds)
      ymin=int(coorf(2,1)/ds)  
      ymax=int(coorf(2,1)/ds)

      do i=1,no
      xmin=min(xmin,int(coorf(1,i)/ds))
      xmax=max(xmax,int(coorf(1,i)/ds))
      ymin=min(ymin,int(coorf(2,i)/ds))
      ymax=max(ymax,int(coorf(2,i)/ds))
      enddo
      
      nx=xmax-xmin+1
      ny=ymax-ymin+1

      if(im.eq.0)write(*,*) "ds",ds,nx,ny

      allocate(xy_count(nx,ny))

      xy_count(:,:)=0
      do i=1,no
      ix=int(coorf(1,i)/ds)-xmin+1
      iy=int(coorf(2,i)/ds)-ymin+1
      xy_count(ix,iy)=xy_count(ix,iy)+1
      enddo

      c_max=0
      do j=1,ny
      do i=1,nx
      c_max=max(c_max,xy_count(ix,iy))
      enddo
      enddo

      if(im.eq.0)write(*,*) "c_max",ds,c_max

      allocate(xy_obs(nx,ny,c_max))

      do j=1,ny
      do i=1,nx
      xy_count(i,j)=0
      do k=1,c_max
      xy_obs(i,j,k)=i
      enddo
      enddo
      enddo

      do i=1,no
      ix=int(coorf(1,i)/ds)-xmin+1  
      iy=int(coorf(2,i)/ds)-ymin+1
      xy_count(ix,iy)=xy_count(ix,iy)+1
      xy_obs(ix,iy,xy_count(ix,iy))=i
      enddo

      inum=0
      no_count(:)=0

      do i=1,isf
      do j=1,2
      tmpcoor(j,1)=coor(j,isfn(i,1))
      tmpcoor(j,2)=coor(j,isfn(i,2))
      tmpcoor(j,3)=coor(j,isfn(i,3))
      enddo
      ix=int(min(tmpcoor(1,1),tmpcoor(1,2),tmpcoor(1,3))/ds)-xmin+1
      iy=int(min(tmpcoor(2,1),tmpcoor(2,2),tmpcoor(2,3))/ds)-ymin+1
      if((ix.ge.1).and.(ix.le.nx).and.
     -   (iy.ge.1).and.(iy.le.ny))then
      xy_c_tmp=xy_count(ix,iy)
      if(xy_c_tmp.ne.0)then
      do it=1,xy_c_tmp
      ic=xy_obs(ix,iy,it)
      cx=coorf(1,ic)
      cy=coorf(2,ic)
            write(*,*) "debug 1716"
            write(*,*) coorf(1,ic),coorf(2,ic),ic,it
            write(*,*) cx,cy
            write(*,*) tmpcoor(1,1),tmpcoor(2,1)
            write(*,*) tmpcoor(1,2),tmpcoor(2,2)
            write(*,*) tmpcoor(1,3),tmpcoor(2,3)
            write(*,*) s,t,im
      call coest(cx,cy,tmpcoor,s,t)
      if((s.ge.0d0).and.(t.ge.0).and.(s+t.le.1))then
        if(no_count(ic).eq.0)then
          no_count(ic)=no_count(ic)+1
          inum=inum+1
          ol(inum)=ic
          wl(1,inum)=   2*s**2 + 2*t**2 + 4*s*t - 3*s - 3*t +1
          wl(2,inum)=   2*s**2                  - 1*s
          wl(3,inum)=            2*t**2               - 1*t 
          wl(4,inum)= - 4*s**2          - 4*s*t + 4*s
          wl(5,inum)=                     4*s*t
          wl(6,inum)=          - 4*t**2 - 4*s*t       + 4*t
          nl(inum)=i
        else
          if((s.eq.0d0).or.(t.eq.0).or.(s+t.eq.1))then
            write(*,*) "edge observation_point"
          else
            write(*,*) "debug 1716"
            write(*,*) coorf(1,ic),coorf(2,ic),ic,it
            write(*,*) cx,cy
            write(*,*) tmpcoor(1,1),tmpcoor(2,1)
            write(*,*) tmpcoor(1,2),tmpcoor(2,2)
            write(*,*) tmpcoor(1,3),tmpcoor(2,3)
            write(*,*) s,t,im
            do j=1,inum
              if(ol(j).eq.ic)then
              write(*,*) coor(1,isfn(nl(j),1)),coor(2,isfn(nl(j),1)),im
              write(*,*) coor(1,isfn(nl(j),2)),coor(2,isfn(nl(j),2)),im
              write(*,*) coor(1,isfn(nl(j),3)),coor(2,isfn(nl(j),3)),im
              endif
            enddo
          endif
        endif
      endif
      enddo
      endif
      endif
      enddo

      nnl=inum

      call MPI_ALLREDUCE(no_count,no_count_all,no,MPI_INTEGER,MPI_SUM
     -                   ,MPI_COMM_WORLD,ierr)
      if(im.eq.0)then
      do i=1,no
      if(no_count_all(i).ne.1)then
      write(*,*) "something wrong 1736",i,no_count_all(i)
!      stop
      endif
      enddo
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) "search_obs end",t2-t1
      t1=t2
      endif

      do i=1,3*no*npc
      tmpu3dl3(i)=0e0
      enddo

      if (nnl>0) then
      do i=1,nnl
      do i_vec=1,npc
      k=nl(i)
      tmp(:)=0
      do ii=1,6
      tmp(1)=tmp(1)+up(3*npc*(isfn(k,ii)-1)+0*npc+i_vec)*wl(ii,i)
      tmp(2)=tmp(2)+up(3*npc*(isfn(k,ii)-1)+1*npc+i_vec)*wl(ii,i)
      tmp(3)=tmp(3)+up(3*npc*(isfn(k,ii)-1)+2*npc+i_vec)*wl(ii,i)
      enddo
      tmpu3dl3(3*(i-1)+1+3*(i_vec-1)*nnl)=tmp(1)
      tmpu3dl3(3*(i-1)+2+3*(i_vec-1)*nnl)=tmp(2)
      tmpu3dl3(3*(i-1)+3+3*(i_vec-1)*nnl)=tmp(3)
      enddo
      enddo

#ifdef UNITFAULT
      if(iunit.eq.1)then
#endif
      write(filename,'(a20,i6.6,a6,i6.6,a4)')
     - './2Doutput/obs/obsid',im,'.disp.',0,'.bin'
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,
     - MPI_MODE_CREATE+MPI_MODE_WRONLY,
     - MPI_INFO_NULL,fh,ierr)
      mpicount=nnl
      mpifsbuf=mpicount*4
      call MPI_FILE_SET_SIZE(fh,mpifsbuf,ierr)
      call MPI_FILE_WRITE(fh,ol,mpicount,MPI_INTEGER,mpistat,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
#ifdef UNITFAULT
      endif
#endif

#ifdef UNITFAULT
      write(filename,'(a15,i6.6,a6,i6.6,a4)')
     - './2Doutput/obs/',im,'.disp.',iunit,'.bin'
#else
      write(filename,'(a15,i6.6,a6,i6.6,a4)')
     - './2Doutput/obs/',im,'.disp.',0,'.bin'
#endif
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,
     - MPI_MODE_CREATE+MPI_MODE_WRONLY,
     - MPI_INFO_NULL,fh,ierr)
      mpicount=nnl*3*npc
      mpifsbuf=mpicount*4
      call MPI_FILE_SET_SIZE(fh,mpifsbuf,ierr)
      call MPI_FILE_WRITE(fh,tmpu3dl3,mpicount,MPI_REAL4,mpistat,ierr)
      call MPI_FILE_CLOSE(fh,ierr)

      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.1)then
      write(*,*) "output surface",t2-t1
      t1=t2
      endif

      end

c_______________________________________________________________________
      subroutine coest(cx,cy,tmpcoor,s,t)
      implicit none
      real*8 tmpcoor(2,3)
      real*8 x21,x31,y21,y31
      real*8 d2121,d2131,d3131,detj,cx,cy,d21c,d31c
      real*8 s,t
	  cx=cx-tmpcoor(1,1)
      cy=cy-tmpcoor(2,1)
      x21=tmpcoor(1,2)-tmpcoor(1,1)
      x31=tmpcoor(1,3)-tmpcoor(1,1)
      y21=tmpcoor(2,2)-tmpcoor(2,1)
      y31=tmpcoor(2,3)-tmpcoor(2,1)
      d2121=x21*x21+y21*y21
      d2131=x21*x31+y21*y31
      d3131=x31*x31+y31*y31
      detj=d2121*d3131-d2131*d2131
      d21c=x21*cx+y21*cy
      d31c=x31*cx+y31*cy
      s=( d3131*d21c-d2131*d31c)/detj
      t=(-d2131*d21c+d2121*d31c)/detj

      end
c_______________________________________________________________________
      subroutine estsurfnode
     -(n,n_size,ne,ne_size,ni,ni_size,nvis,cny,ivis,isf,isfn,im,
     - nad,mpiadi,mpiadlisti,mpli,nadtoi,mpibuf,nadto)
      implicit none
      INCLUDE 'mpif.h'
      integer n,n_size,ni,ni_size,ne,ne_size
	  integer nad,nadto
      integer mpiadi(nad,2),mpiadlisti(2*nadto),mpli(nad+1),nadtoi
      real*8 mpibuf(2*3*npc*(nadto+1)*2)
      integer nvis,isf
      integer cny(10,ne),ivis(nvis),isfn(ne,7)


      integer i,ii,ie,ifl
      integer listmp(4,6),im
      integer flag(n+ni)
	  integer k,ierr
c  
! possible combination of node numbers in surface
      listmp(1,1)=1     
      listmp(1,2)=2    
      listmp(1,3)=4     
      listmp(1,4)=5     
      listmp(1,5)=9     
      listmp(1,6)=8     
      listmp(2,1)=2     
      listmp(2,2)=3     
      listmp(2,3)=4     
      listmp(2,4)=6     
      listmp(2,5)=10     
      listmp(2,6)=9     
      listmp(3,1)=3     
      listmp(3,2)=1     
      listmp(3,3)=4     
      listmp(3,4)=7     
      listmp(3,5)=8     
      listmp(3,6)=10     
      listmp(4,1)=2     
      listmp(4,2)=1     
      listmp(4,3)=3     
      listmp(4,4)=5     
      listmp(4,5)=7     
      listmp(4,6)=6     
c
	  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      flag=0
      do i=1,nvis
      flag(ivis(i))=1
      if(ivis(i).le.0) write(*,*) n,ivis(i)
      if(ivis(i).gt.n) write(*,*) n,ivis(i)
      enddo

      if (im.eq.0) write(*,*) "start sync"

      call sync_int_max
     &  (n_size+ni_size,n+ni,flag,
     & nad,nadtoi,mpiadi,mpiadlisti,mpli,ierr,im,mpibuf)
c
      if (im.eq.0) write(*,*) "end sync"


      isf=0
      do ie=1,ne
      do ii=1,4
      ifl=0
      do i=1,3
      ifl=ifl+flag(cny(listmp(ii,i),ie))
      enddo      
      if(ifl.eq.3)then
      isf=isf+1
	  do i=1,6
! isfn comtains 6 nodes which are located on surface and in the same element
      isfn(isf,i)=cny(listmp(ii,i),ie)
      enddo
      isfn(isf,7)=0!num(ie)
      endif 
      enddo
      enddo
	  
	  
      end
c_______________________________________________________________________
      subroutine sync_int_max
     &  (n_size,n,b,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf)
      implicit none
      INCLUDE 'mpif.h'
      integer n_size,n,nad,nadto
      integer iad,i,j,k,ipc
      integer b(n_size)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      ! temp
      integer mpibuf(2*(nadto+1))
      integer irecvreq(nad),isendreq(nad)
      integer itag,ierr,im,m1,m2,m3,istatus(MPI_STATUS_SIZE,nad)

      do i=1,nadto
      mpibuf(i)=b(mpiadlist(i))
      enddo

      do i=1,nad
        m1=mpl(i)+nadto+1
        m2=mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !m3*10000+im
        call MPI_IRECV(mpibuf(m1),m2,MPI_INTEGER,
     &    m3,itag,MPI_COMM_WORLD,irecvreq(i),ierr)
      enddo

      do i=1,nad
        m1=mpl(i)+1
        m2=mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !im*10000+m3
        call MPI_ISEND(mpibuf(m1),m2,MPI_INTEGER,
     &    m3,itag,MPI_COMM_WORLD,isendreq(i),ierr)
      enddo

      if(nad.ne.0)then
      call MPI_WAITALL(nad,irecvreq(1),istatus,ierr)
      call MPI_WAITALL(nad,isendreq(1),istatus,ierr)
      endif

      do iad=1,nad
        do i = mpl(iad)+1, mpl(iad+1)
          k=mpiadlist(i)
          b(k)=max(b(k),mpibuf(i+nadto))
        enddo
      enddo

      end

