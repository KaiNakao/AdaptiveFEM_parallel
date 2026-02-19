
!---------------------------------------------------------------
      subroutine read_local_setting_hdf(im,kd,
     & n,ne,ib,nab,nad,nadto,fid,itype,dir)
      implicit none
      integer im,kd,i4_0,i4_1,i,itype
      integer n,ne,ib,nab,nad,nadto
      character filename*100
      character dataname*50
      ! input/output to hdf_lib libarary is in 8 byte integers
      integer*8 bufls,buflc,fid,lstart,lsize,ltmp
      integer bufi(1000)
      character(len=*) :: dir
!
      i4_0=0
      i4_1=1
      bufls=1000

      if(itype.eq.0)then
      write(filename,'(2a,i6.6,a8)')
     - trim(dir),'mdata/',im,'.data.h5'
      ltmp=len(trim(filename))
      call set_null_char(filename,ltmp)
      else
!      write(filename,'(a16,i6.6,a8)')
!     - './refine0/mdata/',im,'.data.h5'
!      ltmp=30
      write(filename,'(2a,i6.6,a8)')
     - trim(dir),'cdata/',im,'.data.h5'
      ltmp=len(trim(filename))
      call set_null_char(filename,ltmp)
      endif
      call hdf_open_file(filename,fid)
c
      write(dataname,'(a8)') '/setting'
      ltmp=8
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fid,dataname,bufi,bufls,buflc)
      n=bufi(1)
      ne=bufi(2)
      kd=bufi(3)
c
      write(dataname,'(a3)') '/BC'
      ltmp=3
      call set_null_char(dataname,ltmp)
      lsize=1
      lstart=0
      call hdf_read_int_array_part(fid,dataname,lstart,lsize,bufi)
      ib=bufi(1)
c
      write(dataname,'(a11)') '/ABCsetting'
      ltmp=11
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array_part(fid,dataname,lstart,lsize,bufi)
      nab=bufi(1)
c
      write(dataname,'(a8)') '/MPInode'
      ltmp=8
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array_part(fid,dataname,lstart,lsize,bufi)
      nad=bufi(1)
      if(nad.gt.bufls)then
      write(*,*)'enlarge bufi size',nad,bufls
      stop
      endif
      lsize=2*nad+1
      nadto=0
      call hdf_read_int_array_part(fid,dataname,lstart,lsize,bufi)
      do i=1,nad
      nadto=nadto+bufi(2*(i-1)+2+1)
      enddo
c
      end

!---------------------------------------------------------------
      subroutine read_local_data_hdf
     - (n_size,ne_size,n,ne,ib,nab,nad,nadto,
     - np,bufls,bufds,mpiad,mpiadlist,mpl,dupli,
     - coor,num,cny,ibi,nabi,abcelem,fidm,itype2,itype1)
c     - (n_size,n4_size,ne_size,n,ne,ib,nab,nad,nadto,
c     - n4,ne4,ib4,nab4,nad4,nadto4,map,ncb4,
c     - np,bufls,bufds,
c     - map_ptr,map_ind,map_coes,
c     - mpiad,mpiadlist,mpl,dupli,
c     - mpiad4,mpiadlist4,mpl4,dupli4,
c     - coor,num,cny,ibi,
c     - coor4,num4,cny4,ibi4,
c     - nabc,abc,stress,norm,
c     - nabc4,abc4,
c     - npl,npl4,
c#ifdef USE_CRS4
c     - crb_ind4,crb_ptr4,
c     - liscrs4,lisuval4,liscrsabc4,
c#endif
c     - abcelem,
c     - fidm,fidc)
      implicit none
      integer*4 n_size,nc_size,ne_size,np,i,j,ie,itmp,ii,j1,j2
      integer*4 n,ne,ib,nab,nad,nadto,itype2,itype1
      character dataname*50
      ! input/output to hdf_lib libarary is in 8 byte integers
      integer*8 bufls,bufds,buflc,bufdc,fidm,ltmp
      integer*4 bufi(bufls)
      real*8 bufd(bufds)
      integer*4 nedomain(np)
c--
c      integer*4 map_ptr(n+1),map_ind(map)
c      REAL_4 map_coes(map)
c
      integer*4 mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real*8 dupli(n_size)
c
      real*8 coor(3,n_size)
      integer*4 num(ne_size),cny(itype2,ne_size),ibi(ib)
c
      integer*4 nabi(itype1,nab)
c      real*8 abc(3,6,6,nab),stress(nab,6,6),norm(nab,3)
c
c      integer*4 npl(np+1)
c
#ifdef USE_CRS4
c      integer*4 crb_ind4(ncb4/9),crb_ptr4(n4+1)
c
c      integer*4 liscrs4(4,4,ne4)
c      integer*4 lisuval4(n4)
c      integer*4 liscrsabc4(3,3,nab4)
#endif
c
      integer*4 abcelem(nab)
c--

!--- mapping
c      write(dataname,'(a8)') '/map_ptr'
c      ltmp=8
c      call set_null_char(dataname,ltmp)
c      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
c      do i=1,n+1
c      map_ptr(i)=bufi(i+1)
c      enddo

c      write(dataname,'(a8)') '/map_ind'
c      call set_null_char(dataname,ltmp)
c      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
c      do i=1,map
c      map_ind(i)=bufi(i+1)
c      enddo

c      write(dataname,'(a8)') '/map_coe'
c      call set_null_char(dataname,ltmp)
c      call hdf_read_double_array(fidc,dataname,bufd,bufds,bufdc)
c      do i=1,map
c      map_coes(i)=bufd(i)
c      enddo

!--- MPI mdata
      write(dataname,'(a8)') '/MPInode'
      ltmp=8
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
      do i=1,nad
      mpiad(i,1)=bufi(2*(i-1)+1+1)
      mpiad(i,2)=bufi(2*(i-1)+1+2)
      enddo
      do i=1,nadto
      mpiadlist(i)=bufi(2*nad+1+i)
      enddo
      mpl(1)=0
      do i=1,nad
      mpl(i+1)=mpl(i)+mpiad(i,2)
      enddo
c
      write(dataname,'(a10)') '/duplicate'
      ltmp=10
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
      do i=1,n
      dupli(i)=1./dble(bufi(i))
      enddo

!--- coor mdata
      write(dataname,'(a5)') '/coor'
      ltmp=5
      call set_null_char(dataname,ltmp)
      call hdf_read_double_array(fidm,dataname,bufd,bufds,bufdc)
      do i=1,n
      coor(1,i)=bufd(3*(i-1)+1)
      coor(2,i)=bufd(3*(i-1)+2)
      coor(3,i)=bufd(3*(i-1)+3)
      enddo
!      do i=2,n
!      coor(1,i)=coor(1,i)-coor(1,1)
!      coor(2,i)=coor(2,i)-coor(2,1)
!      enddo
!      coor(1,1)=0.0
!      coor(2,1)=0.0

!--- conn mdata
      write(dataname,'(a5)') '/conn'
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
      do i=1,ne
      do j=1,itype2
      cny(j,i)=bufi((itype2+1)*(i-1)+j)
      enddo
      num(i)=bufi((itype2+1)*(i-1)+(itype2+1))
      enddo

!--- BC mdata
      write(dataname,'(a3)') '/BC'
      ltmp=3
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
      do i=1,ib
      ibi(i)=bufi(i+1)
      enddo

!--- numofelements mdata
c      write(dataname,'(a14)') '/numofelements'
c      ltmp=14
c      call set_null_char(dataname,ltmp)
c      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
c      do ie=1,np
c      nedomain(ie)=bufi(2*(ie-1)+2)
c      enddo
c      npl(1)=0
c      do i=1,np
c      npl(i+1)=npl(i)+nedomain(i)
c      enddo

!--- abc mdata
      write(dataname,'(a8)') '/ABCnode'
      ltmp=8
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
      do i=1,nab
      do j=1,itype1
      nabi(j,i)=bufi(itype1*(i-1)+j)
      enddo
      enddo

#ifdef USECRS4
!--- crs cdata
      write(dataname,'(a8)') '/crs_ind'
      ltmp=8
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
      do i=1,ncb4/9
      crb_ind4(i)=bufi(i+1)
      enddo
c
      write(dataname,'(a8)') '/crs_ptr'
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
      do i=1,n4+1
      crb_ptr4(i)=bufi(i+1)
      enddo
c
      write(dataname,'(a12)') '/crs_connect'
      ltmp=12
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
      itmp=1
      do ie=1,ne4
      do i=1,4
      do j=1,4
      itmp=itmp+1
      liscrs4(i,j,ie)=bufi(itmp)
      enddo
      enddo
      enddo
c
      write(dataname,'(a16)') '/crs_ABC_connect'
      ltmp=16
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
      itmp=1
      do ie=1,nab4
      do i=1,3
      do j=1,3
      itmp=itmp+1
      liscrsabc4(j,i,ie)=bufi(itmp)
      enddo
      enddo
      enddo
#endif

!--- abcelem mdata
      write(dataname,'(a8)') '/abcelem'
      ltmp=8
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidm,dataname,bufi,bufls,buflc)
      do i=1,nab
      abcelem(i)=bufi(i)
      enddo

#ifdef USECRS4
!--- crsblock list cdata
      write(dataname,'(a15)') '/crs_block_list'
      ltmp=15
      call set_null_char(dataname,ltmp)
      call hdf_read_int_array(fidc,dataname,bufi,bufls,buflc)
      do i=1,n4
      lisuval4(i)=bufi(i)
      enddo
#endif
      end
!---------------------------------------------------------------
      subroutine read_outputflag_settings
     - (outputopt,im,n3dvisl,n3dvisa,dir)
      implicit none
      INCLUDE 'mpif.h'
      integer*4 ierr,im
      character filename*100
      integer*4 outputopt
      integer*4 n3dvisl,n3dvisa
      integer*4 i4_1,i4_0
      character(len=*) :: dir
      i4_1=1
      i4_0=0
      if(outputopt.eq.0)then
        write(filename,'(2a,i6.6,a17)')
     - trim(dir),'mdata/',im,'.2Doutputflag.dat'
!        write(filename,'(a8,i6.6,a17)')
!     - './mdata/',im,'.2Doutputflag.dat'
      else
        write(filename,'(2a,i6.6,a17)')
     - trim(dir),'mdata/',im,'.3Doutputflag.dat'
!        write(filename,'(a8,i6.6,a17)')
!     - './mdata/',im,'.3Doutputflag.dat'
      endif
      open(10,file=filename,status='old')
      read(10,*) n3dvisl
      close(10)
      if(im.eq.0)then
      if(outputopt.eq.0)then
        open(10,file='./data/2Doutput.dat',status='old')
      else
        open(10,file='./data/3Doutput.dat',status='old')
      endif
      read(10,*) n3dvisa
      close(10)
      endif
      call MPI_BCAST(n3dvisa,i4_1,MPI_INTEGER,
     - i4_0,MPI_COMM_WORLD,ierr)
      end
c_______________________________________________________________________
      subroutine read_outputflag
     - (outputopt,im,n3dvisl,u3dflagl,dir)
      implicit none
      integer*4 outputopt,im,n3dvisl,u3dflagl(n3dvisl),i
      character filename*100
      character(len=*) :: dir

      if(outputopt.eq.0)then
      write(filename,'(2a,i6.6,a17)')
     - trim(dir),'mdata/',im,'.2Doutputflag.dat'
!        write(filename,'(a8,i6.6,a17)')
!     - './mdata/',im,'.2Doutputflag.dat'
      open(10,file=filename,status='old')
      read(10,*)
      do i=1,n3dvisl
      read(10,*) u3dflagl(i)
      enddo
      close(10)
      else
      write(filename,'(2a,i6.6,a17)')
     - trim(dir),'mdata/',im,'.3Doutputflag.dat'
      open(10,file=filename,status='old')
      read(10,*)
      do i=1,n3dvisl
      read(10,*) u3dflagl(i)
      enddo
      close(10)
      endif
      end
c_______________________________________________________________________
      subroutine read_geometry(n,ne,coor,cny,num)
      real*8 coor(3,n)
      integer cny(10,ne),num(ne)

      open(42,file='./data/geometry.flavia.msh',status='old')
      read(42,*)
      read(42,*)
      read(42,*)
      do i=1,n
      read(42,'(i10,3f16.6)')
     -         itmp,coor(1,itmp),coor(2,itmp),coor(3,itmp)
      enddo
      read(42,*)
      read(42,*)
      read(42,*)
      read(42,*)
      do i=1,ne
      read(42,200) itmp,cny(1,i),cny(2,i),cny(3,i),cny(4,i),
     -                  cny(5,i),cny(6,i),cny(7,i),cny(8,i),
     -                  cny(9,i),cny(10,i),
     -                  num(i)
      enddo
      close(42)
  200 format(12i10)

      end
c_______________________________________________________________________
      subroutine read_geometry4(n,ne,coor,cny,num)
      real*8 coor(3,n)
      integer cny(4,ne),num(ne)

      open(42,file='./data/tet4model3d.flavia.msh',status='old')
      read(42,*)
      read(42,*)
      read(42,*)
      do i=1,n
      read(42,'(i10,3f16.6)')
     -         itmp,coor(1,itmp),coor(2,itmp),coor(3,itmp)
      enddo
      read(42,*)
      read(42,*)
      read(42,*)
      read(42,*)
      do i=1,ne
      read(42,200) itmp,cny(1,i),cny(2,i),cny(3,i),cny(4,i),
     -                  num(i)
      enddo
      close(42)
  200 format(12i10)

      end
c_______________________________________________________________________
      subroutine read_material(kd,rmat)
      real*8 rmat(10,kd)
      open(42,file='./data/material.dat',status='old')
      do j=1,kd
      read(42,*)
      do i=1,10
      read(42,*) rmat(i,j)
      enddo
      enddo
      close(42)
      end

c_______________________________________________________________________
      subroutine read_sig(im,sig)
      use mpi
      implicit none
      integer, intent(in) :: im
      real(8), intent(out) :: sig(nkl)
      integer :: i,itmp,ierr
      integer, parameter :: nkl0 = 0, nkl5 = 5, nkl10 = 10
      if (im == 0) then
        if(nkl == nkl0) then
          continue
        elseif(nkl == nkl5)then
          open(41,file='./data/standarddeviation.dat',status='old')
          do i=1,nkl5
          read(41,*) sig(i)
          enddo
          close(41)
         elseif(nkl == nkl10)then
          open(10, file = 'data/uncertaintysetting.dat', status = 'old')
          read(10, *)
          read(10, *) itmp
          if (itmp /= nkl) then
            write(6, *) "Something is wrong with",
     &                  "data/uncertaintysetting.dat"
            write(6, *) "nkl (arg):", nkl, "nkl (file):", itmp
            stop
          endif
          read(10,*)
          do i = 1, nkl10
            read(10,*) itmp, itmp, sig(i)
          enddo
        else
          write(6,*) "Invalid nkl", nkl
          stop
        endif
        write(*,*) 'sig',sig
      endif
      call MPI_BCAST(sig,nkl,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      end subroutine read_sig

c_______________________________________________________________________
      subroutine read_uvalcoe(im,uvalcoe)
      use mpi
      implicit none
      integer, intent(in) :: im
      real(8), intent(out) :: uvalcoe(npc)
      integer :: i,ierr
      integer, parameter :: npc1=1, npc21=21, npc286=286, npc288=288
      uvalcoe(:) = 1d0
      if (im == 0) then
        if(npc==npc1)then
          continue
        elseif(npc==npc21)then
          uvalcoe(7)=0.5
          uvalcoe(9)=0.5
          uvalcoe(12)=0.5
          uvalcoe(16)=0.5
          uvalcoe(21)=0.5
        elseif(npc==npc288)then
          open(41,file='./data/uvalcoe.dat',status='old')
          read(41,*)
          read(41,*)
          read(41,*)
          do i=1,npc286
          read(41,*) uvalcoe(i)
          enddo
          close(41)
        else
          write(6,*) "Invalid npc", npc
          stop
        endif
        ! write(6,*) "uvalcoe", uvalcoe
      endif
      call MPI_BCAST(uvalcoe,npc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      end subroutine read_uvalcoe

c_______________________________________________________________________
      subroutine read_log_setting(im,log_innerCG,log_innerCG_per_outite)
        use mpi
        implicit none
        integer, intent(in) :: im
        logical, intent(out) :: log_innerCG, log_innerCG_per_outite
        integer :: itmp, ierr
        if(im == 0)then
          open(41, file = 'data/log_setting.dat', status = "old")
          read(41, *)
          read(41, *) itmp
          log_innerCG = itmp == 1
          read(41, *)
          read(41, *) itmp
          log_innerCG_per_outite = itmp == 1
          close(41)
          write(6, *) "output log for inner CG:", log_innerCG
          write(6, *) "output summary for inner CG per outer iteration:"
     &                 ,log_innerCG_per_outite
        endif
        call MPI_Bcast(log_innerCG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(log_innerCG_per_outite,1,
     &                 MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      end subroutine read_log_setting
