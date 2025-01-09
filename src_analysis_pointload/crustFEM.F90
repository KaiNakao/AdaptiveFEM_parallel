
#ifndef DIR_MDATA_CDATA
#define DIR_MDATA_CDATA "./"
#endif
      program mainprogram
      implicit none
      INCLUDE 'mpif.h'
      integer im,ncpu,np,kd,nfa,n2dvisl,n2dvisa
      integer n_size,ne_size,n,ne,ib,nab,nad,nadto
      integer n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4

      ! input/output to hdf_lib libarary is in 8 byte integers
      integer*8 bufls,bufds,fidc,fidm
      integer ierr

      fidm=0
      fidc=1

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, im, ierr)

      if(im.eq.0)then
#ifdef USEIEPENTA
      write(*,*) 'using infinite domain condition'
#else
      write(*,*) 'using dirichlet boundary condition'
#endif
#ifdef PCGE
      write(*,*) 'using PCGE'
#endif
#ifdef TET4GRID
      write(*,*) 'using TET4GRID'
#endif
#ifdef OMPSAME
      write(*,*) 'using OMPSAME'
#else
      write(*,*) 'not using OMPSAME'
#endif
      write(6, *) "read mdata, cdata from ", DIR_MDATA_CDATA
      open(10,file='./data/para_setting.dat',status='old')
      read(10,*)
      read(10,*) np
      if(np.ne.ncpu)then
      write(*,*)'number of mpi processes are inconsistent',np,ncpu
      stop
      endif
      read(10,*)
      read(10,*) np
      close(10)
      open(10,file='./data/setting.dat',status='old')
      read(10,*) !'num of node'
      read(10,*) ! n
      read(10,*) !'num of elem'
      read(10,*) ! ne
      read(10,*) !'num of elem'
      read(10,*)  kd
      close(10)
      call check_macro_parameters
      endif
      call MPI_BCAST(np,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! #ifdef USE_SPLIT_NODE
!       nfa=0
! #else
!       if(im.eq.0)then
!       open(41,file='./data/faultpara.dat',status='old')
!       read(41,*)
!       read(41,*) nfa
!       close(41)
!       endif
!       call MPI_BCAST(nfa,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! #endif

      call read_local_setting_hdf(im,kd, &
            n,ne,ib,nab,nad,nadto,fidm,0,DIR_MDATA_CDATA)
      call compute_nsize(n,n_size)
      call compute_nsize(ne,ne_size)
      call read_local_setting_hdf(im,kd, &
            n4,ne4,ib4,nab4,nad4,nadto4,fidc,1,DIR_MDATA_CDATA)
      call compute_nsize(n4,n4_size)
      call compute_nsize(ne4,ne4_size)
      call read_outputflag_settings &
            (0,im,n2dvisl,n2dvisa,DIR_MDATA_CDATA)
      bufds=n*30
      bufls=n*30

      call main_compute_ni(im,ncpu,np,kd,nfa,n2dvisl, &
            n_size,ne_size,n,ne,ib,nab,nad,nadto,fidm, &
            n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,fidc, &
            bufds,bufls)
      call MPI_FINALIZE(ierr)

      end program mainprogram

      subroutine main_compute_ni(im,ncpu,np,kd,nfa,n2dvisl, &
            n_size,ne_size,n,ne,ib,nab,nad,nadto,fidm, &
            n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,fidc, &
            bufds,bufls)
      implicit none
      INCLUDE 'mpif.h'
      integer im,ncpu,np,nfa,kd,n2dvisl
      real*8 plane(5),ds
      integer ierr,i
      character dataname*50
      ! input/output to hdf_lib library is in 8 byte integers
      integer*8 bufls,bufds,buflc,bufdc,fidc,fidm,ltmp
      integer*8 lsize,lstart

      integer n_size,ne_size,ni_size,nei_size,n,ne,ni,nei
      integer ib,nab,nad,nadto
      real*8 coor(3,n_size)
      integer n4_size,ne4_size,ni4_size,nei4_size,n4,ne4,ni4,nei4
      integer ib4,nab4,nad4,nadto4
      real*8 coor4(3,n4_size)
#ifdef READASCII
      integer cny(10,ne_size),num(ne_size)
      integer cny4(4,ne4_size),num4(ne4_size)
#endif

      if(im.eq.0)then

      open(61,file='./data/modeldomain.dat',status='unknown')
      read(61,*) !'xmin-plane'
      read(61,*) plane(1)
      read(61,*) !'xmax-plane'
      read(61,*) plane(2)
      read(61,*) !'ymin-plane'
      read(61,*) plane(3)
      read(61,*) !'ymin-plane'
      read(61,*) plane(4)
      read(61,*) !'zmin-plane'
      read(61,*) plane(5)
      close(61)

      open(61,file='./data/modeling_setting.dat',status='unknown')
      read(61,*) !'nx, ny'
      read(61,*) !'nx, ny'
      read(61,*) !'ds'
      read(61,*) ds
      close(61)

      endif
      call MPI_BCAST(plane,5,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ds,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      write(dataname,'(a5)') '/coor'
      ltmp=5
      call set_null_char(dataname,ltmp)

      call hdf_read_double_array(fidm,dataname,coor,bufds,bufdc)
#ifdef READASCII
      call read_geometry(n,ne,coor,cny,num)
#endif
#ifdef USEIEPENTA
      call set_ni(n_size,n,coor,plane,ds,ni)
      nei=nab
#else
      ni=0
      nei=0
#endif
      call compute_nsize(ni,ni_size)
      call compute_nsize(nei,nei_size)
      if(im.eq.0)then
      write(*,*) 'ni,nei,ds',ni,nei,ds
      endif

      call hdf_read_double_array(fidc,dataname,coor4,bufds,bufdc)
#ifdef READASCII
      call read_geometry4(n4,ne4,coor4,cny4,num4)
#endif
#ifdef USEIEPENTA
      call set_ni(n4_size,n4,coor4,plane,ds,ni4)
      nei4=nab4
#else
      ni4=0
      nei4=0
#endif
      call compute_nsize(ni4,ni4_size)
      call compute_nsize(nei4,nei4_size)

      if(im.eq.0)then
      write(*,*) 'ni4,nei4',ni4,nei4
      endif

      call main(im,ncpu,np,kd,nfa,plane,ds,n2dvisl, &
            n_size,ne_size,ni_size,nei_size,n,ne,ni,nei, &
            ib,nab,nad,nadto,fidm, &
            n4_size,ne4_size,ni4_size,nei4_size,n4,ne4,ni4,nei4, &
            ib4,nab4,nad4,nadto4,fidc, &
            bufds,bufls)
      end

      subroutine main(im,ncpu,np,kd,nfa,plane,ds,n2dvisl, &
            n_size,ne_size,ni_size,nei_size,n,ne,ni,nei, &
            ib,nab,nad,nadto,fidm, &
            n4_size,ne4_size,ni4_size,nei4_size,n4,ne4,ni4,nei4, &
            ib4,nab4,nad4,nadto4,fidc, &
            bufds,bufls)
      implicit none
      INCLUDE 'mpif.h'
      integer im,ncpu,np,n2dvisl,nfa,kd
      integer ie,itmp,i,ki,iet,ii,it,i1,i2,i3,j

      integer*8 bufls,bufds,buflc,bufdc,fidc,fidm,ltmp
      character dataname*50

      real*8 plane(5),ds
      integer iad

      integer n_size,ne_size,nei_size,ni_size,n,ne,ni,nei
      integer ib,nab,nad,nadto
      real*8 coor(3,n_size),coori(3,n_size+ni_size),dsi
      real*8 uval(6,n_size+ni_size)
      real*8 bcarray(n_size+ni_size)
      integer cny(10,ne_size),num(ne_size)
      integer cnyi(12,nei_size),numi(nei_size),nfla(nei_size)
      real*8 dupli(n_size),duplixyz(3*(n_size+ni_size))
      integer ibi(ib),nabi(6,nab),npl(np+1),abcelem(nab)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      integer mpiadi(nad,2),mpiadlisti(2*nadto),mpli(nad+1),nadtoi
      REAL_4 coors(3,n_size)
!      REAL_4 cooris(3,n_size+ni_size)

! #define CRSSIZE 200
! #define CRSSIZE 1000
#define CRSSIZE 5000
!#define nsta 5

      integer ne4_size,n4_size,ni4_size,nei4_size,n4,ne4,ni4,nei4
      integer ib4,nab4,nad4,nadto4
      real*8 coor4(3,n4_size),coori4(3,n4_size+ni4_size),dsi4
      real*8 uval4(6,n4_size+ni4_size)
      real*8 bcarray4(n4_size+ni4_size)
      integer cny4(4,ne4_size),num4(ne4_size)
      integer cnyi4(6,nei4_size),numi4(nei4_size),nfla4(nei4_size)
      real*8 dupli4(n4_size),duplixyz4(3*(n4_size+ni4_size))
      integer ibi4(ib4),nabi4(3,nab4),npl4(np+1),abcelem4(nab4)
      integer mpiad4(nad4,2),mpiadlist4(nadto4),mpl4(nad4+1)
      integer mpiadi4(nad4,2),mpiadlisti4(2*nadto4)
      integer mpli4(nad4+1),nadtoi4
!      REAL_4 coors4(3,n4_size)
!      REAL_4 cooris4(3,n4_size+ni4_size)

      ! for outer
      integer nird,niwd
      integer pentawptrd(2*nsta+3),pentarptrd((n+ni)*4)
      integer pentawindd((n+ni)*4),pentarindd((n+ni)*CRSSIZE)
      ! for inner fine
      integer nir,niw
      integer pentawptr(2*nsta+3),pentarptr((n+ni)*4)
      integer pentawind((n+ni)*4),pentarind((n+ni)*CRSSIZE)
      ! for inner coarse
      integer nir4,niw4
      integer pentawptr4(2*nsta+3),pentarptr4((n4+ni4)*4)
      integer pentawind4((n4+ni4)*4),pentarind4((n4+ni4)*CRSSIZE)

#if npc == 1 || npc == 21
      real*8 kcrsvald(9,(n+ni)*CRSSIZE)
      real*8 kcrsval(9,(n+ni)*CRSSIZE)
      real*8 kcrsval4(9,(n4+ni4)*CRSSIZE)
#elif npc == 288
      real*8 kcrsvald(9,(n+ni)*CRSSIZE*2)
      real*8 kcrsval(9,(n+ni)*CRSSIZE*2)
      real*8 kcrsval4(9,(n4+ni4)*CRSSIZE*2)
#endif
      integer map4to10(2,n+ni)

!      integer map_f_to_c(n4),map_c_to_f(2,n)
!      real*8 up4(3*(n4+ni4)*npc),rv4(3*(n4+ni4)*npc)

      integer cnysep_read(10,ne),cnysep_write(10,ne)
      integer cnysep_read4(4,ne),cnysep_write4(4,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep_n_ptr4(nsta+2)
      integer nsep,nsepnmap(2*n)
      integer nsep4,nsepnmap4(2*n4)
      integer numsep4(ne),elemmap(ne)
      integer ncolor,color_ind(np*MAXCOLOR)
      integer cnysep_readm(10,ne),cnysep_readm4(4,ne)

      real*8 rmat(10,kd),younglst(2,kd)
      real*8 cort(3),rf(3),ftmp(30),xx(4,3),strike,dip,rake
      real*8 corf(nfa,3),faultp(nfa,4),sig(nkl)
      real*8 up(3*(n+ni)*npc),rv(3*(n+ni)*npc),up2
      real*4 u2dl(3*n2dvisl*npc),utmp(3*(n+ni)),utmp2(3*(n+ni))
      REAL_4 me(10,10),tmps,tmp2s
      integer nr,j1,netmp
      real*8 young,rnyu,cori(3,3),keit(36,36),keit4(18,18)
      real*8 t1,t2
      integer flag(n_size)
      integer iee,id,ierr,ipc,icpu
      integer ibsize,fh,mpicount
      integer u2dflagl(n2dvisl)
      integer mpistat(MPI_STATUS_SIZE)
      real*8 mpibuf(2*3*npc*(nadto+1)*2)
      integer(kind=MPI_OFFSET_KIND) mpifsbuf
      character filename*50
      integer iunit

      integer :: node_old2new(n), node_old2new4(n4)
      real*8 uvalcoe(npc)

      logical log_innerCG,log_innerCG_per_outite

      logical :: use_crs_d, overlap_comm_d
      logical :: overlap_comm_s, persistent_comm_s

#ifdef GFLIB
      ! for GFLIB
      integer :: nload, iload, load_elem_arr(ne), surf_nelem
      integer, allocatable :: surf_cny(:,:)
      real*8, allocatable :: load_arr(:,:)
#endif

      write(*,*) 'here111'

      call IPCG_setting &
            (im,use_crs_d,overlap_comm_d, &
            overlap_comm_s,persistent_comm_s)

      call read_outputflag(0,im,n2dvisl,u2dflagl,DIR_MDATA_CDATA)
      call read_log_setting(im,log_innerCG,log_innerCG_per_outite)

      if(im.eq.0)then
      write(*,*) 'read mdata'
      endif
      call read_local_data_hdf &
            (n_size,ne_size,n,ne,ib,nab,nad,nadto, &
            np,bufls,bufds,mpiad,mpiadlist,mpl,dupli, &
            coor,num,cny,ibi,nabi,abcelem,fidm,10,6)
      if(overlap_comm_s .or. overlap_comm_d)then
      call reorder_nodes &
            (im,n,ne,10,nadto,mpiadlist,dupli,coor,cny,node_old2new)
      call reorder_outflag(im,n,node_old2new,n2dvisl,u2dflagl)
      endif
#ifdef READASCII
      call read_geometry(n,ne,coor,cny,num)
#endif
      bcarray=1
      do i=1,n
      if( (abs(coor(1,i)-plane(1)).le.ds*0.01) .or. &
          (abs(coor(1,i)-plane(2)).le.ds*0.01) .or. &
          (abs(coor(2,i)-plane(3)).le.ds*0.01) .or. &
          (abs(coor(2,i)-plane(4)).le.ds*0.01) .or. &
          (abs(coor(3,i)-plane(5)).le.ds*0.01) )then
      bcarray(i)=0
      endif
      enddo
! length of infinite element in infinite direction
      dsi=5.0*ds
      do i=1,n
      coori(1,i)=coor(1,i)
      coori(2,i)=coor(2,i)
      coori(3,i)=coor(3,i)
      enddo

      if(im.eq.0)then
      write(*,*) 'read cdata'
      endif
      call read_local_data_hdf &
            (n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4, &
            np,bufls,bufds,mpiad4,mpiadlist4,mpl4,dupli4, &
            coor4,num4,cny4,ibi4,nabi4,abcelem4,fidc,4,3)
      if(overlap_comm_s)then
      call reorder_nodes &
            (im,n4,ne4,4,nadto4,mpiadlist4,dupli4,coor4,cny4,node_old2new4)
      endif
#ifdef READASCII
      call read_geometry4(n4,ne4,coor4,cny4,num4)
#endif
      bcarray4=1
      do i=1,n4
      if( (abs(coor4(1,i)-plane(1)).le.ds*0.01) .or. &
          (abs(coor4(1,i)-plane(2)).le.ds*0.01) .or. &
          (abs(coor4(2,i)-plane(3)).le.ds*0.01) .or. &
          (abs(coor4(2,i)-plane(4)).le.ds*0.01) .or. &
          (abs(coor4(3,i)-plane(5)).le.ds*0.01) )then
      bcarray4(i)=0
      endif
      enddo
      dsi4=dsi
      do i=1,n4
      coori4(1,i)=coor4(1,i)
      coori4(2,i)=coor4(2,i)
      coori4(3,i)=coor4(3,i)
      enddo

      call remap_material_id(im,kd,ne,num)
      call remap_material_id(im,kd,ne4,num4)

      if(im.eq.0)then
            print *, "read material"
      call read_material(kd,rmat)
      endif
      call MPI_BCAST(rmat,10*kd,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call makmatlist(kd,rmat,younglst)

#ifdef USEIEPENTA
      call setIEpenta &
            (im,n,ne,cny,num,cnyi,numi,nfla,ni,nad,plane, &
            nei,coor,coori,ds, &
            nadto,mpiad,mpiadlist,mpl, &
            nadtoi,mpiadi,mpiadlisti,mpli,10,12)

      call setIEpenta &
            (im,n4,ne4,cny4,num4,cnyi4,numi4,nfla4,ni4,nad4,plane, &
            nei4,coor4,coori4,ds, &
            nadto4,mpiad4,mpiadlist4,mpl4, &
            nadtoi4,mpiadi4,mpiadlisti4,mpli4,4,6)
#else
      mpli=mpl
      nadtoi=nadto
      mpiadi=mpiad
      do i=1,nadto
      mpiadlisti(i)=mpiadlist(i)
      enddo

      mpli4=mpl4
      nadtoi4=nadto4
      mpiadi4=mpiad4
      do i=1,nadto4
      mpiadlisti4(i)=mpiadlist4(i)
      enddo
#endif
      netmp = 0
      if (im.eq.0) print *, "set_whole_CRS"
      if(use_crs_d) netmp = ne
      call set_whole_CRS &
            (im,kd,n,ni,netmp,nei,10,12,cny,cnyi,num,numi, &
            nfla,younglst,coori,dsi, &
            niwd,pentawptrd,pentawindd, &
            nird,pentarptrd,pentarindd,kcrsvald, &
            nadtoi,mpiadlisti,overlap_comm_d,CRSSIZE)

      call set_whole_CRS &
            (im,kd,n,ni,ne,nei,10,12,cny,cnyi,num,numi, &
            nfla,younglst,coori,dsi, &
            niw,pentawptr,pentawind,nir,pentarptr,pentarind,kcrsval, &
            nadtoi,mpiadlisti,overlap_comm_s,CRSSIZE)

      call set_whole_CRS &
            (im,kd,n4,ni4,ne4,nei4,4,6,cny4,cnyi4,num4,numi4, &
            nfla4,younglst,coori4,dsi4, &
            niw4,pentawptr4,pentawind4, &
            nir4,pentarptr4,pentarind4,kcrsval4, &
            nadtoi4,mpiadlisti4,overlap_comm_s,CRSSIZE)

      if (im.eq.0) print *, "setcnysep"
      call setcnysep(im,n,ne,n4,num,cny,cny4, &
            nsep,nsep4,numsep, &
            nsep_n_ptr,nsep_n_ptr4,nsep_ne_ptr,nsepnmap,nsepnmap4, &
            cnysep_read,cnysep_write,cnysep_read4,cnysep_write4) 
      numsep4=numsep
      if (im.eq.0) print *, "setcolor"
      call setcolor &
            (im,nsep,ne_size,10,4,nsep,ne,np,cnysep_write,numsep, &
            ncolor,color_ind,elemmap,cnysep_write4,numsep4)
      do ie=1,ne
      cnysep_readm(:,ie)=cnysep_read(:,elemmap(ie))
      cnysep_readm4(:,ie)=cnysep_read4(:,elemmap(ie))
      enddo

      if (im.eq.0) print *, "setduplixyz"
      call setduplixyz(n,ni,im,mpibuf, &
            nad,nadtoi,mpiadi,mpiadlisti,mpli,duplixyz)
      call setduplixyz(n4,ni4,im,mpibuf, &
            nad4,nadtoi4,mpiadi4,mpiadlisti4,mpli4,duplixyz4)

      if (im.eq.0) print *, "read_sig"
      call read_sig(im,sig)
      call read_uvalcoe(im,uvalcoe)

#ifdef USE_SPLIT_NODE
      call set_fault_input &
            (im,kd,mpibuf,younglst, &
            n_size,ne_size,n,ne,coor,cny,num, &
            nad,nadtoi,mpiadi,mpiadlisti,mpli, &
            rv &
#ifdef MEASURE_TIME_BARRIER
            ,mp,timerbufd,timerbufi,ind &
#endif
            )
#else
#endif

      uval=0
#ifdef USEIEPENTA
      do ie=1,nei
      nr=numi(ie)
      young=younglst(1,nr)
      rnyu=younglst(2,nr)
      do i1=1,3
      do ii=1,3
      cori(i1,ii)=coori(ii,cnyi(i1,ie))
      enddo
      enddo
      if(nfla(ie).eq.1)then
      call IEpenta2nd1(young,rnyu,cori,dsi,keit)
      endif
      if(nfla(ie).eq.2)then
      call IEpenta2nd2(young,rnyu,cori,dsi,keit)
      endif
      if(nfla(ie).eq.3)then
      call IEpenta2nd3(young,rnyu,cori,dsi,keit)
      endif
      if(nfla(ie).eq.4)then
      call IEpenta2nd4(young,rnyu,cori,dsi,keit)
      endif
      if(nfla(ie).eq.5)then
      call IEpenta2nd5(young,rnyu,cori,dsi,keit)
      endif
      do j1=1,12
      i1=cnyi(j1,ie)
      i2=3*(j1-1)
      uval(1,i1)=uval(1,i1)+keit(i2+1,i2+1)
      uval(2,i1)=uval(2,i1)+keit(i2+2,i2+2)
      uval(3,i1)=uval(3,i1)+keit(i2+3,i2+3)
      uval(4,i1)=uval(4,i1)+keit(i2+1,i2+2)
      uval(5,i1)=uval(5,i1)+keit(i2+1,i2+3)
      uval(6,i1)=uval(6,i1)+keit(i2+2,i2+3)
      enddo
      enddo
#endif
      call compblockdiaginv &
            (n_size,ne_size,n,ne,kd,num,cny,coori,younglst, &
            ni_size,nei_size,ni,nei, &
            nad,nadtoi,mpiadi,mpiadlisti,mpli,im,mpibuf, &
#ifdef MEASURE_TIME_BARRIER
            ,mp,timerbufd,timerbufi,ind, &
#endif
            uval)

      uval4=0
#ifdef USEIEPENTA
      do ie=1,nei4
      nr=numi4(ie)
      young=younglst(1,nr)
      rnyu=younglst(2,nr)
      do i1=1,3
      do ii=1,3
      cori(i1,ii)=coori4(ii,cnyi4(i1,ie))
      enddo
      enddo
      if(nfla4(ie).eq.1)then
      call IEpenta1st1(young,rnyu,cori,dsi4,keit4)
      endif
      if(nfla4(ie).eq.2)then
      call IEpenta1st2(young,rnyu,cori,dsi4,keit4)
      endif
      if(nfla4(ie).eq.3)then
      call IEpenta1st3(young,rnyu,cori,dsi4,keit4)
      endif
      if(nfla4(ie).eq.4)then
      call IEpenta1st4(young,rnyu,cori,dsi4,keit4)
      endif
      if(nfla4(ie).eq.5)then
      call IEpenta1st5(young,rnyu,cori,dsi4,keit4)
      endif
      do j1=1,6
      i1=cnyi4(j1,ie)
      i2=3*(j1-1)
      uval4(1,i1)=uval4(1,i1)+keit4(i2+1,i2+1)
      uval4(2,i1)=uval4(2,i1)+keit4(i2+2,i2+2)
      uval4(3,i1)=uval4(3,i1)+keit4(i2+3,i2+3)
      uval4(4,i1)=uval4(4,i1)+keit4(i2+1,i2+2)
      uval4(5,i1)=uval4(5,i1)+keit4(i2+1,i2+3)
      uval4(6,i1)=uval4(6,i1)+keit4(i2+2,i2+3)
      enddo
      enddo
#endif
      call compblockdiaginv4 &
            (n4_size,ne4_size,n4,ne4,kd,num4,cny4,coori4,younglst, &
            ni4_size,nei4_size,ni4,nei4, &
            nad,nadtoi4,mpiadi4,mpiadlisti4,mpli4,im,mpibuf, &
#ifdef MEASURE_TIME_BARRIER
            ,mp,timerbufd,timerbufi,ind, &
#endif
            uval4)

      call map4to10_init &
            (im,ds, &
            n,ne,ni,nei,coori,cny,cnyi, &
            n4,ne4,ni4,nei4,coori4,cny4,cnyi4, &
            map4to10)

      call print_problem_size(im,n,ne,nadto,ni,n4,nadto4,ni4, &
            niw,nir,niw4,nir4)

#ifdef UNITFAULT
!       do iunit =1,56
       do iunit =1,1
      if(im.eq.0)then
      write(*,*) 'USE POINT SOURCE'
      endif
      call set_fault_point &
            (n,n_size,ne,ne_size,ni,ni_size,im,rv,coor,cny,iunit,ds)
#endif

#ifdef GFLIB
      if (im.eq.0) then
            print *, "read_loads"
            call read_nload(nload)
      endif
      call MPI_BCAST(nload, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      allocate (load_arr(5, nload))
      if (im.eq.0) then
            call read_loads(nload, load_arr)
      endif
      call MPI_BCAST(load_arr, 5*nload, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      call read_surf_nelem(im, surf_nelem)
      if (surf_nelem > 0) then
            allocate (surf_cny(4, surf_nelem))
            call read_surf_cny(im, surf_nelem, surf_cny)
      endif

      load_elem_arr = 0
      do iload = 1, nload
      call pointload(im, n, coor, ne, cny, ni, rv, &
            load_arr, nload, iload, load_elem_arr, &
            surf_nelem, surf_cny)
#endif
      call sync_parallel_block_d &
            (n+ni,n+ni, &
            rv,nad,nadtoi,mpiadi,mpiadlisti,mpli,ierr,im,mpibuf &
#ifdef MEASURE_TIME_BARRIER
            ,mp,timerbufd,timerbufi,ind &
#endif
            )

     
      up=0.0

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t1=MPI_WTIME()
#ifdef USE_PA
      call fapp_start('IPCG', 0, 0)
#endif
      call IPCG &
            (im,np,kd,younglst,sig,mpibuf, &
            n,ne,coori,cny,num,uval,duplixyz,bcarray, &
            ni,nei,cnyi,dsi,numi,nfla, &
            nad,nadtoi,mpiadi,mpiadlisti,mpli, &
            nird,pentarptrd,pentarindd, &
            niwd,pentawptrd,pentawindd,kcrsvald, &
            nir,pentarptr,pentarind, &
            niw,pentawptr,pentawind,kcrsval, &
            nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap, &
            cnysep_readm,cnysep_write, &
            ncolor,color_ind, &
            n4,ne4,coori4,cny4,num4,uval4,duplixyz4,bcarray4, &
            ni4,nei4,cnyi4,dsi4,numi4,nfla4, &
            nad4,nadtoi4,mpiadi4,mpiadlisti4,mpli4, &
            nir4,pentarptr4,pentarind4, &
            niw4,pentawptr4,pentawind4,kcrsval4, &
            map4to10, &
            cnysep_readm4,cnysep_write4,nsep_n_ptr4,nsep4,nsepnmap4, &
            uvalcoe,log_innerCG,log_innerCG_per_outite, &
            use_crs_d,overlap_comm_d,overlap_comm_s,persistent_comm_s, &
            rv,up)
#ifdef USE_PA
      call fapp_stop('IPCG', 0, 0)
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t2=MPI_WTIME()
      if(im.eq.0)then
        write(*,*) 'IPCG took',t2-t1
      endif

!      do i=1,3*n
!      utmp2(i)=up(npc*(i-1)+1)
!      enddo
!      coors=coor
!      call write_tet_vtk_s
!     - (n,ne,n,ne,10,cny,coors,utmp2,im,2)
!      call write_tet_vtk_s
!     - (n,ne,n,ne,10,cny,coors,utmp2-utmp,im,3)
      call inner_product_d_block_sync(im,3*n,up,up,duplixyz,up2)
      if(im.eq.0)then
      write(*,*) 'up2block',up2
      endif
      call inner_product_d_block_sync(im,3*(n+ni),up,up,duplixyz,up2)
      if(im.eq.0)then
      write(*,*) 'up2block_infinite',up2
      endif

      if(n2dvisl.gt.0)then
      do i=1,n2dvisl
      do j=1,npc
      u2dl(3*(i-1)+1+(j-1)*n2dvisl)=up(3*npc*(u2dflagl(i)-1)+0*npc+j)
      u2dl(3*(i-1)+2+(j-1)*n2dvisl)=up(3*npc*(u2dflagl(i)-1)+1*npc+j)
      u2dl(3*(i-1)+3+(j-1)*n2dvisl)=up(3*npc*(u2dflagl(i)-1)+2*npc+j)
      enddo
      enddo
      write(filename,'(a11,i6.6,a6,i6.6,a4)') &
            './2Doutput/',im,'.disp.',0,'.bin'
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, &
            MPI_INFO_NULL,fh,ierr)
      mpicount=n2dvisl*3*npc
      mpifsbuf=mpicount*4
      call MPI_FILE_SET_SIZE(fh,mpifsbuf,ierr)
      call MPI_FILE_WRITE(fh,u2dl,mpicount,MPI_REAL4,mpistat,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      endif
#ifdef UNITFAULT
      enddo
#endif
#ifdef GFLIB
      call output_displacement(im, n, ni, up, iload)
      enddo
      call output_load_elem(im, ne, load_elem_arr)
#endif
      end
