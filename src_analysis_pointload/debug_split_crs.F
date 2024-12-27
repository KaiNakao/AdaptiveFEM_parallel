!#define READASCII
!#define PCGE
!#define USE_SPLIT_NODE
c
      program debug2
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

#ifdef USE_SPLIT_NODE
      nfa=0
#else
      if(im.eq.0)then
      open(41,file='./data/faultpara.dat',status='old')
      read(41,*)
      read(41,*) nfa
      close(41)
      endif
      call MPI_BCAST(nfa,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      call read_local_setting_hdf(im,kd,
     & n,ne,ib,nab,nad,nadto,fidm,0)
      call compute_nsize(n,n_size)
      call compute_nsize(ne,ne_size)
      call read_local_setting_hdf(im,kd,
     & n4,ne4,ib4,nab4,nad4,nadto4,fidc,1)
      call compute_nsize(n4,n4_size)
      call compute_nsize(ne4,ne4_size)
      call read_outputflag_settings(0,im,n2dvisl,n2dvisa)
      bufds=n*30
      bufls=n*30
c
      call debug2_main_compute_ni(im,ncpu,np,kd,nfa,n2dvisl,
     & n_size,ne_size,n,ne,ib,nab,nad,nadto,fidm,
     & n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,fidc,
     & bufds,bufls)
      call MPI_FINALIZE(ierr)
c
      end program debug2
c_______________________________________________________________________
      subroutine debug2_main_compute_ni(im,ncpu,np,kd,nfa,n2dvisl,
     & n_size,ne_size,n,ne,ib,nab,nad,nadto,fidm,
     & n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,fidc,
     & bufds,bufls)
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

      call debug2_main(im,ncpu,np,kd,nfa,plane,ds,n2dvisl,
     & n_size,ne_size,ni_size,nei_size,n,ne,ni,nei,
     & ib,nab,nad,nadto,fidm,
     & n4_size,ne4_size,ni4_size,nei4_size,n4,ne4,ni4,nei4,
     & ib4,nab4,nad4,nadto4,fidc,
     & bufds,bufls)
c
      end
c_______________________________________________________________________
      subroutine debug2_main(im,ncpu,np,kd,nfa,plane,ds,n2dvisl,
     & n_size,ne_size,ni_size,nei_size,n,ne,ni,nei,
     & ib,nab,nad,nadto,fidm,
     & n4_size,ne4_size,ni4_size,nei4_size,n4,ne4,ni4,nei4,
     & ib4,nab4,nad4,nadto4,fidc,
     & bufds,bufls)
      implicit none
      INCLUDE 'mpif.h'
      integer im,ncpu,np,n2dvisl,nfa,kd
      integer ie,itmp,i,ki,iet,ii,it,i1,i2,i3,j

      integer*8 bufls,bufds,buflc,bufdc,fidc,fidm,ltmp
      character dataname*50

      real*8 plane(5),ds
      integer iad
c
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

#define CRSSIZE 200
!#define nsta 5
      integer nir,niw
      integer pentawptr(2*nsta+3),pentarptr((n+ni)*4)
      integer pentawind((n+ni)*4),pentarind((n+ni)*CRSSIZE)
      integer niro,niwo
      integer pentawptro(2*nsta+3),pentarptro((n+ni)*4)
      integer pentawindo((n+ni)*4),pentarindo((n+ni)*CRSSIZE)
      integer niri,niwi
      integer pentawptri(2*nsta+3),pentarptri((n+ni)*4)
      integer pentawindi((n+ni)*4),pentarindi((n+ni)*CRSSIZE)
c
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
      integer nir4,niw4
      integer pentawptr4(2*nsta+3),pentarptr4((n4+ni4)*4)
      integer pentawind4((n4+ni4)*4),pentarind4((n4+ni4)*CRSSIZE)
      integer niro4,niwo4
      integer pentawptro4(2*nsta+3),pentarptro4((n4+ni4)*4)
      integer pentawindo4((n4+ni4)*4),pentarindo4((n4+ni4)*CRSSIZE)

#if npc == 1 || npc == 21
      real*8 kcrsval(9,(n+ni)*CRSSIZE)
      real*8 kcrsval4(9,(n4+ni4)*CRSSIZE)
#elif npc == 288
      real*8 kcrsval(9,(n+ni)*CRSSIZE*2)
      real*8 kcrsvalo(9,(n+ni)*CRSSIZE*2)
      real*8 kcrsvali(9,(n+ni)*CRSSIZE*2)
      real*8 kcrsval4(9,(n4+ni4)*CRSSIZE*2)
      real*8 kcrsvalo4(9,(n4+ni4)*CRSSIZE*2)
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
c
      real*8 rmat(10,kd),younglst(2,kd)
      real*8 cort(3),rf(3),ftmp(30),xx(4,3),strike,dip,rake
      real*8 corf(nfa,3),faultp(nfa,4),sig(nkl)
      real*8 up(3*(n+ni)*npc),rv(3*(n+ni)*npc),up2
      real*4 u2dl(3*n2dvisl*npc),utmp(3*(n+ni)),utmp2(3*(n+ni))
      REAL_4 me(10,10),tmps,tmp2s
      integer nr,j1
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

      real*8 uvalcoe(npc)

      logical log_innerCG,log_innerCG_per_outite

      call read_outputflag(0,im,n2dvisl,u2dflagl)
      call read_log_setting(im,log_innerCG,log_innerCG_per_outite)

      if(im.eq.0)then
      write(*,*) 'read mdata'
      endif
      call read_local_data_hdf
     - (n_size,ne_size,n,ne,ib,nab,nad,nadto,
     - np,bufls,bufds,mpiad,mpiadlist,mpl,dupli,
     - coor,num,cny,ibi,nabi,abcelem,fidm,10,6)
#ifdef READASCII
      call read_geometry(n,ne,coor,cny,num)
#endif
      bcarray=1
      do i=1,n
      if( (abs(coor(1,i)-plane(1)).le.ds*0.01) .or.
     &    (abs(coor(1,i)-plane(2)).le.ds*0.01) .or.
     &    (abs(coor(2,i)-plane(3)).le.ds*0.01) .or.
     &    (abs(coor(2,i)-plane(4)).le.ds*0.01) .or.
     &    (abs(coor(3,i)-plane(5)).le.ds*0.01) )then
      bcarray(i)=0
      endif
      enddo
c length of infinite element in infinite direction
      dsi=5.0*ds
      do i=1,n
      coori(1,i)=coor(1,i)
      coori(2,i)=coor(2,i)
      coori(3,i)=coor(3,i)
      enddo

      if(im.eq.0)then
      write(*,*) 'read cdata'
      endif
      call read_local_data_hdf
     - (n4_size,ne4_size,n4,ne4,ib4,nab4,nad4,nadto4,
     - np,bufls,bufds,mpiad4,mpiadlist4,mpl4,dupli4,
     - coor4,num4,cny4,ibi4,nabi4,abcelem4,fidc,4,3)
#ifdef READASCII
      call read_geometry4(n4,ne4,coor4,cny4,num4)
#endif
      bcarray4=1
      do i=1,n4
      if( (abs(coor4(1,i)-plane(1)).le.ds*0.01) .or.
     &    (abs(coor4(1,i)-plane(2)).le.ds*0.01) .or.
     &    (abs(coor4(2,i)-plane(3)).le.ds*0.01) .or.
     &    (abs(coor4(2,i)-plane(4)).le.ds*0.01) .or.
     &    (abs(coor4(3,i)-plane(5)).le.ds*0.01) )then
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
      call read_material(kd,rmat)
      endif
      call MPI_BCAST(rmat,10*kd,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call makmatlist(kd,rmat,younglst)
c
#ifdef USEIEPENTA
      call setIEpenta
     & (im,n,ne,cny,num,cnyi,numi,nfla,ni,nad,plane,
     & nei,coor,coori,ds,
     & nadto,mpiad,mpiadlist,mpl,
     & nadtoi,mpiadi,mpiadlisti,mpli,10,12)

      call setIEpenta
     & (im,n4,ne4,cny4,num4,cnyi4,numi4,nfla4,ni4,nad4,plane,
     & nei4,coor4,coori4,ds,
     & nadto4,mpiad4,mpiadlist4,mpl4,
     & nadtoi4,mpiadi4,mpiadlisti4,mpli4,4,6)
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
      call set_whole_CRS
     & (im,kd,n,ni,ne,nei,10,12,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,kcrsval,
     & nadtoi,mpiadlisti,.false.,CRSSIZE)
      call set_whole_CRS
     & (im,kd,n,ni,ne,nei,10,12,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niwo,pentawptro,pentawindo,niro,pentarptro,pentarindo,kcrsvalo,
     & nadtoi,mpiadlisti,.true.,CRSSIZE)
      call set_whole_CRS
     & (im,kd,n,ni,0,nei,10,12,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niwi,pentawptri,pentawindi,niri,pentarptri,pentarindi,kcrsvali,
     & nadtoi,mpiadlisti,.false.,CRSSIZE)

      call set_whole_CRS
     & (im,kd,n4,ni4,ne4,nei4,4,6,cny4,cnyi4,num4,numi4,
     & nfla4,younglst,coori4,dsi4,
     & niw4,pentawptr4,pentawind4,
     & nir4,pentarptr4,pentarind4,kcrsval4,
     & nadtoi4,mpiadlisti4,.false.,CRSSIZE)
      call set_whole_CRS
     & (im,kd,n4,ni4,ne4,nei4,4,6,cny4,cnyi4,num4,numi4,
     & nfla4,younglst,coori4,dsi4,
     & niwo4,pentawptro4,pentawindo4,
     & niro4,pentarptro4,pentarindo4,kcrsvalo4,
     & nadtoi4,mpiadlisti4,.true.,CRSSIZE)

      call setcnysep(im,n,ne,n4,num,cny,cny4,
     & nsep,nsep4,numsep,
     & nsep_n_ptr,nsep_n_ptr4,nsep_ne_ptr,nsepnmap,nsepnmap4,
     & cnysep_read,cnysep_write,cnysep_read4,cnysep_write4)
      numsep4=numsep
      call setcolor
     - (im,nsep,ne_size,10,4,nsep,ne,np,cnysep_write,numsep,
     -  ncolor,color_ind,elemmap,cnysep_write4,numsep4)
      do ie=1,ne
      cnysep_readm(:,ie)=cnysep_read(:,elemmap(ie))
      cnysep_readm4(:,ie)=cnysep_read4(:,elemmap(ie))
      enddo

      call setduplixyz(n,ni,im,mpibuf,
     & nad,nadtoi,mpiadi,mpiadlisti,mpli,duplixyz)
      call setduplixyz(n4,ni4,im,mpibuf,
     & nad4,nadtoi4,mpiadi4,mpiadlisti4,mpli4,duplixyz4)

      call read_sig(im,sig)
      call read_uvalcoe(im,uvalcoe)
c
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
      call compblockdiaginv
     & (n_size,ne_size,n,ne,kd,num,cny,coori,younglst,
     & ni_size,nei_size,ni,nei,
     & nad,nadtoi,mpiadi,mpiadlisti,mpli,im,mpibuf,
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind,
#endif
     & uval)

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
      call compblockdiaginv4
     & (n4_size,ne4_size,n4,ne4,kd,num4,cny4,coori4,younglst,
     & ni4_size,nei4_size,ni4,nei4,
     & nad,nadtoi4,mpiadi4,mpiadlisti4,mpli4,im,mpibuf,
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind,
#endif
     & uval4)
c
      call map4to10_init
     & (im,ds,
     & n,ne,ni,nei,coori,cny,cnyi,
     & n4,ne4,ni4,nei4,coori4,cny4,cnyi4,
     & map4to10)

      call debug_compAmat
     &(im,n,ni,ne,kd,num,cny,coor,younglst,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_readm,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kcrsval,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kcrsvalo,
     & niri,pentarptri,pentarindi,
     & niwi,pentawptri,pentawindi,kcrsvali,
     & sig,coori,
     & nad,nadtoi,mpiadi,mpiadlisti,mpli,duplixyz,
     & n4,ni4,
     & nir4,pentarptr4,pentarind4,
     & niw4,pentawptr4,pentawind4,kcrsval4,
     & niro4,pentarptro4,pentarindo4,
     & niwo4,pentawptro4,pentawindo4,kcrsvalo4,
     & nad4,nadto4,mpiad4,mpiadlist4,mpl4,duplixyz4,
     & coori4)

      end

      subroutine debug_compAmat
     &(im,n,ni,ne,kd,num,cny,coor,younglst,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & niri,pentarptri,pentarindi,
     & niwi,pentawptri,pentawindi,kgcrsvali,
     & sig,coori,
     & nad,nadto,mpiad,mpiadlist,mpl,duplixyz,
     & n4,ni4,
     & nir4,pentarptr4,pentarind4,
     & niw4,pentawptr4,pentawind4,kgcrsval4,
     & niro4,pentarptro4,pentarindo4,
     & niwo4,pentawptro4,pentawindo4,kgcrsvalo4,
     & nad4,nadto4,mpiad4,mpiadlist4,mpl4,duplixyz4,
     & coori4)
      implicit none
      integer, intent(in) :: im,n,ni,ne,kd
      integer, intent(in) :: num(ne),cny(10,ne)
      real(8), intent(in) :: coor(3,n),younglst(2,kd)
      integer, intent(in) :: nsep,nsepnmap(nsep)
      integer, intent(in) :: cnysep_read(10,ne),cnysep_write(10,ne)
      integer, intent(in) :: nsep_ne_ptr(nsta+2),numsep(ne)
      integer, intent(in) :: nsep_n_ptr(nsta+2)
      integer, intent(in) :: np,ncolor,color_ind(np*MAXCOLOR)
      integer, intent(in) :: nir,niw
      integer, intent(in) :: pentawptr(2*nsta+3),pentarptr(niw+1)
      integer, intent(in) :: pentawind(niw),pentarind(nir)
      real(8), intent(in) :: kgcrsval(9,nir,2)
      integer, intent(in) :: niro,niwo
      integer, intent(in) :: pentawptro(2*nsta+3),pentarptro(niwo+1)
      integer, intent(in) :: pentawindo(niwo),pentarindo(niro)
      real(8), intent(in) :: kgcrsvalo(9,niro,2)
      integer, intent(in) :: niri,niwi
      integer, intent(in) :: pentawptri(2*nsta+3),pentarptri(niwi+1)
      integer, intent(in) :: pentawindi(niwi),pentarindi(niri)
      real(8), intent(in) :: kgcrsvali(9,niri,2)
      real(8), intent(in) :: sig(nkl)
      real(8), intent(in) :: coori(3,n+ni)
      integer, intent(in) :: n4,ni4
      integer, intent(in) :: nir4,niw4
      integer, intent(in) :: pentawptr4(2*nsta+3),pentarptr4(niw+1)
      integer, intent(in) :: pentawind4(niw4),pentarind4(nir4)
      real(8), intent(in) :: kgcrsval4(9,nir4,2)
      integer, intent(in) :: niro4,niwo4
      integer, intent(in) :: pentawptro4(2*nsta+3),pentarptro4(niwo4+1)
      integer, intent(in) :: pentawindo4(niwo4),pentarindo4(niro4)
      real(8), intent(in) :: kgcrsvalo4(9,niro4,2)
      real(8), intent(in) :: coori4(3,n+ni)

      integer, intent(in) :: nad,nadto
      integer, intent(in) :: mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real(8), intent(in) :: duplixyz(3*(n+ni))
      integer, intent(in) :: nad4,nadto4
      integer, intent(in) :: mpiad4(nad4,2),mpiadlist4(nadto4)
      integer, intent(in) :: mpl4(nad4+1)
      real(8), intent(in) :: duplixyz4(3*(n4+ni4))

      real(8) :: uvalcoe(npc)
      REAL_4 :: kgcrsvals(9,nir,2), kgcrsval4s(9,nir4,2)
      REAL_4 :: kgcrsvalos(9,niro,2), kgcrsvalo4s(9,niro4,2)
      REAL_4 :: sigs(nsta*2), cooris(3,n+ni), coori4s(3,n4+ni4)
      REAL_4 :: uvalcoes(npc),duplixyzs(3*(n+ni)),duplixyz4s(3*(n4+ni4))

      call read_uvalcoe(im,uvalcoe)

      kgcrsvals = kgcrsval
      kgcrsval4s = kgcrsval4
      kgcrsvalos = kgcrsvalo
      kgcrsvalo4s = kgcrsvalo4
      sigs = sig
      cooris = coori
      coori4s = coori4
      uvalcoes = uvalcoe
      duplixyzs = duplixyz
      duplixyz4s = duplixyz4

      if (im == 0) write(6, *) "compAmatd"
      call debug_compAmatd
     &(im,n,ni,ne,kd,num,cny,coor,younglst,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & niri,pentarptri,pentarindi,
     & niwi,pentawptri,pentawindi,kgcrsvali,
     & sig,coori,uvalcoe,
     & nad,nadto,mpiad,mpiadlist,mpl,duplixyz)

      if (im == 0) write(6, *) "compAmats"
      call debug_compAmats
     &(im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsvals,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalos,
     & sigs,cooris,uvalcoes,
     & nad,nadto,mpiad,mpiadlist,mpl,duplixyzs)

      if (im == 0) write(6, *) "compAmat4s"
      call debug_compAmats
     &(im,n4,ni4,
     & nir4,pentarptr4,pentarind4,
     & niw4,pentawptr4,pentawind4,kgcrsval4s,
     & niro4,pentarptro4,pentarindo4,
     & niwo4,pentawptro4,pentawindo4,kgcrsvalo4s,
     & sigs,coori4s,uvalcoes,
     & nad4,nadto4,mpiad4,mpiadlist4,mpl4,duplixyz4s)

      end subroutine debug_compAmat
c_______________________________________________________________________
      subroutine debug_compAmatd
     &(im,n,ni,ne,kd,num,cny,coor,younglst,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & niri,pentarptri,pentarindi,
     & niwi,pentawptri,pentawindi,kgcrsvali,
     & sig,coori,uvalcoe,
     & nad,nadto,mpiad,mpiadlist,mpl,duplixyz)
      implicit none
      integer, intent(in) :: im,n,ni,ne,kd
      integer, intent(in) :: num(ne),cny(10,ne)
      real(8), intent(in) :: coor(3,n),younglst(2,kd)
      integer, intent(in) :: nsep,nsepnmap(nsep)
      integer, intent(in) :: cnysep_read(10,ne),cnysep_write(10,ne)
      integer, intent(in) :: nsep_ne_ptr(nsta+2),numsep(ne)
      integer, intent(in) :: nsep_n_ptr(nsta+2)
      integer, intent(in) :: np,ncolor,color_ind(np*MAXCOLOR)
      integer, intent(in) :: nir,niw
      integer, intent(in) :: pentawptr(2*nsta+3),pentarptr(niw+1)
      integer, intent(in) :: pentawind(niw),pentarind(nir)
      real(8), intent(in) :: kgcrsval(9,nir,2)
      integer, intent(in) :: niro,niwo
      integer, intent(in) :: pentawptro(2*nsta+3),pentarptro(niw+1)
      integer, intent(in) :: pentawindo(niw),pentarindo(nir)
      real(8), intent(in) :: kgcrsvalo(9,nir,2)
      integer, intent(in) :: niri,niwi
      integer, intent(in) :: pentawptri(2*nsta+3),pentarptri(niwi+1)
      integer, intent(in) :: pentawindi(niwi),pentarindi(niri)
      real(8), intent(in) :: kgcrsvali(9,niri,2)
      real(8), intent(in) :: sig(nkl)
      real(8), intent(in) :: coori(3,n+ni)
      real(8), intent(in) :: uvalcoe(npc)
      integer, intent(in) :: nad,nadto
      integer, intent(in) :: mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real(8), intent(in) :: duplixyz(3*(n+ni))

      real(8) :: p(npc,3,n+ni),q(npc,3,n+ni),Aq(npc,3,n+ni)
      real(8) :: q2(npc,3,n+ni),diff(npc,3,n+ni),Aq2_norm(npc,3,n+ni)
      real(8) :: mpibuf(2*3*npc*(nadto+1)),bcarray(n+ni)
      integer :: i,j,ista,ipc,i1,ierr
      real(8) :: p_norm,q_norm,Aq_norm,x,y,z,q_diff
      integer :: npc_
      if(npc==21)then
        npc_ = 21
      elseif(npc==288)then
        npc_ =  286
      endif
      bcarray = 1d0

      p = 0d0
      do i = 1, n+ni
      do j = 1, 3
      do ipc = 1, npc
        x = coori(1,i)
        y = coori(2,i)
        z = coori(3,i)
        p(ipc,j,i) = sin(x*0.11 + y*0.12 + z*0.13+j+ipc*3)
      enddo
      enddo
      enddo

      ! normal matrix-vector multiplication ------------------------------
      ! q = Ap
      q = 0d0
      call compAmatd_CRS
     &(im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & sig,uvalcoe,p,q)
#ifndef USEIEPENTA
      call applydirichletd(n+ni,q,bcarray)
#endif
      call sync_parallel_block_d
     &  (n+ni,n+ni,q,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )

      call inner_product_d_block_sync(im,3*(n+ni),p,p,duplixyz,p_norm)
      call inner_product_d_block_sync(im,3*(n+ni),q,q,duplixyz,q_norm)

      if(im == 0)then
        write(6, *) "p_norm:", p_norm
        write(6, *) "q_norm:", q_norm
      endif
      ! compAmatd_with_sync with use_crs
      call compAmatd_with_sync
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,bcarray,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & .false.,.true.,
     & p,q2
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
      call check_diff_d(im,n+ni,p,q,q2,duplixyz,
     & "compAmatd_with_sync with use_crs")

      ! compAmatd_with_sync with overlap_comm, use_crs
      call compAmatd_with_sync
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,bcarray,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & .true.,.true.,
     & p,q2
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
      call check_diff_d(im,n+ni,p,q,q2,duplixyz,
     & "compAmatd_with_sync with overlap_comm,use_crs")

      ! compAmatd_with_sync with EBE
      call compAmatd_with_sync
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,bcarray,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & niri,pentarptri,pentarindi,
     & niwi,pentawptri,pentawindi,kgcrsvali,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & .false.,.false.,
     & p,q2
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
      call check_diff_d(im,n+ni,p,q,q2,duplixyz,
     & "compAmatd_with_sync with EBE")
      end subroutine debug_compAmatd

      subroutine check_diff_d(im,n,p,q,q2,duplixyz,name)
        implicit none
        integer, intent(in) :: im,n
        real(8), intent(in) :: p(npc,3,n),q(npc,3,n),q2(npc,3,n)
        real(8), intent(in) :: duplixyz(n)
        character(len=*), intent(in) :: name
        real(8) :: diff(npc,3,n), p_norm,q_norm,q_diff
        integer :: ipc,j1,i
        call inner_product_d_block_sync(im,3*n,p,p,duplixyz,p_norm)
        call inner_product_d_block_sync(im,3*n,q2,q2,duplixyz,q_norm)
        diff = q - q2
        call inner_product_d_block_sync
     &  (im,3*n,diff,diff,duplixyz,q_diff)
        if(im == 0)then
          write(6, *) name
          write(6, *) "  p_norm:", p_norm
          write(6, *) "  q_norm:", q_norm
          write(6, *) "  q_diff:", q_diff
        endif
      end subroutine check_diff_d

c_______________________________________________________________________
      subroutine debug_compAmats
     &(im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & sig,coori,uvalcoe,
     & nad,nadto,mpiad,mpiadlist,mpl,duplixyz)
      use mpi
      implicit none
      integer, intent(in) :: im,n,ni
      integer, intent(in) :: nir,niw
      integer, intent(in) :: pentawptr(2*nsta+3),pentarptr(niw+1)
      integer, intent(in) :: pentawind(niw),pentarind(nir)
      REAL_4, intent(in) :: kgcrsval(9,nir,2)
      integer, intent(in) :: niro,niwo
      integer, intent(in) :: pentawptro(2*nsta+3),pentarptro(niw+1)
      integer, intent(in) :: pentawindo(niw),pentarindo(nir)
      REAL_4, intent(in) :: kgcrsvalo(9,nir,2)
      REAL_4, intent(in) :: sig(nkl)
      REAL_4, intent(in) :: coori(3,n+ni)
      REAL_4, intent(in) :: uvalcoe(npc)
      integer, intent(in) :: nad,nadto
      integer, intent(in) :: mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      REAL_4, intent(in) :: duplixyz(3*(n+ni))

      REAL_4 :: p(npc,3,n+ni),q(npc,3,n+ni),Aq(npc,3,n+ni)
      REAL_4 :: q2(npc,3,n+ni),diff(npc,3,n+ni),Aq2_norm(npc,3,n+ni)
      REAL_4 :: mpibuf(2*3*npc*(nadto+1))
      REAL_4 :: bcarray(n+ni)
      integer :: i,j,ista,ipc,i1,ierr
      REAL_4 :: p_norm,q_norm,Aq_norm,x,y,z,q_diff
      integer :: npc_
      integer :: mpireq(nad,2)
      if(npc==21)then
        npc_ = 21
      elseif(npc==288)then
        npc_ =  286
      endif

      p = 0d0
      do i = 1, n+ni
      do j = 1, 3
      do ipc = 1, npc
        x = coori(1,i)
        y = coori(2,i)
        z = coori(3,i)
        p(ipc,j,i) = sin(x*0.11 + y*0.12 + z*0.13+j+ipc*3)*1d-10
      enddo
      enddo
      enddo

      ! normal matrix-vector multiplication ------------------------------
      ! q = Ap
      q = 0d0
      call compAmats_CRS
     &(im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & sig,uvalcoe,p,q)

      call sync_parallel_block_s
     &  (n+ni,n+ni,q,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )

      call inner_product_s_block_sync(im,3*(n+ni),p,p,duplixyz,p_norm)
      call inner_product_s_block_sync(im,3*(n+ni),q,q,duplixyz,q_norm)

      if(im == 0)then
        write(6, *) "normal"
        write(6, *) "  p_norm:", p_norm
        write(6, *) "  q_norm:", q_norm
      endif

      ! normal matrix-vector multiplication using compAmats_with_sync
      bcarray = 1.
      call compAmats_with_sync
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & sig,uvalcoe,bcarray,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & p,q2,.false.,.false.,mpireq
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )

      call inner_product_s_block_sync(im,3*(n+ni),p,p,duplixyz,p_norm)
      call inner_product_s_block_sync(im,3*(n+ni),q2,q2,duplixyz,q_norm)
      diff = q - q2
      call inner_product_s_block_sync
     & (im,3*(n+ni),diff,diff,duplixyz,q_diff)

      if(im == 0)then
        write(6, *) "compAmats_with_sync (overlap_comm=F)"
        write(6, *) "  p_norm:", p_norm
        write(6, *) "  q_norm:", q_norm
        write(6, *) "  q_diff:", q_diff
      endif

      ! matrix-vector multiplication for overlapping communication and computation
      call compAmats_with_sync
     & (im,n,ni,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & sig,uvalcoe,bcarray,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & p,q2,.true.,.false.,mpireq
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )

      call inner_product_s_block_sync(im,3*(n+ni),p,p,duplixyz,p_norm)
      call inner_product_s_block_sync(im,3*(n+ni),q2,q2,duplixyz,q_norm)
      diff = q - q2
      call inner_product_s_block_sync
     & (im,3*(n+ni),diff,diff,duplixyz,q_diff)

      if(im == 0)then
        write(6, *) "compAmats_with_sync"
        write(6, *) "  (overlap_comm=T, persistent_comm=F)"
        write(6, *) "  p_norm:", p_norm
        write(6, *) "  q_norm:", q_norm
        write(6, *) "  q_diff:", q_diff
      endif

      ! matrix-vector multiplication for overlapping communication and computation
      ! with persistent communication
      call sync_parallel_block_s_persistent_comm_init
     & (nad,nadto,mpiad,mpl,im,mpibuf,mpireq)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ! 1 回目
      call compAmats_with_sync
     & (im,n,ni,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & sig,uvalcoe,bcarray,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & p,q2,.true.,.true.,mpireq
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )

      call inner_product_s_block_sync(im,3*(n+ni),p,p,duplixyz,p_norm)
      call inner_product_s_block_sync(im,3*(n+ni),q2,q2,duplixyz,q_norm)
      diff = q - q2
      call inner_product_s_block_sync
     & (im,3*(n+ni),diff,diff,duplixyz,q_diff)

      if(im == 0)then
        write(6, *) "compAmats_with_sync"
        write(6, *) "  (overlap_comm=T, persistent_comm=T) 1"
        write(6, *) "  p_norm:", p_norm
        write(6, *) "  q_norm:", q_norm
        write(6, *) "  q_diff:", q_diff
      endif

      ! 2回目
      call compAmats_with_sync
     & (im,n,ni,
     & niro,pentarptro,pentarindo,
     & niwo,pentawptro,pentawindo,kgcrsvalo,
     & sig,uvalcoe,bcarray,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & p,q2,.true.,.true.,mpireq
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )

      call inner_product_s_block_sync(im,3*(n+ni),p,p,duplixyz,p_norm)
      call inner_product_s_block_sync(im,3*(n+ni),q2,q2,duplixyz,q_norm)
      diff = q - q2
      call inner_product_s_block_sync
     & (im,3*(n+ni),diff,diff,duplixyz,q_diff)

      if(im == 0)then
        write(6, *) "compAmats_with_sync"
        write(6, *) "  (overlap_comm=T, persistent_comm=T) 2"
        write(6, *) "  p_norm:", p_norm
        write(6, *) "  q_norm:", q_norm
        write(6, *) "  q_diff:", q_diff
      endif
      call sync_parallel_persistent_comm_free(nad,mpireq)
      end subroutine debug_compAmats

c_______________________________________________________________________