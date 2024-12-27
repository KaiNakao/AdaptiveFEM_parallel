c_______________________________________________________________________
      subroutine compAmats
     &(n,ni,ne,nei,kd,num,cny,coor,younglst,
     & cnyi,p,q,sig,numi,nfla,dsi,itype,
     & im,nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,qsep,
     & np,ncolor,color_ind)
      implicit none
      integer n,ni,ne,nei,kd,itype,im
      integer np,ncolor,color_ind(np*MAXCOLOR)
      REAL_4 q(3*(n+ni)*npc),p(3*(n+ni)*npc)
      integer cny(10,ne),num(ne),cnyi(12,nei)
      REAL_4 coor(3,(n+ni)),younglst(2,kd)
      REAL_4 young,rnyu,gxy,gyz,gxz,detjmat
      REAL_4 xx(4,3)
      REAL_4 arayi(30),arayo(30),keit(36,36)
      REAL_4 sig(5),arayiie(36),arayoie(36)
      REAL_4 bufi(21*3*12),bufo(21*3*12),bufot(36,21)
      REAL_4 cori(3,3),dsi
      REAL_4 x2,y2,z2,x3,y3,z3,x4,y4,z4
      integer nfla(nei),numi(nei)
      integer ie,in,i,i1,j1,i2,k,ii,ipc,ikd,iu,icolor

      integer nir,niw
      integer pentawptr(nsta+2),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)
      REAL_4 pentacrsval(9,nir)

      integer cnysep_read(10,ne),cnysep_write(10,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep,nsepnmap(nsep)
      REAL_4 qsep(3*nsep*npc)
      REAL_4 tmpq(npc),tmpqq(npc)

c
      q=0.
      qsep=0.
c
      do icolor=1,ncolor
!$OMP PARALLEL DO default(none),
!$OMP& shared(np,icolor,color_ind,numsep,
!$OMP&  younglst,coor,cnysep_read,p,cnysep_write,qsep),
!$OMP& private(iu,ie,in,young,rnyu,gxy,gyz,gxz,i1,j1,xx,i2,bufi,ipc,
!$OMP&  arayi,arayo,detjmat)
      do iu=1,np
      do ie=color_ind((icolor-1)*np+iu)+1,color_ind((icolor-1)*np+iu+1)

      in=numsep(ie)
      young=younglst(1,in)
      rnyu=younglst(2,in)
c
      gxy=young/(2.*(1+rnyu))
      gyz=young/(2.*(1+rnyu))
      gxz=young/(2.*(1+rnyu))
c
      do i1=1,4
      do j1=1,3
      xx(i1,j1)=coor(j1,cnysep_read(i1,ie))
      enddo
      enddo

c----------------------------------------------
      do i1=1,10
      do j1=1,63
      i2=63*(i1-1)+j1
      bufi(i2)=p(63*(cnysep_read(i1,ie)-1)+j1)
      enddo
      enddo
      do ipc=1,npc
      do i1=1,10
      do j1=1,3
      arayi(3*(i1-1)+j1)=bufi(63*(i1-1)+npc*(j1-1)+ipc)
     &                  -bufi(          npc*(j1-1)+ipc)
      enddo
      enddo
      arayo=0.
      call nontet10kus(young,gxy,gyz,gxz,rnyu,xx,arayi,arayo,detjmat)
      do i1=1,10
      do j1=1,3
       qsep(63*(cnysep_write(i1,ie)-1)+(j1-1)*npc+ipc)=
     & qsep(63*(cnysep_write(i1,ie)-1)+(j1-1)*npc+ipc)+
     & arayo(3*(i1-1)+j1)
      enddo
      enddo
      enddo

      enddo ! ie
      enddo ! iu
!$OMP END PARALLEL DO
      enddo ! icolor

      do ikd=1,nsta+1
!$OMP PARALLEL DO default(none),
!$OMP& shared(nsep_n_ptr,ikd,qsep,sig,q,nsepnmap),
!$OMP& private(i,i1,ipc,tmpq,tmpqq)
      do i=nsep_n_ptr(ikd)+1,nsep_n_ptr(ikd+1)
      do i1=1,3
      do ipc=1,npc
      tmpq(ipc)=qsep(3*npc*(i-1)+npc*(i1-1)+ipc)
      enddo
      call compsigs_npc21(ikd,sig,tmpq,tmpqq)
      do ipc=1,npc
       q(3*npc*(nsepnmap(i)-1)+(i1-1)*npc+ipc)=
     & q(3*npc*(nsepnmap(i)-1)+(i1-1)*npc+ipc)+tmpqq(ipc)
      enddo ! ipc
      enddo ! i1
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ikd

#ifdef USEIEPENTA
      call compAmats_iepenta
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,sig,p,q)
#endif
      end

      subroutine compAmats_with_sync
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,bcarray,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & p,q,
     & overlap_comm,persistent_comm,mpireq
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      ! pentawptr は overlap_comm == .false. のときは nsta+2 までしか使わない
      integer,intent(in) :: pentawptr(2*nsta+3),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      REAL_4,intent(in) :: sig(nkl),uvalcoe(npc),p(npc,3,n+ni)
      REAL_4,intent(in) :: bcarray(n+ni)
      integer,intent(in) :: nad,nadto
      integer,intent(in) :: mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
c      REAL_4,intent(inout) :: mpibuf(2*3*npc*(nadto+1))
      integer*2,intent(inout)::
     $     mpibuf(2*(3*npc*(nadto+1)+12*npc*(nadto/DIV_N+nad)))

#define DEBUGKF
#ifdef DEBUGKF
      REAL_4
     $     mpibuf_f(2*(3*npc*(nadto+1)+12*npc*(nadto/DIV_N+nad)))
#endif
      logical,intent(in) :: overlap_comm,persistent_comm
      integer,intent(inout) :: mpireq(nad,2)
#if npc == 1 || npc == 21
      REAL_4,intent(in) :: pentacrsval(9,nir)
#elif npc == 288
      REAL_4,intent(in) :: pentacrsval(9,nir,2)
#endif
      REAL_4,intent(out) :: q(npc,3,n+ni)
      integer :: ierr
#ifdef MEASURE_TIME_BARRIER
      integer, intent(in) :: ind,mp
#include "timerbuf.F"
#endif
      q = 0d0
      call compAmats_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)

#ifdef DEBUGKF
!      write(*,*) 'DEBUGPRINT sync_parallel_block_s'
      call sync_parallel_block_s
     &  (n+ni,n+ni,q,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf_f
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
#else
      call sync_parallel_block_s_part1
     &  (n+ni,n+ni,q,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf,
     &  mpireq,persistent_comm
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
#endif

      if(overlap_comm)then
      call compAmats_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr(nsta+2),pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
      endif

#ifndef DEBUGKF
      call sync_parallel_block_s_part2
     &  (n+ni,n+ni,q,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf,
     &  mpireq
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
#endif

#ifndef USEIEPENTA
      call applydirichlets(n+ni,q,bcarray)
#endif
      end subroutine compAmats_with_sync

      subroutine compAmats_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      REAL_4,intent(in) :: sig(nkl),uvalcoe(npc),p(npc,3,n+ni)
      REAL_4,intent(inout) :: q(npc,3,n+ni)
#if npc == 1
      REAL_4,intent(in) :: pentacrsval(9,nir)
      call compAmats_npc1_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawind,pentacrsval,
     & p,q)
#elif npc == 21
      REAL_4,intent(in) :: pentacrsval(9,nir)
      call compAmats_npc21_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,p,q)
#elif npc == 288
      REAL_4,intent(in) :: pentacrsval(9,nir,2)
      call compAmats_npc288_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
#else
ERROR: invalid npc
#endif
      end subroutine compAmats_CRS

c_______________________________________________________________________
      subroutine compAmats_npc1_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawind,pentacrsval,
     & p,q)
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      REAL_4,intent(in) :: pentacrsval(9,nir)
      REAL_4,intent(in) :: p(npc,3,n+ni)
      REAL_4,intent(out) :: q(npc,3,n+ni)
      integer i,i1,j,i2,ipc
      REAL_4 tmpv(npc,3)

      !q=0.0
!$OMP PARALLEL DO default(none),
!$OMP& shared(niw,pentawind,pentacrsval,
!$OMP&  pentarptr,pentarind,p,q),
!$OMP& private(i,i1,tmpv,j,i2,ipc)
      do i=1,niw
      i1=pentawind(i)
      tmpv=0.0
      do j=pentarptr(i)+1,pentarptr(i+1)
      i2=pentarind(j)
      do ipc=1,npc
      tmpv(ipc,1)=tmpv(ipc,1)
     & + pentacrsval(1,j) * p(ipc,1,i2)
     & + pentacrsval(2,j) * p(ipc,2,i2)
     & + pentacrsval(3,j) * p(ipc,3,i2)
      tmpv(ipc,2)=tmpv(ipc,2)
     & + pentacrsval(4,j) * p(ipc,1,i2)
     & + pentacrsval(5,j) * p(ipc,2,i2)
     & + pentacrsval(6,j) * p(ipc,3,i2)
      tmpv(ipc,3)=tmpv(ipc,3)
     & + pentacrsval(7,j) * p(ipc,1,i2)
     & + pentacrsval(8,j) * p(ipc,2,i2)
     & + pentacrsval(9,j) * p(ipc,3,i2)
      enddo ! ipc
      enddo ! j

      do j=1,3
      do ipc=1,npc
        q(ipc,j,i1) = q(ipc,j,i1) + tmpv(ipc,j)
      enddo ! ipc
      enddo ! j
      enddo ! i
!$OMP END PARALLEL DO
      end

c_______________________________________________________________________
      subroutine compAmats_npc21_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,p,q)
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      REAL_4,intent(in) :: pentacrsval(9,nir)
      REAL_4,intent(in) :: sig(nsta),p(npc,3,n+ni)
      REAL_4,intent(inout) :: q(npc,3,n+ni)
      integer ista,i,i1,j,i2,ipc
      REAL_4 tmpq(npc),tmpqq(npc),tmpv(npc,3)

      !q=0.
      do ista=1,nsta+1
!$OMP PARALLEL DO default(none),
!$OMP& shared(pentawptr,ista,pentawind,sig,pentacrsval,
!$OMP&  pentarptr,pentarind,p,q),
!$OMP& private(i,i1,tmpv,j,i2,ipc,tmpq,tmpqq)
      do i=pentawptr(ista)+1,pentawptr(ista+1)
      i1=pentawind(i)
      tmpv=0.0
      do j=pentarptr(i)+1,pentarptr(i+1)
      i2=pentarind(j)
      do ipc=1,npc
      tmpv(ipc,1)=tmpv(ipc,1)
     & + pentacrsval(1,j) * p(ipc,1,i2)
     & + pentacrsval(2,j) * p(ipc,2,i2)
     & + pentacrsval(3,j) * p(ipc,3,i2)
      tmpv(ipc,2)=tmpv(ipc,2)
     & + pentacrsval(4,j) * p(ipc,1,i2)
     & + pentacrsval(5,j) * p(ipc,2,i2)
     & + pentacrsval(6,j) * p(ipc,3,i2)
      tmpv(ipc,3)=tmpv(ipc,3)
     & + pentacrsval(7,j) * p(ipc,1,i2)
     & + pentacrsval(8,j) * p(ipc,2,i2)
     & + pentacrsval(9,j) * p(ipc,3,i2)
      enddo ! ipc
      enddo ! j

      do j=1,3
      do ipc=1,npc
      tmpq(ipc)=tmpv(ipc,j)
      enddo
      call compsigs_npc21(ista,sig,tmpq,tmpqq)
      do ipc=1,npc
        q(ipc,j,i1) = q(ipc,j,i1) + tmpqq(ipc)
      enddo ! ipc
      enddo ! j
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ista
      end

c_______________________________________________________________________
      subroutine compAmats_npc288_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & sig,uvalcoe,p,q)
      implicit none
#ifdef _OPENMP
      include 'omp_lib.h'
#endif
#ifdef USE_PA
      include 'mpif.h'
      integer ierr
#endif
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      REAL_4,intent(in) :: kgcrsval(9,nir,2)
      REAL_4,intent(in) :: sig(nsta*2),uvalcoe(npc)
      REAL_4,intent(in) :: p(npc,3,n+ni)
      REAL_4,intent(inout) :: q(npc,3,n+ni)

      integer :: ista,i,i1,j,i2,ipc,ikg,iu
      REAL_4 :: tmpv1(npc,3),tmpv2(npc,3)
      REAL_4 :: diag(npc)
      integer nb, nn, tn, nt
      integer is(0:12), ie(0:12)
      integer,parameter:: pf_l1=1
      integer,parameter:: pf_l2=6

#ifdef USE_PA
      call fapp_start('compAmats_npc288_CRS_b1', 0, 0)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call fapp_stop('compAmats_npc288_CRS_b1', 0, 0)
      call fapp_start('compAmats_npc288_CRS', 0, 0)
#endif
c      q=0.
c      diag = 1./uvalcoe

!$OMP PARALLEL default(none),
!$OMP& shared(pentawptr,pentawind,sig,kgcrsval,
!$OMP&  pentarptr,pentarind,p,q,diag,uvalcoe,n,ni),
!$OMP& private(ista,iu,i,i1,tmpv1,tmpv2,j,i2,ipc,
!$OMP&  nb,nn,tn,nt,is,ie)

!$OMP DO
      do ipc=1,npc
         diag(ipc) = 1.0/uvalcoe(ipc)
      enddo
!$OMP ENDDO
         
      do ista=1,nsta+1
         nb = 0
         do i=pentawptr(ista)+1,pentawptr(ista+1)
            nb = nb+pentarptr(i+1)-pentarptr(i)
         enddo
#ifdef _OPENMP
         nt = omp_get_num_threads()
         tn = omp_get_thread_num()
#else
         nt = 1
         tn = 0
#endif
         nb = (nb-1)/nt+1
         nn = 0
         j = 0
         is(:) = pentawptr(ista+1)+1
         ie(:) = 0
         is(j) = pentawptr(ista)+1
         do i=pentawptr(ista)+1,pentawptr(ista+1)
            nn = nn+pentarptr(i+1)-pentarptr(i)
            if(nn .ge. nb) then
               is(j+1) = i+1
               ie(j) = i
               nn = 0
               j = j+1
            endif
         enddo
         ie(j) = pentawptr(ista+1)
         
c!$OMP DO
c      do i=pentawptr(ista)+1,pentawptr(ista+1)
!ocl noprefetch
      do i=is(tn), ie(tn)
        i1=pentawind(i)
c        tmpqq(:,:) = 0.
c        do ikg = 1, 2 ! ikg = 1 for K, 2 for G
          tmpv1(:,:)=0.
          tmpv2(:,:)=0.
          do j=pentarptr(i)+1,pentarptr(i+1)
            i2=pentarind(j)
!ocl prefetch_read(p(1,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(1,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(1,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(65,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(65,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(65,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(129,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(129,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(129,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(193,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(193,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(193,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(257,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(257,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(257,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(1,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(1,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(1,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(65,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(65,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(65,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(129,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(129,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(129,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(193,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(193,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(193,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(257,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(257,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(257,3,pentarind(j+pf_l2)),level=2,strong=1)
            do ipc=1,npc
              tmpv1(ipc,1)=tmpv1(ipc,1)
     &        + kgcrsval(1,j,1) * p(ipc,1,i2)
     &        + kgcrsval(2,j,1) * p(ipc,2,i2)
     &        + kgcrsval(3,j,1) * p(ipc,3,i2)
              tmpv1(ipc,2)=tmpv1(ipc,2)
     &        + kgcrsval(4,j,1) * p(ipc,1,i2)
     &        + kgcrsval(5,j,1) * p(ipc,2,i2)
     &        + kgcrsval(6,j,1) * p(ipc,3,i2)
              tmpv1(ipc,3)=tmpv1(ipc,3)
     &        + kgcrsval(7,j,1) * p(ipc,1,i2)
     &        + kgcrsval(8,j,1) * p(ipc,2,i2)
     &        + kgcrsval(9,j,1) * p(ipc,3,i2)
              tmpv2(ipc,1)=tmpv2(ipc,1)
     &        + kgcrsval(1,j,2) * p(ipc,1,i2)
     &        + kgcrsval(2,j,2) * p(ipc,2,i2)
     &        + kgcrsval(3,j,2) * p(ipc,3,i2)
              tmpv2(ipc,2)=tmpv2(ipc,2)
     &        + kgcrsval(4,j,2) * p(ipc,1,i2)
     &        + kgcrsval(5,j,2) * p(ipc,2,i2)
     &        + kgcrsval(6,j,2) * p(ipc,3,i2)
              tmpv2(ipc,3)=tmpv2(ipc,3)
     &        + kgcrsval(7,j,2) * p(ipc,1,i2)
     &        + kgcrsval(8,j,2) * p(ipc,2,i2)
     &        + kgcrsval(9,j,2) * p(ipc,3,i2)
            enddo ! ipc
          enddo ! j

          if(ista <= nsta) then
             do j=1, 3
                iu = ista
                call compsigs_npc288(iu,sig,q(1,j,i1),tmpv1(1,j),diag)
                iu = nsta+ista
                call compsigs_npc288(iu,sig,q(1,j,i1),tmpv2(1,j),diag)
             enddo
          else
             do j=1,3
                do ipc=1,npc
                   q(ipc,j,i1) = q(ipc,j,i1)
     $                  + (tmpv1(ipc,j)+tmpv2(ipc,j))*diag(j)
                enddo ! ipc
             enddo ! j
          endif
      enddo ! i
!$OMP BARRIER
      enddo ! ista
!$OMP END PARALLEL
#ifdef USE_PA
      call fapp_stop('compAmats_npc288_CRS', 0, 0)
#endif
      end
