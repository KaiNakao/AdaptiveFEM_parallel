c_______________________________________________________________________
      subroutine compAmatd_npc21_EBE
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
      REAL*8 q(3*(n+ni)*npc),p(3*(n+ni)*npc)
      integer cny(10,ne),num(ne),cnyi(12,nei)
      REAL*8 coor(3,(n+ni)),younglst(2,kd)
      REAL*8 young,rnyu,gxy,gyz,gxz,detjmat
      REAL*8 xx(4,3)
      REAL*8 arayi(30),arayo(30),keit(36,36)
      REAL*8 sig(5),arayiie(36),arayoie(36)
      REAL*8 bufi(21*3*12),bufo(21*3*12),bufot(36,21)
      REAL*8 cori(3,3),dsi
      REAL*8 x2,y2,z2,x3,y3,z3,x4,y4,z4
      integer nfla(nei),numi(nei)
      integer ie,in,i,i1,j1,i2,k,ii,ipc,ikd,iu,icolor

      integer nir,niw
      integer pentawptr(nsta+2),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)
      REAL*8 pentacrsval(9,nir)

      integer cnysep_read(10,ne),cnysep_write(10,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep,nsepnmap(nsep)
      REAL*8 qsep(3*nsep*npc)
      REAL*8 tmpq(npc),tmpqq(npc)

c
!$OMP PARALLEL DO default(none),shared(q,n,ni)
      do i=1,3*npc*(n+ni)
      q(i)=0.
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO default(none),shared(qsep,nsep)
      do i=1,3*npc*nsep
      qsep(i)=0.
      enddo
!$OMP END PARALLEL DO
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
      call nontet10kud(young,gxy,gyz,gxz,rnyu,xx,arayi,arayo,detjmat)
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
      call compsigd_npc21(ikd,sig,tmpq,tmpqq)
      do ipc=1,npc
       q(3*npc*(nsepnmap(i)-1)+(i1-1)*npc+ipc)=
     & q(3*npc*(nsepnmap(i)-1)+(i1-1)*npc+ipc)+tmpqq(ipc)
      enddo ! ipc
      enddo ! i1
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ikd

#ifdef USEIEPENTA
      call compAmatd_iepenta
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,sig,p,q)
#endif
      end subroutine compAmatd_npc21_EBE

      subroutine compAmatd_npc288_EBE
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & p,q)
      implicit none
      integer,intent(in) :: im,n,ni,ne,kd
      integer,intent(in) :: num(ne),cny(10,ne)
      real(8),intent(in) :: coor(3,n),younglst(2,kd)
      real(8),intent(in) :: sig(nkl),uvalcoe(npc)
      integer,intent(in) :: nsep,nsepnmap(nsep)
      integer,intent(in) :: cnysep_read(10,ne),cnysep_write(10,ne)
      integer,intent(in) :: nsep_ne_ptr(nsta+2),numsep(ne)
      integer,intent(in) :: nsep_n_ptr(nsta+2)
      integer,intent(in) :: np,ncolor,color_ind(np*MAXCOLOR)
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: pentacrsval(9,nir,2)
      real(8),intent(in) :: p(npc,3,n+ni)
      real(8),intent(out) :: q(npc,3,n+ni)
      integer,parameter :: nkg = 2 ! Kp, Gp を別々に計算
      real(8) :: qsep(npc,3,nsep,nkg)
      real(8) :: young,rnyu,gxy,gyz,gxz,detjmat,rKc,rGc,qtmp
      real(8) :: xx(4,3)
      real(8) :: arayi(3,10),arayo(3,10,nkg)
      real(8) :: bufi(npc,3,12),diag(npc)
      integer :: ie,in,i,i1,j1,i2,k,ii,ipc,ista,iu,icolor,ikg
      real(8) :: tmpq(3,npc),tmpqq(3,npc)

      diag = 1d0/uvalcoe
c
!$OMP PARALLEL DO default(none),shared(q,n,ni)
      do i=1,n+ni
      do j1=1,3
      do ipc=1,npc
        q(ipc,j1,i)=0.d0
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO default(none),shared(qsep,nsep)
      do ikg=1,nkg
      do i=1,nsep
      do j1=1,3
      do ipc=1,npc
        qsep(ipc,j1,i,ikg)=0.d0
      enddo
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
c
      do icolor=1,ncolor
!$OMP PARALLEL DO default(none),
!$OMP& shared(np,icolor,color_ind,numsep,
!$OMP&  younglst,coor,cnysep_read,p,cnysep_write,qsep),
!$OMP& private(iu,ie,in,i,young,rnyu,gxy,gyz,gxz,i1,j1,xx,i2,bufi,ipc,
!$OMP&  arayi,arayo,detjmat,rKc,rGc)
      do iu=1,np
      do ie=color_ind((icolor-1)*np+iu)+1,color_ind((icolor-1)*np+iu+1)

      in=numsep(ie)
      young=younglst(1,in)
      rnyu=younglst(2,in)
c
      rKc=young/3/(1-2*rnyu)
      rGc=young/2/(1+rnyu)
c
      do i1=1,4
      do j1=1,3
        xx(i1,j1)=coor(j1,cnysep_read(i1,ie))
      enddo
      enddo

c----------------------------------------------
      do i1=1,10
      do j1=1,3
      do ipc=1,npc
        bufi(ipc,j1,i1) = p(ipc,j1,cnysep_read(i1,ie))
      enddo
      enddo
      enddo

      do ipc=1,npc
        do i1=1,10
        do j1=1,3
          arayi(j1,i1) = bufi(ipc,j1,i1)
        enddo
        enddo
        arayo=0.d0
        call nontet10kud_KG
     &   (rKc,0d0,xx,arayi,arayo(:,:,1),detjmat)
        call nontet10kud_KG
     &   (0d0,rGc,xx,arayi,arayo(:,:,2),detjmat)
        do ikg=1,nkg
        do i1=1,10
        do j1=1,3
          i = cnysep_write(i1,ie)
          qsep(ipc,j1,i,ikg) = qsep(ipc,j1,i,ikg) + arayo(j1,i1,ikg)
        enddo
        enddo
        enddo
      enddo !ipc

      enddo ! ie
      enddo ! iu
!$OMP END PARALLEL DO
      enddo ! icolor

      do ista=1,nsta+1
!$OMP PARALLEL DO default(none),
!$OMP& shared(nsep_n_ptr,ista,qsep,sig,q,nsepnmap,diag),
!$OMP& private(i,i1,iu,ipc,tmpq,tmpqq,qtmp)
      do i1=nsep_n_ptr(ista)+1,nsep_n_ptr(ista+1)
        tmpqq(:,:) = 0d0
        do ikg=1,nkg
          do j1=1,3
          do ipc=1,npc
            qtmp = qsep(ipc,j1,i1,ikg)
            tmpqq(j1,ipc) = tmpqq(j1,ipc) + qtmp*diag(ipc)
            tmpq(j1,ipc) = qtmp
          enddo
          enddo
          if(ista <= nsta)then
            iu = nsta*(ikg-1)+ista
            call compsigd_npc288(iu,sig,tmpqq,tmpq)
          endif
        enddo ! ikg
        do j1=1,3
        do ipc=1,npc
          i = nsepnmap(i1)
          q(ipc,j1,i) = q(ipc,j1,i) + tmpqq(j1,ipc)
        enddo ! ipc
        enddo
      enddo ! i1
!$OMP END PARALLEL DO
      enddo ! ista

#ifdef USEIEPENTA
      call compAmatd_npc288_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
#endif
      end subroutine compAmatd_npc288_EBE

      subroutine compAmatd_with_sync
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,bcarray,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & overlap_comm,use_crs,
     & p,q
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
      ! overlap_comm = .true. のときでも通信のオーバーラップはされない.
      ! 単精度の行列ベクトル積でoverlap_comm=.true.を使い、
      ! 倍精度の行列ベクトル積でも同じCRSを使うとき用の関数
      implicit none
      integer,intent(in) :: im,n,ni,ne,kd
      integer,intent(in) :: num(ne),cny(10,ne)
      real(8),intent(in) :: coor(3,n),younglst(2,kd)
      real(8),intent(in) :: sig(nkl),uvalcoe(npc),bcarray(n+ni)
      integer,intent(in) :: nsep,nsepnmap(nsep)
      integer,intent(in) :: cnysep_read(10,ne),cnysep_write(10,ne)
      integer,intent(in) :: nsep_ne_ptr(nsta+2),numsep(ne)
      integer,intent(in) :: nsep_n_ptr(nsta+2)
      integer,intent(in) :: np,ncolor,color_ind(np*MAXCOLOR)
      integer,intent(in) :: nir,niw
      ! pentawptr は overlap_comm == .false. のときは nsta+2 までしか使わない
      integer,intent(in) :: pentawptr(2*nsta+3),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: p(npc,3,n+ni)
      integer,intent(in) :: nad,nadto
      integer,intent(in) :: mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real(8),intent(inout) :: mpibuf(2*3*npc*(nadto+1))
      logical,intent(in) :: overlap_comm,use_crs
#if npc == 1 || npc == 21
      real(8),intent(in) :: pentacrsval(9,nir)
#elif npc == 288
      real(8),intent(in) :: pentacrsval(9,nir,2)
#endif
      real(8),intent(out) :: q(npc,3,n+ni)
      integer :: ierr
#ifdef MEASURE_TIME_BARRIER
      integer, intent(in) :: ind,mp
#include "timerbuf.F"
#endif
      q = 0d0
      call compAmatd
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & use_crs,
     & p,q)

      call sync_parallel_block_d
     &  (n+ni,n+ni,q,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )

      if(overlap_comm)then
      if(.not.use_crs)then
        stop "overlap_comm is implemented only for use_crs"
      endif
      call compAmatd
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr(nsta+2),pentawind,pentacrsval,
     & use_crs,
     & p,q)
      endif ! else overlap_comm
#ifndef USEIEPENTA
      call applydirichletd(n+ni,q,bcarray)
#endif
      end subroutine compAmatd_with_sync

      subroutine compAmatd_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: sig(nkl),uvalcoe(npc),p(npc,3,n+ni)
      real(8),intent(out) :: q(npc,3,n+ni)
#if npc == 1
      real(8),intent(in) :: pentacrsval(9,nir)
      call compAmatd_npc1_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawind,pentacrsval,
     & p,q)
#elif npc == 21
      real(8),intent(in) :: pentacrsval(9,nir)
      call compAmatd_npc21_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,p,q)
#elif npc == 288
      real(8),intent(in) :: pentacrsval(9,nir,2)
      call compAmatd_npc288_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
#else
ERROR: invalid npc
#endif
      end subroutine compAmatd_CRS

      subroutine compAmatd
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & use_crs,
     & p,q)
      implicit none
      integer,intent(in) :: im,n,ni,ne,kd
      integer,intent(in) :: num(ne),cny(10,ne)
      real(8),intent(in) :: coor(3,n),younglst(2,kd)
      real(8),intent(in) :: sig(nkl),uvalcoe(npc)
      integer,intent(in) :: nsep,nsepnmap(nsep)
      integer,intent(in) :: cnysep_read(10,ne),cnysep_write(10,ne)
      integer,intent(in) :: nsep_ne_ptr(nsta+2),numsep(ne)
      integer,intent(in) :: nsep_n_ptr(nsta+2)
      integer,intent(in) :: np,ncolor,color_ind(np*MAXCOLOR)
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: p(npc,3,n+ni)
      logical,intent(in) :: use_crs
      real(8),intent(inout) :: q(npc,3,n+ni)
#if npc == 1
      real(8),intent(in) :: pentacrsval(9,nir)
      if(use_crs)then
      call compAmatd_npc1_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawind,pentacrsval,
     & p,q)
      else
        if(im==0)then
          write(6,*) "compAmatd_ncp1_EBE is not implemented"
          stop
        endif
      endif
#elif npc == 21
      real(8),intent(in) :: pentacrsval(9,nir)
      if(use_crs)then
      call compAmatd_npc21_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,p,q)
      else
        if(im==0)then
          write(6,*) "compAmatd_ncp21_EBE is not implemented"
          stop
        endif
      endif
#elif npc == 288
      real(8),intent(in) :: pentacrsval(9,nir,2)
      if(use_crs)then
      call compAmatd_npc288_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,p,q)
      else
      call compAmatd_npc288_EBE
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & p,q)
      endif
#else
ERROR: invalid npc
#endif
      end subroutine compAmatd


c_______________________________________________________________________
      subroutine compAmatd_npc1_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawind,pentacrsval,
     & p,q)
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: pentacrsval(9,nir)
      real(8),intent(in) :: p(npc,3,n+ni)
      real(8),intent(inout) :: q(npc,3,n+ni)
      integer i,i1,j,i2,ipc
      real(8) tmpv(npc,3)

      !q=0d0
!$OMP PARALLEL DO default(none),
!$OMP& shared(niw,pentawind,pentacrsval,
!$OMP&  pentarptr,pentarind,p,q),
!$OMP& private(i,i1,tmpv,j,i2,ipc)
      do i=1,niw
      i1=pentawind(i)
      tmpv=0d0
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
      subroutine compAmatd_npc21_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,p,q)
      implicit none
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: pentacrsval(9,nir)
      real(8),intent(in) :: sig(nsta),p(npc,3,n+ni)
      real(8),intent(out) :: q(npc,3,n+ni)
      integer ista,i,i1,j,i2,ipc
      real(8) tmpq(npc),tmpqq(npc),tmpv(npc,3)

      !q=0d0
      do ista=1,nsta+1
!$OMP PARALLEL DO default(none),
!$OMP& shared(pentawptr,ista,pentawind,sig,pentacrsval,
!$OMP&  pentarptr,pentarind,p,q),
!$OMP& private(i,i1,tmpv,j,i2,ipc,tmpq,tmpqq)
      do i=pentawptr(ista)+1,pentawptr(ista+1)
      i1=pentawind(i)
      tmpv=0d0
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
      call compsigd_npc21(ista,sig,tmpq,tmpqq)
      do ipc=1,npc
        q(ipc,j,i1) = q(ipc,j,i1) + tmpqq(ipc)
      enddo ! ipc
      enddo ! j
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ista
      end

c_______________________________________________________________________
      subroutine compAmatd_npc288_CRS
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,kgcrsval,
     & sig,uvalcoe,p,q)
      implicit none
#ifdef _OPENMP
      include 'omp_lib.h'
#endif
      integer,intent(in) :: im,n,ni
      integer,intent(in) :: nir,niw
      integer,intent(in) :: pentawptr(nsta+2),pentarptr(niw+1)
      integer,intent(in) :: pentawind(niw),pentarind(nir)
      real(8),intent(in) :: kgcrsval(9,nir,2)
      real(8),intent(in) :: sig(nsta*2),uvalcoe(npc)
      real(8),intent(in) :: p(npc,3,n+ni)
      real(8),intent(inout) :: q(npc,3,n+ni)

      integer :: ista,i,i1,j,i2,ipc,ikg,iu
      real(8) :: tmpv1(npc,3),tmpv2(npc,3)
      real(8) :: diag(npc)
      integer nb, nn, tn, nt
      integer is(0:12), ie(0:12)
      integer,parameter:: pf_l1=1
      integer,parameter:: pf_l2=6

      ! call fapp_start('compAmatd_npc288_CRS', 0, 0)
c      q=0d0
c      diag = 1d0/uvalcoe

!$OMP PARALLEL default(none),
!$OMP& shared(pentawptr,pentawind,sig,kgcrsval,
!$OMP&  pentarptr,pentarind,p,q,diag,uvalcoe,n,ni),
!$OMP& private(ista,iu,i,i1,tmpv1,tmpv2,j,i2,ipc,
!$OMP&  nb,nn,tn,nt,is,ie)

!$OMP DO
      do ipc=1,npc
         diag(ipc) = 1d0/uvalcoe(ipc)
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

c!$OMP DO schedule(static,1)
c      do i=pentawptr(ista)+1,pentawptr(ista+1)
!ocl noprefetch
      do i=is(tn), ie(tn)
        i1=pentawind(i)
c        tmpqq(:,:) = 0d0
c        do ikg = 1, 2 ! ikg = 1 for K, 2 for G
          tmpv1(:,:)=0d0
          tmpv2(:,:)=0d0
          do j=pentarptr(i)+1,pentarptr(i+1)
            i2=pentarind(j)
!ocl prefetch_read(p(1,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(1,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(1,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(33,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(33,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(33,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(65,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(65,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(65,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(97,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(97,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(97,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(129,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(129,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(129,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(161,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(161,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(161,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(193,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(193,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(193,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(225,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(225,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(225,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(257,1,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(257,2,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(257,3,pentarind(j+pf_l1)),level=1,strong=1)
!ocl prefetch_read(p(1,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(1,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(1,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(33,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(33,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(33,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(65,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(65,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(65,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(97,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(97,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(97,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(129,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(129,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(129,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(161,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(161,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(161,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(193,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(193,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(193,3,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(225,1,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(225,2,pentarind(j+pf_l2)),level=2,strong=1)
!ocl prefetch_read(p(225,3,pentarind(j+pf_l2)),level=2,strong=1)
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
                call compsigd_npc288(iu,sig,q(1,j,i1),tmpv1(1,j),diag)
                iu = nsta+ista
                call compsigd_npc288(iu,sig,q(1,j,i1),tmpv2(1,j),diag)
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
      enddo ! ista
!$OMP END PARALLEL
      ! call fapp_stop('compAmatd_npc288_CRS', 0, 0)
      end
