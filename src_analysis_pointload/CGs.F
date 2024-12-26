      subroutine CGs
     - (im,np,n,ni,ne,nei,kd,coor,cny,num,cnyi,younglst,b,u,
     & uval,sig,uvalcoe,
     & dsi,numi,nfla,duplixyz,bcarray,eps,itermax,itype,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,icny,icnyi,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,qsep,
     & ncolor,color_ind,output_log,ite,
     & overlap_comm,persistent_comm,mpireq)
      ! itype == 0: compute for K matrix
      ! itype == 1: compute for mapping nodal force to force
      implicit none
#ifdef USE_PA
      include 'mpif.h'
      integer ierr
#endif
      integer n,ni,ne,nei,kd,im,np,itermax,itype,icny,icnyi
      integer ite,i,j,i1,i2,i3,ierr
      REAL_4 coor(3,n+ni),younglst(2,kd)
      integer cny(icny,ne),num(ne),cnyi(icnyi,nei)
      REAL_4 b2,tmp
      REAL_4 alphacg,betacg,r2,rhoa,rhob,pq,zq,err,eps
      REAL_4 b(3*(n+ni)*npc),u(3*(n+ni)*npc),r(3*(n+ni)*npc)
      REAL_4 p(3*(n+ni)*npc),q(3*(n+ni)*npc)
      REAL_4 uval(6,n+ni),z(3*(n+ni)*npc)
      REAL_4 sig(nkl),uvalcoe(npc),coet,dsi
      REAL_4 duplixyz(3*(n+ni)),bcarray(n+ni)
      integer numi(nei),nfla(nei)
      integer nad,nadto
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
c      REAL_4 mpibuf(2*3*npc*(nadto+1))
      integer*2 mpibuf(2*(3*npc*(nadto+1)+12*npc*(nadto/DIV_N+nad)))
c
      integer nir,niw
      integer pentawptr(2*nsta+3),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)

      logical output_log,overlap_comm,persistent_comm
      integer mpireq(nad,2)
#if npc == 1 || npc == 21
      REAL_4 pentacrsval(9,nir)
#elif npc == 288
      REAL_4 pentacrsval(9,nir,2)
#else
ERROR: invalid npc
#endif

c
      integer cnysep_read(icny,ne),cnysep_write(icny,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep,nsepnmap(nsep)
      REAL_4 qsep(3*nsep*npc)

      integer ncolor,color_ind(np*MAXCOLOR)

c
c--------------------------
#ifndef USEIEPENTA
      call applydirichlets(n+ni,b,bcarray)
      call applydirichlets(n+ni,u,bcarray)
#endif
      call inner_product_s_block_sync(im,3*(n+ni),b,b,duplixyz,b2)
!$OMP PARALLEL DO default(none),private(i),shared(n,ni,p)
      do i=1,3*(n+ni)*npc
      p(i)=0.0
      enddo
!$OMP END PARALLEL DO
c___________________________________
      call compAmats_with_sync
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & sig,uvalcoe,bcarray,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & u,r,
     & overlap_comm,persistent_comm,mpireq
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
!$OMP PARALLEL DO default(none),private(i),shared(n,ni,b,r)
      do i=1,3*(n+ni)*npc
      r(i)=b(i)-r(i)
      enddo
!$OMP END PARALLEL DO
c
      call inner_product_s_block_sync(im,3*(n+ni),r,r,duplixyz,r2)
      err=sqrt(abs(r2/b2))
      if(output_log .and. im.eq.0)then
      if(icny.eq.10)then
      write(6,*) 'inner_tet10-loop RE=',err,'init: type',itype
      else
      write(6,*) 'inner_tet4-loop RE=',err,'init: type',itype
      endif
      endif
c
      betacg=0.0
c
      ite=0
      do 101 while (err.gt.eps .and. ite.lt.itermax )
      ite=ite+1
c---------------------------------
c preconditioner
#ifdef USE_PA
      call fapp_start('CGs_OMP_3_b', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('CGs_OMP_3_b', 0, 0)
      call fapp_start('CGs_OMP_3', 0, 0)
#endif
!ocl prefetch_sequential(soft)
!ocl prefetch_strong
!$OMP PARALLEL DO default(none),private(i,j,i1,i2,i3,coet),
!$OMP& shared(n,ni,uvalcoe,z,uval,r)
      do i=1,n+ni
      do j=1,npc
      coet=uvalcoe(j)
      i1=j  + 0*npc  +3*npc*(i-1)
      i2=j  + 1*npc  +3*npc*(i-1)
      i3=j  + 2*npc  +3*npc*(i-1)
      z(i1)=(uval(1,i)*r(i1)+uval(4,i)*r(i2)+uval(5,i)*r(i3))*coet
      z(i2)=(uval(4,i)*r(i1)+uval(2,i)*r(i2)+uval(6,i)*r(i3))*coet
      z(i3)=(uval(5,i)*r(i1)+uval(6,i)*r(i2)+uval(3,i)*r(i3))*coet
      enddo
      enddo
!$OMP END PARALLEL DO
#ifdef USE_PA
      call fapp_stop('CGs_OMP_3', 0, 0)
#endif
c---------------------------------
      if(ite.gt.1)then
      call inner_product_s_block_sync(im,3*(n+ni),z,q,duplixyz,zq)
      betacg=zq/rhoa
      endif
c
#ifdef USE_PA
      call fapp_start('CGs_OMP_4_b', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('CGs_OMP_4_b', 0, 0)
      call fapp_start('CGs_OMP_4', 0, 0)
#endif
!$OMP PARALLEL DO default(none),private(i),shared(n,ni,p,z,betacg)
      do i=1,3*(n+ni)*npc
      p(i)=z(i)+betacg*p(i)
      enddo
!$OMP END PARALLEL DO
#ifdef USE_PA
      call fapp_stop('CGs_OMP_4', 0, 0)
#endif
c

      call compAmats_with_sync
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
c
      call inner_product_s_block_sync(im,3*(n+ni),z,r,duplixyz,rhob)
      call inner_product_s_block_sync(im,3*(n+ni),p,q,duplixyz,pq)
      alphacg=rhob/pq
      rhoa=rhob
c
#ifdef USE_PA
      call fapp_start('CGs_OMP_5_b', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('CGs_OMP_5_b', 0, 0)
      call fapp_start('CGs_OMP_5', 0, 0)
#endif
!$OMP PARALLEL DO default(none),private(i),shared(n,ni,q,alphacg,r,u,p)
      do i=1,3*(n+ni)*npc
      q(i)=-alphacg*q(i)
      r(i)=r(i)+q(i)
      u(i)=u(i)+alphacg*p(i)
      enddo
!$OMP END PARALLEL DO
#ifdef USE_PA
      call fapp_stop('CGs_OMP_5', 0, 0)
#endif
c
      call inner_product_s_block_sync(im,3*(n+ni),r,r,duplixyz,r2)
      err=sqrt(abs(r2/b2))
c
      if(output_log .and. im.eq.0)then
      if(icny.eq.10)then
      write(6,*) 'inner_tet10-loop RE=',err,ite,'-th iteration'
      else
      write(6,*) 'inner_tet4-loop RE=',err,ite,'-th iteration'
      endif
      endif
c
  101 continue
      if(output_log .and. im.eq.0)then
      if(icny.eq.10)then
      write(6,*) 'inner_tet10-loop RE=',err,ite,'-th iteration'
      else
      write(6,*) 'inner_tet4-loop RE=',err,ite,'-th iteration'
      endif
      endif

      end
c_______________________________________________________________________
      subroutine inner_product_d_sync(im,n,p,q,dupli,pq)
      implicit none
      INCLUDE 'mpif.h'
      integer im,n,ierr,i
      real*8 p(n),q(n),dupli(n),pq,tmp
      pq=0.0
#ifndef OMPSAME
!$OMP PARALLEL DO default(none),private(i),shared(n,p,q,dupli)
!$OMP& reduction(+:pq)
#endif
      do i=1,n
      pq=pq+p(i)*q(i)*dupli(i)
      enddo
#ifndef OMPSAME
!$OMP END PARALLEL DO
#endif
      call MPI_ALLREDUCE(pq,tmp,1,
     &  MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      pq=tmp
      end
c_______________________________________________________________________
      subroutine inner_product_d_block_sync(im,n,p,q,dupli,pq)
      implicit none
      INCLUDE 'mpif.h'
      integer im,n,ierr,i,j
      real*8 p(npc,n),q(npc,n),dupli(n),pq,pqtmp,tmp
      pq=0.0
#ifndef OMPSAME
!$OMP PARALLEL DO default(none),private(i,pqtmp),shared(n,p,q,dupli)
!$OMP& reduction(+:pq)
#endif
      do i=1,n
      pqtmp=0.0
      do j=1,npc
      pqtmp=pqtmp+p(j,i)*q(j,i)*dupli(i)
      enddo
      pq=pq+pqtmp
      enddo
#ifndef OMPSAME
!$OMP END PARALLEL DO
#endif
      call MPI_ALLREDUCE(pq,tmp,1,
     &  MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      pq=tmp
      end
c_______________________________________________________________________
      subroutine inner_product_s_sync(im,n,p,q,dupli,pq)
      implicit none
      INCLUDE 'mpif.h'
      integer im,n,ierr,i
      REAL_4 p(n),q(n),dupli(n),pq,tmp
      pq=0.0
#ifndef OMPSAME
!$OMP PARALLEL DO default(none),private(i),shared(n,p,q,dupli)
!$OMP& reduction(+:pq)
#endif
      do i=1,n
      pq=pq+p(i)*q(i)*dupli(i)
      enddo
#ifndef OMPSAME
!$OMP END PARALLEL DO
#endif
      call MPI_ALLREDUCE(pq,tmp,1,
     &  MPI_REAL_4,MPI_SUM,MPI_COMM_WORLD,ierr)
      pq=tmp
      end
c_______________________________________________________________________
      subroutine inner_product_s_block_sync(im,n,p,q,dupli,pq)
      implicit none
      INCLUDE 'mpif.h'
      integer im,n,ierr,i,j
      REAL_4 p(npc,n),q(npc,n),dupli(n),pq,pqtmp,tmp
#ifdef USE_PA
      call fapp_start('inner_product_s_block_sync_b', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('inner_product_s_block_sync_b', 0, 0)
      call fapp_start('inner_product_s_block_sync',0,0)
#endif
      pq=0.0
#ifndef OMPSAME
!$OMP PARALLEL DO default(none),private(i,pqtmp),shared(n,p,q,dupli)
!$OMP& reduction(+:pq)
#endif
      do i=1,n
      pqtmp=0.0
      do j=1,npc
      pqtmp=pqtmp+p(j,i)*q(j,i)*dupli(i)
      enddo
      pq=pq+pqtmp
      enddo
#ifndef OMPSAME
!$OMP END PARALLEL DO
#endif
#ifdef USE_PA
      call fapp_stop('inner_product_s_block_sync',0,0)
      call fapp_start('inner_product_s_block_sync_allr',0,0)
#endif
      call MPI_ALLREDUCE(pq,tmp,1,
     &  MPI_REAL_4,MPI_SUM,MPI_COMM_WORLD,ierr)
      pq=tmp
#ifdef USE_PA
      call fapp_stop('inner_product_s_block_sync_allr',0,0)
#endif
      end
c_______________________________________________________________________
      subroutine applydirichletd(n,u,bcarray)
      implicit none
      integer n,i,j
      real*8 u(3*npc*n),bcarray(n)
!$OMP PARALLEL DO default(none),private(i,j),shared(n,u,bcarray)
      do i=1,n
      do j=1,npc*3
      u(j+npc*3*(i-1))=u(j+npc*3*(i-1))*bcarray(i)
      enddo
      enddo
!$OMP END PARALLEL DO
      end
c_______________________________________________________________________
      subroutine applydirichlets(n,u,bcarray)
      implicit none
      integer n,i,j
      REAL_4 u(3*npc*n),bcarray(n)
!$OMP PARALLEL DO default(none),private(i,j),shared(n,u,bcarray)
      do i=1,n
      do j=1,npc*3
      u(j+npc*3*(i-1))=u(j+npc*3*(i-1))*bcarray(i)
      enddo
      enddo
!$OMP END PARALLEL DO
      end
c_______________________________________________________________________
