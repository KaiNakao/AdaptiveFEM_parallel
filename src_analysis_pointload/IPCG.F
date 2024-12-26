      subroutine IPCG_setting
     & (im,use_crs_d,overlap_comm_d,
     & overlap_comm_s,persistent_comm_s)
        implicit none
        integer, intent(in) :: im
        logical, intent(out) :: use_crs_d, overlap_comm_d
        logical, intent(out) :: overlap_comm_s, persistent_comm_s
        integer :: case_d
        use_crs_d = .true.
        overlap_comm_d = .false.
        overlap_comm_s = .false.
        persistent_comm_s = .false.
!        persistent_comm_s = .false.

        if(im == 0)then
          write(6,*) "use_crs_d:", use_crs_d
          write(6,*) "overlap_comm_d:", overlap_comm_d
          write(6,*) "overlap_comm_s:", overlap_comm_s
          write(6,*) "persistent_comm_s:", persistent_comm_s
        endif
      end subroutine IPCG_setting

      subroutine IPCG
     & (im,np,kd,younglst,sig,mpibuf,
     & n,ne,coor,cny,num,uval,duplixyz,bcarray,
     & ni,nei,cnyi,dsi,numi,nfla,
     & nad,nadto,mpiad,mpiadlist,mpl,
     & nird,pentarptrd,pentarindd,
     & niwd,pentawptrd,pentawindd,pentacrsvald,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & ncolor,color_ind,
     & n4,ne4,coor4,cny4,num4,uval4,duplixyz4,bcarray4,
     & ni4,nei4,cnyi4,dsi4,numi4,nfla4,
     & nad4,nadto4,mpiad4,mpiadlist4,mpl4,
     & nir4,pentarptr4,pentarind4,
     & niw4,pentawptr4,pentawind4,pentacrsval4,
     & map4to10,
     & cnysep_read4,cnysep_write4,nsep_n_ptr4,nsep4,nsepnmap4,
     & uvalcoe,log_innerCG,log_innerCG_per_outite,
     & use_crs_d,overlap_comm_d,overlap_comm_s,persistent_comm_s,
     & b,u)
      implicit none
      INCLUDE 'mpif.h'
      integer kd,im
      integer ite,i,j,i1,i2,i3,ie,ii
      real*8 younglst(2,kd)
      REAL_4 tmps1,tmps2,tmps3
      integer icpu,ncpu

      integer n,ni,ne,nei,np
      integer cny(10,ne),num(ne),cnyi(12,nei)
      real*8 uval(6,n+ni),coor(3,n+ni)
      real*8 duplixyz(3*(n+ni)),bcarray(n+ni),dsi
      integer numi(nei),nfla(nei)
      integer nad,nadto
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      integer nird,niwd
      integer pentawptrd(2*nsta+3),pentarptrd(niwd+1)
      integer pentawindd(niwd),pentarindd(nird)
      integer nir,niw
      integer pentawptr(2*nsta+3),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)

      integer cnysep_read(10,ne),cnysep_write(10,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep,nsepnmap(nsep)
      integer cnysep_read4(4,ne4),cnysep_write4(4,ne4)
      integer nsep_n_ptr4(nsta+2)
      integer nsep4,nsepnmap4(nsep4)
      integer ncolor,color_ind(np*MAXCOLOR)

      integer n4,ni4,ne4,nei4
      integer cny4(4,ne4),num4(ne4),cnyi4(6,nei4)
      real*8 uval4(6,n4+ni4),coor4(3,n4+ni4)
      real*8 duplixyz4(3*(n4+ni4)),bcarray4(n4+ni4),dsi4
      integer numi4(nei4),nfla4(nei4)
      integer nad4,nadto4
      integer mpiad4(nad4,2),mpiadlist4(nadto4),mpl4(nad4+1)
      integer nir4,niw4
      integer pentawptr4(2*nsta+3),pentarptr4(niw4+1)
      integer pentawind4(niw4),pentarind4(nir4)
#if npc == 1 || npc == 21
      real*8 pentacrsvald(9,nird)
      real*8 pentacrsval(9,nir)
      REAL_4 pentacrsvals(9,nir)
      real*8 pentacrsval4(9,nir4)
      REAL_4 pentacrsval4s(9,nir4)
#elif npc == 288
      real*8 pentacrsvald(9,nird,2)
      real*8 pentacrsval(9,nir,2)
      REAL_4 pentacrsvals(9,nir,2)
      real*8 pentacrsval4(9,nir4,2)
      REAL_4 pentacrsval4s(9,nir4,2)
#else
ERROR: invalid npc
#endif

      integer map4to10(2,n+ni)

      real*8 b2,tmp
      real*8 alphacg,betacg,r2,rhoa,rhob,pq,zq,err,eps,z2,p2,q2
      real*8 b(3*(n+ni)*npc),u(3*(n+ni)*npc),r(3*(n+ni)*npc)
      real*8 p(3*(n+ni)*npc),q(3*(n+ni)*npc)
      real*8 qsep(3*nsep*npc)
      REAL_4 qseps(3*nsep*npc)
      real*8 z(3*(n+ni)*npc)
      real*8 sig(nkl),uvalcoe(npc),coet,scalefactor,coe3n
      integer itermax,nprint,ierr
      real*8 mpibuf(2*3*npc*(nadto+1))
      REAL_4 mpibufs(2*3*npc*(nadto+1))
      integer*2
     $     mpibufi4(2*(3*npc*(nadto4+1)+12*npc*(nadto4/DIV_N+nad4)))
      integer*2 mpibufi(2*(3*npc*(nadto+1)+12*npc*(nadto/DIV_N+nad)))
      REAL_4 tmps,sigfs(5),utmp(3*(n+ni))
      REAL_4 younglsts(2,kd),sigs(nkl)
      REAL_4 uvalcoes(npc)
c
      REAL_4 coors(3,n+ni),duplixyzs(3*(n+ni))
      REAL_4 bcarrays(n+ni),epss,dsis
      REAL_4 uvals(6,n+ni)
      REAL_4 rs(3*(n+ni)*npc),zs(3*(n+ni)*npc)

      REAL_4 coor4s(3,n4+ni4),duplixyz4s(3*(n4+ni4))
      REAL_4 bcarray4s(n4+ni4),eps4s,dsi4s
      REAL_4 uval4s(6,n4+ni4)
      REAL_4 z4s(3*(n4+ni4)*npc)
      REAL_4 r4s(3*(n4+ni4)*npc)
      REAL_4 utmps(3*(n+ni))
      real*8 t10,t4,ttmp1,ttmp10,ttmp4

      logical log_innerCG,log_innerCG_per_outite
      integer ite_tet4,ite_tet10
      integer ite_tet4_tot,ite_tet10_tot

      logical use_crs_d,overlap_comm_d
      logical overlap_comm_s,persistent_comm_s
      integer mpireqs(nad,2),mpireq4s(nad,4)
c
      t10=0
      t4=0
c
      ite_tet4_tot = 0
      ite_tet10_tot = 0
c
      sigs=sig
      uvalcoes=uvalcoe
      younglsts=younglst
c
      uvals=uval
      coors=coor
      duplixyzs=duplixyz
      bcarrays=bcarray
      dsis=dsi
      pentacrsvals=pentacrsval

      uval4s=uval4
      coor4s=coor4
      duplixyz4s=duplixyz4
      bcarray4s=bcarray4
      dsi4s=dsi4
      pentacrsval4s=pentacrsval4

      if(persistent_comm_s)then
      call sync_parallel_block_s_persistent_comm_init
     & (nad,nadto,mpiad,mpl,im,mpibufi,mpireqs)
      call sync_parallel_block_s_persistent_comm_init
     & (nad4,nadto4,mpiad4,mpl4,im,mpibufi4,mpireq4s)
      endif

      tmp=0
      do i=1,3*(n+ni)
      tmp=tmp+duplixyz(i)
      enddo
      call MPI_ALLREDUCE(tmp,coe3n,1,
     &  MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(im.eq.0)then
      write(6,*) 'coe3n',coe3n
      endif

      eps=1.0e-8
c--------------------------
      err=1.0
#ifndef USEIEPENTA
      call applydirichletd(n+ni,b,bcarray)
      call applydirichletd(n+ni,u,bcarray)
#endif
      call inner_product_d_block_sync(im,3*(n+ni),b,b,duplixyz,b2)
!$OMP PARALLEL DO default(none),private(i),shared(n,ni,p)
      do i=1,3*(n+ni)*npc
      p(i)=0.0
      enddo
!$OMP END PARALLEL DO

      if(im.eq.0)then
      write(6,*) 'b2',b2
      endif
c___________________________________
      call compAmatd_with_sync
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,bcarray,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nird,pentarptrd,pentarindd,
     & niwd,pentawptrd,pentawindd,pentacrsvald,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & overlap_comm_d,use_crs_d,
     & u,r
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )

      call inner_product_d_block_sync(im,3*(n+ni),r,r,duplixyz,r2)
      if(im.eq.0)then
      write(6,*) 'out-loop r2=',r2
      endif

!$OMP PARALLEL DO default(none),private(i),shared(n,ni,b,r)
      do i=1,3*(n+ni)*npc
      r(i)=b(i)-r(i)
      enddo
!$OMP END PARALLEL DO
      call inner_product_d_block_sync(im,3*(n+ni),r,r,duplixyz,r2)
      err=sqrt(abs(r2/b2))
      if(im.eq.0)then
      write(6,*) 'out-loop RE=',err,'init'
      endif
c
      betacg=0.0
c
      ite=0
      do 101 while (err.gt.eps)
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ttmp1=MPI_WTIME()

      ite=ite+1
c---------------------------------
c preconditioner
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
#ifndef PCGE
      if(im.eq.0 .and. ite == 1)then
      write(*,*) 'use preconditioner'
      endif
      scalefactor=(1.0/sqrt(r2))*(coe3n*npc)
      rs=r*scalefactor
      zs=z*scalefactor
#ifdef TET4GRID
      r4s=0.
      z4s=0.
      do j=1,n+ni
      if(map4to10(1,j) == 0) cycle
      do i=1,3
      do ii=1,npc

       r4s(3*npc*(map4to10(1,j)-1)+ii+(i-1)*npc)=
     & r4s(3*npc*(map4to10(1,j)-1)+ii+(i-1)*npc)+
     & 0.5*rs(3*npc*(j-1)+ii+(i-1)*npc)*duplixyzs(3*(j-1)+i)
       r4s(3*npc*(map4to10(2,j)-1)+ii+(i-1)*npc)=
     & r4s(3*npc*(map4to10(2,j)-1)+ii+(i-1)*npc)+
     & 0.5*rs(3*npc*(j-1)+ii+(i-1)*npc)*duplixyzs(3*(j-1)+i)

       z4s(3*npc*(map4to10(1,j)-1)+ii+(i-1)*npc)=
     & z4s(3*npc*(map4to10(1,j)-1)+ii+(i-1)*npc)+
     & 0.5*zs(3*npc*(j-1)+ii+(i-1)*npc)*duplixyzs(3*(j-1)+i)
       z4s(3*npc*(map4to10(2,j)-1)+ii+(i-1)*npc)=
     & z4s(3*npc*(map4to10(2,j)-1)+ii+(i-1)*npc)+
     &0.5*zs(3*npc*(j-1)+ii+(i-1)*npc)*duplixyzs(3*(j-1)+i)

      enddo
      enddo
      enddo
      call sync_parallel_block_s
     &  (n4+ni4,n4+ni4,r4s,nad4,nadto4,mpiad4,mpiadlist4,mpl4,
     & ierr,im,mpibufs
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      call sync_parallel_block_s
     &  (n4+ni4,n4+ni4,z4s,nad4,nadto4,mpiad4,mpiadlist4,mpl4,
     & ierr,im,mpibufs
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      eps4s=0.075
      itermax=300
!      eps4s=0.125
!      itermax=2000

!      do i=1,3
!      itermax=100
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t4=t4-MPI_WTIME()
      ttmp4=MPI_WTIME()
      call CGs
     & (im,np,n4,ni4,ne4,nei4,kd,
     & coor4s,cny4,num4,cnyi4,younglsts,r4s,z4s,uval4s,sigs,uvalcoes,
     & dsi4s,numi4,nfla4,duplixyz4s,bcarray4s,eps4s,itermax,0,
     & nad4,nadto4,mpiad4,mpiadlist4,mpl4,mpibufi4,4,6,
     & nir4,pentarptr4,pentarind4,
     & niw4,pentawptr4,pentawind4,pentacrsval4s,
     & nsep4,numsep,nsep_n_ptr4,nsep_ne_ptr,nsepnmap4,
     & cnysep_read4,cnysep_write4,qseps,
     & ncolor,color_ind,log_innerCG,ite_tet4,
     & overlap_comm_s,persistent_comm_s,mpireq4s)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t4=t4+MPI_WTIME()
      ttmp4=MPI_WTIME()-ttmp4
!      enddo

      do j=1,n+ni
      if(map4to10(1,j) == 0) cycle
      do i=1,3
      do ii=1,npc
      zs(3*npc*(j-1)+ii+(i-1)*npc)=duplixyzs(3*(j-1)+i)*0.5*(
     &z4s(3*npc*(map4to10(1,j)-1)+ii+(i-1)*npc) +
     &z4s(3*npc*(map4to10(2,j)-1)+ii+(i-1)*npc) )
      enddo
      enddo
      enddo
!~#ifdef DEBUGVTK
!~!      call write_tet_vtk_s
!~!     - (n4,ne4,n4,ne4,4,cny4,coor4s,coor4s,im,1)
!~      z4s=0
!~      do j=1,n+ni
!~      do i=1,3
!~      z4s(3*npc*(map4to10(1,j)-1)+1+(i-1)*npc)=
!~     &z4s(3*npc*(map4to10(1,j)-1)+1+(i-1)*npc)+
!~     & 0.5*coors(i,j)*duplixyzs(3*(j-1)+i)
!~      z4s(3*npc*(map4to10(2,j)-1)+1+(i-1)*npc)=
!~     &z4s(3*npc*(map4to10(2,j)-1)+1+(i-1)*npc)+
!~     & 0.5*coors(i,j)*duplixyzs(3*(j-1)+i)
!~      enddo
!~      enddo
!~      call sync_parallel_block_s
!~     &  (n4+ni4,n4+ni4,z4s,
!~     & nad4,nadto4,mpiad4,mpiadlist4,mpl4,ierr,im,mpibufs
!~#ifdef MEASURE_TIME_BARRIER
!~     - ,mp,timerbufd,timerbufi,ind
!~#endif
!~     - )
!~      do j=1,3*(n4+ni4)
!~      utmps(j)=z4s(npc*(j-1)+1)
!~      enddo
!~      call write_tet_vtk_s
!~     - (n4,ne4,n4,ne4,4,cny4,coor4s,utmps,im,1)
!~      call write_tet4_tri_vtk_s
!~     - (n4+ni4,ne4,nei4,cny4,cnyi4,coor4s,utmps,im,1)
!~
!~      do j=1,3*(n4+ni4)
!~      utmps(j)=duplixyz4s(j)
!~      enddo
!~      call write_tet_vtk_s
!~     - (n4,ne4,n4,ne4,4,cny4,coor4s,utmps,im,2)
!~      call write_tet4_tri_vtk_s
!~     - (n4+ni4,ne4,nei4,cny4,cnyi4,coor4s,utmps,im,2)
!~
!~      do j=1,3*(n+ni)
!~      utmps(j)=duplixyzs(j)
!~      enddo
!~      call write_tet_vtk_s
!~     - (n,ne,n,ne,10,cny,coors,utmps,im,2)
!~      call write_tet_tri_vtk_s
!~     - (n+ni,ne,nei,cny,cnyi,coors,utmps,im,2)
!~
!~      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!~      stop
!~#endif
      call sync_parallel_block_s
     &  (n+ni,n+ni,zs,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibufs
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
#endif
      epss=0.125
      itermax=50
!      itermax=300
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t10=t10-MPI_WTIME()
      ttmp10=MPI_WTIME()
      call CGs
     & (im,np,n,ni,ne,nei,kd,
     & coors,cny,num,cnyi,younglsts,rs,zs,uvals,sigs,uvalcoes,
     & dsis,numi,nfla,duplixyzs,bcarrays,epss,itermax,0,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibufi,10,12,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsvals,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,qseps,
     & ncolor,color_ind,log_innerCG,ite_tet10,
     & overlap_comm_s,persistent_comm_s,mpireqs)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t10=t10+MPI_WTIME()
      ttmp10=MPI_WTIME()-ttmp10

      scalefactor=1.0/scalefactor
      z=zs*scalefactor
      call sync_parallel_block_d
     &  (n+ni,n+ni,z,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      do i=1,3*(n+ni)
      do j=1,npc
      z(npc*(i-1)+j)=z(npc*(i-1)+j)*duplixyz(i)
      enddo
      enddo
c
#endif
c---------------------------------
      if(ite.gt.1)then
      call inner_product_d_block_sync(im,3*(n+ni),z,q,duplixyz,zq)
      betacg=zq/rhoa
      endif
c
      do i=1,3*(n+ni)*npc
      p(i)=z(i)+betacg*p(i)
      enddo
c
      call compAmatd_with_sync
     & (im,n,ni,ne,kd,num,cny,coor,
     & younglst,sig,uvalcoe,bcarray,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,
     & np,ncolor,color_ind,
     & nird,pentarptrd,pentarindd,
     & niwd,pentawptrd,pentawindd,pentacrsvald,
     & nad,nadto,mpiad,mpiadlist,mpl,mpibuf,
     & overlap_comm_d,use_crs_d,
     & p,q
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
c
!      call inner_product_d_block_sync(im,3*(n+ni),z,z,duplixyz,z2)
!      call inner_product_d_block_sync(im,3*(n+ni),q,q,duplixyz,q2)
!      call inner_product_d_block_sync(im,3*(n+ni),p,p,duplixyz,p2)
      call inner_product_d_block_sync(im,3*(n+ni),z,r,duplixyz,rhob)
      call inner_product_d_block_sync(im,3*(n+ni),p,q,duplixyz,pq)
!      if(im.eq.0)then
!      write(6,*) 'z2',z2
!      write(6,*) 'q2',q2
!      write(6,*) 'p2',p2
!      write(6,*) 'rhob',rhob
!      write(6,*) 'pq',pq
!      endif

      alphacg=rhob/pq
      rhoa=rhob
c
      do i=1,3*(n+ni)*npc
      q(i)=-alphacg*q(i)
      r(i)=r(i)+q(i)
      u(i)=u(i)+alphacg*p(i)
      enddo
c
      call inner_product_d_block_sync(im,3*(n+ni),r,r,duplixyz,r2)
      err=sqrt(abs(r2/b2))
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ttmp1=MPI_WTIME()-ttmp1

      if(log_innerCG_per_outite .and. im.eq.0)then
      write(6,'(a,e15.6,i6,2a,2i6)')
     & ' out-loop RE=',err,ite,'-th iteration, ',
     & 'iter (tet4,tet10):',ite_tet4,ite_tet10
      write(6,*) 'tot,ic,if took ',ttmp1,ttmp4,ttmp10
      endif
      ite_tet4_tot = ite_tet4_tot + ite_tet4
      ite_tet10_tot = ite_tet10_tot + ite_tet10
c
  101 continue
      if(im.eq.0)then
      write(6,'(2a,3i6,e15.6)') " total CG iterations",
     & " (outer,tet4,tet10):",ite,ite_tet4_tot,ite_tet10_tot
      write(6,*) 'inner fine,coarse took',t10,t4
      endif

      if(persistent_comm_s)then
        call sync_parallel_persistent_comm_free(nad,mpireqs)
        call sync_parallel_persistent_comm_free(nad4,mpireq4s)
      endif

      end
c_______________________________________________________________________
