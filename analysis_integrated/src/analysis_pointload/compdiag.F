      subroutine compblockdiag(n,ne,kd,num,cny,coor,younglst,uval)
      integer num(ne),cny(10,ne)
      real*8 coor(3,n),younglst(2,kd)
      real*8 uval(6,n)
      real*8 young,rnyu
      real*8 detj,kmat(30,30),xx(4,3),gxy,gyz,gxz

c  1: (1,1), 2: (2,2), 3: (3,3), 4: (1,2), 5: (1,3), 6: (2,3)
c r0=K u
c calAmat
      do ie=1,ne
      in=num(ie)
      young=younglst(1,in)
      rnyu=younglst(2,in)
c
      gxy=young/(2.*(1+rnyu))
      gyz=young/(2.*(1+rnyu))
      gxz=young/(2.*(1+rnyu))
c
      do i1=1,4
      i2=cny(i1,ie)
      xx(i1,1)=coor(1,i2)
      xx(i1,2)=coor(2,i2)
      xx(i1,3)=coor(3,i2)
      enddo
c
      call nontet10(young,gxy,gyz,gxz,rnyu,xx,kmat,detj)
c
      do j1=1,10
      i1=cny(j1,ie)
      i2=3*(j1-1)
      uval(1,i1)=uval(1,i1)+kmat(i2+1,i2+1)
      uval(2,i1)=uval(2,i1)+kmat(i2+2,i2+2)
      uval(3,i1)=uval(3,i1)+kmat(i2+3,i2+3)
      uval(4,i1)=uval(4,i1)+kmat(i2+1,i2+2)
      uval(5,i1)=uval(5,i1)+kmat(i2+1,i2+3)
      uval(6,i1)=uval(6,i1)+kmat(i2+2,i2+3)
      enddo
c
      enddo
      end
c_______________________________________________________________________
      subroutine compblockdiag4(n,ne,kd,num,cny,coor,younglst,uval)
      integer num(ne),cny(4,ne)
      real*8 coor(3,n),younglst(2,kd)
      real*8 uval(6,n)
      real*8 young,rnyu
      real*8 detj,kmat(12,12),xx(4,3),gxy,gyz,gxz

c  1: (1,1), 2: (2,2), 3: (3,3), 4: (1,2), 5: (1,3), 6: (2,3)
c r0=K u
c calAmat
      do ie=1,ne
      in=num(ie)
      young=younglst(1,in)
      rnyu=younglst(2,in)
c
      gxy=young/(2.*(1+rnyu))
      gyz=young/(2.*(1+rnyu))
      gxz=young/(2.*(1+rnyu))
c
      do i1=1,4
      i2=cny(i1,ie)
      xx(i1,1)=coor(1,i2)
      xx(i1,2)=coor(2,i2)
      xx(i1,3)=coor(3,i2)
      enddo
c
      call nontetra4(young,gxy,gyz,gxz,rnyu,xx,kmat)
c
      do j1=1,4
      i1=cny(j1,ie)
      i2=3*(j1-1)
      uval(1,i1)=uval(1,i1)+kmat(i2+1,i2+1)
      uval(2,i1)=uval(2,i1)+kmat(i2+2,i2+2)
      uval(3,i1)=uval(3,i1)+kmat(i2+3,i2+3)
      uval(4,i1)=uval(4,i1)+kmat(i2+1,i2+2)
      uval(5,i1)=uval(5,i1)+kmat(i2+1,i2+3)
      uval(6,i1)=uval(6,i1)+kmat(i2+2,i2+3)
      enddo
c
      enddo
      end
c_______________________________________________________________________
      subroutine sync_uval
     &  (n_size,n,uval,nad,nadto,mpiad,mpiadlist,mpl,ierr,im
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind
#endif
     & )
      use mpi
      implicit none
      integer n_size,n,nad,nadto
      integer iad,i,j,k
      real*8 uval(6,n_size)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      ! temp
      real*8 senduval(6,nadto), recvuval(6,nadto)
      integer irecvreq(nad),isendreq(nad)
      integer itag,ierr,im,m1,m2,m3,istatus(MPI_STATUS_SIZE,nad)
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
!$OMP PARALLEL DO default(none),
!$OMP& shared (nadto,mpiadlist,senduval,uval),
!$OMP& private (i,j,k)
!!OCL PARALLEL
!OCL NORECURRENCE
      do i=1,nadto
      k=mpiadlist(i)
      do j=1,6
      senduval(j,i)=uval(j,k)
      enddo
      enddo
!$OMP END PARALLEL DO
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind)
#endif
      do i=1,nad
        m1=mpl(i)+1
        m2=6*mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !m3*10000+im
        call MPI_IRECV(recvuval(1,m1),m2,MPI_REAL8,
     &    m3,itag,MPI_COMM_WORLD,irecvreq(i),ierr)
      enddo
      do i=1,nad
        m1=mpl(i)+1
        m2=6*mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !im*10000+m3
        call MPI_ISEND(senduval(1,m1),m2,MPI_REAL8,
     &    m3,itag,MPI_COMM_WORLD,isendreq(i),ierr)
      enddo

      if(nad.ne.0)then
      call MPI_WAITALL(nad,irecvreq,istatus,ierr)
      call MPI_WAITALL(nad,isendreq,istatus,ierr)
      endif
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+2)
#endif
!$OMP PARALLEL default(none),
!$OMP& shared (nad,nadto,mpl,mpiadlist,recvuval,uval),
!$OMP& private (iad,i,j,k)
!!OCL SERIAL
      do iad=1,nad
!!OCL PARALLEL
!$OMP DO
!OCL NORECURRENCE
        do i = mpl(iad)+1, mpl(iad+1)
          k = mpiadlist(i)
          do j = 1, 6
            uval(j,k) = uval(j,k) + recvuval(j,i)
          enddo
        enddo
!$OMP END DO
      enddo
!$OMP END PARALLEL
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
#endif
      end
c_______________________________________________________________________
      subroutine compinv(n,uval)
      implicit none
      integer n,i
      real*8 uval(6,n),a(3,3)
      real*8 deno

c  1: (1,1), 2: (2,2), 3: (3,3), 4: (1,2), 5: (1,3), 6: (2,3)

!$OMP PARALLEL DO default(none),private(i,a,deno),shared(uval,n)
      do i=1,n
       a(1,1)=uval(1,i)
       a(2,2)=uval(2,i)
       a(3,3)=uval(3,i)
       a(1,2)=uval(4,i)
       a(1,3)=uval(5,i)
       a(2,3)=uval(6,i)
       a(2,1)=uval(4,i)
       a(3,1)=uval(5,i)
       a(3,2)=uval(6,i)
c
      deno=-(a(1,3)**2*a(2,2)) +
     -  2*a(1,2)*a(1,3)*a(2,3) -
     -  a(1,1)*a(2,3)**2 - a(1,2)**2*a(3,3) +
     -  a(1,1)*a(2,2)*a(3,3)
c
      uval(1,i)=(-a(2,3)**2 + a(2,2)*a(3,3))/deno
      uval(2,i)=(-a(1,3)**2 + a(1,1)*a(3,3))/deno
      uval(3,i)=(-a(1,2)**2 + a(1,1)*a(2,2))/deno
      uval(4,i)=(a(1,3)*a(2,3) - a(1,2)*a(3,3))/deno
      uval(5,i)=(-(a(1,3)*a(2,2)) + a(1,2)*a(2,3))/deno
      uval(6,i)=(a(1,2)*a(1,3) - a(1,1)*a(2,3))/deno
c
      enddo
!$OMP END PARALLEL DO
      end
c_______________________________________________________________________
      subroutine compblockdiaginv
     & (n_size,ne_size,n,ne,kd,num,cny,coor,younglst,
     & ni_size,nei_size,ni,nei,
     & nad,nadto,mpiad,mpiadlist,mpl,im,mpibuf,
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind,
#endif
     & uval)
      implicit none
      integer n_size,ne_size,n,ne,kd,im,nad,nadto
      integer ni_size,nei_size,ni,nei
      integer num(ne_size),cny(10,ne_size)
      real*8 coor(3,n_size),younglst(2,kd),uval(6,n_size+ni_size)
      integer ierr,i,j
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real*8 mpibuf(2*3*npc*(nadto+1))
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
      if(im.eq.0)then
      write(*,*) 'compblockdiaginv'
      endif
      call compblockdiag(n,ne,kd,num,cny,coor,younglst,uval)
      call sync_uval
     & (n_size+ni_size,n+ni,uval,nad,nadto,mpiad,mpiadlist,mpl,
     & ierr,im
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      call compinv(n+ni,uval)
      end
c_______________________________________________________________________
      subroutine compblockdiaginv4
     & (n_size,ne_size,n,ne,kd,num,cny,coor,younglst,
     & ni_size,nei_size,ni,nei,
     & nad,nadto,mpiad,mpiadlist,mpl,im,mpibuf,
#ifdef MEASURE_TIME_BARRIER
     & ,mp,timerbufd,timerbufi,ind,
#endif
     & uval)
      implicit none
      integer n_size,ne_size,n,ne,kd,im,nad,nadto
      integer ni_size,nei_size,ni,nei
      integer num(ne_size),cny(4,ne_size)
      real*8 coor(3,n_size),younglst(2,kd),uval(6,n_size+ni_size)
      integer ierr,i,j
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      real*8 mpibuf(2*3*npc*(nadto+1))
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
      if(im.eq.0)then
      write(*,*) 'compblockdiaginv4'
      endif
      call compblockdiag4(n,ne,kd,num,cny,coor,younglst,uval)
      call sync_uval
     & (n_size+ni_size,n+ni,uval,nad,nadto,mpiad,mpiadlist,mpl,
     & ierr,im
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      call compinv(n+ni,uval)
      end
c_______________________________________________________________________