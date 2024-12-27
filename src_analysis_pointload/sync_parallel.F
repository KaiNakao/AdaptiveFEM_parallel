c_______________________________________________________________________
      subroutine sync_parallel_block_d
     &  (n_size,n,b,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      implicit none
      INCLUDE 'mpif.h'
      integer n_size,n,nad,nadto
      integer iad,i,j,k,ipc
      real*8 b(3*n_size*npc)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      ! temp
      real*8 mpibuf(2*3*npc*(nadto+1))
      integer irecvreq(nad),isendreq(nad)
      integer itag,ierr,im,m1,m2,m3,istatus(MPI_STATUS_SIZE,nad)
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
!$OMP PARALLEL DO default(none),
!$OMP& shared (nadto,mpiadlist,mpibuf,b),
!$OMP& private (i,j,k,ipc)
!!OCL PARALLEL
!OCL NORECURRENCE
      do i=1,nadto
      j=3*npc*(i-1)
      k=3*npc*(mpiadlist(i)-1)
      do ipc=1,3*npc
      mpibuf(j+ipc)=b(k+ipc)
      enddo
      enddo
!$OMP END PARALLEL DO
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind)
#endif
      do i=1,nad
        m1=3*npc*mpl(i)+1+3*nadto*npc
        m2=3*npc*mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !m3*10000+im
        call MPI_IRECV(mpibuf(m1),m2,MPI_REAL8,
     &    m3,itag,MPI_COMM_WORLD,irecvreq(i),ierr)
      enddo
      do i=1,nad
        m1=3*npc*mpl(i)+1
        m2=3*npc*mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !im*10000+m3
        call MPI_ISEND(mpibuf(m1),m2,MPI_REAL8,
     &    m3,itag,MPI_COMM_WORLD,isendreq(i),ierr)
      enddo

      if(nad.ne.0)then
      call MPI_WAITALL(nad,irecvreq(1),istatus,ierr)
      call MPI_WAITALL(nad,isendreq(1),istatus,ierr)
      endif
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+2)
#endif
!$OMP PARALLEL default(none),
!$OMP& shared (nad,nadto,mpl,mpiadlist,mpibuf,b),
!$OMP& private (iad,i,j,k,ipc)
!!OCL SERIAL
      do iad=1,nad
!!OCL PARALLEL
!$OMP DO
!OCL NORECURRENCE
        do i = mpl(iad)+1, mpl(iad+1)
          j=3*(i-1)*npc+3*nadto*npc
          k=3*(mpiadlist(i)-1)*npc
          do ipc=1,3*npc
          b(k+ipc)=b(k+ipc)+mpibuf(j+ipc)
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
      subroutine sync_parallel_block_s
     &  (n_size,n,b,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      implicit none
      INCLUDE 'mpif.h'
      integer n_size,n,nad,nadto
      integer iad,i,j,k,ipc
      REAL_4 b(3*n_size*npc)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      ! temp
      REAL_4 mpibuf(2*3*npc*(nadto+1))
      integer irecvreq(nad),isendreq(nad)
      integer itag,ierr,im,m1,m2,m3,istatus(MPI_STATUS_SIZE,nad)
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
!$OMP PARALLEL DO default(none),
!$OMP& shared (nadto,mpiadlist,mpibuf,b),
!$OMP& private (i,j,k,ipc)
!!OCL PARALLEL
!OCL NORECURRENCE
      do i=1,nadto
      j=3*npc*(i-1)
      k=3*npc*(mpiadlist(i)-1)
      do ipc=1,3*npc
      mpibuf(j+ipc)=b(k+ipc)
      enddo
      enddo
!$OMP END PARALLEL DO
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind)
#endif
      do i=1,nad
        m1=3*npc*mpl(i)+1+3*nadto*npc
        m2=3*npc*mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !m3*10000+im
        call MPI_IRECV(mpibuf(m1),m2,MPI_REAL_4,
     &    m3,itag,MPI_COMM_WORLD,irecvreq(i),ierr)
      enddo
      do i=1,nad
        m1=3*npc*mpl(i)+1
        m2=3*npc*mpiad(i,2)
        m3=mpiad(i,1)
        itag=1 !im*10000+m3
        call MPI_ISEND(mpibuf(m1),m2,MPI_REAL_4,
     &    m3,itag,MPI_COMM_WORLD,isendreq(i),ierr)
      enddo

      if(nad.ne.0)then
      call MPI_WAITALL(nad,irecvreq(1),istatus,ierr)
      call MPI_WAITALL(nad,isendreq(1),istatus,ierr)
      endif
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+2)
#endif
!$OMP PARALLEL default(none),
!$OMP& shared (nad,nadto,mpl,mpiadlist,mpibuf,b),
!$OMP& private (iad,i,j,k,ipc)
!!OCL SERIAL
      do iad=1,nad
!!OCL PARALLEL
!$OMP DO
!OCL NORECURRENCE
        do i = mpl(iad)+1, mpl(iad+1)
          j=3*(i-1)*npc+3*nadto*npc
          k=3*(mpiadlist(i)-1)*npc
          do ipc=1,3*npc
          b(k+ipc)=b(k+ipc)+mpibuf(j+ipc)
          enddo
        enddo
!$OMP END DO
      enddo
!$OMP END PARALLEL
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
#endif
      end

      subroutine sync_parallel_block_s_part1
     &  (n_size,n,b,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf,
     & ireq,persistent_comm
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      use mpi
      implicit none
      integer n_size,n,nad,nadto
      integer iad,i,j,k,ipc
      REAL_4 b(3*n_size*npc)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      logical persistent_comm
      ! temp
      integer*2 mpibuf(2*(3*npc*(nadto+1)+12*npc*(nadto/DIV_N+nad)))
      integer ireq(nad*2)
      integer itag,ierr,im,m1,m2,m3,istatus(MPI_STATUS_SIZE,nad)
      integer ib,nb,nn
      REAL_4 mpibuf_f(npc,3,2,nadto/DIV_N+1)
      REAL_4 max_b(npc,3), min_b(npc,3)
      integer is(nadto), ie(nadto)
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
#ifdef USE_PA
      call fapp_start('sync_parallel_block_s_b1', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('sync_parallel_block_s_b1', 0, 0)
      call fapp_start('sync_parallel_block_s_pack', 0, 0)
#endif
!#define DEBUGINTEGERCOMM
#ifdef DEBUGINTEGERCOMM
      write(*,*) 'DEBUGINTEGERCOMM sync_parallel_block_s_part1'
#endif

!$OMP PARALLEL default(none),
!$OMP& shared (nad,nadto,mpl,mpiadlist,mpibuf,mpibuf_f,b),
!$OMP& private (i,j,k,ipc,iad,nb,ib,is,ie,nn,max_b,min_b)
      do iad=1,nad

         if(iad.eq.1) nn = 0
         nb = (mpl(iad+1)-mpl(iad)-1)/DIV_N+1
         do ib=1,nb
            is(ib) = mpl(iad)+1+DIV_N*(ib-1)
            ie(ib) = min(mpl(iad+1),mpl(iad)+DIV_N*(ib))
         enddo

!$OMP DO 
         do ib=1,nb
         max_b(:,:) = -huge(b(1))
         min_b(:,:) = huge(b(1))
         do i = is(ib), ie(ib)
            k=3*npc*(mpiadlist(i)-1)
            do ipc=1,npc
               max_b(ipc,1) = max(max_b(ipc,1),b(k+npc*0+ipc))
               max_b(ipc,2) = max(max_b(ipc,2),b(k+npc*1+ipc))
               max_b(ipc,3) = max(max_b(ipc,3),b(k+npc*2+ipc))
               min_b(ipc,1) = min(min_b(ipc,1),b(k+npc*0+ipc))
               min_b(ipc,2) = min(min_b(ipc,2),b(k+npc*1+ipc))
               min_b(ipc,3) = min(min_b(ipc,3),b(k+npc*2+ipc))
            enddo
         enddo
         do ipc=1,npc
            mpibuf_f(ipc,1,1,ib)=(max_b(ipc,1)-min_b(ipc,1))/65535.0
            mpibuf_f(ipc,2,1,ib)=(max_b(ipc,2)-min_b(ipc,2))/65535.0
            mpibuf_f(ipc,3,1,ib)=(max_b(ipc,3)-min_b(ipc,3))/65535.0

#ifdef DEBUGINTEGERCOMM
            if(mpibuf_f(ipc,1,1,ib).eq.0.0)then
              mpibuf_f(ipc,1,1,ib)=1.0
            endif
            if(mpibuf_f(ipc,2,1,ib).eq.0.0)then
              mpibuf_f(ipc,2,1,ib)=1.0
            endif
            if(mpibuf_f(ipc,3,1,ib).eq.0.0)then
              mpibuf_f(ipc,3,1,ib)=1.0
            endif
#endif
            mpibuf_f(ipc,1,2,ib) = (max_b(ipc,1)+min_b(ipc,1))/2.0
            mpibuf_f(ipc,2,2,ib) = (max_b(ipc,2)+min_b(ipc,2))/2.0
            mpibuf_f(ipc,3,2,ib) = (max_b(ipc,3)+min_b(ipc,3))/2.0
         enddo
         enddo
      
!$OMP DO
         do ib=1,nb
         do i = is(ib), ie(ib)
            j=3*npc*(i-1)+12*npc*nn
            k=3*npc*(mpiadlist(i)-1)
            do ipc=1,npc
c      mpibuf(j+ipc)=b(k+ipc)
               mpibuf(j+npc*0+ipc)=
     $              int((b(k+npc*0+ipc)-mpibuf_f(ipc,1,2,ib))
     $              /mpibuf_f(ipc,1,1,ib),2)
               mpibuf(j+npc*1+ipc)=
     $              int((b(k+npc*1+ipc)-mpibuf_f(ipc,2,2,ib))
     $              /mpibuf_f(ipc,2,1,ib),2)
               mpibuf(j+npc*2+ipc)=
     $              int((b(k+npc*2+ipc)-mpibuf_f(ipc,3,2,ib))
     $              /mpibuf_f(ipc,3,1,ib),2)
            enddo
         enddo
         enddo
            
         call pack_scalevalue(nb, mpibuf_f(1,1,1,1),
     $        mpibuf(3*npc*mpl(iad+1)+12*npc*nn+1))
         nn = nn+nb
      enddo
!$OMP END PARALLEL
#ifdef USE_PA
      call fapp_stop('sync_parallel_block_s_pack', 0, 0)
#endif
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind)
#endif
#ifdef USE_PA
      call fapp_start('sync_parallel_block_s_b2', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('sync_parallel_block_s_b2', 0, 0)
      call fapp_start('sync_parallel_block_s_sendrecv', 0, 0)
#endif
      if(persistent_comm)then
        if(nad>0)then
          call MPI_STARTALL(nad*2,ireq,ierr)
      !     call FJMPI_PREQUEST_STARTALL(nad*2,ireq,ierr)
        endif
      else
        ! non persistent communication
        nn = 0
        do i=1,nad
          nb = (mpl(i+1)-mpl(i)-1)/DIV_N+1
          m1=3*npc*mpl(i)+1+12*npc*nn
     $         +3*nadto*npc+12*npc*(nadto/DIV_N+nad)
          m2=3*npc*mpiad(i,2)+12*npc*nb
          m3=mpiad(i,1)
          itag=1                !m3*10000+im
          call MPI_IRECV(mpibuf(m1),m2,MPI_INTEGER2,
     &         m3,itag,MPI_COMM_WORLD,ireq(i),ierr)
          nn = nn+nb
        enddo
        nn = 0
        do i=1,nad
          nb = (mpl(i+1)-mpl(i)-1)/DIV_N+1
          m1=3*npc*mpl(i)+1+12*npc*nn
          m2=3*npc*mpiad(i,2)+12*npc*nb
          m3=mpiad(i,1)
          itag=1                !im*10000+m3
          call MPI_ISEND(mpibuf(m1),m2,MPI_INTEGER2,
     &         m3,itag,MPI_COMM_WORLD,ireq(i+nad),ierr)
          nn = nn+nb
        enddo
      endif ! persistent_comm
#ifdef USE_PA
      call fapp_stop('sync_parallel_block_s_sendrecv', 0, 0)
#endif
c      call FJMPI_PROGRESS_START()
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
#endif
      end subroutine sync_parallel_block_s_part1

      subroutine sync_parallel_block_s_part2
     &  (n_size,n,b,nad,nadto,mpiad,mpiadlist,mpl,ierr,im,mpibuf,
     & ireq
#ifdef MEASURE_TIME_BARRIER
     - ,mp,timerbufd,timerbufi,ind
#endif
     - )
      implicit none
      INCLUDE 'mpif.h'
      integer n_size,n,nad,nadto
      integer iad,i,j,k,ipc
      REAL_4 b(3*n_size*npc)
      integer mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      ! temp
      integer*2 mpibuf(2*(3*npc*(nadto+1)+12*npc*(nadto/DIV_N+nad)))
      integer ireq(nad,2)
      integer itag,ierr,im,m1,m2,m3,istatus(MPI_STATUS_SIZE,nad)
      integer ib,nb,nn
      REAL_4 mpibuf_r(npc,3,2,nadto/DIV_N+1)
      integer is(nadto), ie(nadto)
#ifdef MEASURE_TIME_BARRIER
      integer ind,mp
#include "timerbuf.F"
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+1)
#endif
c      call FJMPI_PROGRESS_STOP()
#ifdef USE_PA
      call fapp_start('sync_parallel_block_s_b3', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('sync_parallel_block_s_b3', 0, 0)
      call fapp_start('sync_parallel_block_s_waitall', 0, 0)
#endif
      if(nad.ne.0)then
      call MPI_WAITALL(nad*2,ireq(1,1),MPI_STATUS_IGNORE,ierr)
      endif
#ifdef USE_PA
      call fapp_stop('sync_parallel_block_s_waitall', 0, 0)
#endif
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
      call start_measure_with_barrier
     - (im,timerbufd,timerbufi,ind+2)
#endif
#ifdef USE_PA
      call fapp_start('sync_parallel_block_s_b4', 0, 0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call fapp_stop('sync_parallel_block_s_b4', 0, 0)
      call fapp_start('sync_parallel_block_s_unpack', 0, 0)
#endif
!$OMP PARALLEL default(none),
!$OMP& shared (nad,nadto,mpl,mpiadlist,mpibuf,mpibuf_r,b),
!$OMP& private (iad,i,j,k,ipc,ib,nb,is,ie,nn)
      do iad=1,nad
        if(iad.eq.1) nn=0
        nb = (mpl(iad+1)-mpl(iad)-1)/DIV_N+1
        do ib=1,nb
           is(ib) = mpl(iad)+1+DIV_N*(ib-1)
           ie(ib) = min(mpl(iad+1),mpl(iad)+DIV_N*(ib))
        enddo
        call unpack_scalevalue(nb,
     $       mpibuf(3*npc*mpl(iad+1)+12*npc*nn+1
     $       +3*nadto*npc+12*npc*(nadto/DIV_N+nad)),
     $       mpibuf_r(1,1,1,1))
!$OMP DO
!OCL NORECURRENCE
        do ib=1,nb
        do i = is(ib), ie(ib)
          j=3*(i-1)*npc+12*npc*nn
     $          +3*nadto*npc+12*npc*(nadto/DIV_N+nad)
          k=3*(mpiadlist(i)-1)*npc
          do ipc=1,npc
c          b(k+ipc)=b(k+ipc)+mpibuf(j+ipc)
             b(k+npc*0+ipc) = b(k+npc*0+ipc)
     $            +real(mpibuf(j+npc*0+ipc))*mpibuf_r(ipc,1,1,ib)
     $            +mpibuf_r(ipc,1,2,ib)
             b(k+npc*1+ipc) = b(k+npc*1+ipc)
     $            +real(mpibuf(j+npc*1+ipc))*mpibuf_r(ipc,2,1,ib)
     $            +mpibuf_r(ipc,2,2,ib)
             b(k+npc*2+ipc) = b(k+npc*2+ipc)
     $            +real(mpibuf(j+npc*2+ipc))*mpibuf_r(ipc,3,1,ib)
     $            +mpibuf_r(ipc,3,2,ib)
          enddo
        enddo
        enddo
!$OMP END DO
        nn = nn+nb
      enddo
!$OMP END PARALLEL
#ifdef USE_PA
      call fapp_stop('sync_parallel_block_s_unpack', 0, 0)
#endif
#ifdef MEASURE_TIME_BARRIER
      call stop_measure(im,mp,timerbufd,timerbufi)
#endif
      end subroutine sync_parallel_block_s_part2

      subroutine sync_parallel_block_s_persistent_comm_init
     & (nad,nadto,mpiad,mpl,im,mpibuf,mpireq)
      use mpi
      implicit none
      integer, intent(in) :: im,nad,nadto,mpiad(nad,2),mpl(nad+1)
      integer*2, intent(inout) :: mpibuf(2*(3*npc*(nadto+1)
     $     +12*npc*(nadto/DIV_N+nad)))
      integer, intent(out) :: mpireq(nad,2)
      integer :: i,m1,m2,m3,itag,ierr
      integer is(nadto),ie(nadto),nn,nb,ib
      integer,save:: tagbase=1

      nn = 0
      do i=1,nad
        nb = (mpl(i+1)-mpl(i)-1)/DIV_N+1
        m1=3*npc*mpl(i)+1+12*npc*nn
     $        +3*nadto*npc+12*npc*(nadto/DIV_N+nad)
        m2=3*npc*mpiad(i,2)+12*npc*nb
        m3=mpiad(i,1)
c        itag=1 !m3*10000+im
        itag = tagbase
        call MPI_RECV_INIT
      !   call FJMPI_PREQUEST_RECV_INIT
     &   (mpibuf(m1),m2,MPI_INTEGER2,
     &    m3,itag,MPI_COMM_WORLD,mpireq(i,1),ierr)
        nn = nn+nb
      enddo
      nn = 0
      do i=1,nad
        nb = (mpl(i+1)-mpl(i)-1)/DIV_N+1
        m1=3*npc*mpl(i)+1+12*npc*nn
        m2=3*npc*mpiad(i,2)+12*npc*nb
        m3=mpiad(i,1)
c        itag=1 !im*10000+m3
        itag = tagbase
        call MPI_SEND_INIT
      !   call FJMPI_PREQUEST_SEND_INIT
     &   (mpibuf(m1),m2,MPI_INTEGER2,
     &    m3,itag,MPI_COMM_WORLD,mpireq(i,2),ierr)
        nn = nn+nb
      enddo
      tagbase = tagbase+1
      end subroutine sync_parallel_block_s_persistent_comm_init

      subroutine sync_parallel_persistent_comm_free(nad,mpireq)
      use mpi
      implicit none
      integer,intent(in) :: nad
      integer,intent(inout) :: mpireq(nad*2)
      integer :: i,ierr
      do i = 1, nad*2
        call MPI_REQUEST_FREE(mpireq(i),ierr)
      enddo
      end subroutine sync_parallel_persistent_comm_free
