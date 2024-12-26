      subroutine set_whole_CRS
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & nadto,mpiadlist,
     & overlap_comm,crssize)
      implicit none
      integer,intent(in) :: im,kd,n,ni,ne,nei,cnysize,cnysizei
      integer,intent(in) :: cny(cnysize,ne),num(ne)
      integer,intent(in) :: cnyi(cnysizei,nei),numi(nei)
      real(8),intent(in) :: younglst(2,kd),dsi,coori(3,n+ni)
      integer,intent(in) :: nfla(nei)
      integer,intent(in) :: crssize
      integer,intent(in) :: nadto,mpiadlist(nadto)
      logical,intent(in) :: overlap_comm
      integer,intent(out) :: pentawptr(2*nsta+3),pentarptr((ni+n)*4)
      integer,intent(out) :: pentawind((n+ni)*4)
      integer,intent(out) :: pentarind((ni+n)*crssize)
      integer,intent(out) :: niw,nir
#if npc == 1 || npc == 21
      real(8),intent(out) :: crsval(9,(n+ni)*crssize)
      if(overlap_comm)then
        if(im==0)then
          write(6,*) "overlap_comm is not implemented ",
     &               "for set_whole_CRS_KplusG"
          stop
        endif
      endif
      call set_whole_CRS_KplusG
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & crssize)
#elif npc == 288
      real(8),intent(out) :: crsval(9,(n+ni)*crssize*2)
      if(overlap_comm)then
      call set_whole_CRS_KandG_overlap_comm
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & nadto,mpiadlist,
     & crssize)
      else
      call set_whole_CRS_KandG
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & crssize)
      endif
#else
ERROR: invalid nkl
#endif
      end subroutine set_whole_CRS

      subroutine set_whole_CRS_KplusG
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & crssize)
      implicit none
      integer,intent(in) :: im,kd,n,ni,ne,nei,cnysize,cnysizei
      integer,intent(in) :: cny(cnysize,ne),num(ne)
      integer,intent(in) :: cnyi(cnysizei,nei),numi(nei)
      real(8),intent(in) :: younglst(2,kd),dsi,coori(3,n+ni)
      integer,intent(in) :: nfla(nei)
      integer,intent(in) :: crssize
      integer,intent(out) :: pentawptr(nsta+2),pentarptr((ni+n)*4)
      integer,intent(out) :: pentawind((n+ni)*4)
      integer,intent(out) :: pentarind((ni+n)*crssize)
      integer,intent(out) :: niw,nir
      real(8),intent(out) :: crsval(9,(n+ni)*crssize)
c
      integer numtmp(ne),numitmp(nei)
      integer pentasize(n+ni)
      integer pentaind(crssize,n+ni)
      integer i,j,i1,j1,i2,j2,k,ie,ista,nmax,ii
      real(8) ket(3*cnysize,3*cnysize),young,rnyu,cori(3,3),cor(4,3)
      real(8) gxy,gyz,gxz,detj
      real(8) keit(3*cnysizei,3*cnysizei)

!     initialization, preparation
      do ie=1,ne
        if(num(ie).gt.nsta)then
          numtmp(ie)=nsta+1
        else
          numtmp(ie)=num(ie)
        endif
      enddo ! ie

      do ie=1,nei
        if(numi(ie).gt.nsta)then
          numitmp(ie)=nsta+1
        else
          numitmp(ie)=numi(ie)
        endif
      enddo ! ie

      nir=0
      niw=0
      pentawptr=0
      pentarptr=0
      nmax=0

      do ista=1,nsta+1

!      write(*,*) 'crs table construction',ista,nsta+1

      pentasize=0
      pentaind=-1
!     1-1. pentaind, pentasize for normal elements
      do ie=1,ne
        if(numtmp(ie).eq.ista)then
          do j1=1,cnysize
            i1=cny(j1,ie)
            do j2=1,cnysize
              i2=cny(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 315
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE 1',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
315           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

!     1-2. pentaind, pentasize for IEpenta elements
      do ie=1,nei
        if(numitmp(ie).eq.ista)then
          do j1=1,cnysizei
            i1=cnyi(j1,ie)
            do j2=1,cnysizei
              i2=cnyi(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 314
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE 2',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
314           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

!     2. pentarind, pentarptr, pentawind, pntwptr
      do i=1,n+ni
        if(pentasize(i).ne.0)then
          niw=niw+1
          pentawind(niw)=i
          call sort_bubble(pentasize(i),pentaind(1,i))
          do j=1,pentasize(i)
            nir=nir+1
            pentarind(nir)=pentaind(j,i)
          enddo
          pentarptr(niw+1)=nir
        endif
      enddo
      pentawptr(ista+1)=niw

      do i=1,n+ni
      if(pentasize(i).gt.nmax)then
      nmax=pentasize(i)
      endif
      enddo

      enddo ! ista


      do i=1,nir
      do j=1,9
      crsval(j,i)=0.0
      enddo
      enddo

      ! check crs table
      do ista=1,nsta+1
!     3-1. crsval for normal elements
      do ie=1,ne
        if(numtmp(ie).eq.ista)then

          young=younglst(1,num(ie))
          rnyu=younglst(2,num(ie))
          gxy=young/(2.*(1+rnyu))
          gyz=young/(2.*(1+rnyu))
          gxz=young/(2.*(1+rnyu))
          do i1=1,4
          do ii=1,3
          cor(i1,ii)=coori(ii,cny(i1,ie))
          enddo
          enddo
          if(cnysize.eq.10)then
            call nontet10(young,gxy,gyz,gxz,rnyu,cor,ket,detj)
          else ! tet4 elem
            call nontetra4(young,gxy,gyz,gxz,rnyu,cor,ket)
          endif

          do j1=1,cnysize
            i1=cny(j1,ie)
            do i=pentawptr(ista)+1,pentawptr(ista+1)
              if(pentawind(i).eq.i1)then
                do j2=1,cnysize
                  i2=cny(j2,ie)
                  do k=pentarptr(i)+1,pentarptr(i+1)
                    if(pentarind(k).eq.i2)then
                    crsval(1,k)=crsval(1,k)+ket(3*(j1-1)+1,3*(j2-1)+1)
                    crsval(2,k)=crsval(2,k)+ket(3*(j1-1)+1,3*(j2-1)+2)
                    crsval(3,k)=crsval(3,k)+ket(3*(j1-1)+1,3*(j2-1)+3)
                    crsval(4,k)=crsval(4,k)+ket(3*(j1-1)+2,3*(j2-1)+1)
                    crsval(5,k)=crsval(5,k)+ket(3*(j1-1)+2,3*(j2-1)+2)
                    crsval(6,k)=crsval(6,k)+ket(3*(j1-1)+2,3*(j2-1)+3)
                    crsval(7,k)=crsval(7,k)+ket(3*(j1-1)+3,3*(j2-1)+1)
                    crsval(8,k)=crsval(8,k)+ket(3*(j1-1)+3,3*(j2-1)+2)
                    crsval(9,k)=crsval(9,k)+ket(3*(j1-1)+3,3*(j2-1)+3)
                    goto 144
                    endif
                  enddo ! k
                  write(*,*) 'pentaw not found 1',i1,i2,ie
                  stop
144               k=0
                enddo ! j2
                goto 143
              endif
            enddo ! i
            write(*,*) 'pentaw not found 1',i1
            stop
143         i=0
          enddo ! j1
        endif
      enddo ! ie

!     3-2. crsval for normal IEpenta elements
      do ie=1,nei
        if(numitmp(ie).eq.ista)then

          young=younglst(1,numi(ie))
          rnyu=younglst(2,numi(ie))
          do i1=1,3
          do ii=1,3
          cori(i1,ii)=coori(ii,cnyi(i1,ie))
          enddo
          enddo
          if(cnysizei.eq.12)then
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
          else ! tet4 elem
            if(nfla(ie).eq.1)then
            call IEpenta1st1(young,rnyu,cori,dsi,keit)
            endif
            if(nfla(ie).eq.2)then
            call IEpenta1st2(young,rnyu,cori,dsi,keit)
            endif
            if(nfla(ie).eq.3)then
            call IEpenta1st3(young,rnyu,cori,dsi,keit)
            endif
            if(nfla(ie).eq.4)then
            call IEpenta1st4(young,rnyu,cori,dsi,keit)
            endif
            if(nfla(ie).eq.5)then
            call IEpenta1st5(young,rnyu,cori,dsi,keit)
            endif
          endif

          do j1=1,cnysizei
            i1=cnyi(j1,ie)
            do i=pentawptr(ista)+1,pentawptr(ista+1)
              if(pentawind(i).eq.i1)then
                do j2=1,cnysizei
                  i2=cnyi(j2,ie)
                  do k=pentarptr(i)+1,pentarptr(i+1)
                    if(pentarind(k).eq.i2)then
                    crsval(1,k)=crsval(1,k)+keit(3*(j1-1)+1,3*(j2-1)+1)
                    crsval(2,k)=crsval(2,k)+keit(3*(j1-1)+1,3*(j2-1)+2)
                    crsval(3,k)=crsval(3,k)+keit(3*(j1-1)+1,3*(j2-1)+3)
                    crsval(4,k)=crsval(4,k)+keit(3*(j1-1)+2,3*(j2-1)+1)
                    crsval(5,k)=crsval(5,k)+keit(3*(j1-1)+2,3*(j2-1)+2)
                    crsval(6,k)=crsval(6,k)+keit(3*(j1-1)+2,3*(j2-1)+3)
                    crsval(7,k)=crsval(7,k)+keit(3*(j1-1)+3,3*(j2-1)+1)
                    crsval(8,k)=crsval(8,k)+keit(3*(j1-1)+3,3*(j2-1)+2)
                    crsval(9,k)=crsval(9,k)+keit(3*(j1-1)+3,3*(j2-1)+3)
                    goto 142
                    endif
                  enddo ! k
                  write(*,*) 'pentaw not found 2',i1,ie
                  stop
142               k=0
                enddo ! j2
                goto 141
              endif
            enddo ! i
            write(*,*) 'pentaw not found 2',i1
            stop
141         i=0
          enddo ! j1
        endif
      enddo ! ie
      enddo ! ista
!      write(*,*) 'check crs table succeed'

      end subroutine set_whole_CRS_KplusG

c_______________________________________________________________________
      subroutine set_whole_CRS_KandG
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & crssize)
      implicit none
      integer,intent(in) :: im,kd,n,ni,ne,nei,cnysize,cnysizei
      integer,intent(in) :: cny(cnysize,ne),num(ne)
      integer,intent(in) :: cnyi(cnysizei,nei),numi(nei)
      real(8),intent(in) :: younglst(2,kd),dsi,coori(3,n+ni)
      integer,intent(in) :: nfla(nei)
      integer,intent(in) :: crssize
      integer,intent(out) :: pentawptr(nsta+2),pentarptr((ni+n)*4)
      integer,intent(out) :: pentawind((n+ni)*4)
      integer,intent(out) :: pentarind((ni+n)*crssize)
      integer,intent(out) :: niw,nir
      real(8),intent(out) :: crsval(9,(n+ni)*crssize*2)
c
      integer numtmp(ne),numitmp(nei)
      integer pentasize(n+ni)
      integer pentaind(crssize,n+ni)
      integer i,j,i1,j1,i2,j2,i3,k,ie,ista,nmax,ii,ikg
      real(8) ket(3*cnysize,3*cnysize,2),young,rnyu,cori(3,3),cor(4,3)
      real(8) gxy,gyz,gxz,rKc,rGc
      real(8) keit(3*cnysizei,3*cnysizei,2)

!     initialization, preparation
      do ie=1,ne
        if(num(ie).gt.nsta)then
          numtmp(ie)=nsta+1
        else
          numtmp(ie)=num(ie)
        endif
      enddo ! ie

      do ie=1,nei
        if(numi(ie).gt.nsta)then
          numitmp(ie)=nsta+1
        else
          numitmp(ie)=numi(ie)
        endif
      enddo ! ie

      nir=0
      niw=0
      pentawptr=0
      pentarptr=0
      nmax=0

      do ista=1,nsta+1

!      write(*,*) 'crs table construction',ista,nsta+1

      pentasize=0
      pentaind=-1
!     1-1. pentaind, pentasize for normal elements
      do ie=1,ne
        if(numtmp(ie).eq.ista)then
          do j1=1,cnysize
            i1=cny(j1,ie)
            do j2=1,cnysize
              i2=cny(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 315
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE 1',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
315           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

!     1-2. pentaind, pentasize for IEpenta elements
      do ie=1,nei
        if(numitmp(ie).eq.ista)then
          do j1=1,cnysizei
            i1=cnyi(j1,ie)
            do j2=1,cnysizei
              i2=cnyi(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 314
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE 2',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
314           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

!     2. pentarind, pentarptr, pentawind, pntwptr
      do i=1,n+ni
        if(pentasize(i).ne.0)then
          niw=niw+1
          pentawind(niw)=i
          call sort_bubble(pentasize(i),pentaind(1,i))
          do j=1,pentasize(i)
            nir=nir+1
            pentarind(nir)=pentaind(j,i)
          enddo
          pentarptr(niw+1)=nir
        endif
      enddo
      pentawptr(ista+1)=niw

      do i=1,n+ni
      if(pentasize(i).gt.nmax)then
      nmax=pentasize(i)
      endif
      enddo

      enddo ! ista


      do i=1,nir*2
      do j=1,9
      crsval(j,i)=0.0
      enddo
      enddo

      ! check crs table
      do ista=1,nsta+1
!     3-1. crsval for normal elements
      do ie=1,ne
        if(numtmp(ie).eq.ista)then

          young=younglst(1,num(ie))
          rnyu=younglst(2,num(ie))
          gxy=young/(2.*(1+rnyu))
          gyz=young/(2.*(1+rnyu))
          gxz=young/(2.*(1+rnyu))
          rKc=young/3/(1-2*rnyu)
          rGc=young/2/(1+rnyu)
          do i1=1,4
          do ii=1,3
          cor(i1,ii)=coori(ii,cny(i1,ie))
          enddo
          enddo
          if(cnysize.eq.10)then
            call nontet10_KG(rKc,0d0,cor,ket(:,:,1))
            call nontet10_KG(0d0,rGc,cor,ket(:,:,2))
          else ! tet4 elem
            call nontetra4_KG(rKc,0d0,cor,ket(:,:,1))
            call nontetra4_KG(0d0,rGc,cor,ket(:,:,2))
          endif

          do j1=1,cnysize
            i1=cny(j1,ie)
            do i=pentawptr(ista)+1,pentawptr(ista+1)
              if(pentawind(i).eq.i1)then
                do j2=1,cnysize
                  i2=cny(j2,ie)
                  do k=pentarptr(i)+1,pentarptr(i+1)
                    if(pentarind(k).eq.i2)then
                      do ikg = 1, 2
                        i3 = (ikg-1)*nir
                        crsval(1,k+i3)=crsval(1,k+i3)
     &                    +ket(3*(j1-1)+1,3*(j2-1)+1,ikg)
                        crsval(2,k+i3)=crsval(2,k+i3)
     &                    +ket(3*(j1-1)+1,3*(j2-1)+2,ikg)
                        crsval(3,k+i3)=crsval(3,k+i3)
     &                    +ket(3*(j1-1)+1,3*(j2-1)+3,ikg)
                        crsval(4,k+i3)=crsval(4,k+i3)
     &                    +ket(3*(j1-1)+2,3*(j2-1)+1,ikg)
                        crsval(5,k+i3)=crsval(5,k+i3)
     &                    +ket(3*(j1-1)+2,3*(j2-1)+2,ikg)
                        crsval(6,k+i3)=crsval(6,k+i3)
     &                    +ket(3*(j1-1)+2,3*(j2-1)+3,ikg)
                        crsval(7,k+i3)=crsval(7,k+i3)
     &                    +ket(3*(j1-1)+3,3*(j2-1)+1,ikg)
                        crsval(8,k+i3)=crsval(8,k+i3)
     &                    +ket(3*(j1-1)+3,3*(j2-1)+2,ikg)
                        crsval(9,k+i3)=crsval(9,k+i3)
     &                    +ket(3*(j1-1)+3,3*(j2-1)+3,ikg)
                      enddo
                      goto 144
                    endif
                  enddo ! k
                  write(*,*) 'pentaw not found 1',i1,i2,ie
                  stop
144               k=0
                enddo ! j2
                goto 143
              endif
            enddo ! i
            write(*,*) 'pentaw not found 1',i1
            stop
143         i=0
          enddo ! j1
        endif
      enddo ! ie

!     3-2. crsval for normal IEpenta elements
      do ie=1,nei
        if(numitmp(ie).eq.ista)then

          young=younglst(1,numi(ie))
          rnyu=younglst(2,numi(ie))
          rKc=young/3/(1-2*rnyu)
          rGc=young/2/(1+rnyu)
          do i1=1,3
          do ii=1,3
          cori(i1,ii)=coori(ii,cnyi(i1,ie))
          enddo
          enddo
          if(cnysizei.eq.12)then
            if(nfla(ie).eq.1)then
            call IEpenta2nd1_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd1_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.2)then
            call IEpenta2nd2_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd2_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.3)then
            call IEpenta2nd3_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd3_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.4)then
            call IEpenta2nd4_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd4_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.5)then
            call IEpenta2nd5_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd5_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
          else ! tet4 elem
            if(nfla(ie).eq.1)then
            call IEpenta1st1_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st1_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.2)then
            call IEpenta1st2_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st2_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.3)then
            call IEpenta1st3_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st3_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.4)then
            call IEpenta1st4_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st4_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.5)then
            call IEpenta1st5_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st5_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
          endif

          do j1=1,cnysizei
            i1=cnyi(j1,ie)
            do i=pentawptr(ista)+1,pentawptr(ista+1)
              if(pentawind(i).eq.i1)then
                do j2=1,cnysizei
                  i2=cnyi(j2,ie)
                  do k=pentarptr(i)+1,pentarptr(i+1)
                    if(pentarind(k).eq.i2)then
                      do ikg = 1, 2
                        i3 = nir*(ikg-1)
                        crsval(1,k+i3)=crsval(1,k+i3)
     &                    + keit(3*(j1-1)+1,3*(j2-1)+1,ikg)
                        crsval(2,k+i3)=crsval(2,k+i3)
     &                    + keit(3*(j1-1)+1,3*(j2-1)+2,ikg)
                        crsval(3,k+i3)=crsval(3,k+i3)
     &                    + keit(3*(j1-1)+1,3*(j2-1)+3,ikg)
                        crsval(4,k+i3)=crsval(4,k+i3)
     &                    + keit(3*(j1-1)+2,3*(j2-1)+1,ikg)
                        crsval(5,k+i3)=crsval(5,k+i3)
     &                    + keit(3*(j1-1)+2,3*(j2-1)+2,ikg)
                        crsval(6,k+i3)=crsval(6,k+i3)
     &                    + keit(3*(j1-1)+2,3*(j2-1)+3,ikg)
                        crsval(7,k+i3)=crsval(7,k+i3)
     &                    + keit(3*(j1-1)+3,3*(j2-1)+1,ikg)
                        crsval(8,k+i3)=crsval(8,k+i3)
     &                    + keit(3*(j1-1)+3,3*(j2-1)+2,ikg)
                        crsval(9,k+i3)=crsval(9,k+i3)
     &                    + keit(3*(j1-1)+3,3*(j2-1)+3,ikg)
                      enddo
                      goto 142
                    endif
                  enddo ! k
                  write(*,*) 'pentaw not found 2',i1,ie
                  stop
142               k=0
                enddo ! j2
                goto 141
              endif
            enddo ! i
            write(*,*) 'pentaw not found 2',i1
            stop
141         i=0
          enddo ! j1
        endif
      enddo ! ie
      enddo ! ista
!      write(*,*) 'check crs table succeed'
      end subroutine set_whole_CRS_KandG


c_______________________________________________________________________
      subroutine set_whole_CRS_KandG_overlap_comm
     & (im,kd,n,ni,ne,nei,cnysize,cnysizei,cny,cnyi,num,numi,
     & nfla,younglst,coori,dsi,
     & niwo,pentawptro,pentawindo,niro,pentarptro,pentarindo,crsval,
     & nadto,mpiadlist,
     & crssize)
      implicit none
      integer,intent(in) :: im,kd,n,ni,ne,nei,cnysize,cnysizei
      integer,intent(in) :: cny(cnysize,ne),num(ne)
      integer,intent(in) :: cnyi(cnysizei,nei),numi(nei)
      real(8),intent(in) :: younglst(2,kd),dsi,coori(3,n+ni)
      integer,intent(in) :: nfla(nei)
      integer,intent(in) :: nadto,mpiadlist(nadto)
      integer,intent(in) :: crssize
      integer,intent(out) :: pentawptro(2*nsta+3),pentarptro((ni+n)*4)
      integer,intent(out) :: pentawindo((n+ni)*4)
      integer,intent(out) :: pentarindo((ni+n)*crssize)
      integer,intent(out) :: niwo,niro
      real(8),intent(out) :: crsval(9,(n+ni)*crssize*2)
c
      integer :: niw,nir
      integer :: pentawptr(nsta+2),pentarptr((ni+n)*4)
      integer :: pentawind((n+ni)*4)
      integer :: pentarind((ni+n)*crssize)
c
      integer numtmp(ne),numitmp(nei)
      integer pentasize(n+ni)
      integer pentaind(crssize,n+ni)
      integer i,j,i1,j1,i2,j2,i3,k,ie,ista,nmax,ii,ikg
      real(8) ket(3*cnysize,3*cnysize,2),young,rnyu,cori(3,3),cor(4,3)
      real(8) gxy,gyz,gxz,rKc,rGc
      real(8) keit(3*cnysizei,3*cnysizei,2)

      integer icomm, istao
      integer node_flag(n+ni)
      integer, parameter :: NOCOMM = 2, COMM = 1

!     initialization, preparation
      do ie=1,ne
        if(num(ie).gt.nsta)then
          numtmp(ie)=nsta+1
        else
          numtmp(ie)=num(ie)
        endif
      enddo ! ie

      do ie=1,nei
        if(numi(ie).gt.nsta)then
          numitmp(ie)=nsta+1
        else
          numitmp(ie)=numi(ie)
        endif
      enddo ! ie

      nir=0
      niw=0
      pentawptr=0
      pentarptr=0
      nmax=0

      do ista=1,nsta+1

!      write(*,*) 'crs table construction',ista,nsta+1

      pentasize=0
      pentaind=-1
!     1-1. pentaind, pentasize for normal elements
      do ie=1,ne
        if(numtmp(ie).eq.ista)then
          do j1=1,cnysize
            i1=cny(j1,ie)
            do j2=1,cnysize
              i2=cny(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 315
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE 1',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
315           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

!     1-2. pentaind, pentasize for IEpenta elements
      do ie=1,nei
        if(numitmp(ie).eq.ista)then
          do j1=1,cnysizei
            i1=cnyi(j1,ie)
            do j2=1,cnysizei
              i2=cnyi(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 314
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE 2',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
314           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

!     2-1. pentarind, pentarptr, pentawind, pntwptr
      do i=1,n+ni
        if(pentasize(i).ne.0)then
          niw=niw+1
          pentawind(niw)=i
          call sort_bubble(pentasize(i),pentaind(1,i))
          do j=1,pentasize(i)
            nir=nir+1
            pentarind(nir)=pentaind(j,i)
          enddo
          pentarptr(niw+1)=nir
        endif
      enddo
      pentawptr(ista+1)=niw

      do i=1,n+ni
      if(pentasize(i).gt.nmax)then
      nmax=pentasize(i)
      endif
      enddo

      enddo ! ista

!     2-2. reorder pentarind, pentarptr, pentawind, pntwptr
      ! NOCOMM = 2, COMM = 1
      node_flag(:) = NOCOMM
      do j = 1, nadto
        i = mpiadlist(j)
        node_flag(i) = COMM
      enddo

      niwo = 0
      niro = 0
      pentawptro(1) = 0
      pentarptro(1) = 0
      istao = 0
      do icomm = 1, 2 ! icomm = NOCOMM or COMM
      do ista = 1, nsta + 1
        do j1 = pentawptr(ista)+1, pentawptr(ista+1)
          i1 = pentawind(j1)
          if (node_flag(i1) /= icomm) cycle

          do j2 = pentarptr(j1)+1, pentarptr(j1+1)
            niro = niro + 1
            pentarindo(niro) = pentarind(j2)
          enddo
          niwo = niwo + 1
          pentawindo(niwo) = i1
          pentarptro(niwo+1) = niro
        enddo ! j1
        istao = istao + 1
        pentawptro(istao+1) = niwo

      enddo ! ista
      enddo ! icomm
c
      do i=1,nir*2
      do j=1,9
      crsval(j,i)=0.0
      enddo
      enddo

      ! check crs table
      do ista=1,nsta+1
!     3-1. crsval for normal elements
      do ie=1,ne
        if(numtmp(ie).eq.ista)then

          young=younglst(1,num(ie))
          rnyu=younglst(2,num(ie))
          gxy=young/(2.*(1+rnyu))
          gyz=young/(2.*(1+rnyu))
          gxz=young/(2.*(1+rnyu))
          rKc=young/3/(1-2*rnyu)
          rGc=young/2/(1+rnyu)
          do i1=1,4
          do ii=1,3
          cor(i1,ii)=coori(ii,cny(i1,ie))
          enddo
          enddo
          if(cnysize.eq.10)then
            call nontet10_KG(rKc,0d0,cor,ket(:,:,1))
            call nontet10_KG(0d0,rGc,cor,ket(:,:,2))
          else ! tet4 elem
            call nontetra4_KG(rKc,0d0,cor,ket(:,:,1))
            call nontetra4_KG(0d0,rGc,cor,ket(:,:,2))
          endif

          do j1=1,cnysize
            i1=cny(j1,ie)
            icomm = node_flag(i1)
            istao = ista + (icomm-1)*(nsta+1)
            do i=pentawptro(istao)+1,pentawptro(istao+1)
              if(pentawindo(i).eq.i1)then
                do j2=1,cnysize
                  i2=cny(j2,ie)
                  do k=pentarptro(i)+1,pentarptro(i+1)
                    if(pentarindo(k).eq.i2)then
                      do ikg = 1, 2
                        i3 = (ikg-1)*niro
                        crsval(1,k+i3)=crsval(1,k+i3)
     &                    +ket(3*(j1-1)+1,3*(j2-1)+1,ikg)
                        crsval(2,k+i3)=crsval(2,k+i3)
     &                    +ket(3*(j1-1)+1,3*(j2-1)+2,ikg)
                        crsval(3,k+i3)=crsval(3,k+i3)
     &                    +ket(3*(j1-1)+1,3*(j2-1)+3,ikg)
                        crsval(4,k+i3)=crsval(4,k+i3)
     &                    +ket(3*(j1-1)+2,3*(j2-1)+1,ikg)
                        crsval(5,k+i3)=crsval(5,k+i3)
     &                    +ket(3*(j1-1)+2,3*(j2-1)+2,ikg)
                        crsval(6,k+i3)=crsval(6,k+i3)
     &                    +ket(3*(j1-1)+2,3*(j2-1)+3,ikg)
                        crsval(7,k+i3)=crsval(7,k+i3)
     &                    +ket(3*(j1-1)+3,3*(j2-1)+1,ikg)
                        crsval(8,k+i3)=crsval(8,k+i3)
     &                    +ket(3*(j1-1)+3,3*(j2-1)+2,ikg)
                        crsval(9,k+i3)=crsval(9,k+i3)
     &                    +ket(3*(j1-1)+3,3*(j2-1)+3,ikg)
                      enddo
                      goto 144
                    endif
                  enddo ! k
                  write(*,*) 'pentaw not found 1',i1,i2,ie
                  stop
144               k=0
                enddo ! j2
                goto 143
              endif
            enddo ! i
            write(*,*) 'pentaw not found 1',i1
            stop
143         i=0
          enddo ! j1
        endif
      enddo ! ie

!     3-2. crsval for normal IEpenta elements
      do ie=1,nei
        if(numitmp(ie).eq.ista)then

          young=younglst(1,numi(ie))
          rnyu=younglst(2,numi(ie))
          rKc=young/3/(1-2*rnyu)
          rGc=young/2/(1+rnyu)
          do i1=1,3
          do ii=1,3
          cori(i1,ii)=coori(ii,cnyi(i1,ie))
          enddo
          enddo
          if(cnysizei.eq.12)then
            if(nfla(ie).eq.1)then
            call IEpenta2nd1_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd1_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.2)then
            call IEpenta2nd2_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd2_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.3)then
            call IEpenta2nd3_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd3_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.4)then
            call IEpenta2nd4_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd4_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.5)then
            call IEpenta2nd5_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta2nd5_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
          else ! tet4 elem
            if(nfla(ie).eq.1)then
            call IEpenta1st1_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st1_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.2)then
            call IEpenta1st2_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st2_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.3)then
            call IEpenta1st3_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st3_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.4)then
            call IEpenta1st4_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st4_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
            if(nfla(ie).eq.5)then
            call IEpenta1st5_KG(rKc,0d0,cori,dsi,keit(:,:,1))
            call IEpenta1st5_KG(0d0,rGc,cori,dsi,keit(:,:,2))
            endif
          endif

          do j1=1,cnysizei
            i1=cnyi(j1,ie)
            icomm = node_flag(i1)
            istao = ista + (icomm-1)*(nsta+1)
            do i=pentawptro(istao)+1,pentawptro(istao+1)
              if(pentawindo(i).eq.i1)then
                do j2=1,cnysizei
                  i2=cnyi(j2,ie)
                  do k=pentarptro(i)+1,pentarptro(i+1)
                    if(pentarindo(k).eq.i2)then
                      do ikg = 1, 2
                        i3 = niro*(ikg-1)
                        crsval(1,k+i3)=crsval(1,k+i3)
     &                    + keit(3*(j1-1)+1,3*(j2-1)+1,ikg)
                        crsval(2,k+i3)=crsval(2,k+i3)
     &                    + keit(3*(j1-1)+1,3*(j2-1)+2,ikg)
                        crsval(3,k+i3)=crsval(3,k+i3)
     &                    + keit(3*(j1-1)+1,3*(j2-1)+3,ikg)
                        crsval(4,k+i3)=crsval(4,k+i3)
     &                    + keit(3*(j1-1)+2,3*(j2-1)+1,ikg)
                        crsval(5,k+i3)=crsval(5,k+i3)
     &                    + keit(3*(j1-1)+2,3*(j2-1)+2,ikg)
                        crsval(6,k+i3)=crsval(6,k+i3)
     &                    + keit(3*(j1-1)+2,3*(j2-1)+3,ikg)
                        crsval(7,k+i3)=crsval(7,k+i3)
     &                    + keit(3*(j1-1)+3,3*(j2-1)+1,ikg)
                        crsval(8,k+i3)=crsval(8,k+i3)
     &                    + keit(3*(j1-1)+3,3*(j2-1)+2,ikg)
                        crsval(9,k+i3)=crsval(9,k+i3)
     &                    + keit(3*(j1-1)+3,3*(j2-1)+3,ikg)
                      enddo
                      goto 142
                    endif
                  enddo ! k
                  write(*,*) 'pentaw not found 2',i1,ie
                  stop
142               k=0
                enddo ! j2
                goto 141
              endif
            enddo ! i
            write(*,*) 'pentaw not found 2',i1
            stop
141         i=0
          enddo ! j1
        endif
      enddo ! ie
      enddo ! ista
!      write(*,*) 'check crs table succeed'
      end subroutine set_whole_CRS_KandG_overlap_comm
