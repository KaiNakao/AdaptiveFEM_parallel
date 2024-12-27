      subroutine set_ni(n_size,n,coor,plane,ds,ni)
      implicit none
      integer n_size,n,ni
      real*8 coor(3,n_size),ds,plane(5)
      integer i

      ni=0
      do i=1,n
      if(abs(coor(1,i)-plane(1)).le.ds*0.01)then
      ni=ni+1
      endif
      enddo
      do i=1,n
      if(abs(coor(1,i)-plane(2)).le.ds*0.01)then
      ni=ni+1
      endif
      enddo
      do i=1,n
      if(abs(coor(2,i)-plane(3)).le.ds*0.01)then
      ni=ni+1
      endif
      enddo
      do i=1,n
      if(abs(coor(2,i)-plane(4)).le.ds*0.01)then
      ni=ni+1
      endif
      enddo
      do i=1,n
      if(abs(coor(3,i)-plane(5)).le.ds*0.01)then
      ni=ni+1
      endif
      enddo
      end
c_______________________________________________________________________
      subroutine setIEpenta
     & (im,n,ne,cny,num,cnyi,numi,nfla,ni,nad,plane,
     & nei,coor,coori,ds,
     & nadto,mpiad,mpiadlist,mpl,
     & nadto1,mpiad1,mpiadlist1,mpl1,itype,itypei)
      implicit none
      integer im,n,ne,nad,itype,itypei
      integer cny(itype,ne),num(ne),cnyi(itypei,nei),numi(nei)
      integer nfla(nei)
      integer ni,nei
      real*8 plane(5),ds
      real*8 coor(3,n),coori(3,n+ni)
      integer nadto,mpiad(nad,2),mpiadlist(nadto),mpl(nad+1)
      integer nadto1,mpiad1(nad,2),mpiadlist1(2*nadto),mpl1(nad+1)

      integer nadto2,mpiad2(nad,2),mpiadlist2(2*nadto),mpl2(nad+1)
      integer idcoor(5)
      integer inum,iecount,i,id,iad,ie
      integer flag(n)

      idcoor(1)=1
      idcoor(2)=1
      idcoor(3)=2
      idcoor(4)=2
      idcoor(5)=3
      inum=n
      iecount=0
      nfla=0
      mpiadlist1=-1
      mpiadlist2=-1
      mpiad1=mpiad
      do i=1,nadto
      mpiadlist1(i)=mpiadlist(i)
      enddo
      mpl1=mpl
c
      do id=1,5
c
      flag=0
      do i=1,n
      if(abs(coor(idcoor(id),i)-plane(id)).le.ds*0.01)then
      inum=inum+1
      flag(i)=inum
      coori(1:3, inum) = coor(1:3,i)
      endif
      enddo
      mpiad2=mpiad1
      mpl2(1)=mpl1(1)
      do iad=1,nad
      mpl2(iad+1)=mpl2(iad)
      ! copy nodes in mpl
      do i=mpl1(iad)+1,mpl1(iad+1)
        mpl2(iad+1)=mpl2(iad+1)+1 ! shift pointer
        mpiadlist2(mpl2(iad+1))=mpiadlist1(i) ! set index
      enddo
      ! insert bc nodes
      do i=mpl(iad)+1,mpl(iad+1)
        if(flag(mpiadlist(i)).ne.0)then
          mpiad2(iad,2)=mpiad2(iad,2)+1 ! enlarge size
          mpl2(iad+1)=mpl2(iad+1)+1 ! shift pointer
          mpiadlist2(mpl2(iad+1))=flag(mpiadlist(i)) ! set index
        endif
      enddo
      enddo
      mpiad1=mpiad2
      mpiadlist1=mpiadlist2
      mpl1=mpl2
      if(itypei.eq.12)then
      call tetIE
     & (n,ne,cny,num,flag,id,cnyi,numi,nfla,ni,
     & nei,iecount,coor,coori)
      else
      call tetIE4
     & (n,ne,cny,num,flag,id,cnyi,numi,nfla,ni,
     & nei,iecount,coor,coori)
      endif
c
      enddo
      nadto1=mpl1(nad+1)
      do ie=1,nei
      if(nfla(ie).eq.0)then
      write(*,*) 'something wrong with nfla'
      stop
      endif
      enddo
      end
c_______________________________________________________________________
      subroutine tetIE
     & (n,ne,cny,num,flag,id,cnyi,numi,nfla,ni,
     & nei,iecount,coor,coori)
      implicit none
      integer cny(10,ne),num(ne)
      integer flag(n),cnynl(4),id,iecount
      integer cnyi(12,nei),numi(nei),nfla(nei)
      integer n,ne,ni,nei
      integer iee,ie,ii,i1,i2,i3,i4,i5,i6,j
      real*8 coor(3,n),coori(3,n+ni)
c
#ifndef USEIEPENTA
      return
#endif

      do ie=1,ne
c
      cnynl=0
      do ii=1,4
      if(flag(cny(ii,ie)).gt.0)then
      cnynl(ii)=1
      endif
      enddo
c

      i1=1
      i2=2
      i3=4
      i4=5
      i5=9
      i6=8
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=cny(i4,ie)
      cnyi(5,iee)=cny(i5,ie)
      cnyi(6,iee)=cny(i6,ie)
      cnyi(7,iee)=flag(cny(i1,ie))
      cnyi(8,iee)=flag(cny(i2,ie))
      cnyi(9,iee)=flag(cny(i3,ie))
      cnyi(10,iee)=flag(cny(i4,ie))
      cnyi(11,iee)=flag(cny(i5,ie))
      cnyi(12,iee)=flag(cny(i6,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! coori(j,flag(cny(i4,ie)))=coor(j,cny(i4,ie))
      ! coori(j,flag(cny(i5,ie)))=coor(j,cny(i5,ie))
      ! coori(j,flag(cny(i6,ie)))=coor(j,cny(i6,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif
      i1=2
      i2=3
      i3=4
      i4=6
      i5=10
      i6=9
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=cny(i4,ie)
      cnyi(5,iee)=cny(i5,ie)
      cnyi(6,iee)=cny(i6,ie)
      cnyi(7,iee)=flag(cny(i1,ie))
      cnyi(8,iee)=flag(cny(i2,ie))
      cnyi(9,iee)=flag(cny(i3,ie))
      cnyi(10,iee)=flag(cny(i4,ie))
      cnyi(11,iee)=flag(cny(i5,ie))
      cnyi(12,iee)=flag(cny(i6,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! coori(j,flag(cny(i4,ie)))=coor(j,cny(i4,ie))
      ! coori(j,flag(cny(i5,ie)))=coor(j,cny(i5,ie))
      ! coori(j,flag(cny(i6,ie)))=coor(j,cny(i6,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif
      i1=3
      i2=1
      i3=4
      i4=7
      i5=8
      i6=10
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=cny(i4,ie)
      cnyi(5,iee)=cny(i5,ie)
      cnyi(6,iee)=cny(i6,ie)
      cnyi(7,iee)=flag(cny(i1,ie))
      cnyi(8,iee)=flag(cny(i2,ie))
      cnyi(9,iee)=flag(cny(i3,ie))
      cnyi(10,iee)=flag(cny(i4,ie))
      cnyi(11,iee)=flag(cny(i5,ie))
      cnyi(12,iee)=flag(cny(i6,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! coori(j,flag(cny(i4,ie)))=coor(j,cny(i4,ie))
      ! coori(j,flag(cny(i5,ie)))=coor(j,cny(i5,ie))
      ! coori(j,flag(cny(i6,ie)))=coor(j,cny(i6,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif
      i1=2
      i2=1
      i3=3
      i4=5
      i5=7
      i6=6
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=cny(i4,ie)
      cnyi(5,iee)=cny(i5,ie)
      cnyi(6,iee)=cny(i6,ie)
      cnyi(7,iee)=flag(cny(i1,ie))
      cnyi(8,iee)=flag(cny(i2,ie))
      cnyi(9,iee)=flag(cny(i3,ie))
      cnyi(10,iee)=flag(cny(i4,ie))
      cnyi(11,iee)=flag(cny(i5,ie))
      cnyi(12,iee)=flag(cny(i6,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! coori(j,flag(cny(i4,ie)))=coor(j,cny(i4,ie))
      ! coori(j,flag(cny(i5,ie)))=coor(j,cny(i5,ie))
      ! coori(j,flag(cny(i6,ie)))=coor(j,cny(i6,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif
c
      enddo
c
      end
c_______________________________________________________________________
      subroutine tetIE4
     & (n,ne,cny,num,flag,id,cnyi,numi,nfla,ni,
     & nei,iecount,coor,coori)
      implicit none
      integer cny(4,ne),num(ne)
      integer flag(n),cnynl(4),id,iecount
      integer cnyi(6,nei),numi(nei),nfla(nei)
      integer n,ne,ni,nei
      integer iee,ie,ii,i1,i2,i3,j
      real*8 coor(3,n),coori(3,n+ni)
c
      do ie=1,ne
c
      cnynl=0
      do ii=1,4
      if(flag(cny(ii,ie)).gt.0)then
      cnynl(ii)=1
      endif
      enddo
c
      i1=1
      i2=2
      i3=4
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=flag(cny(i1,ie))
      cnyi(5,iee)=flag(cny(i2,ie))
      cnyi(6,iee)=flag(cny(i3,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif

      i1=2
      i2=3
      i3=4
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=flag(cny(i1,ie))
      cnyi(5,iee)=flag(cny(i2,ie))
      cnyi(6,iee)=flag(cny(i3,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif

      i1=3
      i2=1
      i3=4
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=flag(cny(i1,ie))
      cnyi(5,iee)=flag(cny(i2,ie))
      cnyi(6,iee)=flag(cny(i3,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif

      i1=2
      i2=1
      i3=3
      if(cnynl(i1)+ cnynl(i2)+ cnynl(i3).eq.3)then
      iecount=iecount+1
      iee=iecount
      cnyi(1,iee)=cny(i1,ie)
      cnyi(2,iee)=cny(i2,ie)
      cnyi(3,iee)=cny(i3,ie)
      cnyi(4,iee)=flag(cny(i1,ie))
      cnyi(5,iee)=flag(cny(i2,ie))
      cnyi(6,iee)=flag(cny(i3,ie))
      ! do j=1,3
      ! coori(j,flag(cny(i1,ie)))=coor(j,cny(i1,ie))
      ! coori(j,flag(cny(i2,ie)))=coor(j,cny(i2,ie))
      ! coori(j,flag(cny(i3,ie)))=coor(j,cny(i3,ie))
      ! enddo
      numi(iee)=num(ie)
      nfla(iee)=id
      endif
c
      enddo
      end
c_______________________________________________________________________
      subroutine setIEpentaCRS
     & (im,kd,n,ni,nei,cnysize,cnyi,numi,nfla,younglst,coori,dsi,
     & niw,pentawptr,pentawind,nir,pentarptr,pentarind,crsval,
     & crssize)
      implicit none
      integer im,kd,n,ni,nei,cnysize,cnyi(cnysize,nei),numi(nei)
      real*8 younglst(2,kd),dsi,coori(3,n+ni)
      integer nfla(nei)
      integer pentawptr(nsta+2),pentarptr(ni*4)
      integer crssize
      integer pentawind(ni*4),pentarind(ni*crssize)
      integer niw,nir
      real*8 crsval(9,ni*crssize)
c
      integer numitmp(nei)
      integer pentasize(n+ni)
      integer pentaind(crssize,n+ni)
      integer i,j,i1,j1,i2,j2,k,ie,ista,nmax,ii
      real*8 keit(3*cnysize,3*cnysize),young,rnyu,cori(3,3)

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
      do ie=1,nei
        if(numitmp(ie).eq.ista)then
          do j1=1,cnysize
            i1=cnyi(j1,ie)
            do j2=1,cnysize
              i2=cnyi(j2,ie)
              do k=1,pentasize(i1)
                if(pentaind(k,i1).eq.i2)then
                  goto 314
                endif
              enddo ! k
              pentasize(i1)=pentasize(i1)+1
              if(pentasize(i1).eq.crssize)then
                write(*,*) 'enlarge CRSSIZE',im,crssize,pentaind(:,i1)
                stop
              endif
              pentaind(pentasize(i1),i1)=i2
314           k=0
            enddo ! j2
          enddo ! j1
        endif
      enddo ! ie

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
      do ie=1,nei
        if(numitmp(ie).eq.ista)then

          young=younglst(1,numi(ie))
          rnyu=younglst(2,numi(ie))
          do i1=1,3
          do ii=1,3
          cori(i1,ii)=coori(ii,cnyi(i1,ie))
          enddo
          enddo
          if(cnysize.eq.12)then
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

          do j1=1,cnysize
            i1=cnyi(j1,ie)
            do i=pentawptr(ista)+1,pentawptr(ista+1)
              if(pentawind(i).eq.i1)then
                do j2=1,cnysize
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
                  write(*,*) 'pentaw not found',i1,ie
                  stop
142               k=0
                enddo ! j2
                goto 141
              endif
            enddo ! i
            write(*,*) 'pentaw not found',i1
            stop
141         i=0
          enddo ! j1
        endif
      enddo
      enddo
!      write(*,*) 'check crs table succeed'

      end
c
c_______________________________________________________________________
      subroutine compAmatd_iepenta
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,sig,p,q)
      implicit none
      integer im,n,ni
      integer nir,niw
      integer pentawptr(nsta+2),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)
      real*8 pentacrsval(9,nir)
      real*8 sig(nsta),p(3*npc*(n+ni)),q(3*npc*(n+ni))
      integer ista,i,i1,j,i2,ipc
      real*8 tmpq(npc),tmpqq(npc),tmpv(npc,3)

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
     & + pentacrsval(1,j) * p(3*npc*(i2-1)+npc*0+ipc)
     & + pentacrsval(2,j) * p(3*npc*(i2-1)+npc*1+ipc)
     & + pentacrsval(3,j) * p(3*npc*(i2-1)+npc*2+ipc)
      tmpv(ipc,2)=tmpv(ipc,2)
     & + pentacrsval(4,j) * p(3*npc*(i2-1)+npc*0+ipc)
     & + pentacrsval(5,j) * p(3*npc*(i2-1)+npc*1+ipc)
     & + pentacrsval(6,j) * p(3*npc*(i2-1)+npc*2+ipc)
      tmpv(ipc,3)=tmpv(ipc,3)
     & + pentacrsval(7,j) * p(3*npc*(i2-1)+npc*0+ipc)
     & + pentacrsval(8,j) * p(3*npc*(i2-1)+npc*1+ipc)
     & + pentacrsval(9,j) * p(3*npc*(i2-1)+npc*2+ipc)
      enddo ! ipc
      enddo ! j

      do j=1,3
      do ipc=1,npc
      tmpq(ipc)=tmpv(ipc,j)
      enddo
      call compsigd_npc21(ista,sig,tmpq,tmpqq)
      do ipc=1,npc
        q(3*npc*(i1-1)+npc*(j-1)+ipc)=
     &  q(3*npc*(i1-1)+npc*(j-1)+ipc)+tmpqq(ipc)
      enddo ! ipc
      enddo ! j
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ista

      end
c_______________________________________________________________________
      subroutine compAmats_iepenta
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,sig,p,q)
      implicit none
      integer im,n,ni
      integer nir,niw
      integer pentawptr(nsta+2),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)
      REAL_4 pentacrsval(9,nir)
      REAL_4 sig(nsta),p(3*npc*(n+ni)),q(3*npc*(n+ni))
      integer ista,i,i1,j,i2,ipc
      REAL_4 tmpq(npc),tmpqq(npc),tmpv(npc,3)

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
     & + pentacrsval(1,j) * p(3*npc*(i2-1)+npc*0+ipc)
     & + pentacrsval(2,j) * p(3*npc*(i2-1)+npc*1+ipc)
     & + pentacrsval(3,j) * p(3*npc*(i2-1)+npc*2+ipc)
      tmpv(ipc,2)=tmpv(ipc,2)
     & + pentacrsval(4,j) * p(3*npc*(i2-1)+npc*0+ipc)
     & + pentacrsval(5,j) * p(3*npc*(i2-1)+npc*1+ipc)
     & + pentacrsval(6,j) * p(3*npc*(i2-1)+npc*2+ipc)
      tmpv(ipc,3)=tmpv(ipc,3)
     & + pentacrsval(7,j) * p(3*npc*(i2-1)+npc*0+ipc)
     & + pentacrsval(8,j) * p(3*npc*(i2-1)+npc*1+ipc)
     & + pentacrsval(9,j) * p(3*npc*(i2-1)+npc*2+ipc)
      enddo ! ipc
      enddo ! j

      do j=1,3
      do ipc=1,npc
      tmpq(ipc)=tmpv(ipc,j)
      enddo
      call compsigs_npc21(ista,sig,tmpq,tmpqq)
      do ipc=1,npc
        q(3*npc*(i1-1)+npc*(j-1)+ipc)=
     &  q(3*npc*(i1-1)+npc*(j-1)+ipc)+tmpqq(ipc)
      enddo ! ipc
      enddo ! j
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ista

      end
c_______________________________________________________________________