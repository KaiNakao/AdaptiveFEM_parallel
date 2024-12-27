!#define PRINT_COLOR_LOG
#ifdef MAXCOLOR
      subroutine setcolor
     - (im,n_size,ne_size,elem_size,elem4_size,
     - n,ne,np,cny,num,
     -  ncolor,color_ind,elemmap,cny4,num4)
      implicit none
      ! elem_size: 10 for tet10, 4 for tet4, 8 for vox8, 27 for vox
      integer n_size,ne_size,n,ne,np,elem_size,elem4_size,im
      integer ncolor,color_ind(np*MAXCOLOR)
      integer color_balance(np*MAXCOLOR)
      integer cny(elem_size,ne_size),num(ne_size)
      integer cny4(elem4_size,ne_size),num4(ne_size)
      integer cny4tmp(elem4_size,ne_size)
      integer elemmap(ne_size)

      integer eptr(ne+1),eind(ne*elem_size)
      integer vwgt(ne),vsize(ne)
      integer ncommon,nparts
      real tpwgts(np)
      integer options(0:39),objval
      integer epart(ne),npart(n)
      integer epart2(ne),coloring(ne)
      integer*8 itmp

      integer numcommon(n),commonlist(np,n)
      real dummy
      integer ie,i,j,k,ii,ip,err
      integer ntmp,netmp,ntmp1,netmp1,loop
      integer nmapping(n),nemapping(ne),nemapping1(ne)

      integer METIS_PartMeshDual
      external METIS_PartMeshDual
!      write(*,*) 'coloring: implement elemmap later'

      coloring=0
      do ie=1,ne
      nemapping(ie)=ie
      nemapping1(ie)=ie
      enddo

      vwgt=1
      vsize=1
      ncommon=2
      nparts=np
      tpwgts=1.0/np
      call METIS_SetDefaultOptions(options)

      ntmp=n
      netmp=ne

      ! set mesh connectivity
      do ie=1,netmp
      eptr(ie)=elem_size*(ie-1)
      do j=1,elem_size
      eind(elem_size*(ie-1)+j)=cny(j,ie)-1 ! to c style
      enddo
      enddo
      eptr(netmp+1)=elem_size*netmp
#ifdef PRINT_COLOR_LOG
      write(*,*) 'start coloring',ne,n,'into',np,'parts'
#endif
      do loop=1,10

      if(nparts.gt.1)then
      ! partition mesh using METIS_PartMeshDual
      err = METIS_PartMeshDual
     - (netmp,ntmp,eptr,eind,vwgt,vsize,ncommon,nparts,
     - tpwgts,options,objval,epart,npart)
      if(err.ne.1)then
        write(*,*) 'error in METIS_PartMeshDual',err
        stop
      endif
      else
      do ie=1,netmp
      epart(ie)=0
      enddo
      endif
      do ie=1,netmp
      epart(ie)=epart(ie)+1 ! to fortran style
      enddo

      ! extract mesh that is nonoverlapping with other partitions
      numcommon=0
      do ie=1,netmp
      ip=epart(ie)
      do j=1,elem_size
        i=eind(elem_size*(ie-1)+j)+1 ! to fortran style
        do k=1,numcommon(i)
          if(commonlist(k,i).eq.ip)then
            goto 294
          endif
        enddo
        numcommon(i)=numcommon(i)+1
        commonlist(numcommon(i),i)=ip
294     dummy=0.0
      enddo
      enddo

      epart2=0
      do ie=1,netmp
      do j=1,elem_size
        i=eind(elem_size*(ie-1)+j)+1 ! to fortran style
        if(loop.eq.1)then
          if(numcommon(i).ne.1)then
          goto 295
          endif
        else
          ip=epart(ie)
          do ii=1,numcommon(i)
            if(commonlist(ii,i).gt.ip)then
            goto 295
            endif
          enddo
        endif
      enddo
      epart2(ie)=epart(ie)
      coloring(nemapping(ie))=epart(ie)+(loop-1)*np
295   dummy=0.0
      enddo

      ! extract residual mesh
      netmp1=0
      ntmp1=0
      nmapping=0
      do ie=1,netmp
      if(epart2(ie).eq.0)then
        netmp1=netmp1+1
        nemapping1(netmp1)=nemapping(ie)
        do j=1,elem_size
          i=eind(elem_size*(ie-1)+j)+1 ! to fortran style
          if(nmapping(i).eq.0)then
            ntmp1=ntmp1+1
            nmapping(i)=ntmp1
          endif
          eind(elem_size*(netmp1-1)+j)=nmapping(i)-1 ! to c style
        enddo
      endif
      enddo
#ifdef PRINT_COLOR_LOG
      write(*,*) 'loop,ne_left,n_left',loop,netmp1,ntmp1
#endif

#ifdef DEBUG
      netmp1=0
      do ie=1,ne
      if(coloring(ie).eq.0)then
        netmp1=netmp1+1
      endif
      enddo
      write(*,*) 'loop,ne_left,n_left',loop,netmp1,ntmp1
#endif

      netmp=netmp1
      ntmp=ntmp1
      ncolor=loop
      nemapping=nemapping1

      if(netmp.lt.100)then
        goto 255
      endif

      enddo ! loop

255   dummy=0.0

      netmp=0
      do ie=1,ne
        if(coloring(ie).eq.0)then
          netmp=netmp+1
          coloring(ie)=ncolor*np+1
        endif
      enddo
      if(netmp.ne.0)then
        ncolor=ncolor+1
!        write(*,*) 'assign residual elements to proc 0',netmp
      endif

      color_balance=0
      do ie=1,ne
      ip=coloring(ie)
      color_balance(ip)=color_balance(ip)+1
      enddo

      color_ind=0
      do i=1,ncolor*np
      color_ind(i+1)=color_ind(i)+color_balance(i)
      enddo

#ifdef PRINT_COLOR_LOG
      do i=1,ncolor
      do j=1,np
      write(*,*) 'color balance',i,j,color_balance(np*(i-1)+j)
      enddo
      enddo
#endif
      color_balance=0
      itmp=0
      do ie=1,ne
      ip=coloring(ie)
      color_balance(ip)=color_balance(ip)+1
      ii=color_ind(ip)+color_balance(ip)
      epart(ii)=num(ie)

      elemmap(ii)=ie

      do j=1,elem_size
      eind(elem_size*(ii-1)+j)=cny(j,ie)
      itmp=itmp+cny(j,ie)
      enddo
      do j=1,elem4_size
      cny4tmp(j,ii)=cny4(j,ie)
      enddo
      enddo

!      elemmap=nemapping
!      do ie=1,nab
!      abcelem(ie)=elemmap(abcelem(ie))
!      enddo

      nemapping=0
      do ie=1,ne
      nemapping(elemmap(ie))=nemapping(elemmap(ie))+1
      enddo
      do ie=1,ne
      if(nemapping(ie).ne.1)then
      write(*,*) 'error in nemapping',im,ie,nemapping(ie)
      endif
      enddo


#ifdef PRINT_COLOR_LOG
      do i=1,ncolor
      do j=1,np
      write(*,*) 'color balance',i,j,
     - color_balance(np*(i-1)+j),
     - color_ind(np*(i-1)+j)
      enddo
      enddo
      i=ncolor+1
      j=1
      write(*,*) 'color balance',i,j,
     - color_balance(np*(i-1)+j),
     - color_ind(np*(i-1)+j)
#endif

      do ie=1,ne
      num(ie)=epart(ie)
      num4(ie)=epart(ie)
      do j=1,elem_size
      cny(j,ie)=eind(elem_size*(ie-1)+j)
      enddo
      do j=1,elem4_size
      cny4(j,ie)=cny4tmp(j,ie)
      enddo
      enddo

      ! check coloring
      do ie=1,ne
      do j=1,elem_size
      itmp=itmp-cny(j,ie)
      enddo
      enddo
      if(itmp.ne.0)then
        write(*,*) 'ERROR in coloring',itmp
        stop
      endif

      ! check coloring
      do i=1,ncolor
      npart=0
      do j=1,np
      do ie=color_ind(np*(i-1)+j)+1,color_ind(np*(i-1)+j+1)
      do k=1,elem_size
        ii=cny(k,ie)
        if(npart(ii).eq.0)then
          npart(ii)=j
        endif
        if(npart(ii).ne.j)then
          write(*,*) 'ERROR in coloring',ii,j,npart(ii),k,ie
          stop
        endif
      enddo
      enddo
      enddo
      enddo

#ifdef PRINT_COLOR_LOG
      write(*,*) 'coloring finished'
#endif
      end
#endif
      subroutine dummy_func(i,j)
      implicit none
      integer i,j
      i=j
      end
