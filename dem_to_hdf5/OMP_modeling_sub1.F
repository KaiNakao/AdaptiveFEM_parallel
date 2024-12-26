c make subroutine nfem2        -> need to use
c make subroutine voxelnode2   -> no need to use
c make subroutine pyramid2     -> no need to use
c modify subroutine valest
c modify subroutine deldouble
c modify subroutine nfem
c_______________________________________________________________________
      subroutine 
     -modeling(nx,ny,nz,ns,sur,oris,ds,freq,now,rmat,ix,iy,im,jm,nk
     - ,coeout,nout,noutdep)
!$      use omp_lib
c nt ��4�̔{��
      parameter(nt=400,net=400) ! original
!      parameter(nt=50,net=400) ! made nt smaller for reducing memory usage
      integer*2 numori(nx,ny,nz)
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
      integer nf(nx+1,ny+1,nz+1)
      real*8 ds,c2(ns),freq,now,coeout
      real*8 sur(nx+1,ny+1,ns),ori(3),rmat(ns,10),oris(3)  
      integer nsur(nx+1,ny+1,ns)
      integer cnyn(nx*ny*nt,5),numsur(nx+1,ny+1)
      character filename*50
      integer nfemelem
c      INTEGER FUNCTION OMP_GET_THREAD_NUM()
c
      write(*,*)'using version 20150307'
      write(*,*)'nt,net',nt,net
      write(*,*)'nx,ny,nz',nx,ny,nz
      write(*,*)'ns,ds',ns,ds
!$omp parallel
      np=omp_get_num_threads()
!$omp end parallel

      write(6,*) 'num of thread =', np

      do is=1,ns
      c2(is)=rmat(is,2)*rmat(is,8)
      enddo
cc 
c surface input
      do is=1,ns
      do j=1,ny+1
      do i=1,nx+1
      sur(i,j,is)=sur(i,j,is)-oris(3)
      enddo
      enddo
      enddo
c
       id=1
c
c      write(6,*) id,'-thread begin for',im,jm,' region'
c
c      open(31,file='./data/model.flavia.msh',status='unknown')
      itmp=ix*100000000+iy*10000+im*100+jm
!      write(filename,'(a29)')
!     - './data/tet4model3d.flavia.msh'
!      open(31+3*id,file=filename,status='unknown')
!      write(filename,'(a29)')
!     - './data/tet4model2d.flavia.msh'
!      open(33+3*id,file=filename,status='unknown')
c      open(42,file='./data/coor.dat',status='unknown')
c      open(43,file='./data/cnynfem.dat',status='unknown')
c      open(44,file='./data/cnyvfem.dat',status='unknown')
c      open(51,file='./data/surface.dat',status='unknown')
      write(filename,'(a26)')
     - './data/tet4vox8setting.dat'
!     - './data/tet4setting.dat'
      open(32+3*id,file=filename,status='unknown')
c
!      open(65,file='./data/tet4setting2d.dat',status='unknown')
      open(65,file='./data/tet4vox8setting2d.dat',status='unknown')
      write(65,*) 
      write(65,*) 'num of node'
      write(65,*) (nx+1)*(ny+1)
      write(65,*) 'num of elem'
      write(65,*) nx*ny	
      close(65)
c
c surface input
      ori(1)=0.
      ori(2)=0.
      ori(3)=0.
cc
c surface modification
c      call surfacemod(nx,ny,nz,ns,ori,ds,sur)
      call surfacemod80(nx,ny,nz,ns,ori,ds,sur)
c surface modification

ccc______________________________________________________
c num data set. can be compacted! 
c num data 
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
      numori(i,j,k)=0
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
      do is=1,ns
      nzt=int((sur(i,j,is)-ori(3))/ds) 
      if(abs(sur(i,j,is)-nzt*ds+ori(3)).gt.0.99*ds)then
      nzt=nzt+1
      endif
      do k=1,nzt
      numori(i,j,k)=is
      enddo
      enddo
      enddo
      enddo
c num data
cc


c num data modification 
      do is=1,ns
      do j=1,ny
      do i=1,nx
       if(
     -((abs(sur(i,j,is)-sur(i+1,j,is))+
     -  abs(sur(i,j,is)-sur(i+1,j+1,is))+
     -  abs(sur(i,j,is)-sur(i,j+1,is))).lt.0.01*ds).and.
     -  (mod(sur(i,j,is),ds).lt.0.01*ds))then
       nzt=int((sur(i,j,is)-ori(3))/ds)
       numori(i,j,nzt)=101
       numori(i,j,nzt+1)=101
      else
       do j1=0,1
       do i1=0,1
       if(abs(mod(sur(i+i1,j+j1,is)-ori(3),ds)).ge.0.01*ds)then
       nzt=int((sur(i+i1,j+j1,is)-ori(3))/ds)
       if(sur(i+i1,j+j1,is)-ori(3)-nzt*ds.gt.0.01*ds)then
       nzt=nzt+1 
       endif
       numori(i,j,nzt)=101
       endif
       enddo
       enddo
       endif
      enddo
      enddo
      enddo
c num data modification 
c
      do k=1,2*nz-1
      do j=1,2*ny-1
      do i=1,2*nx-1
      num(i,j,k)=0
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny
      do i=1,nx
      num(2*i-1,2*j-1,2*k-1)=numori(i,j,k)
      enddo
      enddo
      enddo
ccc______________________________________________________
c nk
c      call calcnk(nx,ny,nz,nk)        
      call octree(nx,ny,nz,nk,ns,num,ds,freq,c2,now,coeout,nout,noutdep)

      write(6,*) 'octree'
c
c node estimation in VFEM
      call voxelnode(nx,ny,nz,num,nk,nf,n) 
c      call voxelnode2(nx,ny,nz,num,nk,nf,n) 

      write(6,*) 'voxel'

c node estimation in VFEM
cc
      call pyramid(nx,ny,nz,nk,nt,net,nf,num,cnyn,nfemelem,np)
c      call pyramid2(nx,ny,nz,nk,nt,net,nf,num,cnyn,nfemelem)
c
      write(6,*) 'pyramid'

c      call nfem(nx,ny,nz,ns,nt,net,
c     - sur,nf,numori,ori,ds,n,nsur,nfemelem,cnyn,numsur)
      call nfem2(nx,ny,nz,ns,nt,net,
     - sur,nf,numori,ori,ds,n,nsur,nfemelem,cnyn,numsur,np)

      write(6,*) 'nfem'

      call output
     -(nx,ny,nz,ns,nk,nt,ds,ori,nf,nsur,n,sur,nfemelem,cnyn,num,oris,id)

!      write(33+3*id,300) '"bar"', 3, ' Tetrahedra', 4
!  300 format('MESH',a10,'  dimension',i5,'  ElemType',a11,'  Nnode',i5)
c Node
!      write(33+3*id,'(a)') 'Coordinates'
!      write(33+3*id,'(a)')'# node_number     coor_x    coor_y    coor_z'
      itmp=0
      do j=1,ny+1
      do i=1,nx+1
      itmp=itmp+1
!      write(33+3*id,'(i10,3f16.6)') numsur(i,j),
!     - ds*(i-1)+ori(1)+oris(1),
!     - ds*(j-1)+ori(2)+oris(2),sur(i,j,1)+oris(3)
      enddo
      enddo
!      write(33+3*id,'(a)') 'end coordinates'
!      write(33+3*id,*)
c
!      write(33+3*id,*)''
!      write(33+3*id,*) 'MESH     "bar"  dimension    3  
!     - ElemType   Quadrilateral   Nnode    4'
!      write(33+3*id,*)'Elements'
!      write(33+3*id,*)'# element node_1 - node_4 material_number'
       inum2=0
       do j=1,ny     
       do i=1,nx
       inum2=inum2+1
!    	 write(33+3*id,200) inum2,numsur(i,j),numsur(i+1,j),
!     - numsur(i+1,j+1),numsur(i,j+1),1
      enddo
      enddo
!  200 format(11i10)
!      write(33+3*id,'(a)') 'end elements'
!      close(33+3*id)

      end
c_______________________________________________________________________
      subroutine pyramid(nx,ny,nz,nk,nt,net,nf,num,cnyn,nfemelem,np)
!$      use omp_lib
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
	integer nf(nx+1,ny+1,nz+1)
	integer cnyn(nx*ny*nt,5)
	real*8 coortmp(nt,3)
	integer cnytmp(net,4),numtmp(nt)
!	integer cnyntmp(np,nx*ny*nt/4,5),nfemtmp(np),in,ine
	integer cnyntmp(np,nx*ny*nt*4/np,5),nfemtmp(np),in,ine
      integer nfemelem,cnyntmpsize
      cnyntmpsize=nx*ny*nt*4/np

      nfemtmp=0
c
      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)

!$OMP parallel private(nex,ney,nez,itet,ii,jj,kk,
!$OMP& coortmp,numtmp,cnytmp,in,ine,i,j,k,ip,i3,i2)
!$OMP do   
        do j=1,ny/nd
c        
        ip=omp_get_thread_num()+1
c        write(6,*) ip,np       

c        
        do k=1,nz/nd
        do i=1,nx/nd
         nex=nd+2*nd*(i-1)
         ney=nd+2*nd*(j-1)
         nez=nd+2*nd*(k-1)
          if((num(nex,ney,nez).gt.0)
     -       .and.(num(nex,ney,nez).lt.100))then
            itet=0
             do kk=nez/2-nd2+1,nez/2+nd2+1
             do jj=ney/2-nd2+1,ney/2+nd2+1
             do ii=nex/2-nd2+1,nex/2+nd2+1
c            do kk=nez/2-nd2+1,nez/2+nd2+1,nd2
c            do jj=ney/2-nd2+1,ney/2+nd2+1,nd2
c            do ii=nex/2-nd2+1,nex/2+nd2+1,nd2
             if(nf(ii,jj,kk).gt.0)then
              itet=itet+1
              if(itet.eq.(nt-1))then
              write(*,*) 'enlarge nt'
              stop
              endif
	      coortmp(itet,1)=ii
	      coortmp(itet,2)=jj
	      coortmp(itet,3)=kk
	      numtmp(itet)=nf(ii,jj,kk)
             endif
            enddo
            enddo
            enddo

           if(itet.gt.8)then
            call deldouble(nt,itet,coortmp,numtmp)
c        write(6,*) ip,np,itet       
            call Delaunay(nt,net,itet,coortmp,cnytmp,in,ine)
c        write(6,*) ip,np       
             if(ine.gt.0)then
	      do i3=1,ine
            nfemtmp(ip)=nfemtmp(ip)+1
!            if(nfemtmp(ip).eq.(nx*ny*nt/4-1))then
            if(nfemtmp(ip).eq.(cnyntmpsize-1))then
            write(*,*) 'enlarge cnytmpsize nt'
            stop
            endif
              do i2=1,4
            cnyntmp(ip,nfemtmp(ip),i2)=numtmp(cnytmp(i3,i2))
              enddo
            cnyntmp(ip,nfemtmp(ip),5)=num(nex,ney,nez)
              enddo
             endif
            num(nex,ney,nez)=0
            endif !itet.eq.8
           endif !
       enddo !i
       enddo !j
       enddo !k
!$OMP end parallel
      enddo !ik

c      write(6,*) 123

      nfemelem=0
      do ip=1,np
      if(nfemtmp(ip).gt.0)then
      do i=1,nfemtmp(ip)
      nfemelem=nfemelem+1
      if(nfemelem.eq.(nx*ny*nt-1))then
      write(*,*) 'enlarge nt'
      stop
      endif
      do i2=1,5
      cnyn(nfemelem,i2)=cnyntmp(ip,i,i2)
      enddo 
      enddo
      endif
      enddo

c      write(6,*) 124,nfemelem


      end
c_______________________________________________________________________
!      subroutine pyramid2(nx,ny,nz,nk,nt,net,nf,num,cnyn,nfemelem)
!      integer num(2*nx-1,2*ny-1,2*nz-1)
!	integer nf(nx+1,ny+1,nz+1)
!	integer cnyn(nx*ny*nt,5)
!	real*8 coortmp(nt,3)
!	integer cnytmp(net,4),numtmp(nt)
!c
!      nfemelem=0
!c
!      do ik=2,nk
!       nd=2**(ik-1)
!       nd2=2**(ik-2)
!        do k=1,nz/nd
!        do j=1,ny/nd
!        do i=1,nx/nd
!         nex=nd+2*nd*(i-1)
!         ney=nd+2*nd*(j-1)
!         nez=nd+2*nd*(k-1)
!          if((num(nex,ney,nez).gt.0)
!     -       .and.(num(nex,ney,nez).lt.100))then
!           itet=0
!            do kk=nez/2-nd2+1,nez/2+nd2+1,nd2
!            do jj=ney/2-nd2+1,ney/2+nd2+1,nd2
!            do ii=nex/2-nd2+1,nex/2+nd2+1,nd2
!             if(nf(ii,jj,kk).gt.0)then
!              itet=itet+1  
!      if(itet.eq.(nt-1))then
!      write(*,*) 'enlarge nt'
!      stop
!      endif
!	      coortmp(itet,1)=ii
!	      coortmp(itet,2)=jj
!	      coortmp(itet,3)=kk
!	      numtmp(itet)=nf(ii,jj,kk)
!             endif
!            enddo
!            enddo
!            enddo
!
!           if(itet.gt.8)then
!c_____add_____
!	     itet=itet+1
!	      ii=nex/2-nd2+1 + nd2
!	      jj=ney/2-nd2+1 + nd2
!	      kk=nez/2-nd2+1 + nd2
!	     coortmp(itet,1)=ii
!	     coortmp(itet,2)=jj
!	     coortmp(itet,3)=kk
!	      numtmp(itet)=nf(ii,jj,kk)
!c_____add_____
!            call deldouble(nt,itet,coortmp,numtmp)
!            call Delaunay(nt,net,itet,coortmp,cnytmp,in,ine)
!             if(ine.gt.0)then
!	      do i3=1,ine
!              nfemelem=nfemelem+1
!              if(nfemelem.eq.(nx*ny*nt-1))then
!              write(*,*) 'enlarge nt'
!              stop
!              endif
!              do i2=1,4
!	      cnyn(nfemelem,i2)=numtmp(cnytmp(i3,i2))
!              enddo
!	      cnyn(nfemelem,5)=num(nex,ney,nez)
!              enddo
!             endif
!            num(nex,ney,nez)=0
!            endif !itet.eq.8
!           endif !
!       enddo !i
!       enddo !j
!       enddo !k
!      enddo !ik
!
!      end
c_______________________________________________________________________
!      subroutine nfem(nx,ny,nz,ns,nt,
!     -  net,sur,nf,num,ori,ds,n,nsur,nfemelem,cnyn,numsur)
!      integer num(nx,ny,nz)
!	integer nf(nx+1,ny+1,nz+1)
!      real*8 ds,coortmp(nt,3)
!	real*8 sur(nx+1,ny+1,ns),ori(3),ep,df
!      integer nsur(nx+1,ny+1,ns)
!	integer numtmp(nt)
!	integer cnyn(nx*ny*nt,5),numsur(nx+1,ny+1)
!	real*8 cortmp(nt,3),cl(3),dfx,dfy,dfz
!	integer cnytmp(net,4),nutmp(nt),lu(0:101,2),rlist(4),ndest(2,2)
!	real*8 eptmp
!c
!c
!      ep=ds*0.01
!      eptmp=ds*0.1
!c
!c--- setting of nsur(i,j,is)
!      inum=n
!	do is=1,ns
!	do j=1,ny+1
!	do i=1,nx+1
!c
!      k=int((sur(i,j,is)-ori(3))/ds)
!      if(sur(i,j,is)-ori(3)-k*ds.gt.0.5*ds)then
!      k=k+1
!      endif
!c
!	df=mod(sur(i,j,is)-ori(3),ds)
!c
!	if(is.ge.2)then
!	if(abs(sur(i,j,is-1)-sur(i,j,is)).lt.ep)then
!      nsur(i,j,is)=nsur(i,j,is-1)      
!	else
!c____mod trial__________
!c_____org
!c      if((abs(df).ge.ep).or.(nf(i,j,k+1).eq.0))then
!c_____mod
!      if(
!     &  ((abs(df).ge.ep).and.(abs(df).le.ds-ep))
!     &  .or.
!     &  (nf(i,j,k+1).eq.0)
!     &  )then
!c____mod trial__________
!c
!      inum=inum+1
!      nsur(i,j,is)=inum
!	else
!      nsur(i,j,is)=nf(i,j,k+1)
!c
!      endif
!	endif
!      else
!c____mod trial__________
!c_____org
!      if((abs(df).ge.ep).or.(nf(i,j,k+1).eq.0))then
!c_____mod
!c      if(
!c     &  ((abs(df).ge.ep).and.(abs(df).le.ds-ep))
!c     &  .or.
!c     &  (nf(i,j,k+1).eq.0)
!c     &  )then
!c____mod trial__________
!      inum=inum+1
!      nsur(i,j,is)=inum
!	else
!      nsur(i,j,is)=nf(i,j,k+1)
!      endif
!	endif
!      enddo
!      enddo
!      enddo
!c
!      do j=1,ny+1
!      do i=1,ny+1
!      numsur(i,j)=nsur(i,j,1)
!      enddo
!      enddo
!
!c
!      do j=1,ny
!      do i=1,nx
!      inum=0
!c     node definision (i,j is fixed)
!      do k=1,nz
!       if(num(i,j,k).gt.100)then
!        if(inum.eq.0)then
!         do k1=0,1
!         do j1=0,1
!         do i1=0,1
!          inum=inum+1
!          numtmp(inum)=nf(i+i1,j+j1,k+k1)
!          coortmp(inum,1)=(i+i1-1)*ds
!          coortmp(inum,2)=(j+j1-1)*ds
!          coortmp(inum,3)=(k+k1-1)*ds
!         enddo
!         enddo
!         enddo
!        else
!         do k1=0,1
!         do j1=0,1
!         do i1=0,1
!          do it=1,inum
!           if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
!	    goto 101
!	   endif
!          enddo
!          inum=inum+1
!          numtmp(inum)=nf(i+i1,j+j1,k+k1)
!          coortmp(inum,1)=(i+i1-1)*ds
!          coortmp(inum,2)=(j+j1-1)*ds
!          coortmp(inum,3)=(k+k1-1)*ds
!  101 damy=0 
!         enddo
!         enddo
!         enddo
!        endif
!       endif
!      enddo !k
!c
!      do is=1,ns
!      do j1=0,1
!      do i1=0,1
!      if(nsur(i+i1,j+j1,is).gt.n)then
!      inum=inum+1
!      numtmp(inum)=nsur(i+i1,j+j1,is)
!      coortmp(inum,1)=(i+i1-1)*ds
!      coortmp(inum,2)=(j+j1-1)*ds
!      coortmp(inum,3)=sur(i+i1,j+j1,is)
!      endif
!      enddo
!      enddo
!      enddo
!c
!cc
!c
!      do l2=1,2
!      do l1=1,101
!	lu(l1,l2)=-1
!      enddo
!      enddo
!      do k=1,nz
!	lu(num(i,j,k),2)=k
!      enddo
!      do k=nz,1,-1
!	lu(num(i,j,k),1)=k-1
!      enddo
!c
!      if(ns.ge.2)then
!	 do is=1,ns-1
!c_________________________________
!c there is VOXEL elem in taget-layer
!       if(lu(is,1).gt.0)then
!c
!       itet=0
!	 do i1=1,inum
!       ix=int((coortmp(i1,1)-ori(1))/ds)
!       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
!	 ix=ix+1
!       endif
!       iy=int((coortmp(i1,2)-ori(2))/ds)
!       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
!	 iy=iy+1
!       endif
!c    upper boundary of target-layer
!       if((coortmp(i1,3).le.sur(ix+1,iy+1,is)+ep).and.
!     -    (numtmp(i1).ne.0)
!     - .and.(coortmp(i1,3).ge.lu(is,2)*ds-ep+ori(3)))then
!	 itet=itet+1
!       do ii=1,3
!	 cortmp(itet,ii)=coortmp(i1,ii)
!	 enddo
!       nutmp(itet)=numtmp(i1)
!       endif
!       enddo
!       call deldouble(nt,itet,cortmp,nutmp)
!       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)    
!
!      if(ine.gt.0)then
!       do i3=1,ine
!cc
!       call valest(nt,net,cortmp,cnytmp,i3,ds,key)
!        if(key.gt.0)then
!	 nfemelem=nfemelem+1
!      if(nfemelem.eq.(nx*ny*nt-1))then
!      write(*,*) 'enlarge nt'
!      stop
!      endif
!          do i2=1,4
!	   cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
!          enddo
!	 cnyn(nfemelem,5)=is
!        endif
!       enddo
!      endif !ine.gt.0
!c
!       itet=0
!	 do i1=1,inum
!       ix=int((coortmp(i1,1)-ori(1))/ds)
!       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
!	 ix=ix+1
!       endif
!       iy=int((coortmp(i1,2)-ori(2))/ds)
!       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
!	 iy=iy+1
!       endif
!c    lower boundary of target-layer
!       if((coortmp(i1,3).ge.sur(ix+1,iy+1,is+1)-ep).and.
!     -    (numtmp(i1).ne.0)
!     - .and.(coortmp(i1,3).le.lu(is,1)*ds+ep+ori(3)))then
!	 itet=itet+1
!       do ii=1,3
!	 cortmp(itet,ii)=coortmp(i1,ii)
!	 enddo
!       nutmp(itet)=numtmp(i1)
!       endif
!       enddo
!
!       call deldouble(nt,itet,cortmp,nutmp)
!       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)    
!      if(ine.gt.0)then
!	 do i3=1,ine
!       call valest(nt,net,cortmp,cnytmp,i3,ds,key)
!      if(key.gt.0)then
!	nfemelem=nfemelem+1
!      if(nfemelem.eq.(nx*ny*nt-1))then
!      write(*,*) 'enlarge nt'
!      stop
!      endif
!       do i2=1,4
!	 cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
!       enddo
!	 cnyn(nfemelem,5)=is
!       endif
!       enddo
!      endif
!c
!c_________________________________
!c there is NO VOXEL elem in taget-layer
!	 else
!c
!       itet=0
!	 do i1=1,inum
!c
!       ix=int((coortmp(i1,1)-ori(1))/ds)
!       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
!	 ix=ix+1
!       endif
!       iy=int((coortmp(i1,2)-ori(2))/ds)
!       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
!	 iy=iy+1
!       endif
!c
!       if((coortmp(i1,3).le.sur(ix+1,iy+1,is)+ep).and.
!     -    (numtmp(i1).ne.0)
!     - .and.(coortmp(i1,3).ge.sur(ix+1,iy+1,is+1)-ep))then
!	 itet=itet+1
!       do ii=1,3
!	 cortmp(itet,ii)=coortmp(i1,ii)
!	 enddo
!       nutmp(itet)=numtmp(i1)
!       endif
!       enddo
!c    
!      ndest(1,1)=
!     -abs(nsur(i  ,j  ,is)-nsur(i  ,j  ,is+1))+
!     -abs(nsur(i+1,j  ,is)-nsur(i+1,j  ,is+1))+
!     -abs(nsur(i  ,j+1,is)-nsur(i  ,j+1,is+1))
!      ndest(2,1)=
!     -abs(nsur(i  ,j  ,is)-nsur(i  ,j  ,is+1))+
!     -abs(nsur(i+1,j  ,is)-nsur(i+1,j  ,is+1))+
!     -abs(nsur(i+1,j+1,is)-nsur(i+1,j+1,is+1))
!      ndest(2,2)=
!     -abs(nsur(i+1,j+1,is)-nsur(i+1,j+1,is+1))+
!     -abs(nsur(i+1,j  ,is)-nsur(i+1,j  ,is+1))+
!     -abs(nsur(i  ,j+1,is)-nsur(i  ,j+1,is+1))
!      ndest(1,2)=
!     -abs(nsur(i  ,j  ,is)-nsur(i  ,j  ,is+1))+
!     -abs(nsur(i+1,j+1,is)-nsur(i+1,j+1,is+1))+
!     -abs(nsur(i  ,j+1,is)-nsur(i  ,j+1,is+1))
!      il=0
!	do j2=1,2
!	do i2=1,2
!	if(ndest(i2,j2).eq.0)then
!      il=il+1
!	rlist(il)=nsur(i-1+i2,j-1+j2,is)
!	endif
!      enddo
!      enddo
!      if(il.gt.0)then
!      do i2=1,il
!      do i1=1,itet 
!      if(rlist(i2).eq.nutmp(i1))then
!	nutmp(i1)=0
!	endif
!      enddo
!      enddo
!	endif
!       call deldouble(nt,itet,cortmp,nutmp)
!       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)    
!       if(ine.gt.0)then
!        do i3=1,ine
!         call valest(nt,net,cortmp,cnytmp,i3,ds,key)
!         if(key.eq.1)then
!	  nfemelem=nfemelem+1
!      if(nfemelem.eq.(nx*ny*nt-1))then
!      write(*,*) 'enlarge nt'
!      stop
!      endif
!           do i2=1,4
!	    cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
!           enddo
!	  cnyn(nfemelem,5)=is
!         endif
!        enddo
!       endif
!
!      endif
!       enddo
!      endif
!c_________________________________
!c
!cc
!c lowest layer
!       itet=0
!	 do i1=1,inum
!       ix=int((coortmp(i1,1)-ori(1))/ds)
!       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
!	 ix=ix+1
!       endif
!       iy=int((coortmp(i1,2)-ori(2))/ds)
!       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
!	 iy=iy+1
!       endif
!       if((coortmp(i1,3).le.sur(ix+1,iy+1,ns)+ep).and.
!     -    (numtmp(i1).ne.0))then
!	 itet=itet+1
!       do ii=1,3
!	 cortmp(itet,ii)=coortmp(i1,ii)
!	 enddo
!       nutmp(itet)=numtmp(i1)
!       endif
!       enddo
!
!       call deldouble(nt,itet,cortmp,nutmp)
!       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)     
!      if(ine.gt.0)then
!	 do i3=1,ine
!       call valest(nt,net,cortmp,cnytmp,i3,ds,key)
!       if(key.gt.0)then
!       nfemelem=nfemelem+1
!      if(nfemelem.eq.(nx*ny*nt-1))then
!      write(*,*) 'enlarge nt'
!      stop
!      endif
!       do i2=1,4
!	 cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
!       enddo
!	 cnyn(nfemelem,5)=ns
!       endif
!       enddo
!       endif
!
!      enddo
!      enddo
!c
!      end
c_______________________________________________________________________
      subroutine nfem2(nx,ny,nz,ns,nt,
     -  net,sur,nf,num,ori,ds,n,nsur,nfemelem,cnyn,numsur,np)
!$      use omp_lib
      integer*2 num(nx,ny,nz)
	integer nf(nx+1,ny+1,nz+1)
      real*8 ds,coortmp(nt,3)
	real*8 sur(nx+1,ny+1,ns),ori(3),ep,df
      integer nsur(nx+1,ny+1,ns)
	integer numtmp(nt)
	integer cnyn(nx*ny*nt,5),numsur(nx+1,ny+1)
	real*8 cortmp(nt,3),cl(3),dfx,dfy,dfz
	integer cnytmp(net,4),nutmp(nt),lu(0:101,2),rlist(4),ndest(2,2)
	real*8 eptmp
      integer nfemelem
c
!	integer cnyntmp(np,nx*ny*nt/4,5),nfemtmp(np),in,ine
	integer cnyntmp(np,nx*ny*nt*4/np,5),nfemtmp(np),in,ine
c
      integer cnyntmpsize
      cnyntmpsize=nx*ny*nt*4/np

      ep=ds*0.01
c
c*** setting of nsur(i,j,is)
c
      inum=n
      do is=1,ns
      do j=1,ny+1
      do i=1,nx+1
       k=int((sur(i,j,is)-ori(3))/ds)
        if(sur(i,j,is)-ori(3)-k*ds.gt.0.5*ds)then
         k=k+1
        endif
       df=mod(sur(i,j,is)-ori(3),ds)
c
	if(is.ge.2)then
	if(abs(sur(i,j,is-1)-sur(i,j,is)).lt.ep)then
      nsur(i,j,is)=nsur(i,j,is-1)      
	else
c____mod trial__________
c_____org
c      if((abs(df).ge.ep).or.(nf(i,j,k+1).eq.0))then
c_____mod
      if(
     &  ((abs(df).ge.ep).and.(abs(df).le.ds-ep))
     &  .or.
     &  (nf(i,j,k+1).eq.0)
     &  )then
c____mod trial__________
c
      inum=inum+1
      nsur(i,j,is)=inum
	else
      nsur(i,j,is)=nf(i,j,k+1)
c
      endif
	endif
      else
c____mod trial__________
c_____org
      if((abs(df).ge.ep).or.(nf(i,j,k+1).eq.0))then
c_____mod
c      if(
c     &  ((abs(df).ge.ep).and.(abs(df).le.ds-ep))
c     &  .or.
c     &  (nf(i,j,k+1).eq.0)
c     &  )then
c____mod trial__________
      inum=inum+1
      nsur(i,j,is)=inum
	else
      nsur(i,j,is)=nf(i,j,k+1)
      endif
	endif
      enddo
      enddo
      enddo
c
      do j=1,ny+1
      do i=1,nx+1
      numsur(i,j)=nsur(i,j,1)
      enddo
      enddo    
     
      nfemtmp=0
c
c*** i j LOOP start
c
!$OMP parallel private(ip,i,j,k,itri,inum,k1,i1,j1,numtmp,coortmp,
!$OMP& it,is,l1,l2,lu,ix,iy,itet,cortmp,nutmp,cnytmp,
!$OMP& in,ine,i3,key,ndest,il,i2,j2,rlist)
!$OMP do   
      do j=1,ny

        ip=omp_get_thread_num()+1
        
      do i=1,nx
      do itri=0,1
c     node definision (i,j is fixed)
      inum=0
      do k=1,nz
       if(num(i,j,k).gt.100)then
        if(inum.eq.0)then
         do k1=0,1
          if(itri.eq.0)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            i1=0
            j1=0
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            i1=1
            j1=0
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            i1=0
            j1=1
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
          else
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            i1=1
            j1=0
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            i1=0
            j1=1
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            i1=1
            j1=1
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
          endif
         enddo
        else
         do k1=0,1
          if(itri.eq.0)then
            i1=0
            j1=0
               do it=1,inum
               if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
	       goto 201
	       endif
               enddo
            inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
  201 damy=0 
            i1=1
            j1=0
               do it=1,inum
               if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
	       goto 202
	       endif
               enddo
            inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
  202 damy=0 
            i1=0
            j1=1
               do it=1,inum
               if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
	       goto 203
	       endif
               enddo
            inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
  203 damy=0 
          else
            i1=1
            j1=0
               do it=1,inum
               if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
	       goto 204
	       endif
               enddo
            inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
  204 damy=0 
            i1=0
            j1=1
               do it=1,inum
               if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
	       goto 205
	       endif
               enddo
            inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
  205 damy=0 
            i1=1
            j1=1
               do it=1,inum
               if(numtmp(it).eq.nf(i+i1,j+j1,k+k1))then
	       goto 206
	       endif
               enddo
            inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nf(i+i1,j+j1,k+k1)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=(k+k1-1)*ds
  206 damy=0 
          endif
         enddo
        endif
       endif
      enddo !k
c
      do is=1,ns
          if(itri.eq.0)then
            i1=0
            j1=0
              if(nsur(i+i1,j+j1,is).gt.n)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nsur(i+i1,j+j1,is)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=sur(i+i1,j+j1,is)
              endif
            i1=1
            j1=0
              if(nsur(i+i1,j+j1,is).gt.n)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nsur(i+i1,j+j1,is)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=sur(i+i1,j+j1,is)
              endif
            i1=0
            j1=1
              if(nsur(i+i1,j+j1,is).gt.n)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nsur(i+i1,j+j1,is)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=sur(i+i1,j+j1,is)
              endif
          else
            i1=1
            j1=0
              if(nsur(i+i1,j+j1,is).gt.n)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nsur(i+i1,j+j1,is)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=sur(i+i1,j+j1,is)
              endif
            i1=0
            j1=1
              if(nsur(i+i1,j+j1,is).gt.n)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nsur(i+i1,j+j1,is)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=sur(i+i1,j+j1,is)
              endif
            i1=1
            j1=1
              if(nsur(i+i1,j+j1,is).gt.n)then
          inum=inum+1
	    if(inum.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            numtmp(inum)=nsur(i+i1,j+j1,is)
            coortmp(inum,1)=(i+i1-1)*ds
            coortmp(inum,2)=(j+j1-1)*ds
            coortmp(inum,3)=sur(i+i1,j+j1,is)
              endif
          endif
      enddo
c
cc
c
      do l2=1,2
      do l1=1,101
	lu(l1,l2)=-1
      enddo
      enddo
      do k=1,nz
	lu(num(i,j,k),2)=k
      enddo
      do k=nz,1,-1
	lu(num(i,j,k),1)=k-1
      enddo
c
      if(ns.ge.2)then
	 do is=1,ns-1
c_________________________________
c there is VOXEL elem in taget-layer
       if(lu(is,1).gt.0)then
c
       itet=0
	 do i1=1,inum
       ix=int((coortmp(i1,1)-ori(1))/ds)
       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
	 ix=ix+1
       endif
       iy=int((coortmp(i1,2)-ori(2))/ds)
       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
	 iy=iy+1
       endif
c    upper boundary of target-layer
       if((coortmp(i1,3).le.sur(ix+1,iy+1,is)+ep).and.
     -    (numtmp(i1).ne.0)
     - .and.(coortmp(i1,3).ge.lu(is,2)*ds-ep+ori(3)))then
	 itet=itet+1
	    if(itet.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
       do ii=1,3
	 cortmp(itet,ii)=coortmp(i1,ii)
	 enddo
       nutmp(itet)=numtmp(i1)
       endif
       enddo
       
       ine=0
       if(itet.gt.0)then
       call deldouble(nt,itet,cortmp,nutmp)
       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)    
       endif

      if(ine.gt.0)then
       do i3=1,ine
cc
       call valest(nt,net,cortmp,cnytmp,i3,ds,key)
        if(key.gt.0)then
ccccc
c	 nfemelem=nfemelem+1
c          do i2=1,4
c	   cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
c          enddo
c	 cnyn(nfemelem,5)=is
c
            nfemtmp(ip)=nfemtmp(ip)+1
!	    if(nfemtmp(ip).eq.(nx*ny*nt/4-1))then
            if(nfemtmp(ip).eq.(cnyntmpsize-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
            do i2=1,4
            cnyntmp(ip,nfemtmp(ip),i2)=nutmp(cnytmp(i3,i2))
            enddo
            cnyntmp(ip,nfemtmp(ip),5)=is
ccccc
        endif
       enddo
      endif !ine.gt.0
c
       itet=0
	 do i1=1,inum
       ix=int((coortmp(i1,1)-ori(1))/ds)
       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
	 ix=ix+1
       endif
       iy=int((coortmp(i1,2)-ori(2))/ds)
       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
	 iy=iy+1
       endif
c    lower boundary of target-layer
       if((coortmp(i1,3).ge.sur(ix+1,iy+1,is+1)-ep).and.
     -    (numtmp(i1).ne.0)
     - .and.(coortmp(i1,3).le.lu(is,1)*ds+ep+ori(3)))then
	 itet=itet+1
	    if(itet.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
       do ii=1,3
	 cortmp(itet,ii)=coortmp(i1,ii)
	 enddo
       nutmp(itet)=numtmp(i1)
       endif
       enddo

       ine=0
       if(itet.gt.0)then
       call deldouble(nt,itet,cortmp,nutmp)
       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)    
       endif

      if(ine.gt.0)then
	 do i3=1,ine
       call valest(nt,net,cortmp,cnytmp,i3,ds,key)
      if(key.gt.0)then
cccccc
c	nfemelem=nfemelem+1
c       do i2=1,4
c	 cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
c       enddo
c	 cnyn(nfemelem,5)=is
            nfemtmp(ip)=nfemtmp(ip)+1
!          if(nfemtmp(ip).eq.(nx*ny*nt/4-1))then
          if(nfemtmp(ip).eq.(cnyntmpsize-1))then
          write(*,*) 'enlarge cnytmpsize nt'
          stop
          endif
            do i2=1,4
            cnyntmp(ip,nfemtmp(ip),i2)=nutmp(cnytmp(i3,i2))
            enddo
            cnyntmp(ip,nfemtmp(ip),5)=is
cccccc
       endif
       enddo
      endif
c
c_________________________________
c there is NO VOXEL elem in taget-layer
	 else
c
       itet=0
	 do i1=1,inum
c
       ix=int((coortmp(i1,1)-ori(1))/ds)
       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
	 ix=ix+1
       endif
       iy=int((coortmp(i1,2)-ori(2))/ds)
       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
	 iy=iy+1
       endif
c
       if((coortmp(i1,3).le.sur(ix+1,iy+1,is)+ep).and.
     -    (numtmp(i1).ne.0)
     - .and.(coortmp(i1,3).ge.sur(ix+1,iy+1,is+1)-ep))then
	 itet=itet+1
	    if(itet.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
       do ii=1,3
	 cortmp(itet,ii)=coortmp(i1,ii)
	 enddo
       nutmp(itet)=numtmp(i1)
       endif
       enddo
c    
      ndest(1,1)=
     -abs(nsur(i  ,j  ,is)-nsur(i  ,j  ,is+1))+
     -abs(nsur(i+1,j  ,is)-nsur(i+1,j  ,is+1))+
     -abs(nsur(i  ,j+1,is)-nsur(i  ,j+1,is+1))
      ndest(2,1)=
     -abs(nsur(i  ,j  ,is)-nsur(i  ,j  ,is+1))+
     -abs(nsur(i+1,j  ,is)-nsur(i+1,j  ,is+1))+
     -abs(nsur(i+1,j+1,is)-nsur(i+1,j+1,is+1))
      ndest(2,2)=
     -abs(nsur(i+1,j+1,is)-nsur(i+1,j+1,is+1))+
     -abs(nsur(i+1,j  ,is)-nsur(i+1,j  ,is+1))+
     -abs(nsur(i  ,j+1,is)-nsur(i  ,j+1,is+1))
      ndest(1,2)=
     -abs(nsur(i  ,j  ,is)-nsur(i  ,j  ,is+1))+
     -abs(nsur(i+1,j+1,is)-nsur(i+1,j+1,is+1))+
     -abs(nsur(i  ,j+1,is)-nsur(i  ,j+1,is+1))
      il=0
	do j2=1,2
	do i2=1,2
	if(ndest(i2,j2).eq.0)then
      il=il+1
	rlist(il)=nsur(i-1+i2,j-1+j2,is)
	endif
      enddo
      enddo
      if(il.gt.0)then
      do i2=1,il
      do i1=1,itet 
      if(rlist(i2).eq.nutmp(i1))then
!      nutmp(i1)=0
	endif
      enddo
      enddo
	endif

       ine=0
       if(itet.gt.0)then
       call deldouble(nt,itet,cortmp,nutmp)
       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine)    
       endif
       
       if(ine.gt.0)then
        do i3=1,ine
         call valest(nt,net,cortmp,cnytmp,i3,ds,key)
         if(key.eq.1)then
ccccc
c	  nfemelem=nfemelem+1
c           do i2=1,4
c	    cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
c           enddo
c	  cnyn(nfemelem,5)=is
            nfemtmp(ip)=nfemtmp(ip)+1
!          if(nfemtmp(ip).eq.(nx*ny*nt/4-1))then
          if(nfemtmp(ip).eq.(cnyntmpsize-1))then
          write(*,*) 'enlarge cnytmpsize nt'
          stop
          endif
            do i2=1,4
            cnyntmp(ip,nfemtmp(ip),i2)=nutmp(cnytmp(i3,i2))
            enddo
            cnyntmp(ip,nfemtmp(ip),5)=is
ccccc
         endif
        enddo
       endif

      endif
       enddo
      endif
c_________________________________
c
cc
c lowest layer
       itet=0
	 do i1=1,inum
       ix=int((coortmp(i1,1)-ori(1))/ds)
       if(coortmp(i1,1)-ori(1)-ix*ds.gt.0.5*ds)then
	 ix=ix+1
       endif
       iy=int((coortmp(i1,2)-ori(2))/ds)
       if(coortmp(i1,2)-ori(2)-iy*ds.gt.0.5*ds)then
	 iy=iy+1
       endif
       if((coortmp(i1,3).le.sur(ix+1,iy+1,ns)+ep).and.
     -    (numtmp(i1).ne.0))then
	 itet=itet+1
	    if(itet.eq.(nt-1))then
          write(*,*) 'enlarge nt'
          stop
          endif
       do ii=1,3
	 cortmp(itet,ii)=coortmp(i1,ii)
	 enddo
       nutmp(itet)=numtmp(i1)
       endif
       enddo

       ine=0
       if(itet.gt.0)then
       call deldouble(nt,itet,cortmp,nutmp)
       call Delaunay(nt,net,itet,cortmp,cnytmp,in,ine) 
       endif
           
      if(ine.gt.0)then
	 do i3=1,ine
       call valest(nt,net,cortmp,cnytmp,i3,ds,key)
       if(key.gt.0)then
ccccc
c       nfemelem=nfemelem+1
c       do i2=1,4
c	 cnyn(nfemelem,i2)=nutmp(cnytmp(i3,i2))
c       enddo
c	 cnyn(nfemelem,5)=ns
            nfemtmp(ip)=nfemtmp(ip)+1
!          if(nfemtmp(ip).eq.(nx*ny*nt/4-1))then
            if(nfemtmp(ip).eq.(cnyntmpsize-1))then
          write(*,*) 'enlarge cnytmpsize nt'
          stop
          endif
            do i2=1,4
            cnyntmp(ip,nfemtmp(ip),i2)=nutmp(cnytmp(i3,i2))
            enddo
            cnyntmp(ip,nfemtmp(ip),5)=is
ccccc
       endif
       enddo
       endif

      enddo ! itri
      enddo ! i
      enddo ! j
c
!$OMP end parallel

      do ip=1,np
      if(nfemtmp(ip).gt.0)then
      do i=1,nfemtmp(ip)
      nfemelem=nfemelem+1
      if(nfemelem.eq.(nx*ny*nt-1))then
      write(*,*) 'enlarge nt'
      stop
      endif
      do i2=1,5
      cnyn(nfemelem,i2)=cnyntmp(ip,i,i2)
      enddo 
      enddo
      endif
      enddo
      
      end
c_______________________________________________________________________
      subroutine valest(nt,net,coor,cny,ie,ds,key)
	real*8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,vol
      real*8 coor(nt,3),ds
	integer cny(net,4) 

      nod1=cny(ie,1)
      nod2=cny(ie,2)
      nod3=cny(ie,3)
      nod4=cny(ie,4)
      x1=coor(nod1,1)
      y1=coor(nod1,2)
      z1=coor(nod1,3)
      x2=coor(nod2,1)
      y2=coor(nod2,2)
      z2=coor(nod2,3)
      x3=coor(nod3,1)
      y3=coor(nod3,2)
      z3=coor(nod3,3)
      x4=coor(nod4,1)
      y4=coor(nod4,2)
      z4=coor(nod4,3)

	x2=x2-x1
	y2=y2-y1
	z2=z2-z1
	x3=x3-x1
	y3=y3-y1
	z3=z3-z1
	x4=x4-x1
	y4=y4-y1
	z4=z4-z1
c
      vol= (-(x4*y3*z2) + x3*y4*z2 + x4*y2*z3 - x2*y4*z3 - 
     -  x3*y2*z4 + x2*y3*z4) 
c_____add_____
	vol=vol/6.
c_____add_____

c
      if(vol.lt.ds*ds*ds*0.001)then
	key=0
	else
	key=1
	endif
      end
c_______________________________________________________________________
      subroutine deldouble(nt,itet,cortmp,nutmp)
	real*8 cortmp(nt,3),cortmpd(nt,3)
	integer nutmp(nt),nutmpd(nt)
c____modifyed_________________
c--- comment out-----
c      nutmpd(1)=nutmp(1)
c      cortmpd(1,1)=cortmp(1,1)
c      cortmpd(1,2)=cortmp(1,2)
c      cortmpd(1,3)=cortmp(1,3)
c--- comment out-----
c--- added-----------
	ist=0
      do i=1,itet
	ist=ist+1
	if(nutmp(i).ne.0)then
	 nutmpd(1)=nutmp(i)
	 cortmpd(1,1)=cortmp(i,1)
	 cortmpd(1,2)=cortmp(i,2)
	 cortmpd(1,3)=cortmp(i,3)
	 exit
	endif
      enddo
c--- added-----------
c____modifyed_________________
c
      i0=1
c____modifyed_________________
c      do i=1,itet
c      do j=1,i0
      do i=ist,itet
      do j=1,i0
c____modifyed_________________
	if((nutmp(i).eq.nutmpd(j)).or.(nutmp(i).eq.0))then
	goto 101
	endif
	enddo 
      i0=i0+1
	nutmpd(i0)=nutmp(i)
      cortmpd(i0,1)=cortmp(i,1)
      cortmpd(i0,2)=cortmp(i,2)
      cortmpd(i0,3)=cortmp(i,3)   
  101 damy=0
	enddo

	do i=1,i0
	nutmp(i)=nutmpd(i)
      cortmp(i,1)=cortmpd(i,1)
      cortmp(i,2)=cortmpd(i,2)
      cortmp(i,3)=cortmpd(i,3)   
	enddo
      itet=i0

      end
c_______________________________________________________________________
      subroutine orgdeldouble(nt,itet,cortmp,nutmp)
	real*8 cortmp(nt,3),cortmpd(nt,3)
	integer nutmp(nt),nutmpd(nt)

      nutmpd(1)=nutmp(1)
      cortmpd(1,1)=cortmp(1,1)
      cortmpd(1,2)=cortmp(1,2)
      cortmpd(1,3)=cortmp(1,3)

      i0=1
      do i=1,itet
      do j=1,i0
	if((nutmp(i).eq.nutmpd(j)).or.(nutmp(i).eq.0))then
	goto 101
	endif
	enddo 
      i0=i0+1
	nutmpd(i0)=nutmp(i)
      cortmpd(i0,1)=cortmp(i,1)
      cortmpd(i0,2)=cortmp(i,2)
      cortmpd(i0,3)=cortmp(i,3)   
  101 damy=0
	enddo

	do i=1,i0
	nutmp(i)=nutmpd(i)
      cortmp(i,1)=cortmpd(i,1)
      cortmp(i,2)=cortmpd(i,2)
      cortmp(i,3)=cortmpd(i,3)   
	enddo
      itet=i0

      end
c_______________________________________________________________________
      subroutine output
     -(nx,ny,nz,ns,nk,nt,ds,ori,nf,nsur,n,sur,nfemelem,cnyn,num,oris,id)
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
	integer nf(nx+1,ny+1,nz+1)
      integer nsur(nx+1,ny+1,ns)
      real*8 ds,ori(3),sur(nx+1,ny+1,ns),oris(3)
	integer cnyn(nx*ny*nt,5)
	integer notmp(8)
      integer nfemelem
      integer fid
      character filename*50
c
c      write(41,'(a)') 'ds'
c      write(41,*) ds
c      write(41,'(a)') 'nx, ny, nz'
c      write(41,*) nx,ny,nz
c      write(41,'(a)') 'ori(1), ori(2), ori(3)'
c      write(41,*) ori(1),ori(2),ori(3)
c
!      write(31+3*id,300) '"bar"', 3, ' Tetrahedra', 4
!  300 format('MESH',a10,'  dimension',i5,'  ElemType',a11,'  Nnode',i5)
c Node
!      write(31+3*id,'(a)') 'Coordinates'
!      write(31+3*id,'(a)')'# node_number     coor_x    coor_y    coor_z'
c      write(42,'(a)') 'Coordinates'
c      write(42,'(a)') '# node_number     coor_x    coor_y    coor_z'
      itmp=0
!DEC$ NOVECTOR
      do k=1,nz+1
!DEC$ NOVECTOR
      do j=1,ny+1
!DEC$ NOVECTOR
      do i=1,nx+1 
	if(nf(i,j,k).gt.0)then
      itmp=itmp+1
!      write(31+3*id,'(i10,3f16.6)') nf(i,j,k),
!     - ds*(i-1)+ori(1)+oris(1),
!     - ds*(j-1)+ori(2)+oris(2),ds*(k-1)+ori(3)+oris(3)
c	write(42,'(i10,3f16.6)') nf(i,j,k),
c     - ds*(i-1)+ori(1)+oris(1),
c     - ds*(j-1)+ori(2)+oris(2),ds*(k-1)+ori(3)+oris(3)
      endif
	enddo
      enddo
      enddo

      if(ns.ge.2)then
!DEC$ NOVECTOR
	do is=1,ns-1
!DEC$ NOVECTOR
	do j=1,ny+1
!DEC$ NOVECTOR
	do i=1,nx+1
      if(nsur(i,j,is).eq.nsur(i,j,is+1))then
      nsur(i,j,is)=0
	endif
      enddo
      enddo
      enddo
	endif

!DEC$ NOVECTOR
	do is=1,ns
!DEC$ NOVECTOR
	do j=1,ny+1
!DEC$ NOVECTOR
	do i=1,nx+1
      if(nsur(i,j,is).gt.n)then
      itmp=itmp+1
!      write(31+3*id,'(i10,3f16.6)') nsur(i,j,is),
!     - ds*(i-1)+ori(1)+oris(1),
!     - ds*(j-1)+ori(2)+oris(2),sur(i,j,is)+oris(3)
c      write(42,'(i10,3f16.6)') nsur(i,j,is),
c     - ds*(i-1)+ori(1)+oris(1),
c     - ds*(j-1)+ori(2)+oris(2),sur(i,j,is)+oris(3)
      endif
      enddo
      enddo
      enddo
!      write(31+3*id,'(a)') 'end coordinates'
!      write(31+3*id,*)
c
      write(32+3*id,'(a)') 'num of node'
      write(32+3*id,*) itmp
cc
!      write(31+3*id,'(a)') 'Elements'
!      write(31+3*id,'(a)') '# element node_1 node_2 node_3 
!     -node_4 material_number'
c      write(43,'(a)') 'Elements'
c      write(43,'(a)') '# element node_1 node_2 node_3 
c     -node_4 material_number'
!!DEC$ NOVECTOR
!      do i=1,nfemelem
!      write(31+3*id,802)
!     -  i,cnyn(i,1),cnyn(i,2),cnyn(i,3),cnyn(i,4),cnyn(i,5)
c     -  i,cnyn(i,1),cnyn(i,3),cnyn(i,2),cnyn(i,4),cnyn(i,5)
c       write(43,802) i,cnyn(i,1),cnyn(i,2),cnyn(i,3),cnyn(i,4),cnyn(i,5)
!  802 format(10i10)
!      enddo
c      write(31+3*id,'(a)') 'end elements'

c      write(32+3*id,'(a)') 'num of nfem elem'
c      write(32+3*id,*) nfemelem


cc
c VFEM
      inum2=nfemelem
      inumvox=0
c	write(31+3*id,*)''
c	write(31+3*id,*) 'MESH     "bar"  dimension    3  ElemType   Hexahedra
c     -   Nnode    8'
c	write(31+3*id,*)'Elements'
c	write(31+3*id,*)'# element node_1 - node_8 material_number'
c	write(44,*)'Elements'
c	write(44,*)'# element node_1 - node_8 material_number'
c lowest level VFEM 
!DEC$ NOVECTOR
       do k=1,2*nz-1,2     
!DEC$ NOVECTOR
       do j=1,2*ny-1,2   
!DEC$ NOVECTOR
       do i=1,2*nx-1,2
	 if((num(i,j,k).ne.0).and.(num(i,j,k).lt.100))then
c	 inum2=inum2+1
c    	 write(31+3*id,200) inum2,
c     - nf((i+1)/2  ,(j+1)/2  ,(k+1)/2  ),
c     - nf((i+1)/2+1,(j+1)/2  ,(k+1)/2  ),
c     - nf((i+1)/2+1,(j+1)/2+1,(k+1)/2  ),
c    - nf((i+1)/2  ,(j+1)/2+1,(k+1)/2  ),
c     - nf((i+1)/2  ,(j+1)/2  ,(k+1)/2+1),
c     - nf((i+1)/2+1,(j+1)/2  ,(k+1)/2+1),
c     - nf((i+1)/2+1,(j+1)/2+1,(k+1)/2+1),
c     - nf((i+1)/2  ,(j+1)/2+1,(k+1)/2+1),(num(i,j,k)-1)*nk+1

        notmp(1)=nf((i+1)/2  ,(j+1)/2  ,(k+1)/2  )
        notmp(2)=nf((i+1)/2+1,(j+1)/2  ,(k+1)/2  )
        notmp(3)=nf((i+1)/2+1,(j+1)/2+1,(k+1)/2  )
        notmp(4)=nf((i+1)/2  ,(j+1)/2+1,(k+1)/2  )
        notmp(5)=nf((i+1)/2  ,(j+1)/2  ,(k+1)/2+1)
        notmp(6)=nf((i+1)/2+1,(j+1)/2  ,(k+1)/2+1)
        notmp(7)=nf((i+1)/2+1,(j+1)/2+1,(k+1)/2+1)
        notmp(8)=nf((i+1)/2  ,(j+1)/2+1,(k+1)/2+1)

	 inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(7),notmp(5),notmp(6),num(i,j,k)
       inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(3),notmp(1),notmp(5),notmp(4),num(i,j,k)
       inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(2),notmp(5),notmp(1),num(i,j,k)
       inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(6),notmp(5),notmp(2),num(i,j,k)
       inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(4),notmp(5),notmp(7),num(i,j,k)
       inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(4),notmp(5),notmp(7),notmp(8),num(i,j,k)
       inumvox=inumvox+1
       endif	        
       enddo
       enddo
       enddo
c multi level VFEM 
!DEC$ NOVECTOR
      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)
!DEC$ NOVECTOR
       do k=1,nz/nd     
!DEC$ NOVECTOR
       do j=1,ny/nd     
!DEC$ NOVECTOR
       do i=1,nx/nd
       nex=nd+2*nd*(i-1)
       ney=nd+2*nd*(j-1)
       nez=nd+2*nd*(k-1)            
	 if(num(nex,ney,nez).ne.0)then
c	 inum2=inum2+1
c    	 write(31+3*id,200) inum2,
c     - nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2),
c     - nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2),
c     - nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2+nd),
c     - (num(nex,ney,nez)-1)*nk+ik
c
        notmp(1)=nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2)
        notmp(2)=nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2)
        notmp(3)=nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2)
        notmp(4)=nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2)
        notmp(5)=nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2+nd)
        notmp(6)=nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2+nd)
        notmp(7)=nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2+nd)
        notmp(8)=nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2+nd)

	 inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(7),notmp(5),notmp(6),num(nex,ney,nez)
	 inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(3),notmp(1),notmp(5),notmp(4),num(nex,ney,nez)
	 inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(2),notmp(5),notmp(1),num(nex,ney,nez)
	 inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(6),notmp(5),notmp(2),num(nex,ney,nez)
	 inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(4),notmp(5),notmp(7),num(nex,ney,nez)
	 inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(4),notmp(5),notmp(7),notmp(8),num(nex,ney,nez)
       inumvox=inumvox+1

       endif
       enddo
       enddo
       enddo    
      enddo    

  200 format(11i10)
!      write(31+3*id,'(a)') 'end elements'
      write(32+3*id,'(a)') 'num of tet elem'
      write(32+3*id,*) nfemelem+6*inumvox
      write(32+3*id,'(a)') 'num of vox elem'
      write(32+3*id,*) 0
      write(32+3*id,'(a)') 'num of material properties'
      write(32+3*id,*) ns


c      write(32+3*id,'(a)') 'num of octree level'
c      write(32+3*id,*) nk
c      write(32+3*id,'(a)') 'num of materials'
c      write(32+3*id,*) ns
c      write(32+3*id,'(a)') 'base ds'
c      write(32+3*id,*) ds


!      close(31+3*id)
      close(32+3*id)
c      close(42)
c      close(43)
c      close(44)  
      write(*,*)'nodes,tet4elems,vox4elems'
      write(*,*) itmp,nfemelem,inumvox

!      write(filename,'(a17)') 'data/tet4model.h5'
      write(filename,'(a21)') 'data/tet4vox8model.h5'
      call set_null_char(filename,21)
      fid=0
      call hdf_create_file(filename,fid)
      call output_hdf5_nodes
     -(nx,ny,nz,ns,nk,ds,ori,nf,nsur,n,sur,oris,itmp,fid)
      call output_hdf5_elements
     -(nx,ny,nz,nk,nt,nf,nfemelem,cnyn,num,inum2,fid)
!      call output_hdf5_elements_vox
!     -(nx,ny,nz,nk,nt,nf,nfemelem,cnyn,num,inumvox,fid)
      call hdf_close_file(fid)
      end
c_______________________________________________________________________
      subroutine output_hdf5_nodes
     -(nx,ny,nz,ns,nk,ds,ori,nf,nsur,n,sur,oris,nnodes,fid)
!$      use omp_lib
      implicit none
      integer nx,ny,nz,ns,nk,n,nnodes,itmp,i,j,k,l,is,fid
	integer nf(nx+1,ny+1,nz+1)
      integer nsur(nx+1,ny+1,ns)
      real*8 ds,ori(3),sur(nx+1,ny+1,ns),oris(3)
      character dataname*50
c
      real*8 t1
      real*8 coordbuf(nnodes*3)
      integer*8 coordbufsize,filenamesize
c
      itmp=0
!DEC$ NOVECTOR
      do k=1,nz+1
!DEC$ NOVECTOR
      do j=1,ny+1
!DEC$ NOVECTOR
      do i=1,nx+1 
	if(nf(i,j,k).gt.0)then
      itmp=itmp+1
!      write(31+3*id,'(i10,3f16.6)') nf(i,j,k),
!     - ds*(i-1)+ori(1)+oris(1),
!     - ds*(j-1)+ori(2)+oris(2),ds*(k-1)+ori(3)+oris(3)
      l=nf(i,j,k)
!      if(l.ne.itmp) write(*,*)'1 l.ne.itmp',l,itmp
      coordbuf(3*(l-1)+1)=ds*(i-1)+ori(1)+oris(1)
      coordbuf(3*(l-1)+2)=ds*(j-1)+ori(2)+oris(2)
      coordbuf(3*(l-1)+3)=ds*(k-1)+ori(3)+oris(3)
      endif
	enddo
      enddo
      enddo

!DEC$ NOVECTOR
	do is=1,ns
!DEC$ NOVECTOR
	do j=1,ny+1
!DEC$ NOVECTOR
	do i=1,nx+1
      if(nsur(i,j,is).gt.n)then
      itmp=itmp+1
!      write(31+3*id,'(i10,3f16.6)') nsur(i,j,is),
!     - ds*(i-1)+ori(1)+oris(1),
!     - ds*(j-1)+ori(2)+oris(2),sur(i,j,is)+oris(3)
      l=nsur(i,j,is)
!      if(l.ne.itmp) write(*,*)'2 l.ne.itmp',l,itmp
      coordbuf(3*(l-1)+1)=ds*(i-1)+ori(1)+oris(1)
      coordbuf(3*(l-1)+2)=ds*(j-1)+ori(2)+oris(2)
      coordbuf(3*(l-1)+3)=sur(i,j,is)+oris(3)
      endif
      enddo
      enddo
      enddo
      
!      t1=omp_get_wtime()
!      write(filename,'(a27)')
!     - './data/tet4model3d.node.msh'
!      open(31,file=filename,status='unknown')
!      do i=1,nnodes
!        write(31,'(i10,3f16.6)') i,coordbuf(3*(i-1)+1),
!     -  coordbuf(3*(i-1)+2),coordbuf(3*(i-1)+3)
!      enddo
!      close(31)
!      t1=omp_get_wtime()-t1
!      write(*,*)'write node ASCII took',t1
      
      write(dataname,'(a5)') '/coor'
      call set_null_char(dataname,5)
      coordbufsize=nnodes*3
#ifdef NSQUASH
      do i=1,nnodes
      coordbuf(3*(i-1)+3)=coordbuf(3*(i-1)+3)*NSQUASH
      enddo
#endif
      call hdf_write_double_array
     - (fid,dataname,coordbuf,coordbufsize)
     
      end
c_______________________________________________________________________
      subroutine output_hdf5_elements
     -(nx,ny,nz,nk,nt,nf,nfemelem,cnyn,num,nelems,fid)
!$      use omp_lib
      implicit none
      integer i,j,k,ik,nx,ny,nz,nk,nt,nd,nd2,inum2,nex,ney,nez,nelems
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
	integer nf(nx+1,ny+1,nz+1)
	integer cnyn(nx*ny*nt,5)
	integer notmp(8)
      integer nfemelem
      character dataname*50
!
      real*8 t1
      integer*8 connbuf(5*nelems)
      integer*8 connbufsize,filenamesize
      integer*8 fid
!

!DEC$ NOVECTOR
      do i=1,nfemelem
!      write(31+3*id,802)
!     -  i,cnyn(i,1),cnyn(i,2),cnyn(i,3),cnyn(i,4),cnyn(i,5)
      connbuf(5*(i-1)+1)=cnyn(i,1)
      connbuf(5*(i-1)+2)=cnyn(i,2)
      connbuf(5*(i-1)+3)=cnyn(i,3)
      connbuf(5*(i-1)+4)=cnyn(i,4)
      connbuf(5*(i-1)+5)=cnyn(i,5)
      enddo
c
c VFEM
      inum2=nfemelem
c lowest level VFEM 
!DEC$ NOVECTOR
       do k=1,2*nz-1,2     
!DEC$ NOVECTOR
       do j=1,2*ny-1,2   
!DEC$ NOVECTOR
       do i=1,2*nx-1,2
	 if((num(i,j,k).ne.0).and.(num(i,j,k).lt.100))then
c	 inum2=inum2+1
c    	 write(31+3*id,200) inum2,
c     - nf((i+1)/2  ,(j+1)/2  ,(k+1)/2  ),
c     - nf((i+1)/2+1,(j+1)/2  ,(k+1)/2  ),
c     - nf((i+1)/2+1,(j+1)/2+1,(k+1)/2  ),
c    - nf((i+1)/2  ,(j+1)/2+1,(k+1)/2  ),
c     - nf((i+1)/2  ,(j+1)/2  ,(k+1)/2+1),
c     - nf((i+1)/2+1,(j+1)/2  ,(k+1)/2+1),
c     - nf((i+1)/2+1,(j+1)/2+1,(k+1)/2+1),
c     - nf((i+1)/2  ,(j+1)/2+1,(k+1)/2+1),(num(i,j,k)-1)*nk+1

        notmp(1)=nf((i+1)/2  ,(j+1)/2  ,(k+1)/2  )
        notmp(2)=nf((i+1)/2+1,(j+1)/2  ,(k+1)/2  )
        notmp(4)=nf((i+1)/2+1,(j+1)/2+1,(k+1)/2  )
        notmp(3)=nf((i+1)/2  ,(j+1)/2+1,(k+1)/2  )
        notmp(5)=nf((i+1)/2  ,(j+1)/2  ,(k+1)/2+1)
        notmp(6)=nf((i+1)/2+1,(j+1)/2  ,(k+1)/2+1)
        notmp(8)=nf((i+1)/2+1,(j+1)/2+1,(k+1)/2+1)
        notmp(7)=nf((i+1)/2  ,(j+1)/2+1,(k+1)/2+1)

       inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(7),notmp(5),notmp(6),num(i,j,k)
       connbuf(5*(inum2-1)+1)=notmp(2)!3)
       connbuf(5*(inum2-1)+2)=notmp(7)!7)
       connbuf(5*(inum2-1)+3)=notmp(3)!5)
       connbuf(5*(inum2-1)+4)=notmp(5)!6)
       connbuf(5*(inum2-1)+5)=num(i,j,k)
       
       inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(3),notmp(1),notmp(5),notmp(4),num(i,j,k)
       connbuf(5*(inum2-1)+1)=notmp(4)!3)
       connbuf(5*(inum2-1)+2)=notmp(6)!1)
       connbuf(5*(inum2-1)+3)=notmp(8)!5)
       connbuf(5*(inum2-1)+4)=notmp(7)!4)
       connbuf(5*(inum2-1)+5)=num(i,j,k)
       
       inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(2),notmp(5),notmp(1),num(i,j,k)
       connbuf(5*(inum2-1)+1)=notmp(2)!3)
       connbuf(5*(inum2-1)+2)=notmp(7)!2)
       connbuf(5*(inum2-1)+3)=notmp(4)!5)
       connbuf(5*(inum2-1)+4)=notmp(3)!1)
       connbuf(5*(inum2-1)+5)=num(i,j,k)
       
       inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(6),notmp(5),notmp(2),num(i,j,k)
       connbuf(5*(inum2-1)+1)=notmp(3)!3)
       connbuf(5*(inum2-1)+2)=notmp(2)!6)
       connbuf(5*(inum2-1)+3)=notmp(5)!5)
       connbuf(5*(inum2-1)+4)=notmp(1)!2)
       connbuf(5*(inum2-1)+5)=num(i,j,k)
       
       inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(4),notmp(5),notmp(7),num(i,j,k)
       connbuf(5*(inum2-1)+1)=notmp(4)!3)
       connbuf(5*(inum2-1)+2)=notmp(6)!4)
       connbuf(5*(inum2-1)+3)=notmp(7)!5)
       connbuf(5*(inum2-1)+4)=notmp(2)!7)
       connbuf(5*(inum2-1)+5)=num(i,j,k)
       
       inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(4),notmp(5),notmp(7),notmp(8),num(i,j,k)
       connbuf(5*(inum2-1)+1)=notmp(6)!4)
       connbuf(5*(inum2-1)+2)=notmp(7)!5)
       connbuf(5*(inum2-1)+3)=notmp(2)!7)
       connbuf(5*(inum2-1)+4)=notmp(5)!8)
       connbuf(5*(inum2-1)+5)=num(i,j,k)

       endif	        
       enddo
       enddo
       enddo
c multi level VFEM 
!DEC$ NOVECTOR
      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)
!DEC$ NOVECTOR
       do k=1,nz/nd     
!DEC$ NOVECTOR
       do j=1,ny/nd     
!DEC$ NOVECTOR
       do i=1,nx/nd
       nex=nd+2*nd*(i-1)
       ney=nd+2*nd*(j-1)
       nez=nd+2*nd*(k-1)            
	 if(num(nex,ney,nez).ne.0)then
c	 inum2=inum2+1
c    	 write(31+3*id,200) inum2,
c     - nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2),
c     - nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2),
c     - nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2+nd),
c     - (num(nex,ney,nez)-1)*nk+ik
c
        notmp(1)=nf(nex/2+1-nd2   ,ney/2+1-nd2   ,nez/2+1-nd2)
        notmp(2)=nf(nex/2+1-nd2+nd,ney/2+1-nd2   ,nez/2+1-nd2)
        notmp(4)=nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2)
        notmp(3)=nf(nex/2+1-nd2   ,ney/2+1-nd2+nd,nez/2+1-nd2)
        notmp(5)=nf(nex/2+1-nd2   ,ney/2+1-nd2   ,nez/2+1-nd2+nd)
        notmp(6)=nf(nex/2+1-nd2+nd,ney/2+1-nd2   ,nez/2+1-nd2+nd)
        notmp(8)=nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2+nd)
        notmp(7)=nf(nex/2+1-nd2   ,ney/2+1-nd2+nd,nez/2+1-nd2+nd)
!       notmp=1

       inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(7),notmp(5),notmp(6),num(nex,ney,nez)
       connbuf(5*(inum2-1)+1)=notmp(2)
       connbuf(5*(inum2-1)+2)=notmp(7)
       connbuf(5*(inum2-1)+3)=notmp(3)
       connbuf(5*(inum2-1)+4)=notmp(5)
       connbuf(5*(inum2-1)+5)=num(nex,ney,nez)
       
       inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(3),notmp(1),notmp(5),notmp(4),num(nex,ney,nez)
       connbuf(5*(inum2-1)+1)=notmp(4)
       connbuf(5*(inum2-1)+2)=notmp(6)
       connbuf(5*(inum2-1)+3)=notmp(8)
       connbuf(5*(inum2-1)+4)=notmp(7)
       connbuf(5*(inum2-1)+5)=num(nex,ney,nez)
       
       inum2=inum2+1
!    	 write(31+3*id,200) 
!     - inum2,notmp(3),notmp(2),notmp(5),notmp(1),num(nex,ney,nez)
       connbuf(5*(inum2-1)+1)=notmp(2)
       connbuf(5*(inum2-1)+2)=notmp(7)
       connbuf(5*(inum2-1)+3)=notmp(4)
       connbuf(5*(inum2-1)+4)=notmp(3)
       connbuf(5*(inum2-1)+5)=num(nex,ney,nez)
       
       inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(6),notmp(5),notmp(2),num(nex,ney,nez)
       connbuf(5*(inum2-1)+1)=notmp(3)
       connbuf(5*(inum2-1)+2)=notmp(2)
       connbuf(5*(inum2-1)+3)=notmp(5)
       connbuf(5*(inum2-1)+4)=notmp(1)
       connbuf(5*(inum2-1)+5)=num(nex,ney,nez)
       
        inum2=inum2+1
!      write(31+3*id,200) 
!     - inum2,notmp(3),notmp(4),notmp(5),notmp(7),num(nex,ney,nez)
       connbuf(5*(inum2-1)+1)=notmp(4)
       connbuf(5*(inum2-1)+2)=notmp(6)
       connbuf(5*(inum2-1)+3)=notmp(7)
       connbuf(5*(inum2-1)+4)=notmp(2)
       connbuf(5*(inum2-1)+5)=num(nex,ney,nez)
       
       inum2=inum2+1
!    	 write(31+3*id,200)
!     - inum2,notmp(4),notmp(5),notmp(7),notmp(8),num(nex,ney,nez)
       connbuf(5*(inum2-1)+1)=notmp(6)
       connbuf(5*(inum2-1)+2)=notmp(7)
       connbuf(5*(inum2-1)+3)=notmp(2)
       connbuf(5*(inum2-1)+4)=notmp(5)
       connbuf(5*(inum2-1)+5)=num(nex,ney,nez)

       endif
       enddo
       enddo
       enddo    
      enddo    

!      t1=omp_get_wtime()
!      write(filename,'(a27)')
!     - './data/tet4model3d.elem.msh'
!      open(31,file=filename,status='unknown')
!      do i=1,nelems
!        write(31,200) i,connbuf(5*(i-1)+1),
!     -  connbuf(5*(i-1)+2),connbuf(5*(i-1)+3),
!     -  connbuf(5*(i-1)+4),connbuf(5*(i-1)+5)
!      enddo
!  200 format(11i10)
!      close(31)
!      t1=omp_get_wtime()-t1
!      write(*,*)'write elem ASCII took',t1
      
!      write(dataname,'(a5)') '/conn'
!      call set_null_char(dataname,5)
!      connbufsize=nelems*5
!      call hdf_write_long_array
!     - (fid,dataname,connbuf,connbufsize)
      write(dataname,'(a10)') '/conn_tet4'
      call set_null_char(dataname,10)
      connbufsize=nelems*5
      call hdf_write_long_array
     - (fid,dataname,connbuf,connbufsize)
      write(dataname,'(a10)') '/conn_vox8'
      call set_null_char(dataname,10)
      connbufsize=0
      call hdf_write_long_array
     - (fid,dataname,connbuf,connbufsize)
      end
c_______________________________________________________________________
      subroutine output_hdf5_elements_vox
     -(nx,ny,nz,nk,nt,nf,nfemelem,cnyn,num,nvoxelems,fid)
!$      use omp_lib
      implicit none
      integer i,j,k,ik,nx,ny,nz,nk,nt,nd,nd2
      integer inum2,nex,ney,nez,nvoxelems
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
      integer nf(nx+1,ny+1,nz+1)
      integer cnyn(nx*ny*nt,5)
      integer notmp(8)
      integer nfemelem
      character dataname*50
!
      real*8 t1
      integer*8 connbuf(5*nfemelem)
      integer*8 connbufvox(9*nvoxelems)
      integer*8 connbufsize,filenamesize
      integer*8 fid,jj
!

!DEC$ NOVECTOR
      do i=1,nfemelem
!      write(31+3*id,802)
!     -  i,cnyn(i,1),cnyn(i,2),cnyn(i,3),cnyn(i,4),cnyn(i,5)
      connbuf(5*(i-1)+1)=cnyn(i,1)
      connbuf(5*(i-1)+2)=cnyn(i,2)
      connbuf(5*(i-1)+3)=cnyn(i,3)
      connbuf(5*(i-1)+4)=cnyn(i,4)
      connbuf(5*(i-1)+5)=cnyn(i,5)
      enddo
c
c VFEM
      inum2=0
c lowest level VFEM 
!DEC$ NOVECTOR
       do k=1,2*nz-1,2     
!DEC$ NOVECTOR
       do j=1,2*ny-1,2   
!DEC$ NOVECTOR
       do i=1,2*nx-1,2
       if((num(i,j,k).ne.0).and.(num(i,j,k).lt.100))then
c	 inum2=inum2+1
c    	 write(31+3*id,200) inum2,
c     - nf((i+1)/2  ,(j+1)/2  ,(k+1)/2  ),
c     - nf((i+1)/2+1,(j+1)/2  ,(k+1)/2  ),
c     - nf((i+1)/2+1,(j+1)/2+1,(k+1)/2  ),
c    - nf((i+1)/2  ,(j+1)/2+1,(k+1)/2  ),
c     - nf((i+1)/2  ,(j+1)/2  ,(k+1)/2+1),
c     - nf((i+1)/2+1,(j+1)/2  ,(k+1)/2+1),
c     - nf((i+1)/2+1,(j+1)/2+1,(k+1)/2+1),
c     - nf((i+1)/2  ,(j+1)/2+1,(k+1)/2+1),(num(i,j,k)-1)*nk+1

        notmp(1)=nf((i+1)/2  ,(j+1)/2  ,(k+1)/2  )
        notmp(2)=nf((i+1)/2+1,(j+1)/2  ,(k+1)/2  )
        notmp(3)=nf((i+1)/2+1,(j+1)/2+1,(k+1)/2  )
        notmp(4)=nf((i+1)/2  ,(j+1)/2+1,(k+1)/2  )
        notmp(5)=nf((i+1)/2  ,(j+1)/2  ,(k+1)/2+1)
        notmp(6)=nf((i+1)/2+1,(j+1)/2  ,(k+1)/2+1)
        notmp(7)=nf((i+1)/2+1,(j+1)/2+1,(k+1)/2+1)
        notmp(8)=nf((i+1)/2  ,(j+1)/2+1,(k+1)/2+1)
       inum2=inum2+1
       do jj=1,8
       connbufvox(9*(inum2-1)+jj)=notmp(jj)
       enddo
       connbufvox(9*(inum2-1)+9)=num(i,j,k)
       endif 
       enddo
       enddo
       enddo
c multi level VFEM 
!DEC$ NOVECTOR
      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)
!DEC$ NOVECTOR
       do k=1,nz/nd     
!DEC$ NOVECTOR
       do j=1,ny/nd     
!DEC$ NOVECTOR
       do i=1,nx/nd
       nex=nd+2*nd*(i-1)
       ney=nd+2*nd*(j-1)
       nez=nd+2*nd*(k-1)            
	 if(num(nex,ney,nez).ne.0)then
c	 inum2=inum2+1
c    	 write(31+3*id,200) inum2,
c     - nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2),
c     - nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2),
c     - nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2+nd),
c     - nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2+nd),
c     - (num(nex,ney,nez)-1)*nk+ik
c
        notmp(1)=nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2)
        notmp(2)=nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2)
        notmp(3)=nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2)
        notmp(4)=nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2)
        notmp(5)=nf(nex/2+1-nd2,ney/2+1-nd2,nez/2+1-nd2+nd)
        notmp(6)=nf(nex/2+1-nd2+nd,ney/2+1-nd2,nez/2+1-nd2+nd)
        notmp(7)=nf(nex/2+1-nd2+nd,ney/2+1-nd2+nd,nez/2+1-nd2+nd)
        notmp(8)=nf(nex/2+1-nd2,ney/2+1-nd2+nd,nez/2+1-nd2+nd)

       inum2=inum2+1
       do jj=1,8
       connbufvox(9*(inum2-1)+jj)=notmp(jj)
       enddo
       connbufvox(9*(inum2-1)+9)=num(nex,ney,nez)
       endif
       enddo
       enddo
       enddo    
      enddo    

!      t1=omp_get_wtime()
!      write(filename,'(a27)')
!     - './data/tet4model3d.elem.msh'
!      open(31,file=filename,status='unknown')
!      do i=1,nelems
!        write(31,200) i,connbuf(5*(i-1)+1),
!     -  connbuf(5*(i-1)+2),connbuf(5*(i-1)+3),
!     -  connbuf(5*(i-1)+4),connbuf(5*(i-1)+5)
!      enddo
!  200 format(11i10)
!      close(31)
!      t1=omp_get_wtime()-t1
!      write(*,*)'write elem ASCII took',t1
      
      write(dataname,'(a10)') '/conn_tet4'
      call set_null_char(dataname,10)
      connbufsize=nfemelem*5
      call hdf_write_long_array
     - (fid,dataname,connbuf,connbufsize)
      write(dataname,'(a10)') '/conn_vox8'
      call set_null_char(dataname,10)
      connbufsize=nvoxelems*9
      call hdf_write_long_array
     - (fid,dataname,connbufvox,connbufsize)
      end
c_______________________________________________________________________
c_______________________________________________________________________
      subroutine voxelnode(nx,ny,nz,num,nk,nf,n) 
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
	integer nf(nx+1,ny+1,nz+1)

      do k=1,nz+1
      do j=1,ny+1
      do i=1,nx+1
      nf(i,j,k)=0
      enddo
      enddo
      enddo

      do k=1,2*nz-1,2
      do j=1,2*ny-1,2
      do i=1,2*nx-1,2
       if((num(i,j,k).ne.0).and.(num(i,j,k).lt.100))then
        nex=(i-1)/2+1
        ney=(j-1)/2+1
        nez=(k-1)/2+1
        do k1=0,1
        do j1=0,1
        do i1=0,1
         nf(nex+i1,ney+j1,nez+k1)=1
        enddo
        enddo
        enddo
       endif
      enddo
      enddo
      enddo

      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)
        do k=1,nz/nd     
        do j=1,ny/nd     
        do i=1,nx/nd
         nex=nd+2*nd*(i-1)
         ney=nd+2*nd*(j-1)
         nez=nd+2*nd*(k-1)      
         if(num(nex,ney,nez).ne.0)then
          do k1=0,1
          do j1=0,1
          do i1=0,1
       nf(nex/2+1-nd2+i1*nd,ney/2+1-nd2+j1*nd,nez/2+1-nd2+k1*nd)=1
          enddo
          enddo
          enddo
         endif       
        enddo
        enddo
        enddo
      enddo
c
      n=0
      do k=1,nz+1
      do j=1,ny+1
      do i=1,nx+1
       if(nf(i,j,k).gt.0)then
        n=n+1
        nf(i,j,k)=n
       endif
      enddo
      enddo
      enddo

      end
c_______________________________________________________________________
      subroutine voxelnode2(nx,ny,nz,num,nk,nf,n) 
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
	integer nf(nx+1,ny+1,nz+1)

      do k=1,nz+1
      do j=1,ny+1
      do i=1,nx+1
      nf(i,j,k)=0
      enddo
      enddo
      enddo

      do k=1,2*nz-1,2
      do j=1,2*ny-1,2
      do i=1,2*nx-1,2
       if((num(i,j,k).ne.0).and.(num(i,j,k).lt.100))then
        nex=(i-1)/2+1
        ney=(j-1)/2+1
        nez=(k-1)/2+1
        do k1=0,1
        do j1=0,1
        do i1=0,1
         nf(nex+i1,ney+j1,nez+k1)=1
        enddo
        enddo
        enddo
       endif
      enddo
      enddo
      enddo

      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)
        do k=1,nz/nd     
        do j=1,ny/nd     
        do i=1,nx/nd
         nex=nd+2*nd*(i-1)
         ney=nd+2*nd*(j-1)
         nez=nd+2*nd*(k-1)      
         if(num(nex,ney,nez).ne.0)then
          do k1=0,1
          do j1=0,1
          do i1=0,1

       nf(nex/2+1-nd2+i1*nd,ney/2+1-nd2+j1*nd,nez/2+1-nd2+k1*nd)=1
          enddo
          enddo
          enddo
         endif       
        enddo
        enddo
        enddo
      enddo
c
c_____add_____
      do ik=2,nk
       nd=2**(ik-1)
       nd2=2**(ik-2)
        do k=1,nz/nd
        do j=1,ny/nd
        do i=1,nx/nd
         nex=nd+2*nd*(i-1)
         ney=nd+2*nd*(j-1)
         nez=nd+2*nd*(k-1)
          if((num(nex,ney,nez).gt.0)
     -       .and.(num(nex,ney,nez).lt.100))then
           itet=0
            do kk=nez/2-nd2+1,nez/2+nd2+1,nd2
            do jj=ney/2-nd2+1,ney/2+nd2+1,nd2
            do ii=nex/2-nd2+1,nex/2+nd2+1,nd2
             if(nf(ii,jj,kk).gt.0)then
              itet=itet+1
             endif
            enddo
            enddo
            enddo

            if(itet.gt.8)then
	     ii=nex/2-nd2+1 + nd2
	     jj=ney/2-nd2+1 + nd2
	     kk=nez/2-nd2+1 + nd2
	     nf(ii,jj,kk)=1
            endif !itet.eq.8
           endif !
       enddo !i
       enddo !j
       enddo !k
      enddo !ik
c_____add_____
      n=0
      do k=1,nz+1
      do j=1,ny+1
      do i=1,nx+1
       if(nf(i,j,k).gt.0)then
        n=n+1
        nf(i,j,k)=n
       endif
      enddo
      enddo
      enddo

c      end
      end
c_______________________________________________________________________
      subroutine surfacemod(nx,ny,nz,ns,ori,ds,sur)
      real*8 ds,dif,sur(nx+1,ny+1,ns),ori(3)
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,1)-ori(3),ds)
      if(dif.le.0.25*ds)then
      sur(i,j,1)=sur(i,j,1)-dif
	endif
      enddo
      enddo
	if(ns.ge.2)then
	do is=2,ns
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,is)-ori(3),ds)
      if(dif.le.0.25*ds)then
      sur(i,j,is)=sur(i,j,is)-dif
	endif
	if(dif.ge.0.75*ds)then
      sur(i,j,is)=sur(i,j,is)-dif+ds
	endif
      if(sur(i,j,is-1)-sur(i,j,is).le.0.25*ds)then
      sur(i,j,is)=sur(i,j,is-1)
      endif
      enddo
      enddo
	enddo
	endif
      end
c_______________________________________________________________________
      subroutine 
     - octree(nx,ny,nz,nk,ns,num,ds,freq,c2,now,coeout,nout,noutdep)
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
      integer nl(8,3)
	integer isame(0:100+ns)
      real*8 ds,c2(ns),freq,gnow(0:ns,2),now
      real*8 coeout
c
      coeout=1
      gnow(0,1)=0
      do is=1,ns
      gnow(is,1)=c2(is)/freq/(now*1.0)/ds      
      enddo
      do is=0,ns
      gnow(is,2)=gnow(is,1)*coeout      
      enddo      
c
      do ik=2,nk
      nd=2**(ik-1)
      nd2=2**(ik-2)
       do k=1,nz/nd     
       do j=1,ny/nd     
       do i=1,nx/nd
        nex=nd+2*nd*(i-1)
        ney=nd+2*nd*(j-1)
        nez=nd+2*nd*(k-1)      
c
       iout=1
c
        do is=0,100+ns
        isame(is)=0
        enddo
        do k1=nez-nd2,nez+nd2,2*nd2    
        do j1=ney-nd2,ney+nd2,2*nd2    
        do i1=nex-nd2,nex+nd2,2*nd2
         isame(num(i1,j1,k1))=isame(num(i1,j1,k1))+1
        enddo
        enddo
        enddo
c
        do is=0,ns
        if((isame(is).eq.8).and.(2*nd2.le.gnow(is,iout)))then
c

        if(ik.ge.3)then
        call twotoone(nx,ny,nz,num,i,j,k,ik,key) 
         if(key.gt.0)then
          goto 102
         endif
         endif
c
         if(ik.eq.2)then
          k1=nez-3*nd2
           if(k1.ge.1)then 	 
	    do j1=ney-nd2,ney+nd2,2*nd2    
             do i1=nex-nd2,nex+nd2,2*nd2
              if(num(i1,j1,k1).gt.100)then
               goto 102
	      endif
             enddo
            enddo
           endif
          k1=nez+3*nd2
	   if(k1.le.2*nz-1)then 	 
	    do j1=ney-nd2,ney+nd2,2*nd2    
             do i1=nex-nd2,nex+nd2,2*nd2
              if(num(i1,j1,k1).gt.100)then
               goto 102
	      endif
             enddo
            enddo
	   endif
	 endif
c
         do k1=nez-nd2,nez+nd2,2*nd2    
         do j1=ney-nd2,ney+nd2,2*nd2    
         do i1=nex-nd2,nex+nd2,2*nd2
         num(i1,j1,k1)=0
         enddo
         enddo
         enddo
         num(nex,ney,nez)=is
         goto 102

        endif
        enddo
  102  damy=0
c
       enddo !i
       enddo !j
       enddo !k
      enddo  !ik

      end
c_______________________________________________________________________
      subroutine twotoone(nx,ny,nz,num,i,j,k,ik,key) 
      integer*2 num(2*nx-1,2*ny-1,2*nz-1)
      key=0
c
      nd=2**(ik-1)
      nd2=2**(ik-2)
      nd3=2**(ik-3)
      nex=nd+2*nd*(i-1)
      ney=nd+2*nd*(j-1)
      nez=nd+2*nd*(k-1)      
c
c i
      if(i.ne.1)then
      do k1=1,4
      do j1=1,4
      nex1=nex-nd3-2*nd2
      ney1=ney-nd3-nd2+(j1-1)*nd2
      nez1=nez-nd3-nd2+(k1-1)*nd2
      key=key+num(nex1,ney1,nez1)
      enddo
      enddo
	endif
      if(i.ne.nx/nd)then
      do k1=1,4
      do j1=1,4
      nex1=nex+nd3+2*nd2
      ney1=ney-nd3-nd2+(j1-1)*nd2
      nez1=nez-nd3-nd2+(k1-1)*nd2
      key=key+num(nex1,ney1,nez1)
      enddo
      enddo
	endif
c j
      if(j.ne.1)then
      do k1=1,4
      do i1=1,4
      nex1=nex-nd3-nd2+(i1-1)*nd2
      ney1=ney-nd3-2*nd2
      nez1=nez-nd3-nd2+(k1-1)*nd2
      key=key+num(nex1,ney1,nez1)
      enddo
      enddo
	endif
      if(j.ne.ny/nd)then
      do k1=1,4
      do i1=1,4
      nex1=nex-nd3-nd2+(i1-1)*nd2
      ney1=ney+nd3+2*nd2
      nez1=nez-nd3-nd2+(k1-1)*nd2
      key=key+num(nex1,ney1,nez1)
      enddo
      enddo
      endif
c k
      if(k.ne.1)then
      do j1=1,4
      do i1=1,4
      nex1=nex-nd3-nd2+(i1-1)*nd2
      ney1=ney-nd3-nd2+(j1-1)*nd2
      nez1=nez-nd3-2*nd2
      key=key+num(nex1,ney1,nez1)
      enddo
      enddo
	endif
      if(k.ne.nz/nd)then
      do j1=1,4
      do i1=1,4
      nex1=nex-nd3-nd2+(i1-1)*nd2
      ney1=ney-nd3-nd2+(j1-1)*nd2
      nez1=nez+nd3+2*nd2
      key=key+num(nex1,ney1,nez1)
      enddo
      enddo
      endif
      end
c_______________________________________________________________________
      subroutine surfacemod2(nx,ny,nz,ns,ori,ds,sur)
      real*8 ds,dif,sur(nx+1,ny+1,ns),ori(3)
	integer isur(nx+1,ny+1)

      do j=1,ny+1
      do i=1,nx+1
	isur(i,j)=0
      enddo
      enddo

c-- for 1st layer
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,1)-ori(3),ds)
      if(dif.le.0.25*ds)then
      sur(i,j,1)=sur(i,j,1)-dif
	endif
      enddo
      enddo
c-- for 2nd layer
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,2)-ori(3),ds)
      if(sur(i,j,2).lt.49999.)then
	if(dif.le.0.25*ds)then
	sur(i,j,2)=sur(i,j,2)-dif
	endif
	if(dif.ge.0.75*ds)then
	sur(i,j,2)=sur(i,j,2)-dif+ds*0.75
	isur(i,j)=1
	endif
      endif
      enddo
      enddo
c
      do j=1,ny+1
      do i=1,nx+1
      if(isur(i,j).eq.0)then
	dif=mod(sur(i,j,2)-ori(3),ds)
	if(dif.le.0.25*ds)then
	sur(i,j,2)=sur(i,j,2)-dif
	endif
	if(dif.ge.0.75*ds)then
	sur(i,j,2)=sur(i,j,2)-dif+ds
	endif
	if(sur(i,j,2-1)-sur(i,j,2).le.0.25*ds)then
	sur(i,j,2)=sur(i,j,2-1)
	endif
c      else
c	dif=mod(sur(i,j,2)-ori(3),ds)
c	if(dif.ge.0.7199*ds)then
c	sur(i,j,2)=sur(i,j,2)-dif+ds*0.75
c	endif
      endif
      enddo
      enddo
c
c
c-- for 3rd- layer
	if(ns.ge.2)then
	do is=3,ns
c	do is=2,ns
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,is)-ori(3),ds)
      if(dif.le.0.25*ds)then
      sur(i,j,is)=sur(i,j,is)-dif
	endif
	if(dif.ge.0.75*ds)then
      sur(i,j,is)=sur(i,j,is)-dif+ds
	endif
      if(sur(i,j,is-1)-sur(i,j,is).le.0.25*ds)then
      sur(i,j,is)=sur(i,j,is-1)
      endif
      enddo
      enddo
	enddo
	endif
      end
c_______________________________________________________________________
#ifdef NSQUASH
      subroutine surfacemod80(nx,ny,nz,ns,ori,ds,sur)
      real*8 ds,dif,sur(nx+1,ny+1,ns),ori(3)
      real*8 tmp
      tmp=1.0/NSQUASH

      do j=1,ny+1
      do i=1,nx+1
      dif=mod(sur(i,j,1)-ori(3),ds)
      if(dif.le.0.125*tmp*ds)then
      sur(i,j,1)=sur(i,j,1)-dif
      else if((dif.gt.0.125*tmp*ds).and.(dif.le.0.25*tmp*ds))then
      sur(i,j,1)=sur(i,j,1)-dif+0.25*tmp*ds
      endif
      enddo
      enddo
      if(ns.ge.2)then
      do is=2,ns
      do j=1,ny+1
      do i=1,nx+1
      dif=mod(sur(i,j,is)-ori(3),ds)
      if((dif.gt.0.).and.(dif.le.0.125*tmp*ds))then
      sur(i,j,is)=sur(i,j,is)-dif
      else if((dif.gt.0.125*tmp*ds).and.(dif.le.0.25*tmp*ds))then
      sur(i,j,is)=sur(i,j,is)-dif+0.25*tmp*ds
      endif

!      if((dif.gt.0.75*ds).and.(dif.le.0.875*ds))then
!      sur(i,j,is)=sur(i,j,is)-dif+0.75*ds
!      else if(dif.gt.0.875*ds)then
!      sur(i,j,is)=sur(i,j,is)-dif+ds
      if((dif.gt.(1.0-0.25*tmp)*ds).and.(dif.le.(1.0-0.125*tmp)*ds))then
      sur(i,j,is)=sur(i,j,is)-dif+(1.0-0.25*tmp)*ds
      else if(dif.gt.(1.0-0.125*tmp)*ds)then
      sur(i,j,is)=sur(i,j,is)-dif+ds

      endif
      if(is.eq.2)then
      if(sur(i,j,is-1)-sur(i,j,is).lt.0.25*tmp*ds)then
      sur(i,j,is)=sur(i,j,is-1)
      endif
      else
      if(sur(i,j,is-1)-sur(i,j,is).le.0.25*tmp*ds)then
      sur(i,j,is)=sur(i,j,is-1)
      endif
      endif
      enddo
      enddo
      enddo
      endif
      end
#else
c_______________________________________________________________________
      subroutine surfacemod80(nx,ny,nz,ns,ori,ds,sur)
      real*8 ds,dif,sur(nx+1,ny+1,ns),ori(3)
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,1)-ori(3),ds)
c      if(dif.le.0.25*ds)then
c      if((dif.gt.0.).and.(dif.le.0.125*ds))then
      if(dif.le.0.125*ds)then
      sur(i,j,1)=sur(i,j,1)-dif
      else if((dif.gt.0.125*ds).and.(dif.le.0.25*ds))then
      sur(i,j,1)=sur(i,j,1)-dif+0.25*ds
	endif
      enddo
      enddo
	if(ns.ge.2)then
	do is=2,ns
      do j=1,ny+1
      do i=1,nx+1
	dif=mod(sur(i,j,is)-ori(3),ds)
c      if(dif.le.0.25*ds)then
c      sur(i,j,is)=sur(i,j,is)-dif
      if((dif.gt.0.).and.(dif.le.0.125*ds))then
      sur(i,j,is)=sur(i,j,is)-dif
      else if((dif.gt.0.125*ds).and.(dif.le.0.25*ds))then
      sur(i,j,is)=sur(i,j,is)-dif+0.25*ds
	endif
c	if(dif.ge.0.75*ds)then
c      sur(i,j,is)=sur(i,j,is)-dif+ds
      if((dif.gt.0.75*ds).and.(dif.le.0.875*ds))then
      sur(i,j,is)=sur(i,j,is)-dif+0.75*ds
      else if(dif.gt.0.875*ds)then
      sur(i,j,is)=sur(i,j,is)-dif+ds
	endif
c      if(sur(i,j,is-1)-sur(i,j,is).le.0.25*ds)then
	if(is.eq.2)then
      if(sur(i,j,is-1)-sur(i,j,is).lt.0.25*ds)then
      sur(i,j,is)=sur(i,j,is-1)
      endif
	else
      if(sur(i,j,is-1)-sur(i,j,is).le.0.25*ds)then
      sur(i,j,is)=sur(i,j,is-1)
      endif
	endif
      enddo
      enddo
	enddo
	endif
      end
#endif
c_______________________________________________________________________
