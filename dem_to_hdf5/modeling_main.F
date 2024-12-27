      implicit none
      real*8 ds,freq,coeout,now
      integer nxs,nys,ns,nk,nout,noutdep

      open(31,file='./data/modeling_setting.dat',status='old')
      read(31,*) !'nx, ny'
      read(31,*) nxs, nys
      read(31,*) !'ds'
      read(31,*) ds
      read(31,*) !'num of layer'
      read(31,*) ns
      read(31,*) !'nk 2**nk: nx and ny = 2**nk * alpha'
      read(31,*) nk
      read(31,*) !'freq'
      read(31,*) freq
      read(31,*) !'now'
      read(31,*) now
      close(31) 
      write(*,*) 'use squash factor',NSQUASH
      if(NSQUASH.ne.1)then
      write(*,*) 'equivalent ds',ds*NSQUASH
      write(*,*) 'equivalent now',now/NSQUASH
      endif
      call mainmodeling
     - (nxs,nys,ds,ns,nk,freq,now,coeout,nout,noutdep)

      end
c_______________________________________________________________________
      subroutine mainmodeling
     - (nxs,nys,ds,ns,nk,freq,now,coeout,nout,noutdep)
c modeling
c npx*npy=num of cpu
      parameter(npx=1,npy=1)	!no touch
c domain for one cpu: how many unit is incuded
      parameter(nxm=1,nym=1)
      real*8 rmat(ns,10)  
      real*8 ds,ori(3),surf(nxs+1,nys+1,ns),surmax
      real*8 freq,xt,yt,pl(5),tmp,coeout
      real*8 minx(3),now,rnsquash
      character*64 filename

      ori(1)=0.
      ori(2)=0.
      ori(3)=0.

      open(61,file='./data/modeldomain.dat',status='unknown')
      write(61,*) 'xmin-plane'
      write(61,*) ori(1)
      write(61,*) 'xmax-plane'
      write(61,*) ori(1)+ds*nxs
      write(61,*) 'ymin-plane'
      write(61,*) ori(2)
      write(61,*) 'ymin-plane'
      write(61,*) ori(2)+ds*nys
      write(61,*) 'zmin-plane'
      write(61,*) ori(3)
      close(61)

      rnsquash=1.0/NSQUASH
      do is=1,ns
      write(filename,'(a10,i4.4,a4)')
     - './data/sur',is,'.dat'
      open(31,file=filename,status='old')
      read(31,*) !'nx, ny'
      read(31,*) !nx, ny
      read(31,*) !'DEM data'
      do j=1,nys+1
      do i=1,nxs+1
      read(31,*) surf(i,j,is)
      surf(i,j,is)=surf(i,j,is)*rnsquash
      enddo
      enddo
      close(31)
      enddo


c_______________________________________________________________________
c_______________________________________________________________________
cc
c___________________________________
cc
c read material data
      call readmaterial(ns,rmat)
c read material data

c level set for octree
cc
c domain decomposition
      myrank=0
      ix=mod(myrank,npx)+1
      iy=(myrank-ix+1)/npx+1   
      ijm=1
      im=mod(ijm-1,nxm)+1
      jm=int((ijm-im)/nxm)+1

c modeling
      call searchsurfmax(nxs,nys,ns,surf,surmax)
      call calcnz(surmax,ori,ds,nk,nzs)
      call modeling
     -   (nxs,nys,nzs,ns,surf,ori,ds,freq,now,rmat,ix,iy,im,jm,nk,
     -  coeout,nout,noutdep)
c modeling

      stop
      end
c_______________________________________________________________________
      subroutine calcnz(surmax,ori,ds,nk,nz)
      real*8 surmax,ori(3),ds
      nz=int((surmax-ori(3))/ds)
      if(surmax-ori(3)-ds*nz.gt.0.1*ds)then
      nz=nz+1
      endif     
      nz=(int(nz/(2**nk))+1)*2**nk
      end
c_______________________________________________________________________
      subroutine searchsurfmax(nx,ny,ns,surf,surmax)
      real*8 surmax,surf(nx+1,ny+1,ns)
      surmax=surf(1,1,1)
      do j=1,ny+1
      do i=1,nx+1
      if(surf(i,j,1).gt.surmax)then
      surmax=surf(i,j,1)
      endif
      enddo
      enddo
      end
c_______________________________________________________________________
      subroutine readmaterial(ns,rmat)
      real*8 rmat(ns,10)  
      open(51,file='./data/material.dat',status='old')
      do is=1,ns
      read(51,*)
      do i=1,10
      read(51,*) rmat(is,i)
      enddo
      enddo
      close(51)
      end
