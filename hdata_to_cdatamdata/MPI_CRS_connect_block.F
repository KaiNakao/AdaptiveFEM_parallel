c_______________________________________________________________________
      subroutine crs_connect_block_lib(n,ne,ib,fid)      
      implicit none
      integer*4 n,ne,ib
      integer*4 ntm,nab
      integer*8 fid,one,zero,intbufsize,realbufsize
      character dataname*50

      one=1
      zero=0
      
      write(dataname,'(a10)') '/crs_nodes'
      call set_null_char(dataname,10)
      call hdf_read_int_array_part(fid,dataname,zero,one,ntm)
c
      write(dataname,'(a11)') '/ABCsetting'
      call set_null_char(dataname,11)
      call hdf_read_int_array_part(fid,dataname,zero,one,nab)
c
      intbufsize=n*100 !n*100
      realbufsize=n*30 !n*30
      
      call main_CRS_connect_block
     - (n,ne,ntm,ib,nab,fid,intbufsize,realbufsize)
      
      end
c_______________________________________________________________________
!~      INCLUDE 'mpif.h'
!~      integer*4 myrank,ncpu,ierr
!~	integer im,mp
!~      real*8 dt
!~
!~      call MPI_INIT(ierr)
!~      call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
!~      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
!~	
!~      open(10,file='./data/para_setting.dat',status='unknown')
!~      read(10,*) !'num of MPI thread'
!~      read(10,*) mp
!~c      read(10,*) !'num of OpenMP thread'
!~c      read(10,*) np
!~      close(10)
!~      
!~      if(0.NE.MOD(mp,ncpu))then
!~      write(*,*) 'mp is not divisible by ncpu'
!~	call exit()
!~	endif
!~      do im=myrank*(mp/ncpu)+1,(myrank+1)*(mp/ncpu)
!~      write(*,*) myrank,im 
!~!      do im=1,mp 

!~c      open(10,file='./data/setting.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a12)')
!~     - './mdata/',im-1,'.setting.dat'
!~      open(10,file=filename,status='unknown')
!~      read(10,*) !'num of node' 
!~      read(10,*)  n
!~      read(10,*) !'num of elem' 
!~      read(10,*)  ne
!~      read(10,*) !'num of materials' 
!~      read(10,*)  kd
!~      close(10)
c
!      open(10,file='./data/time_setting.dat',status='unknown')
!      read(10,*) !'num of increment'
!      read(10,*) nt
!      read(10,*) !'num of increment'
!      read(10,*) dt
!      close(10)
c
!~c      open(10,file='./data/crs_nodes.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a14)')
!~     - './mdata/',im-1,'.crs_nodes.dat'
!~      open(10,file=filename,status='unknown')
!~      read(10,*) 
!~      read(10,*) ntm
!~      read(10,*) 
!~      read(10,*) mc
!~      close(10)    
c
!~c      open(10,file='./data/BC.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a7)')
!~     - './mdata/',im-1,'.BC.dat'
!~      open(10,file=filename,status='old')
!~      read(10,*) 
!~      read(10,*) ib
!~      close(10)
!~c
!~      write(filename,'(a8,i6.6,a15)')
!~     - './mdata/',im-1,'.ABCsetting.dat'
!~      open(10,file=filename,status='unknown')
!~c      open(10,file='./data/ABCsetting.dat',status='unknown')
!~      read(10,*) 
!~      read(10,*) nab
!~      close(10)
c
!      call main_crs_connect_block(n,ne,mc,ntm,ib,nab,im,fid)
      
!~      enddo
!~c
!~      call MPI_FINALIZE(ierr)
!      end
c_______________________________________________________________________
      subroutine main_crs_connect_block
     - (n,ne,ntm,ib,nab,fid,intbufsize,realbufsize)
      implicit none
      real*8 coor(n,3) !,rmat(kd,10),dt,fmin,fmax,alp,bet,alpha,beta
      integer*4 num(ne),cny(ne,10)
!      real*8 se(3*n),younglst(kd,2),ke(30,30),xx(4,3),me(10,10)
!      real*8 rho,young,rnyu,tmp1,tmp2
!	real*8 abc(nab,3,6,6)
      integer*4 nabc(nab,6)
      integer*4 ibc(ib),n,ne,ib,i,j,ie,ntm,nab
      integer*4 i1,i2,i3,i4,i5,itmp,j1,j2,i_val
      real*8 damy
      integer*8 intbufsize,realbufsize,bufcount,fid
      integer*4 intbuf(intbufsize)
!      real*8 realbuf(realbufsize)
!      character filename*50
      integer*4 numcrs(n,ntm),liscrs(n) !,numtmp(ntm),numtmpo(ntm)
!      integer*4 cr_ind(3*n*3*ntm),cr_ptr(3*n)
      integer*4 liscrssum(n)
      character dataname*50

!      call read_material(kd,rmat)
!      call read_geometry(n,ne,coor,cny,num,im)
      call read_tet10_geometry_hdf(n,ne,coor,cny,num,fid)
!      call makmatlist(kd,rmat,younglst)

!      call read_targetfreq(fmin,fmax)
!      call calc_alp_bet(fmin,fmax,alp,bet)
c
!      call read_BC(ib,ibc,fid)
      call MPI_ABC_read_BC(ib,ibc,fid)
c
c-------------------------------------------------------------------
!~c      open(10,file='./data/ABCnode.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a12)')
!~     - './mdata/',im-1,'.ABCnode.dat'
!~      open(10,file=filename,status='old')
!~      do i=1,nab
!~      do j=1,6
!~      read(10,*) nabc(i,j)
!~      enddo
!~      enddo
!~      close(10)
!~c      open(10,file='./data/ABCcoe.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a11)')
!~     - './mdata/',im-1,'.ABCcoe.dat'
!~      open(10,file=filename,status='old')
!~	do i=1,nab
!~	do ii=1,3
!~	do j2=1,6
!~	do j1=1,6
!~      read(10,*) abc(i,ii,j1,j2)
!~	enddo
!~	enddo
!~	enddo
!~	enddo
!~      close(10)
c-------------------------------------------------------------------
!~      numcrs=0
!~c      open(10,file='./data/crs_nodes.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a14)')
!~     - './mdata/',im-1,'.crs_nodes.dat'
!~      open(10,file=filename,status='old')
!~      do i=1,5
!~      read(10,*) 
!~      enddo
!~      do i=1,n
!~      read(10,*) liscrs(i)
!~      enddo
!~      read(10,*) 
!~      do i=1,n
!~      do j=1,liscrs(i)
!~      read(10,*) numcrs(i,j)
!~      enddo
!~      enddo
!~      close(10)
!~c
!~      liscrssum=0
!~      do i=1,n-1
!~      liscrssum(i+1)=liscrssum(i)+liscrs(i)
!~      enddo
!~      
!~c      open(31,file='./data/crs_connect.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a16)')
!~     - './mdata/',im-1,'.crs_connect.dat'
!~      open(31,file=filename,status='unknown')
!~      write(31,*) 'num of comp num in crs'
!~      write(31,*) 10*10*ne
!~      write(31,*) 'comp num in crs'
!~c      open(41,file='./data/crs_ABC_connect.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a20)')
!~     - './mdata/',im-1,'.crs_ABC_connect.dat'
!~      open(41,file=filename,status='unknown')
!~      write(41,*) 'num of comp num in crs'
!~      write(41,*) 6*6*nab
!~      write(41,*) 'comp num in crs'


      write(dataname,'(a8)') '/ABCnode'
      call set_null_char(dataname,8)
      call hdf_read_int_array(fid,dataname,
     - intbuf,intbufsize,bufcount)
      do i=1,nab
      do j=1,6
      nabc(i,j)=intbuf(6*(i-1)+j)
      enddo
      enddo
      
!      write(dataname,'(a7)') '/ABCcoe'
!      call set_null_char(dataname,7)
!      call hdf_read_double_array(fid,dataname,
!     - realbuf,realbufsize,bufcount)
!      bufcount=0
!      do i=1,nab
!      do ii=1,3
!      do j2=1,6
!      do j1=1,6
!      bufcount=bufcount+1
!      abc(i,ii,j1,j2)=realbuf(bufcount)
!      enddo
!      enddo
!      enddo
!      enddo

      numcrs=0
      write(dataname,'(a10)') '/crs_nodes'
      call set_null_char(dataname,10)
      call hdf_read_int_array(fid,dataname,
     - intbuf,intbufsize,bufcount)
      do i=1,n
      liscrs(i)=intbuf(i+2)
      enddo
      bufcount=2+n
      do i=1,n
      do j=1,liscrs(i)
      bufcount=bufcount+1
      numcrs(i,j)=intbuf(bufcount)
      enddo
      enddo
c
      liscrssum=0
      do i=1,n-1
      liscrssum(i+1)=liscrssum(i)+liscrs(i)
      enddo
      bufcount=1
      intbuf(1)=10*10*ne


c
!      crsval=0.
!      se=0.
c----------------------------------------------
      do ie=1,ne
c
!      in=num(ie)
!      rho=rmat(in,3)
!      young=younglst(in,1) 
!	rnyu=younglst(in,2)
!      alpha=rmat(in,4)*alp 
!      beta=rmat(in,4)*bet 
!	tmp1=4/dt/dt+2*alpha/dt
!	tmp2=1+2*beta/dt
	
!      do j1=1,3
!      do i1=1,4
!      xx(i1,j1)=coor(cny(ie,i1),j1)
!      enddo
!      enddo
c      call tetra102(young,rnyu,xx,ke)
c      call tetra10me2(rho,xx,me)
c
      do i1=1,10
      i2=cny(ie,i1)
      do i3=1,10
      i4=cny(ie,i3)
c
      do i5=1,liscrs(i2)
      if(numcrs(i2,i5).eq.i4)then
      itmp=i5
      goto 252
      endif      
      enddo    
  252 damy=0
c
!      write(31,*) liscrssum(i2)+itmp
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=liscrssum(i2)+itmp
c
c      do ii=1,3
c      do jj=1,3
c      crsval(3*(i2-1)+ii,3*(itmp-1)+jj)=
c     -  crsval(3*(i2-1)+ii,3*(itmp-1)+jj)  
c     - +ke(3*(i1-1)+ii,3*(i3-1)+jj)*tmp2
c      enddo
c      enddo
c      do ii=1,3
c      crsval(3*(i2-1)+ii,3*(itmp-1)+ii)=
c     -  crsval(3*(i2-1)+ii,3*(itmp-1)+ii)  
c     - +me(i1,i3)*tmp1
c      enddo
c
      enddo
      enddo
c      
      enddo
      write(dataname,'(a12)') '/crs_connect'
      call set_null_char(dataname,12)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      bufcount=1
      intbuf(1)=6*6*nab
      
c----------------------------------------------
      do i=1,nab
      do j2=1,6
      i2=nabc(i,j2)
      do j1=1,6
      i4=nabc(i,j1)
      do i5=1,liscrs(i2)
      if(numcrs(i2,i5).eq.i4)then
      itmp=i5
      goto 251
      endif      
      enddo    
  251 damy=0     
!      write(41,*) liscrssum(i2)+itmp
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=liscrssum(i2)+itmp
c      do ii=1,3
c      crsval(3*(i2-1)+ii,3*(itmp-1)+ii)=
c     - crsval(3*(i2-1)+ii,3*(itmp-1)+ii)+abc(i,ii,j2,j1)*2./dt
c	enddo
      enddo
      enddo
      enddo

!      close(31)
!      close(41)

!~c------------------------------------------------------------
!~c      open(21,file='./data/crs_ptr.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a12)')
!~     - './mdata/',im-1,'.crs_ptr.dat'
!~      open(21,file=filename,status='unknown')
!~      write(21,*) 'num of components'
!~      write(21,*) n+1
!~      write(21,*) 'components'
!~      i_val=0
!~      write(21,*) 1
!~      do i=1,n
!~      i_val=i_val+liscrs(i)
!~      write(21,*) i_val+1
!~      enddo
!~      close(21)   
!~
!~c      open(21,file='./data/ori_crs_val_block.dat',status='unknown')
!~c      write(21,*) 'num of components'
!~c      write(21,*) 9*i_val
!~c      write(21,*) 'components'
!~c      do i1=1,n
!~c      do i2=1,liscrs(i1)
!~c      do i4=1,3
!~c      do i3=1,3
!~c      write(21,*) crsval(3*(i1-1)+i3,3*(i2-1)+i4)
!~c      enddo
!~c      enddo
!~c      enddo
!~c      enddo
!~c      close(21)    
!~c   
!~c      open(21,file='./data/crs_ind.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a12)')
!~     - './mdata/',im-1,'.crs_ind.dat'
!~      open(21,file=filename,status='unknown')
!~      write(21,*) 'num of components'
!~      write(21,*) i_val
!~      write(21,*) 'components'
!~      do i=1,n
!~      do j=1,liscrs(i)
!~      write(21,*) numcrs(i,j)
!~      enddo
!~      enddo
!~      close(21)
!~c
!~c      open(21,file='./data/crs_block_list.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a19)')
!~     - './mdata/',im-1,'.crs_block_list.dat'
!~      open(21,file=filename,status='unknown')
!~      do i=1,n
!~      do j=1,liscrs(i)
!~      if(i.eq.numcrs(i,j))then
!~      write(21,*) liscrssum(i)+j      
!~      endif 
!~      enddo
!~      enddo     
!~      close(21)
!~c
      write(dataname,'(a16)') '/crs_ABC_connect'
      call set_null_char(dataname,16)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      intbuf(1)=n+1
      i_val=0
      intbuf(2)=1
      do i=1,n
      i_val=i_val+liscrs(i)
      intbuf(2+i)=i_val+1
      enddo
      
      bufcount=2+n-1
      call check_increment(bufcount,intbufsize)
      write(dataname,'(a8)') '/crs_ptr'
      call set_null_char(dataname,8)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)
   
      intbuf(1)=i_val
      bufcount=1
      do i=1,n
      do j=1,liscrs(i)
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=numcrs(i,j)
      enddo
      enddo
      write(dataname,'(a8)') '/crs_ind'
      call set_null_char(dataname,8)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      bufcount=0
      do i=1,n
      do j=1,liscrs(i)
      if(i.eq.numcrs(i,j))then
      call check_increment(bufcount,intbufsize)
      intbuf(bufcount)=liscrssum(i)+j
      endif 
      enddo
      enddo     
      write(dataname,'(a15)') '/crs_block_list'
      call set_null_char(dataname,15)
      call hdf_write_int_array(fid,dataname,intbuf,bufcount)

      end	   
!~c_______________________________________________________________________
!~      subroutine sort(list,list2,n1,n2)
!~      integer list(n1),list2(n1)
!~     
!~      do j=0,n2-2
!~      do i=2,n2-j
!~      i2=list(i)
!~      i1=list(i-1)      
!~      if(i1.gt.i2)then
!~      list(i)=i1
!~      list(i-1)=i2
!~      k1=list2(i)
!~      k2=list2(i-1)     
!~      list2(i)=k2
!~      list2(i-1)=k1
!~      endif
!~      enddo
!~      enddo
!~
!~      end
!~c_______________________________________________________________________
!~      subroutine makmatlist(kd,rmat,younglst)
!~      real*8 rmat(kd,10)
!~	real*8 younglst(kd,2)
!~	real*8 c1tmp,c2tmp,rhtmp,ygtmp,rntmp
!~      do it=1,kd
!~      c1tmp=rmat(it,1)
!~	c2tmp=rmat(it,2)
!~	rhtmp=rmat(it,3)
!~      call mak_E_nyu(c1tmp,c2tmp,rhtmp,ygtmp,rntmp)
!~      younglst(it,1)=ygtmp
!~      younglst(it,2)=rntmp
!~      enddo
!~      end
!~c_______________________________________________________________________
!~      subroutine mak_E_nyu(c1,c2,rho,E,nyu)
!~      real*8 c1,c2,rho
!~	real*8 E,nyu
!~      E=((3*c1**2*c2**2-4*c2**4)*rho)/(c1**2-c2**2)
!~	nyu=(c1**2-2*c2**2)/(2*c1**2-2*c2**2)
!~c      E=c1**2*rho-4/3.*c2**2*rho
!~c	nyu=c2**2*rho
!~	end
!~c_______________________________________________________________________
!~      subroutine read_material(kd,rmat)
!~	real*8 rmat(kd,10)
!~      open(42,file='./data/material.dat',status='unknown')
!~      do j=1,kd
!~	read(42,*)
!~      do i=1,10
!~      read(42,*) rmat(j,i)
!~	enddo
!~	enddo
!~	close(42)
!~      end    
!~c_______________________________________________________________________
!~      subroutine read_geometry(n,ne,coor,cny,num,im)
!~      character filename*50
!~	real*8 coor(n,3)
!~	integer cny(ne,10),num(ne)
!~
!~c      open(42,file='./data/geometry.flavia.msh',status='unknown')
!~      write(filename,'(a8,i6.6,a20)')
!~     - './mdata/',im-1,'.geometry.flavia.msh'
!~      open(42,file=filename,status='unknown')
!~      read(42,*) 
!~      read(42,*) 
!~      read(42,*) 
!~      do i=1,n
!~      read(42,'(i10,3f16.6)') 
!~     -         itmp,coor(itmp,1),coor(itmp,2),coor(itmp,3)
!~      enddo
!~      read(42,*) 
!~      read(42,*) 
!~      read(42,*) 
!~      read(42,*) 
!~      do i=1,ne
!~      read(42,200) itmp,cny(i,1),cny(i,2),cny(i,3),cny(i,4),
!~     -                  cny(i,5),cny(i,6),cny(i,7),cny(i,8),
!~     -                  cny(i,9),cny(i,10),
!~     -                  num(i)  
!~      enddo
!~      close(42)
!~  200 format(12i10)
!~      end
!~c_______________________________________________________________________
!~      subroutine read_BC(ib,ibc,im)
!~      character filename*50
!~	integer ibc(ib)
!~c      open(42,file='./data/BC.dat',status='unknown')
!~      write(filename,'(a8,i6.6,a7)')
!~     - './mdata/',im-1,'.BC.dat'
!~      open(42,file=filename,status='unknown')
!~     	read(42,*)
!~	read(42,*)
!~	read(42,*)
!~      do i=1,ib
!~      read(42,*) ibc(i)
!~	enddo
!~	close(42)
!~      end
!~c_______________________________________________________________________
!~      subroutine calc_alp_bet(fmin,fmax,alp,bet)
!~      real*8 alp,bet,fmin,fmax,pi
!~      pi=atan(1.0)*4.0
!~c C= alp M + beta K
!~c Integrate[(1/2/Q - (alp/2/w + bet w/2))^2, {w, wmin, wmax}]
!~c wmin=2 Pi fmin, wmax=2 Pi fmax
!~c alp & bet is computed for h=1.0
!~      alp=(-2*fmax*fmin*Pi*(3*(fmax**2 - fmin**2) + 
!~     -      2*(fmax**2 + fmax*fmin + fmin**2)*Log(fmin/fmax)))/
!~     -  (fmax - fmin)**3
!~	bet=(3*(fmax**2 - fmin**2 + 2*fmax*fmin*Log(fmin/fmax)))/
!~     -  (2.*(fmax - fmin)**3*Pi) 
!~c
!~      end
!~c_______________________________________________________________________
!~      subroutine read_targetfreq(fmin,fmax)
!~      real*8 fmin,fmax
!~      open(10,file='./data/freq_setting.dat',status='unknown')
!~	read(10,*)
!~	read(10,*) fmin
!~	read(10,*)
!~	read(10,*) fmax
!~      close(10)
!~      end
