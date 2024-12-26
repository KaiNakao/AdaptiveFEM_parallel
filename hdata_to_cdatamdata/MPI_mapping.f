c_______________________________________________________________________
! map tet4 to tet10 mesh
! requires that elements in tet4 mesh is in same order of tet10 mesh,
! also, each second order node in tet10 mesh should be in the order 
! designated by cnyele
! time for mapping is O(N) where N is model size
      subroutine mapping_pairs_lib(n,ne,n4,ne4,fidm,fidc)      
      implicit none
      integer*4 n,ne,n4,ne4,ie,i,iptr,ii
      integer*4 num(ne),cny(ne,10),num4(ne4),cny4(ne4,4)
      integer*4 pairs(n,4),pairsize(n),pair1,pair2
      integer*4 cnyele(5:10,2)
      real*8 coor(n,3),coor4(n4,3),tmp,one,half
      integer*4 ptrbuf(4*n),indbuf(4*n)
      real*8 coebuf(4*n)
      integer*8 fidm,fidc,ptrcount,indcount,coecount,bufsize,isize
      character dataname*50
c
      bufsize=4*n
      cnyele(5,1)=1
      cnyele(5,2)=2
      cnyele(6,1)=2
      cnyele(6,2)=3
      cnyele(7,1)=1
      cnyele(7,2)=3
      cnyele(8,1)=1
      cnyele(8,2)=4
      cnyele(9,1)=2
      cnyele(9,2)=4
      cnyele(10,1)=4
      cnyele(10,2)=3
      one=1.0
      half=0.5
      call read_tet10_geometry_hdf(n,ne,coor,cny,num,fidm)
      call read_tet4_geometry_hdf(n4,ne4,coor4,cny4,num4,fidc)

      do i=1,n
      pairsize(i)=-1
      pairs(i,1)=-2
      pairs(i,2)=-2
      pairs(i,3)=-2
      pairs(i,4)=-2
      enddo
c
      do ie=1,ne

      do i=1,4

! error check
      if(pairsize(cny(ie,i)).eq.1)then
        if( pairs(cny(ie,i),1).ne.cny4(ie,i) ) then
          write(*,*)'error:2',ie,i
          stop
        endif
        if( pairs(cny(ie,i),2).ne.-1 ) then
          write(*,*)'error:3',ie,i
          stop
        endif
      endif
      if(pairsize(cny(ie,i)).eq.2)then
        write(*,*)'error:4',ie,i
        stop
      endif
! end error check

      pairsize(cny(ie,i))=1
      pairs(cny(ie,i),1)=cny4(ie,i)
      pairs(cny(ie,i),2)=-1
      enddo

      do i=5,10
      pair1=cny4(ie,cnyele(i,1))
      pair2=cny4(ie,cnyele(i,2))
      if(pair1.gt.pair2)then
      pair1=cny4(ie,cnyele(i,2))
      pair2=cny4(ie,cnyele(i,1))
      endif

! error check
!      if(pairsize(cny(ie,i)).eq.2)then
!        if( pairs(cny(ie,i),1).ne.pair1 ) then
!          write(*,*) pairs(cny(ie,i),1),cny4(ie,cnyele(i,1))
!          write(*,*) pairs(cny(ie,i),2),cny4(ie,cnyele(i,2))
!          write(*,*)'error:5',ie,i
!          stop
!        endif
!        if( pairs(cny(ie,i),2).ne.pair2 ) then
!          write(*,*)'error:6',ie,i
!          write(*,*) pairs(cny(ie,i),1),cny4(ie,cnyele(i,1))
!          write(*,*) pairs(cny(ie,i),2),cny4(ie,cnyele(i,2))
!          stop
!        endif
!      endif
!      if( pairsize(cny(ie,i)).eq.1 ) then
!        write(*,*)'error:7',ie,i
!        stop
!      endif
! end error check
      if(pairsize(cny(ie,i)).eq.2)then
      if( (pairs(cny(ie,i),1).ne.pair1).or.
     -    (pairs(cny(ie,i),2).ne.pair2) )then
      pairsize(cny(ie,i))=4
      pairs(cny(ie,i),3)=pair1
      pairs(cny(ie,i),4)=pair2
      endif
      endif
      if(pairsize(cny(ie,i)).eq.-1)then
      pairsize(cny(ie,i))=2
      pairs(cny(ie,i),1)=pair1
      pairs(cny(ie,i),2)=pair2
      endif

      enddo ! i
      enddo ! ie

! error check
      do i=1,n
      if( pairs(i,1).lt.1 ) then
        write(*,*)'error:01',i,pairsize(i),pairs(i,1),pairs(i,2)
        stop
      endif
      if( (pairsize(i).eq.1).and.(pairs(i,2).ne.-1))then
        write(*,*)'error:02',i,pairsize(i),pairs(i,1),pairs(i,2)
        stop
      endif
      if( (pairsize(i).eq.2).and.(pairs(i,2).lt.1))then
        write(*,*)'error:03',i,pairsize(i),pairs(i,1),pairs(i,2)
        stop
      endif
      enddo
      do i=1,n
      do ii=1,3
        if(pairsize(i).eq.1)then
          tmp=coor4(pairs(i,1),ii)
        elseif(pairsize(i).eq.2)then
          tmp=(coor4(pairs(i,1),ii)+coor4(pairs(i,2),ii))*0.5
        elseif(pairsize(i).eq.4)then
          tmp=(
     -    coor4(pairs(i,1),ii)+coor4(pairs(i,2),ii)+
     -    coor4(pairs(i,3),ii)+coor4(pairs(i,4),ii))*0.25
        else
          write(*,*)'error:1',i,ii,pairsize(i),pairs(i,1),pairs(i,2)
          stop
        endif
        if( abs(coor(i,ii)-tmp).gt.0.001 )then
          write(*,*)'error:2',i,ii,pairsize(i),pairs(i,1),pairs(i,2)
          stop
        endif
      enddo
      enddo
! end error check
      iptr=1

      ptrbuf(1)=n+1
      ptrbuf(2)=iptr
      ptrcount=2
      
      do i=1,n
      iptr=iptr+pairsize(i)
      enddo
      
      indbuf(1)=iptr-1
      indcount=1
      coecount=0

      iptr=1
      do i=1,n
      if(pairsize(i).eq.1)then
      iptr=iptr+1
      call check_increment(ptrcount,bufsize)
      ptrbuf(ptrcount)=iptr
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,1)
      call check_increment(coecount,bufsize)
      coebuf(coecount)=one
      endif
      if(pairsize(i).eq.2)then
      iptr=iptr+2
      call check_increment(ptrcount,bufsize)
      ptrbuf(ptrcount)=iptr
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,1)
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,2)
      call check_increment(coecount,bufsize)
      coebuf(coecount)=half
      call check_increment(coecount,bufsize)
      coebuf(coecount)=half
      endif
      if(pairsize(i).eq.4)then
      iptr=iptr+4
      call check_increment(ptrcount,bufsize)
      ptrbuf(ptrcount)=iptr
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,1)
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,2)
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,3)
      call check_increment(indcount,bufsize)
      indbuf(indcount)=pairs(i,4)
      call check_increment(coecount,bufsize)
      coebuf(coecount)=0.25
      call check_increment(coecount,bufsize)
      coebuf(coecount)=0.25
      call check_increment(coecount,bufsize)
      coebuf(coecount)=0.25
      call check_increment(coecount,bufsize)
      coebuf(coecount)=0.25
      endif
      if(pairsize(i).eq.-1)then
      write(*,*) 'something wrong here',i,pairsize(i)
      stop
      endif
      enddo

      write(dataname,'(a8)') '/map_ptr'
      isize=8
      call set_null_char(dataname,isize)
      call hdf_write_int_array(fidc,dataname,ptrbuf,ptrcount)
      write(dataname,'(a8)') '/map_ind'
      isize=8
      call set_null_char(dataname,isize)
      call hdf_write_int_array(fidc,dataname,indbuf,indcount)
      write(dataname,'(a8)') '/map_coe'
      isize=8
      call set_null_char(dataname,isize)
      call hdf_write_double_array(fidc,dataname,coebuf,coecount)
      
! error check
!      write(filename,'(a8,i6.6,a12)')
!     - './ccata/',im-1,'.map_ind.dat'
!      open(13,file=filename,status='unknown')
!      read(13,*)
!      read(13,*)
!      read(13,*)
!
!      do i=1,n
!      if(pairsize(i).eq.1)then
!        read(13,*) pair1
!        if(pair1.ne.pairs(i,1))then
!          write(*,*)'error:91',i,pairs(i,1),pair1
!          stop
!        endif
!      elseif(pairsize(i).eq.2)then
!        read(13,*) pair1
!        read(13,*) pair2
!        if(pair1.gt.pair2)then
!          iptr=pair1
!          pair1=pair2
!          pair2=iptr
!        endif
!        if(pair1.ne.pairs(i,1))then
!          write(*,*)'error:92',i,pairs(i,1),pairs(i,2),pair1,pair2
!        endif
!        if(pair2.ne.pairs(i,2))then
!          write(*,*)'error:93',i,pairs(i,1),pairs(i,2),pair2,pair2
!        endif
!      endif
!      enddo
!      close(13)
!
      end
c_______________________________________________________________________
!! map tet4 to tet10 mesh
!! tet4 mesh is mapped to tet10 mesh by based on position of nodes
!! element order can be non-matching between tet4 and tet10 mesh
!! time for mapping is O(N^2) where N is model size
!!      subroutine main_midpoint(n,ne,n4,ne4,kd,im,fidm,fidc)
!      subroutine mapping_midpoint_lib(n,ne,n4,ne4,kd,im,fidm,fidc)      
!      implicit none
!      real*8 coor(n,3),coor4(n4,3)
!      integer num(ne),cny(ne,10),num4(ne4),cny4(ne4,4)
!      character filename*50
!      real*8 vol1,vol2,vol3,vol4,vols
!      real*8 coe(5,3),co(3),cl(3,3),vole,cc(4)
!      real*8 coelist(4*n)
!      integer indlist(4*n)
!      integer combination(6,2),k1,k2
!      integer im,iptr,id,i,ii,kk,k,ne,n,n4,ne4,kd
!      character dataname*50
!      integer fidm,fidc,bufcount,bufsize,intbuf(2*n)
!      real*8 realbuf(2*n)
!c
!      bufsize=2*n
!c
!      combination(1,1)=1
!      combination(1,2)=2
!      combination(2,1)=1
!      combination(2,2)=3
!      combination(3,1)=1
!      combination(3,2)=4
!      combination(4,1)=2
!      combination(4,2)=3
!      combination(5,1)=2
!      combination(5,2)=4
!      combination(6,1)=3
!      combination(6,2)=4
!
!      call read_tet10_geometry_hdf(n,ne,coor,cny,num,fidm)
!      call read_tet4_geometry_hdf(n4,ne4,coor4,cny4,num4,fidc)
!ccc
!
!      coelist=0
!      indlist=0
!      iptr=1
!      write(dataname,'(a8)') '/map_ptr'
!      call set_null_char(dataname,8)
!      intbuf(1)=n+1
!      intbuf(2)=iptr
!      bufcount=2
!c
!	do id=1,n
!	
!      if((mod(id,1000).eq.0).and.(im.eq.0))then
!      write(*,*) 'processing',id,'out of',n
!      endif
!	  
!	!check if tet10 point exists in tet4 point
!      do i=1,n4
!      if( (abs(coor(id,1)-coor4(i,1)).lt.0.0001).and.
!     -    (abs(coor(id,2)-coor4(i,2)).lt.0.0001).and.
!     -    (abs(coor(id,3)-coor4(i,3)).lt.0.0001) )then
!	coelist(iptr)=1.0
!      indlist(iptr)=i
!      iptr=iptr+1
!	goto 132
!      endif
!      enddo !i
!	 
!	!check if tet10 point is on edge of tet4
!      do k=1,ne4
!      do kk=1,4
!      do ii=1,3
!      coe(kk,ii)=coor4(cny4(k,kk),ii)
!      enddo !ii
!      enddo !kk
!      do kk=1,6
!      k1=combination(kk,1)
!      k2=combination(kk,2)
!      if( (abs(coor(id,1)-(coe(k1,1)+coe(k2,1))*0.5).lt.0.0001).and.
!     -    (abs(coor(id,2)-(coe(k1,2)+coe(k2,2))*0.5).lt.0.0001).and.
!     -    (abs(coor(id,3)-(coe(k1,3)+coe(k2,3))*0.5).lt.0.0001) )then
!      coelist(iptr)=0.5
!      indlist(iptr)=cny4(k,k1)
!      iptr=iptr+1
!      coelist(iptr)=0.5
!      indlist(iptr)=cny4(k,k2)
!      iptr=iptr+1
!	goto 132
!      endif
!      enddo !kk
!      enddo !k
!	
!      write(6,*) 'something wrong?'
!      stop
!132   i=0
!      call check_increment(bufcount,bufsize)
!      intbuf(bufcount)=iptr
!      enddo ! id
!      call hdf_write_long_array(fidc,dataname,intbuf,bufcount)
!c
!      if(iptr.gt.bufsize)then
!      write(*,*) 'enlarge intbufsize',iptr,bufsize
!      endif
!      intbuf(1)=iptr-1
!      do i=1,iptr-1
!      intbuf(i+1)=indlist(i)
!      enddo
!      write(dataname,'(a8)') '/map_ind'
!      call set_null_char(dataname,8)
!      bufcount=iptr
!      call hdf_write_long_array(fidc,dataname,intbuf,bufcount)
!      write(dataname,'(a8)') '/map_coe'
!      call set_null_char(dataname,8)
!      bufcount=iptr-1
!      call hdf_write_double_array(fidc,dataname,coelist,bufcount)
!      
!      end	
!c_______________________________________________________________________
!      subroutine main(n,ne,n4,ne4,kd,im)
!      real*8 coor(n,3),coor4(n4,3)
!      integer num(ne),cny(ne,10),num4(ne4),cny4(ne4,4)
!      character filename*50
!      real*8 vol1,vol2,vol3,vol4,vols
!      real*8 coe(5,3),co(3),cl(3,3),vole,cc(4)
!      real*8 coelist(4*n)
!      integer indlist(4*n)
!      integer count,ktmp
!c
!      call read_geometry(n,ne,coor,cny,num,im)
!c
!      call read_geometry4(n4,ne4,coor4,cny4,num4,im)
!ccc
!
!      iptr=1
!      write(filename,'(a8,i6.6,a12)')
!     - './cdata/',im-1,'.map_ptr.dat'
!      open(10,file=filename,status='unknown')
!      write(filename,'(a8,i6.6,a12)')
!     - './cdata/',im-1,'.map_ind.dat'
!      open(11,file=filename,status='unknown')
!      write(filename,'(a8,i6.6,a12)')
!     - './cdata/',im-1,'.map_coe.dat'
!      open(12,file=filename,status='unknown')
!
!      coelist=0
!      indlist=0
!
!      write(10,*) 'num of components'
!      write(10,*) n+1
!      write(10,*) 'components'
!      write(10,*) iptr
!c
!      do id=1,n
!c  
!      do ii=1,3
!      co(ii)=coor(id,ii)
!      enddo ! ii
!c
!      do k=1,ne4
!      do kk=1,4
!      do ii=1,3
!      coe(kk,ii)=coor4(cny4(k,kk),ii)
!      enddo !ii
!      enddo !kk
!      do ii=1,3
!      coe(5,ii)=co(ii)
!      enddo !ii 
!c----------------------
!      do ii=1,3
!      cl(1,ii)=coe(2,ii)-coe(1,ii)
!      cl(2,ii)=coe(3,ii)-coe(1,ii)
!      cl(3,ii)=coe(4,ii)-coe(1,ii)
!      enddo !ii
!      call tetvol(cl,vole)
!c----------------------
!      do ii=1,3
!      cl(1,ii)=coe(2,ii)-coe(1,ii)
!      cl(2,ii)=coe(3,ii)-coe(1,ii)
!      cl(3,ii)=coe(5,ii)-coe(1,ii)
!      enddo !ii
!      call tetvol(cl,vol3)
!      do ii=1,3
!      cl(1,ii)=coe(4,ii)-coe(1,ii)
!      cl(2,ii)=coe(2,ii)-coe(1,ii)
!      cl(3,ii)=coe(5,ii)-coe(1,ii)
!      enddo !ii
!      call tetvol(cl,vol2)
!      do ii=1,3
!      cl(1,ii)=coe(3,ii)-coe(1,ii)
!      cl(2,ii)=coe(4,ii)-coe(1,ii)
!      cl(3,ii)=coe(5,ii)-coe(1,ii)
!      enddo !ii
!      call tetvol(cl,vol1)
!      do ii=1,3
!      cl(1,ii)=coe(4,ii)-coe(2,ii)
!      cl(2,ii)=coe(3,ii)-coe(2,ii)
!      cl(3,ii)=coe(5,ii)-coe(2,ii)
!      enddo !ii
!      call tetvol(cl,vol4)
!c
!      vol1=vol1/vole
!      vol2=vol2/vole
!      vol3=vol3/vole
!      vol4=vol4/vole
!c
!      if(abs(vol1+vol2+vol3+vol4-1).le.1/1000.)then
!      count=1
!      call test_vol(vol1,count)
!      call test_vol(vol2,count)
!      call test_vol(vol3,count)
!      call test_vol(vol4,count)
!      if(count.eq.1)then
!      ktmp=k
!      goto 123
!      endif
!      endif
!
!      enddo !k
!
!      write(6,*) 'something wrong?'
!      stop
!
!123   damy=0.
!      cc(1)=1-vol1-vol2-vol3
!      cc(2)=vol1
!      cc(3)=vol2
!      cc(4)=vol3
!c
!      count=0
!      do ii=1,4
!!      if(cc(ii).ne.0)then
!      if(abs(cc(ii)).gt.0.0001)then
!      coelist(iptr)=cc(ii)
!      indlist(iptr)=cny4(ktmp,ii)
!c      write(12,*) cc(ii)
!c      write(11,*) cny4(k,ii)
!      iptr=iptr+1
!      count=count+1
!      endif
!      enddo !ii
!      if((count.ne.1).and.(count.ne.2))then
!      write(*,*) 'count id',count,id,ktmp
!      write(*,*) 'cny4'
!      write(*,*) cny4(ktmp,1),cny4(ktmp,2),cny4(ktmp,3),cny4(ktmp,4)
!      write(*,*) 'cc'
!      write(*,*) cc(1),cc(2),cc(3),cc(4)
!      write(*,*) 'coe'
!      do ii=1,5
!      write(*,*) coe(ii,1),coe(ii,2),coe(ii,3)
!      enddo
!      stop
!      endif
!
!      write(10,*) iptr       
!c
!      enddo ! id
!      close(10)
!c
!      write(11,*) 'num of components'
!      write(11,*) iptr-1
!      write(11,*) 'components'
!      write(12,*) 'num of components'
!      write(12,*) iptr-1 
!      write(12,*) 'components'
!      do i=1,iptr-1
!      write(11,*) indlist(i)
!      write(12,*) coelist(i)
!      enddo 
!      close(11)
!      close(12)
!
!      end
c_______________________________________________________________________
      subroutine read_tet4_geometry_hdf(n,ne,coor,cny,num,fid)
      implicit none
      integer*4 n,ne,i,j
      real*8 coor(n,3)
      integer*4 cny(ne,4),num(ne)
      character dataname*50
      integer*4 connbuf(ne*5)
      integer*8 fid,isize
      real*8 coorbuf(n*3)
      integer*8 bufsize,bufcount
 !
      write(dataname,'(a5)') '/coor'
      isize=5
      call set_null_char(dataname,isize)
      bufsize=n*3
      call hdf_read_double_array
     -(fid,dataname,coorbuf,bufsize,bufcount)
      
      write(dataname,'(a5)') '/conn'
      isize=5
      call set_null_char(dataname,isize)
      bufsize=ne*5
      call hdf_read_int_array
     -(fid,dataname,connbuf,bufsize,bufcount)

      do i=1,n
      do j=1,3
      coor(i,j)=coorbuf(3*(i-1)+j)
      enddo
      enddo
            
      do i=1,ne
      do j=1,4
      cny(i,j)=connbuf(5*(i-1)+j)
      enddo
      num(i)=connbuf(5*(i-1)+5)
      enddo
      
      end
c_______________________________________________________________________ 
!      subroutine tetvol(cl,vol)
!      real*8 vol,cl(3,3)
!      vol= abs((cl(1,3)*(-(cl(2,2)*cl(3,1)) + cl(2,1)*cl(3,2)) + 
!     -    cl(1,2)*(cl(2,3)*cl(3,1) - cl(2,1)*cl(3,3)) + 
!     -    cl(1,1)*(-(cl(2,3)*cl(3,2)) + cl(2,2)*cl(3,3)))/6.)
!      end      
c_______________________________________________________________________ 
!      subroutine test_vol(vol,count)
!      real*8 vol
!      integer* count,counttmp
!      counttmp=0
!      if(abs(vol-1.0).lt.0.0001)counttmp=1
!      if(abs(vol-0.5).lt.0.0001)counttmp=1
!      if(abs(vol-0.0).lt.0.0001)counttmp=1
!      if(counttmp.eq.0)count=0
!      end
c_______________________________________________________________________ 
