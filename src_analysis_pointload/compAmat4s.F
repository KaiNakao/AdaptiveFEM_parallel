c_______________________________________________________________________
      subroutine compAmat4s
     &(n,ni,ne,nei,kd,num,cny,coor,younglst,
     & cnyi,p,q,sig,numi,nfla,dsi,itype,
     & im,nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,
     & nsep,numsep,nsep_n_ptr,nsep_ne_ptr,nsepnmap,
     & cnysep_read,cnysep_write,qsep,
     & np,ncolor,color_ind)
      implicit none
      integer n,ni,ne,nei,kd,itype,im
      integer np,ncolor,color_ind(np*MAXCOLOR)
      REAL_4 q(3*(n+ni)*npc),p(3*(n+ni)*npc)
      integer cny(4,ne),num(ne),cnyi(6,nei)
      REAL_4 coor(3,(n+ni)),younglst(2,kd)
      REAL_4 young,rnyu,gxy,gyz,gxz,detjmat
      REAL_4 xx(4,3)
      REAL_4 arayi(12),arayo(12),keit(18,18)
      REAL_4 sig(5),arayiie(18),arayoie(18)
      REAL_4 bufi(21*3*6),bufo(21*3*6),bufot(18,21)
      REAL_4 cori(3,3),dsi
      REAL_4 me(4,4)
      REAL_4 x2,y2,z2,x3,y3,z3,x4,y4,z4
      integer nfla(nei),numi(nei)
      integer ie,in,i,i1,j1,i2,k,ii,ipc,ikd,iu,icolor

      integer nir,niw
      integer pentawptr(nsta+2),pentarptr(niw+1)
      integer pentawind(niw),pentarind(nir)
      REAL_4 pentacrsval(9,nir)

      integer cnysep_read(4,ne),cnysep_write(4,ne)
      integer nsep_ne_ptr(nsta+2),numsep(ne)
      integer nsep_n_ptr(nsta+2)
      integer nsep,nsepnmap(nsep)
      REAL_4 qsep(3*nsep*npc)
      REAL_4 tmpq(npc),tmpqq(npc)

c
!$OMP PARALLEL DO default(none),shared(q,n,ni)
      do i=1,3*npc*(n+ni)
      q(i)=0.
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO default(none),shared(qsep,nsep)
      do i=1,3*npc*nsep
      qsep(i)=0.
      enddo
!$OMP END PARALLEL DO
c
      do icolor=1,ncolor
!$OMP PARALLEL DO default(none),
!$OMP& shared(np,icolor,color_ind,numsep,
!$OMP&  younglst,coor,cnysep_read,p,cnysep_write,qsep),
!$OMP& private(iu,ie,in,young,rnyu,gxy,gyz,gxz,i1,j1,xx,i2,bufi,ipc,
!$OMP&  arayi,arayo,detjmat)
      do iu=1,np
      do ie=color_ind((icolor-1)*np+iu)+1,color_ind((icolor-1)*np+iu+1)
      in=numsep(ie)
      young=younglst(1,in)
      rnyu=younglst(2,in)
c
      gxy=young/(2.*(1+rnyu))
      gyz=young/(2.*(1+rnyu))
      gxz=young/(2.*(1+rnyu))
c
      do i1=1,4
      do j1=1,3
      xx(i1,j1)=coor(j1,cnysep_read(i1,ie))
      enddo
      enddo

c----------------------------------------------
      do i1=1,4
      do j1=1,63
      i2=63*(i1-1)+j1
      bufi(i2)=p(63*(cnysep_read(i1,ie)-1)+j1)
      enddo
      enddo
      do ipc=1,npc
      do i1=1,4
      do j1=1,3
      arayi(3*(i1-1)+j1)=bufi(63*(i1-1)+npc*(j1-1)+ipc)
     &                  -bufi(          npc*(j1-1)+ipc)
      enddo
      enddo
      arayo=0.
      call tet4kus(young,gxy,gyz,gxz,rnyu,xx,arayi,arayo,detjmat)
      do i1=1,4
      do j1=1,3
       qsep(63*(cnysep_write(i1,ie)-1)+(j1-1)*npc+ipc)=
     & qsep(63*(cnysep_write(i1,ie)-1)+(j1-1)*npc+ipc)+
     & arayo(3*(i1-1)+j1)
      enddo
      enddo
      enddo

      enddo ! ie
      enddo ! iu
!$OMP END PARALLEL DO
      enddo ! icolor

      do ikd=1,nsta+1
!$OMP PARALLEL DO default(none),
!$OMP& shared(nsep_n_ptr,ikd,qsep,sig,q,nsepnmap),
!$OMP& private(i,i1,ipc,tmpq,tmpqq)
      do i=nsep_n_ptr(ikd)+1,nsep_n_ptr(ikd+1)
      do i1=1,3
      do ipc=1,npc
      tmpq(ipc)=qsep(3*npc*(i-1)+npc*(i1-1)+ipc)
      enddo
      call compsigs_npc21(ikd,sig,tmpq,tmpqq)
      do ipc=1,npc
       q(3*npc*(nsepnmap(i)-1)+(i1-1)*npc+ipc)=
     & q(3*npc*(nsepnmap(i)-1)+(i1-1)*npc+ipc)+tmpqq(ipc)
      enddo ! ipc
      enddo ! i1
      enddo ! i
!$OMP END PARALLEL DO
      enddo ! ikd

#ifdef USEIEPENTA
      call compAmats_iepenta
     & (im,n,ni,
     & nir,pentarptr,pentarind,
     & niw,pentawptr,pentawind,pentacrsval,sig,p,q)
#endif
      end
