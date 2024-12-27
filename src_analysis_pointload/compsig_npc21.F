      subroutine compsigd_npc21(ikd,sig,tmpq,tmpqq)
      implicit none
      integer, parameter :: nkl5 = 5
      REAL*8 sig(nkl5),tmpq(npc),tmpqq(npc)
      integer ikd
      integer j

      do j=1,npc
      tmpqq(j)=tmpq(j)
      enddo
      tmpqq(7)=tmpq(7)*2.
      tmpqq(9)=tmpq(9)*2.
      tmpqq(12)=tmpq(12)*2.
      tmpqq(16)=tmpq(16)*2.
      tmpqq(21)=tmpq(21)*2.

      if(ikd.eq.1)then
      tmpqq(1)=tmpqq(1)+tmpq(2)*sig(1)
      tmpqq(2)=tmpqq(2)+(tmpq(1)+tmpq(7)*2.)*sig(1)
      tmpqq(3)=tmpqq(3)+tmpq(8)*sig(1)
      tmpqq(4)=tmpqq(4)+tmpq(10)*sig(1)
      tmpqq(5)=tmpqq(5)+tmpq(13)*sig(1)
      tmpqq(6)=tmpqq(6)+tmpq(17)*sig(1)
      tmpqq(7)=tmpqq(7)+tmpq(2)*2.*sig(1)
      tmpqq(8)=tmpqq(8)+tmpq(3)*sig(1)
      tmpqq(10)=tmpqq(10)+tmpq(4)*sig(1)
      tmpqq(13)=tmpqq(13)+tmpq(5)*sig(1)
      tmpqq(17)=tmpqq(17)+tmpq(6)*sig(1)
      endif

      if(ikd.eq.2)then
      tmpqq(1)=tmpqq(1)+tmpq(3)*sig(2)
      tmpqq(2)=tmpqq(2)+tmpq(8)*sig(2)
      tmpqq(3)=tmpqq(3)+(tmpq(1)+2.*tmpq(9))*sig(2)
      tmpqq(4)=tmpqq(4)+tmpq(11)*sig(2)
      tmpqq(5)=tmpqq(5)+tmpq(14)*sig(2)
      tmpqq(6)=tmpqq(6)+tmpq(18)*sig(2)
      tmpqq(8)=tmpqq(8)+tmpq(2)*sig(2)
      tmpqq(9)=tmpqq(9)+tmpq(3)*2.*sig(2)
      tmpqq(11)=tmpqq(11)+tmpq(4)*sig(2)
      tmpqq(14)=tmpqq(14)+tmpq(5)*sig(2)
      tmpqq(18)=tmpqq(18)+tmpq(6)*sig(2)
      endif

      if(ikd.eq.3)then
      tmpqq(1)=tmpqq(1)+tmpq(4)*sig(3)
      tmpqq(2)=tmpqq(2)+tmpq(10)*sig(3)
      tmpqq(3)=tmpqq(3)+tmpq(11)*sig(3)
      tmpqq(4)=tmpqq(4)+(tmpq(1)+2.*tmpq(12))*sig(3)
      tmpqq(5)=tmpqq(5)+tmpq(15)*sig(3)
      tmpqq(6)=tmpqq(6)+tmpq(19)*sig(3)
      tmpqq(10)=tmpqq(10)+tmpq(2)*sig(3)
      tmpqq(11)=tmpqq(11)+tmpq(3)*sig(3)
      tmpqq(12)=tmpqq(12)+tmpq(4)*2.*sig(3)
      tmpqq(15)=tmpqq(15)+tmpq(5)*sig(3)
      tmpqq(19)=tmpqq(19)+tmpq(6)*sig(3)
      endif

      if(ikd.eq.4)then
      tmpqq(1)=tmpqq(1)+tmpq(5)*sig(4)
      tmpqq(2)=tmpqq(2)+tmpq(13)*sig(4)
      tmpqq(3)=tmpqq(3)+tmpq(14)*sig(4)
      tmpqq(4)=tmpqq(4)+tmpq(15)*sig(4)
      tmpqq(5)=tmpqq(5)+(tmpq(1)+2.*tmpq(16))*sig(4)
      tmpqq(6)=tmpqq(6)+tmpq(20)*sig(4)
      tmpqq(13)=tmpqq(13)+tmpq(2)*sig(4)
      tmpqq(14)=tmpqq(14)+tmpq(3)*sig(4)
      tmpqq(15)=tmpqq(15)+tmpq(4)*sig(4)
      tmpqq(16)=tmpqq(16)+tmpq(5)*2.*sig(4)
      tmpqq(20)=tmpqq(20)+tmpq(6)*sig(4)
      endif

      if(ikd.eq.5)then
      tmpqq(1)=tmpqq(1)+tmpq(6)*sig(5)
      tmpqq(2)=tmpqq(2)+tmpq(17)*sig(5)
      tmpqq(3)=tmpqq(3)+tmpq(18)*sig(5)
      tmpqq(4)=tmpqq(4)+tmpq(19)*sig(5)
      tmpqq(5)=tmpqq(5)+tmpq(20)*sig(5)
      tmpqq(6)=tmpqq(6)+(tmpq(1)+2.*tmpq(21))*sig(5)
      tmpqq(17)=tmpqq(17)+tmpq(2)*sig(5)
      tmpqq(18)=tmpqq(18)+tmpq(3)*sig(5)
      tmpqq(19)=tmpqq(19)+tmpq(4)*sig(5)
      tmpqq(20)=tmpqq(20)+tmpq(5)*sig(5)
      tmpqq(21)=tmpqq(21)+tmpq(6)*2.*sig(5)
      endif
      end
c_______________________________________________________________________
      subroutine compsigs_npc21(ikd,sig,tmpq,tmpqq)
      implicit none
      integer, parameter :: nkl5 = 5
      REAL_4 sig(nkl5),tmpq(npc),tmpqq(npc)
      integer ikd
      integer j

      do j=1,npc
      tmpqq(j)=tmpq(j)
      enddo
      tmpqq(7)=tmpq(7)*2.
      tmpqq(9)=tmpq(9)*2.
      tmpqq(12)=tmpq(12)*2.
      tmpqq(16)=tmpq(16)*2.
      tmpqq(21)=tmpq(21)*2.

      if(ikd.eq.1)then
      tmpqq(1)=tmpqq(1)+tmpq(2)*sig(1)
      tmpqq(2)=tmpqq(2)+(tmpq(1)+tmpq(7)*2.)*sig(1)
      tmpqq(3)=tmpqq(3)+tmpq(8)*sig(1)
      tmpqq(4)=tmpqq(4)+tmpq(10)*sig(1)
      tmpqq(5)=tmpqq(5)+tmpq(13)*sig(1)
      tmpqq(6)=tmpqq(6)+tmpq(17)*sig(1)
      tmpqq(7)=tmpqq(7)+tmpq(2)*2.*sig(1)
      tmpqq(8)=tmpqq(8)+tmpq(3)*sig(1)
      tmpqq(10)=tmpqq(10)+tmpq(4)*sig(1)
      tmpqq(13)=tmpqq(13)+tmpq(5)*sig(1)
      tmpqq(17)=tmpqq(17)+tmpq(6)*sig(1)
      endif

      if(ikd.eq.2)then
      tmpqq(1)=tmpqq(1)+tmpq(3)*sig(2)
      tmpqq(2)=tmpqq(2)+tmpq(8)*sig(2)
      tmpqq(3)=tmpqq(3)+(tmpq(1)+2.*tmpq(9))*sig(2)
      tmpqq(4)=tmpqq(4)+tmpq(11)*sig(2)
      tmpqq(5)=tmpqq(5)+tmpq(14)*sig(2)
      tmpqq(6)=tmpqq(6)+tmpq(18)*sig(2)
      tmpqq(8)=tmpqq(8)+tmpq(2)*sig(2)
      tmpqq(9)=tmpqq(9)+tmpq(3)*2.*sig(2)
      tmpqq(11)=tmpqq(11)+tmpq(4)*sig(2)
      tmpqq(14)=tmpqq(14)+tmpq(5)*sig(2)
      tmpqq(18)=tmpqq(18)+tmpq(6)*sig(2)
      endif

      if(ikd.eq.3)then
      tmpqq(1)=tmpqq(1)+tmpq(4)*sig(3)
      tmpqq(2)=tmpqq(2)+tmpq(10)*sig(3)
      tmpqq(3)=tmpqq(3)+tmpq(11)*sig(3)
      tmpqq(4)=tmpqq(4)+(tmpq(1)+2.*tmpq(12))*sig(3)
      tmpqq(5)=tmpqq(5)+tmpq(15)*sig(3)
      tmpqq(6)=tmpqq(6)+tmpq(19)*sig(3)
      tmpqq(10)=tmpqq(10)+tmpq(2)*sig(3)
      tmpqq(11)=tmpqq(11)+tmpq(3)*sig(3)
      tmpqq(12)=tmpqq(12)+tmpq(4)*2.*sig(3)
      tmpqq(15)=tmpqq(15)+tmpq(5)*sig(3)
      tmpqq(19)=tmpqq(19)+tmpq(6)*sig(3)
      endif

      if(ikd.eq.4)then
      tmpqq(1)=tmpqq(1)+tmpq(5)*sig(4)
      tmpqq(2)=tmpqq(2)+tmpq(13)*sig(4)
      tmpqq(3)=tmpqq(3)+tmpq(14)*sig(4)
      tmpqq(4)=tmpqq(4)+tmpq(15)*sig(4)
      tmpqq(5)=tmpqq(5)+(tmpq(1)+2.*tmpq(16))*sig(4)
      tmpqq(6)=tmpqq(6)+tmpq(20)*sig(4)
      tmpqq(13)=tmpqq(13)+tmpq(2)*sig(4)
      tmpqq(14)=tmpqq(14)+tmpq(3)*sig(4)
      tmpqq(15)=tmpqq(15)+tmpq(4)*sig(4)
      tmpqq(16)=tmpqq(16)+tmpq(5)*2.*sig(4)
      tmpqq(20)=tmpqq(20)+tmpq(6)*sig(4)
      endif

      if(ikd.eq.5)then
      tmpqq(1)=tmpqq(1)+tmpq(6)*sig(5)
      tmpqq(2)=tmpqq(2)+tmpq(17)*sig(5)
      tmpqq(3)=tmpqq(3)+tmpq(18)*sig(5)
      tmpqq(4)=tmpqq(4)+tmpq(19)*sig(5)
      tmpqq(5)=tmpqq(5)+tmpq(20)*sig(5)
      tmpqq(6)=tmpqq(6)+(tmpq(1)+2.*tmpq(21))*sig(5)
      tmpqq(17)=tmpqq(17)+tmpq(2)*sig(5)
      tmpqq(18)=tmpqq(18)+tmpq(3)*sig(5)
      tmpqq(19)=tmpqq(19)+tmpq(4)*sig(5)
      tmpqq(20)=tmpqq(20)+tmpq(5)*sig(5)
      tmpqq(21)=tmpqq(21)+tmpq(6)*2.*sig(5)
      endif
      end
c_______________________________________________________________________
