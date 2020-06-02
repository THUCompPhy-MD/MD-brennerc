c 
c move simulated diamond tip   
c 
c*********************************************************************
      subroutine move_tip 
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
C
C move tip down
C
      do i=1,np
           if((ktype(i).le.2).and.(itr(i).eq.2)) then
                     r0(i,2) = r0(i,2)  - 0.0005
           endif
      enddo
c
c add velocity to carbon tip for first step
      if(ttime.le.delta) then
           do i=1,np
                if((ktype(i).le.2).and.(itr(i).ne.2)) then
                     r1(i,2) = -0.0005d0
                endif
           enddo
      endif
      return 
      end 
