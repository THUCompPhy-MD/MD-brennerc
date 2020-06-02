c 
c calculated load on simulated diamond tip   
c 
c*********************************************************************
      subroutine load  
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
C
      xload = 0.0d0 
      do i=1,np
           if((ktype(i).le.2).and.(itr(i).eq.2)) then
                xload = xload + rnp(i,2)  
           endif
      enddo
      write(6,100) lstep, xload * 1.6022 
      return 
100   format(i8,f12.6) 
      end 
c
