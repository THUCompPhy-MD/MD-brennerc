c 
c reset lists of thermostated atoms according to their position; 
c all others converted to moving atoms 
c 
c*********************************************************************
      subroutine reset_itr 
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      parameter( tpo = 20.0d0) 

      nma = 0
      nta = 0

      do i=1,np
           nma = nma + 1 
           mlist(nma) = i
           itr(i) = 0
           if(r0(i,2).lt.tpo) then
                nta = nta + 1
                nlist(nta) = i
                itr(i) = 1
           endif 
      enddo 
      return 
      end 
