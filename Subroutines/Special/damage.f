c 
c count # carbon-carbon bonds before and after 
C simulation; use any differences to characterize
C if damage via C-C bond breaking has occurred.  
c 
c*********************************************************************
      subroutine damage 
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      common/dama/ncc1,ncc2 
      dist = 1.8d0*1.80d0 

      ncc1 = 0 
      do i=1,np-1
           do j=i+1,np 
                if((ktype(i).eq.1).and.(ktype(j).eq.1)) then 
                     rsq = 0.0d0
                     do k=1,3 
                          rr = r0(i,k)-r0(j,k)
                          rr = rr - CUBE(k)*ANINT(RR/CUBE(k))
                          rsq = rsq + rr*rr 
                     endif 
                     if(rsq.lt.dist) ncc1 = ncc1 + 1 
                endif 
           enddo
      enddo 
      if(lstep.eq.1) then
           ncc2 = ncc1
           return
      else 
           if(ncc2.ne.ncc1) then
                call xmol 
                include 'close.inc'
                stop
           endif
      endif 
      return
      end
