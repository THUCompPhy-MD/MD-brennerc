c 
c subroutine for finding atom with maximum ke  
c 
c*********************************************************************
      subroutine max_ke 
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      ccm = 1.0d3/(2.0d0*1.602189d0*6.023d0)  
      emmax = 0.0d0  
      do i=1,np 
           emv = 0.0d0 
           do n=1,3 
                emv = emv + ((r1(i,n)/delta)**2)*xmass(ktype(i))*ccm 
           enddo
           if(emv.gt.emmax) then
                 nemax = i
                 emmax = emv 
           endif 
      enddo 
      write(55,100) nemax,emmax,(r0(nemax,n),n=1,3) 
      call flush(55) 
      return
100   format(i6,f7.2,3x,3f7.2) 
      end
