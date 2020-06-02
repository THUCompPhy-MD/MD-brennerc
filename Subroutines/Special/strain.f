c 
c    strain system 
c 
c*********************************************************************
      subroutine strain  
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      parameter(idstr = 2) 
      dimension csav(3), rsav(npmax,3) 
c 
c store original data
c
      if(ttime.lt.delta) then
           do j=1,3
              csav(j) = cube(j)
              do i=1,np
                 rsav(i,j) = r0(i,j)
              enddo
           enddo
      endif 
C
C  apply strain 
C
      strain = 0.01*(lstep - kvc/2) 
      do j=1,3
           cube(j) = csav(j)*strain
           cube2(j) = cube(j)/2.0d0
           do i=1,np
                r0(i,j) = rsav(i,j)*strain
           enddo
      enddo 
      return 
      end 
