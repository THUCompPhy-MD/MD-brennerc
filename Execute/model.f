      SUBROUTINE MODEL
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
c
c set forces to zero
C
      DO J=1,3
           DO I=1,NP
                rnp(i,j) = 0.0d0
           ENDDO
      ENDDO

C set potential energy to zero
      tote = 0.0d0

       sss = noa(1)+noa(2)+noa(3)+noa(4)
       if((ipot.eq.1).and.(sss.ne.0)) CALL CAGUTS
       if((ipot.eq.2).and.((noa(1)+noa(2)+noa(3)).ne.0)) call tight_bind
       if(ILJ.eq.1) call ljguts
c      call ljcont
      return
      end
