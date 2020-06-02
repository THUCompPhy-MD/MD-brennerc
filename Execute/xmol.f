      SUBROUTINE XMOL
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      idum=3
C
      WRITE(1,1800) HEAD
      WRITE(1,100) NP,IDUM,NRA,NLA
      WRITE(1,300) TTIME,DELTA
      WRITE(1,300) (CUBE(N),N=1,3)
C
      DO 11 I=1,NP
           WRITE(1,350) I,KT2(KTYPE(I)),(R0(I,N),N=1,3)
11    CONTINUE
C
      call flush(1)
c write one-half pair energy for each atom
        do i=1,np
             write(85,*) i,ktype(i),eatom(i)
        enddo
      call flush(85)
  100 FORMAT(4I6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(20A2)
      return
      end
