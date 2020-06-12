      subroutine setpc
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'

C
C NORDSIECK PREDICTOR-CORRECTOR COEFFICIENTS
C
      F02=1.0d0/6.0d0
      F12=5.0d0/6.0d0
      F32=1.0d0/3.0d0
c
      return 
      end 
c 
      SUBROUTINE CPRED
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'
      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
           IF(NMA.NE.0) CALL PRED
           TTIME=TTIME+DELTA
           time = time + delta
      ENDIF
      return
      end 

      SUBROUTINE ccorr
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'
      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
              IF(NMA.NE.0) CALL CORR
      endif
      return
      end 

      SUBROUTINE PRED
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'
C
C THIRD ORDER PREDICTOR
C
      DO 31 J=1,3
           DO 30 II=1,NMA
                I = mlist(II)
                R0(I,J)=R0(I,J)+R1(I,J)+R2(I,J)+R3(I,J)
                R1(I,J)=R1(I,J) + 2.0d0*R2(I,J) + 3.0d0*R3(I,J)
                R2(I,J)=R2(I,J) + 3.0d0*R3(I,J)
                R0(I,J)=R0(I,J)-CUBE(J)*ANINT(R0(I,J)/CUBE(J))
   30      CONTINUE
   31 CONTINUE
      RETURN
      END
C
C*****
C
      SUBROUTINE CORR
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
C
      DIMENSION COM(3)
c
C  THIRD ORDER CORRECTOR
C
      DE=DELTSQ/ECONV
      DO 1 J=1,3
      DO 2 II=1,NMA
           I = mlist(II)
           XXM=XMASS(KTYPE(I))
           RI=R2(I,J) - DE*RNP(I,J)/XXM
           R0(I,J)=R0(I,J) - RI*F02
           R1(I,J)=R1(I,J) - RI*F12
           R2(I,J)=R2(I,J) - RI
           R3(I,J)=R3(I,J) - RI*F32
C
           R0(I,J)=R0(I,J)-CUBE(J)*ANINT(R0(I,J)/CUBE(J))


C
   2  CONTINUE
   1  CONTINUE

      RETURN
      END

