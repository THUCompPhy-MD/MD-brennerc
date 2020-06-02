      SUBROUTINE HCHANG
C
C  PROGRAM TO CHANGE STEP SIZE
C
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
c
      IF(KFLAG.LT.0) THEN
           IF(NLA.GE.NP) RETURN
           IST=NLA+1
      ELSE
           IST=NRA+1
      ENDIF
      XMAX=0.0d0
      ALPH=1.0d0
C
      DLMAX=5.00d0
      DLMIN=0.25d0
C
      DO 2 J=1,3
           DO 1 I=IST,NP
                DR=ABS(R3(I,J)-R4(I,J))
                XMAX=DMAX1(XMAX,DR)
1          CONTINUE
2     CONTINUE
      IF(XMAX.EQ.0.0d0) RETURN
C
      ALPH=1.2d0*((EPAR1/XMAX)**EPAR2)
C****
c      ALPH=DLMIN/DELTA
C****
      DOLD=DELTA
      DELTA=DELTA*ALPH
c
      IF(DELTA.GT.DLMIN) GO TO 3
      ALPH=DLMIN/DOLD
      DELTA=DLMIN
      GO TO 4
3     CONTINUE
      IF(DELTA.LT.DLMAX) GO TO 4
      ALPH=DLMAX/DOLD
      DELTA=DLMAX
4     CONTINUE
C
      DELTSQ=DELTA*DELTA/2.0d0
      BET=WD*PI/6.0d0/DELTA
      GSIG=SQRT(2.0d0*TR*BET)
      DNLA=3.0d0*DELTA*DELTA*FLOAT(NTA)
C
      ALP2=ALPH*ALPH
      ALP3=ALP2*ALPH
      DO 6 J=1,3
      DO 5 I=1,NP
           R1(I,J)=R1(I,J)*(ALPH)
           R2(I,J)=R2(I,J)*(ALP2)
           R3(I,J)=R3(I,J)*(ALP3)
           R4(I,J)=R3(I,J)
5     CONTINUE
6     CONTINUE
C
      RETURN
      END

