      subroutine read_data
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
C
C READ INPUT DATA
C
      READ(13,*) KVC,MAXKB,KFLAG,nxmol
      READ(13,*) PSEED,RLL,TEM
      READ(13,*) IPOT

      ktmax = 4
      ilj = 0

1     continue

      READ(13,*,end=2) natom, xma, epst , sigt
      if(natom.le.0) go to 1
      ilj = 1
      if(kt(natom).eq.0) then
           ktmax = ktmax + 1
           if(ktmax.gt.ntypes) then
                write(*,*) 'Maximum ntypes of ',ntypes,' exceeded'
                write(*,*) 'Change NTYPES in common_n.inc and recompile'
                include 'close.inc'
                stop
           endif
           kt(natom) = ktmax
           kt2(ktmax) = natom
      endif
      xmass(kt(natom)) = xma
      sig(kt(natom),kt(natom)) = sigt
      eps(kt(natom),kt(natom)) = epst
      go to 1

2     continue

C
C RESTART INFORMATION FROM PREVIOUS RUN
C

      READ(11,1800) HEAD
      READ(11,*) NP
      READ(11,*) TTIME,DELTA
      READ(11,*) (CUBE(N),N=1,3)

      if(np.gt.npmax) then
           write(*,*) 'np= ',np,' greater than npmax= ',npmax
           write(*,*) 'increase npmax and recompile'
           include 'close.inc'
           stop
      endif

      nma = 0
      nta = 0

      DO I=1,NP
           READ(11,*) K,NATOM,(R0(K,N),N=1,3),itr(k)
           do n=1,3
                if(R0(K,N).gt.cube(n)/2.0d0) R0(K,N)=R0(K,N)-cube(n)
           enddo
           KTYPE(K)=KT(NATOM)

           if(ktype(k).eq.0) then
                write(6,*) 'unknown atom type for atom ',k
                include 'close.inc'
                stop
           endif
           noa(ktype(k)) = noa(ktype(k)) + 1

C
c set up list of nonrigid and thermostated atoms
           if(itr(k).ne.2) then
                nma = nma + 1
                mlist(nma) = k
                if(itr(k).eq.1) then
                     nta = nta + 1
                     nlist(nta) = k
                endif
           endif
c
      ENDDO
C
      DO I=1,NP
           READ(11,*) K,(R1(K,N),N=1,3)
           DO N=1,3
                R1(k,N)=R1(k,N)*DELTA
           ENDDO
      ENDDO
C
      DO I=1,NP
           READ(11,*) K,(R2(K,N),N=1,3)
      ENDDO
C
      DO I=1,NP
           READ(11,*) K,(R3(K,N),N=1,3)
      ENDDO
C
      DO I=1,NP
           READ(11,*) K,(R4(K,N),N=1,3)
      ENDDO
C
C ESTABLISH INITIAL POSITIONS FOR NEIGHBOR LIST UPDATE
C
      DO 8 I=1,3
           DO 7 J=1,NP
                R0L(J,I)=R0(J,I)
7          CONTINUE
8     CONTINUE
C
      VOL=CUBE(1)*CUBE(2)*CUBE(3)
C
      DO 6 I=1,3
           CUBE2(I)=CUBE(I)/2.0D0
6     CONTINUE
C
      DELTSQ=DELTA*DELTA/2.0D0
      TTCONV=2.0d0/3.0d0/FLOAT(NP)

      return
  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(20A2)

      end 

      subroutine write_data1
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      WRITE(9,500) HEAD
      WRITE(9,600) KVC*MAXKB,MAXKB
      WRITE(9,700) 1,1,1,ECONV
      WRITE(9,900) TEM,PSEED
      call flush(9)
      return
  500 FORMAT('*CLASSICAL DYNAMICS SIMULATION OF ',20A2)
  600 FORMAT(/,'TOTAL STEPS= ',I6,' DATA WRITTEN EVERY ',I4,
C    &' STEPS WITH ',F7.4,' fs/STEP')
     &' STEPS')
  700 FORMAT(/,'UNITS OF LENGTH, MASS, TIME AND ENERGY:',I2,' A ',I2,
     &' AMU ',I2,' fs ',F8.4,' eV')
  800 FORMAT(/,'NEIGHBOR LIST PARAMETERS: ',F12.3,' A ')
  900 FORMAT(/,'LANGEVIN PARAMETERS: ',F12.3,' k PSEED= ',F9.6,/)
      end 

      subroutine write_data2
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'

c
C  CALCULATE KINETIC ENERGY
C
      XX=0.0d0
      DO J=1,3
           DO I=1,NP
                XX=XX+(R1(I,J)**2)*XMASS(KTYPE(I))
           enddo
      enddo

      EETOT=TOTE+XX/(4.0d0*DELTSQ)*ECONV
C
      WRITE(9,1500) LSTEP,EETOT/FLOAT(NP),
     &   ENPR*XX/6.0d0/DELTSQ*ECONV,DELTA,TIME
      call flush(9)
      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
           WRITE(6,*) 'TOTE: ',TOTE
      ENDIF

      return
 1500 FORMAT(I6,F14.8,F10.3,2F9.3)
      end 

      subroutine write_data3  
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      REWIND 11
      WRITE(11,1800) HEAD
      WRITE(11,100) NP,IDUM,NRA,NLA
      WRITE(11,300) TTIME,DELTA
      WRITE(11,300) (CUBE(N),N=1,3)
C
      DO 11 I=1,NP
           WRITE(11,350) I,KT2(KTYPE(I)),(R0(I,N),N=1,3),itr(i)
11    CONTINUE
C
      DO 12 I=1,NP
           WRITE(11,360) I,((R1(I,N)/DELTA),N=1,3)
12    CONTINUE
C
      DO 13 I=1,NP
           WRITE(11,360) I,((R2(I,N)),N=1,3)
13    CONTINUE
C
      DO 14 I=1,NP
           WRITE(11,360) I,((R3(I,N)),N=1,3)
14    CONTINUE
C
      DO 15 I=1,NP
           WRITE(11,360) I,((R4(I,N)),N=1,3)
15    CONTINUE

      return 

  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
  500 FORMAT('*CLASSICAL DYNAMICS SIMULATION OF ',20A2)
  600 FORMAT(/,'TOTAL STEPS= ',I6,' DATA WRITTEN EVERY ',I4,
C    &' STEPS WITH ',F7.4,' fs/STEP')
     &' STEPS')
  700 FORMAT(/,'UNITS OF LENGTH, MASS, TIME AND ENERGY:',I2,' A ',I2,
     &' AMU ',I2,' fs ',F8.4,' eV')
  800 FORMAT(/,'NEIGHBOR LIST PARAMETERS: ',F12.3,' A ')
  900 FORMAT(/,'LANGEVIN PARAMETERS: ',F12.3,' k PSEED= ',F9.6,/)
 1200 FORMAT('NEIGHBOR LIST UPDATES: ',/)
 1300 FORMAT(10I4,/)
 1400 FORMAT(8X,'ENERGY(eV)',5X,'T',6X,'TSTEP(fs)  TIME(fs)',/)
 1500 FORMAT(I6,F14.8,F10.3,2F9.3)
 1800 FORMAT(20A2)
      end 

