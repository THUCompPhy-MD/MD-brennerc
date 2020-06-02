c 
c subroutine for adding incoming ion 
c 
c*********************************************************************
      subroutine bombard
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
      dimension rbom(3),ein(3),ven(3) 
c
c add incoming atom
c
      open(23,file='bombard.d',status='old')
c 
c read position and atomic # for incoming atom
c 
      READ(23,*) K,NATOM,(Rbom(i),i=1,3),iit 
  
c
c no new read for negative k 
c 
      if(k.lt.0) return 
C
C read energies in eV for each direction
C
      READ(23,*) (ein(n),n=1,3)

      close(23) 

      np = np + 1 
      if(np.gt.npmax) then
           write(*,*) 'np= ',np,' greater than npmax= ',npmax
           write(*,*) 'increase npmax and recompile'
      endif

      do n=1,3
           if(Rbom(n).gt.cube(n)/2.0d0) Rbom(n)=Rbom(N)-cube(n)
           r0(np,n) = Rbom(n) 
      enddo

      KTYPE(np)=KT(NATOM)
      if(KTYPE(np).eq.0) then
            write(*,*) 'unknown atom type for incoming atom'
            stop
      endif 

      itr(np) = 0 

c 
c convert to velocities in A/fs
c 
      do n=1,3 
           r1(np,n)= sign(1.0d0,ein(n))*sqrt(2.0d0*
     &               abs(ein(n))/xmass(ktype(np))
     &               *1.602189d-19*6.023d+23*1.0d+3)
     &               *1.0d10/1.0d15*delta 
           R0L(np,n)=R0(np,n)
           r2(np,n) = 0.0d0
           r3(np,n) = 0.0d0
           r4(np,n) = 0.0d0
      enddo
c
c write out information back out with negative k; this 
c keeps it from being added again for restart
c
      open(23,file='bombard.d',status='old')
      write(23,350) -K,NATOM,(Rbom(i),i=1,3),iit
      write(23,300) (ein(n),n=1,3)
      close(23)

      return
  300 format(3f16.6) 
  350 FORMAT(2I5,3E20.11,I3)
      end
