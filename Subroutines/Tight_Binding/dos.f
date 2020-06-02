      subroutine dos
C     calculates global and local densities of electronic states
      implicit real*8(a-h,o-z)
      include 'common_files.inc'

C******************************
C write eigenvectors
c
      do i=1,nmax
           do j=1,nmax
                write(85,900) i,j,VR(i,j)
           enddo
      enddo

C******************************
c write eigenenergies and calculate total density of 
C states

      emax = e(2*neigen)
      emintb = e(1)
      ndiv = 1000 
      del = (emax-emintb)/float(ndiv) 
      do i=1,ndiv
          idos(i) = 0.0d0
      enddo 

      do i=1,2*neigen
          write(86,*) i,e(i)
          n = nint((e(i)-emintb)/del) +1 
          idos(n) = idos(n) + 1 
      enddo 

c write out total density of states 

      do i=1,ndiv+1
          write(87,101) del*float(i-1)+emintb,
     &                  float(idos(i))/float(np)
      enddo 

      ido = 0
      if(ido.eq.0) then
          return
      else 
C*******************************************
c calculate local density of states for atom #=iatom
c
10    continue 
      write(*,*) 'enter iatom'
      read(*,*) iatom
      if(iatom.eq.0) return 

      do i=1,ndiv
          xdos(i) = 0.0d0
      enddo

c set index for finding correct atomic orbitals for atom
      ii = (iatom-1)*4 + 1

      do i=1,2*neigen
          wt = vr(ii,i)**2+vr(ii+1,i)**2+
     &         vr(ii+2,i)**2+vr(ii+3,i)**2
          n = nint((e(i)-emintb)/del) +1
          xdos(n) = xdos(n) + wt 
      enddo

      write(88,100) ndiv+1,e(neigen)

      do i=1,ndiv+1
          write(88,101) del*float(i-1)+emintb,xdos(i)
      enddo

      go to 10 
      endif 

900   format(2i6,f12.6)
100   format(i5,f12.6)
101   format(2f12.6)
       end

