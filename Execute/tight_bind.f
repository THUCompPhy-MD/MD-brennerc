      subroutine tight_bind
C
C----------------------------------------------------------------------------
C
C                  TIGHT BINDING ELECTRONIC STRUCTURE
C
C       Original Version by K. Laasonen and R. Virkkunen
C
C
C       The program reads the atomic coordinates and calculates the
C       tight-binding matrix H, the core-core energy, and the total
C       energy.
C
C       The eigenvalues e are calculated from equation  H x = e x
C       (where H is hermitian matrix) using attached eispack routines.
C
C       Units:   eV = 1    A = 1
C
C       Only the gamma point in used in this version
C
C
C----------------------------------------------------------------------------
C
C
       implicit real*8(a-h,o-z)
      include 'common_files.inc'
       external E_corr
C
C      naapurit = neighbours, paikat = places (coordinates) of the atoms,
C

C----------------------------------------------------------------------------
C
C
C         TB parametrization for the Hamiltonian is from
C         Xu, Wang, Chan, and Ho, J. Phys.: Condens. Matter 4 (1992) 6047.
C

C      neigen=number of occupied eigenstates
C      Emintb,Emax: Minimum and maximum eigenenergies
C      basis=basis vectors of the lattice
C      rcut=distance determining if a neighbor atom is a nearest neighbor or notC
C
      if(np.ne.natx) then
            write(*,*) 'set natx=',np,' in common_tb.inc and recompile'
            include 'close.inc'
            stop
      endif
      do i=1,3
        do j=1,3
           basis(i,j) = 0.0d0
        enddo
      enddo
      do i=1,3
           basis(i,i) = cube(i)
      enddo

       rcut = 2.60d0
c      rcut = 2.0d0
       neigen=2*Np
c
C      k-points (only the gamma-point is used in this version)
c       read(8,*) Nk
c       do 40 i=1,Nk
c         read(8,*) (xkt(i,j),j=1,3)
c       do 40 j=1,3
c         xkt(i,j) = 2*pi/scale*xkt(i,j)
c40     continue

C----------------------------------------------------------------------------

C energy for separated atoms

c classical contribution:
       e_cl = np*f0
c quantum contribution:
       e_qu = np*(2.0d0*tbpar(1) + 2.0d0*tbpar(2))
       e_separ = e_cl + e_qu
       ip = 1
       n  = 4*Np
       iik = 0
       NE_yli = 0

C      determines the nearest neighbors
       CALL neighbor

C      calculates Hamiltonian matrix
       CALL laske_H


       ifail = 1
C eispack call
       ierr = 0
       matz = 1


c real Hermitian Matrix
       call rs(nmax,nmax,HR,E,matz,VR,fv1,fv2,ierr)

c
c write out all eigenvectors 
c
c      do i=1,256
c           do j=1,256
c                write(36,900) i,j,VR(i,j)
c           enddo
c      enddo
c900   format(2i6,f12.6)

c
c calculate pair energy 
       ecl = E_corr(i)

C calculate forces
      call force

C
C      calculate the potential energy
C
       summa = 0.0d0
       do i = 1,neigen
           summa = summa + 2.*E(i)
       enddo

       tote = summa + ecl  - E_separ
       write(*,*) 'band energy,bond energy= ',summa,summa - E_qu

C  write final global and local densities of states, eigenenergies 
C and eigenvectors 
      if(lstep.eq.kvc) call dos

      return

94     format('  BS energy ',F13.6,'  Core energy ',F13.6, /,
     .           'total PE/atom ',f13.6)
       end

