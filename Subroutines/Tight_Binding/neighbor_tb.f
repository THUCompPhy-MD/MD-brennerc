       subroutine neighbor
C   --------------------------------------------------------------------------
C
C      determines the nearest neighbors of the atoms
C      also handles the periodic boundary conditions
C
C    Variables:
C      paikat   -  coordinates of the atoms
C      naapurit -  neighbours of the atoms  (real coordinates)
C      Nntb       -  number of the neighbours of the atoms
C      Nident   -  identification number of the atom in the cell
C      rcut     -  nearest neighbour cut-off
C

       implicit real*8(a-h,o-z)
      include 'common_files.inc'


       rcut2 = rcut*rcut
       do 30 i = 1,Np
         n = 0
         do 20 k = 0,26
           n3 = k/9 - 1
           n2 = (k-9*(n3+1))/3 - 1
           n1 = mod(k,3) - 1
           vec1(1) = n1*basis(1,1) + n2*basis(2,1) + n3*basis(3,1)
           vec1(2) = n1*basis(1,2) + n2*basis(2,2) + n3*basis(3,2)
           vec1(3) = n1*basis(1,3) + n2*basis(2,3) + n3*basis(3,3)
           do 10 j = 1,Np
             vec2(1) = r0(j,1) + vec1(1)
             vec2(2) = r0(j,2) + vec1(2)
             vec2(3) = r0(j,3) + vec1(3)


             dist = (r0(i,1)-vec2(1))**2 +
     .              (r0(i,2)-vec2(2))**2 + (r0(i,3)-vec2(3))**2


             if ((dist.lt.rcut2).and.(dist.gt.0.01)) then
               if (dist.lt.4.84) iik = iik + 1
               n = n + 1
               Nident(i,n) = j
               xnaapurit(i,n,1) = vec2(1)
               xnaapurit(i,n,2) = vec2(2)
               xnaapurit(i,n,3) = vec2(3)
               if (iprint.gt.2)
     .           write(*,'(i4,3f9.4,i5)') i,vec2(1),vec2(2),vec2(3),j
             endif
10         continue
20       continue
         Nntb(i) = n
30     continue
       return
       end


