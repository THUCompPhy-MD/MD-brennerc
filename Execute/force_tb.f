       subroutine FORCE

       implicit real*8(a-h,o-z)
      include 'common_files.inc'

C
C     Calculates the forces acting on the atoms: the Hellman-Feynman
C     force and the force due to the classical repulsion potential
C
C     The Hellman-Feynman force is
C     F_n = sum_j <j| e^ik.R  dH/dR_n |j> + ik <j| e^ik.R H |j>
C
C     Variables:
C        DH     --  Hamiltonin derivative
C        VR     --  Eigenvector (real part) VR(*,i) (i:th eig.vec.)
C        HF     --  dH |j> assistant vector
C        FHF    --  Hellman-Feynman force
C        F(n,3) --  the force acting on particle n
C


       p56 = tbpar(5)-tbpar(6)

       do 2 i=1,4*Np
       do 2   j=1,neigen
           VATB(j,i)=VR(i,j)
2      continue

       do j=1,3
         do k=1,Np
             FHF(k,j) = 0.0d0
         enddo
       enddo


       do 50 i  = 1,Np
         do 40 k  = 1,Nntb(i)
           j  = Nident(i,k) - 1
           x  = r0(i,1)-xnaapurit(i,k,1)
           y  = r0(i,2)-xnaapurit(i,k,2)
           z  = r0(i,3)-xnaapurit(i,k,3)
           rsq = x*x+y*y+z*z
           r  = sqrt(rsq)
           apu= smooth(r)/r
           dsm= Dsmooth(r)
C          dsm = derivative of the function smooth
           x = x/r
           y = y/r
           z = z/r
           capu2 =  dsm - 2.*apu
           capu3 =  dsm - 3.*apu
           capu4 =  dsm - 4.*apu

           DH(1,1,1) = capu2*tbpar(3)*x
           DH(1,1,2) = capu2*tbpar(3)*y
           DH(1,1,3) = capu2*tbpar(3)*z

           DH(1,2,1) = (capu3*x*x + apu)*tbpar(4)
           DH(1,2,2) = capu3*y*x*tbpar(4)
           DH(1,2,3) = capu3*z*x*tbpar(4)
           DH(1,3,1) = capu3*x*y*tbpar(4)
           DH(1,3,2) = (capu3*y*y + apu)*tbpar(4)
           DH(1,3,3) = capu3*tbpar(4)*z*y
           DH(1,4,1) = capu3*tbpar(4)*x*z
           DH(1,4,2) = capu3*tbpar(4)*y*z
           DH(1,4,3) = (capu3*z*z + apu)*tbpar(4)

           DH(2,2,1) = (p56*(capu4*x*x+2.*apu) + capu2*tbpar(6))*x
           DH(2,2,2) = (p56*capu4*x*x + capu2*tbpar(6))*y
           DH(2,2,3) = (p56*capu4*x*x + capu2*tbpar(6))*z
           DH(3,3,1) = (p56*capu4*y*y + capu2*tbpar(6))*x
           DH(3,3,2) = (p56*(capu4*y*y+2.*apu) + capu2*tbpar(6))*y
           DH(3,3,3) = (p56*capu4*y*y + capu2*tbpar(6))*z
           DH(4,4,1) = (p56*capu4*z*z + capu2*tbpar(6))*x
           DH(4,4,2) = (p56*capu4*z*z + capu2*tbpar(6))*y
           DH(4,4,3) = (p56*(capu4*z*z+2.*apu) + capu2*tbpar(6))*z

           DH(2,3,1) = p56*(capu4*x*x + apu)*y
           DH(2,3,2) = p56*(capu4*y*y + apu)*x
           DH(2,3,3) = capu4*p56*x*y*z
           DH(2,4,1) = p56*(capu4*x*x + apu)*z
           DH(2,4,2) = capu4*p56*x*y*z
           DH(2,4,3) = p56*(capu4*z*z + apu)*x
           DH(3,4,1) = capu4*p56*x*y*z
           DH(3,4,2) = p56*(capu4*y*y + apu)*z
           DH(3,4,3) = p56*(capu4*z*z + apu)*y


C           if ((i.eq.0).and.(j.eq.1).and.(iprint.gt.2)) write(*,'(10f7.3)')
C     .     real(DH(1,1,2)),real(DH(1,2,2)),real(DH(2,2,2)),real(DH(1,3,2)),
C     .     real(DH(2,3,2)),real(DH(3,3,2)),real(DH(1,4,2)),real(DH(2,4,2)),
C     .     real(DH(3,4,2)),real(DH(4,4,2))

           do 20 m = 1,3
             DH(2,1,m) = -DH(1,2,m)
             DH(3,1,m) = -DH(1,3,m)
             DH(4,1,m) = -DH(1,4,m)
             DH(3,2,m) =  DH(2,3,m)
             DH(4,2,m) =  DH(2,4,m)
             DH(4,3,m) =  DH(3,4,m)
             DH(1,1,m) =  0.5*DH(1,1,m)
             DH(2,2,m) =  0.5*DH(2,2,m)
             DH(3,3,m) =  0.5*DH(3,3,m)
             DH(4,4,m) =  0.5*DH(4,4,m)
20         continue

c           do 30 m = 1,3
c           do 30 l = 1,4
c           do 30 o = l+1,4
c           do 30 i2 = 1,neigen
c             FHF(i,m) = FHF(i,m) + 4.*VATB(i2,4*(i-1)+l)*DH(l,o,m)*
c     .       VATB(i2,4*j+o)+ 4.*VATB(i2,4*j+l)*DH(o,l,m)*VATB(i2,4*(i-1)+o)
c30         continue
           do 30 m = 1,3
           do 30 l = 1,4
           do 30 o = l,4
           do 30 i2 = 1,neigen
             FHF(i,m) = FHF(i,m) + 4.*VATB(i2,4*(i-1)+l)*DH(l,o,m)*
     .       VATB(i2,4*j+o)+ 4.*VATB(i2,4*j+l)*DH(o,l,m)
     .         *VATB(i2,4*(i-1)+o)

30         continue

40       continue
C                  END  k = 1,neigen


C      Hellman-Feynman forces are calculated,
C      Next the forces due to the classical potential are calculated

C      H-F force is in principle real
C
C       F(i,1) = -real(FHF(i,1))
C       F(i,2) = -real(FHF(i,2))
C       F(i,3) = -real(FHF(i,3))

       rnp(i,1) = rnp(i,1)-FHF(i,1)
       rnp(i,2) = rnp(i,2)-FHF(i,2)
       rnp(i,3) = rnp(i,3)-FHF(i,3)

       do 490 k = 1,Nntb(i)
         x = r0(i,1)-xnaapurit(i,k,1)
         y = r0(i,2)-xnaapurit(i,k,2)
         z = r0(i,3)-xnaapurit(i,k,3)
         r = sqrt(x*x+y*y+z*z)
         apu2 = 2.0d0*dwsum(i)*Dkorj(r)/r
         rnp(i,1) = rnp(i,1) - apu2*x
         rnp(i,2) = rnp(i,2) - apu2*y
         rnp(i,3) = rnp(i,3) - apu2*z
490    continue

50     continue
C                  END  i = 1,Np


       return
       end

