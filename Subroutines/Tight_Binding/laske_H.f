       subroutine laske_H


C      Tight-binding matrix is from
C
C      Slater and Koster  Phys. Rev.  94, 1498  (1954)
C
C      H(i,j) = sum_R <f(i,0) | H | f(j,R)> e^(ik.R)


       implicit real*8(a-h,o-z)
      include 'common_files.inc'

C      dimensions of H are (numb.of atoms)*(numb.of basis)
C      for C 4*(numb.of atoms) because minimal basis (s,px,py,pz) is used
       do i=1,4*Np
            do j=1,4*Np
                 HR(i,j) = 0.0d0
            enddo
C      diagonal elements Es and Ep
            if (mod(i,4).eq.1) then
               HR(i,i) = tbpar(1)
            else
               HR(i,i) = tbpar(2)
            endif
       enddo

       do 40 i=0,Np-2
       do 40 k=1,Nntb(i+1)
         j  = Nident(i+1,k) - 1
         x  = r0(i+1,1)-xnaapurit(i+1,k,1)
         y  = r0(i+1,2)-xnaapurit(i+1,k,2)
         z  = r0(i+1,3)-xnaapurit(i+1,k,3)
         rsq = x*x+y*y+z*z
         r  = sqrt(rsq)
c divide by r**2 done in function smooth 
         apu= smooth(r)
         x = x/r
         y = y/r
         z = z/r

         HR(j*4+1,i*4+1) = apu*tbpar(3)
         HR(j*4+2,i*4+1) = apu*tbpar(4)*x
         HR(j*4+3,i*4+1) = apu*tbpar(4)*y
         HR(j*4+4,i*4+1) = apu*tbpar(4)*z
         HR(j*4+2,i*4+2) = apu*(tbpar(5)*x*x+tbpar(6)*(1-x*x))
         HR(j*4+3,i*4+3) = apu*(tbpar(5)*y*y+tbpar(6)*(1-y*y))
         HR(j*4+4,i*4+4) = apu*(tbpar(5)*z*z+tbpar(6)*(1-z*z))

         HR(j*4+3,i*4+2) = apu*(tbpar(5)-tbpar(6))*x*y
         HR(j*4+4,i*4+2) = apu*(tbpar(5)-tbpar(6))*x*z
         HR(j*4+4,i*4+3) = apu*(tbpar(5)-tbpar(6))*y*z

         HR(j*4+1,i*4+2) = -HR(j*4+2,i*4+1)
         HR(j*4+1,i*4+3) = -HR(j*4+3,i*4+1)
         HR(j*4+1,i*4+4) = -HR(j*4+4,i*4+1)

         HR(j*4+2,i*4+3) =  HR(j*4+3,i*4+2)
         HR(j*4+2,i*4+4) =  HR(j*4+4,i*4+2)
         HR(j*4+3,i*4+4) =  HR(j*4+4,i*4+3)

40     continue

       do 50 i=1,4*Np
       do 50 j=i,4*Np
50        HR(i,j) =  HR(j,i)

       return
       end

