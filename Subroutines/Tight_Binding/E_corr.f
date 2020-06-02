
               function E_corr(m)
C   --------------------------------------------------------------------------

       implicit real*8(a-h,o-z)
      include 'common_files.inc'

       external xkorjaus

       E_corr = 0.0d0
       do i = 1,Np
         s = 0.0d0
         do k = 1,Nntb(i)
           dist = sqrt( (r0(i,1)-xnaapurit(i,k,1))**2 +
     .                  (r0(i,2)-xnaapurit(i,k,2))**2 +
     .                  (r0(i,3)-xnaapurit(i,k,3))**2 )
           s = s + xkorjaus(dist)
         enddo
         wsum(i) = (f0 + s*(f1 + s*(f2 + s*(f3 + s*f4))))
         dwsum(i) = f1 + s*(2.0d0*f2 + s*(3.0d0*f3 + s*4.0d0*f4))
         E_corr = E_corr + wsum(i)
       enddo
10     continue
       return
       end

