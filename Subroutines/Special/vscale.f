c*********************************************************************
      subroutine vscale
      IMPLICIT REAL*8(A-H,O-Z)
      dimension rtemp(2000,3) ,all(3),ey(3) ,ckeep(3)
      include 'common_files.inc'
c scale volume according to potential energy
c
c      go to 10
      do j=1,3
           ckeep(j)=cube(j)
           do i=1,np
                rtemp(i,j) = r0(i,j)
           enddo
      enddo
c
      kkeep = kflag
      mkeep = maxkb
      maxkb = 1
      kflag = 8
      call model
      escale=tote
      scale = -0.0010d0
      write(*,*) 'enter scaling direction'
      read(*,*) isc
78    continue
c change here to control scaling direction
      j = isc
c      do j=2,2 
           cube(j) = ckeep(j)*(scale + 1.0d0)
           CUBE2(j)=CUBE(j)/2.0D0
           do i=1,np
                r0(i,j) = rtemp(i,j)*(scale + 1.0d0)
           enddo
c      enddo

      call model
      write(6,*) 'scale,escale,tote: ',scale,escale,tote
      if(tote. gt.escale) then
           if(scale.eq.0.-0010d0) then
                scale = -scale
           else
                scale = -scale/2.0d0
                if(abs(scale).lt.0.0002d0) goto 79
           endif
           goto 78
      else
           escale=tote
           do j=1,3
           ckeep(j)=cube(j)
           do i=1,np
           rtemp(i,j) = r0(i,j)
           enddo
           enddo
           goto 78
      endif
79    continue
      do j=1,3
      cube(j)=ckeep(j)
      CUBE2(j)=CUBE(j)/2.0D0
      do i=1,np
      r0(i,j) = rtemp(i,j)
      enddo
      enddo
      kflag = kkeep
      maxkb = mkeep
c      write(42,*) cube(1)*cube(2)*cube(3) , escale
      call flush(42)
      return
      end
