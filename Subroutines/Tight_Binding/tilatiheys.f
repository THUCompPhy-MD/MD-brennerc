       subroutine tilatiheys
C      density of states=tilatiheys

       implicit real*8(a-h,o-z)
      include 'common_files.inc'


       write(6,*) 'nmax,nelen,Emax,Emintb: ',nmax,nelen,Emax,Emintb
       do 10 i=1,nmax
          l = int(NElen*(E(i)-Emintb)/(Emax-Emintb))
          if ((l .gt. 0) .and. (l .le. NElen)) then
             idos(l)=idos(l)+1
          else
             NE_yli=NE_yli+1
          endif
10     continue

       return
       end

