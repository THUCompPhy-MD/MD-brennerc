                         function xkorjaus(x)
C   --------------------------------------------------------------------------
C       Repulsive pair-potential
C
C        E_rep = A*(ro/x)^m*exp(m(aa-(x/rctb)^nc))
C

       implicit real*8(a-h,o-z)
      include 'common_files.inc'


      if(x.lt.2.60d0) then
           if(x.lt.2.57d0) then
                xkorjaus=Am*((do/x)**xmtb)*exp(xmtb*(aa1-(x/dc)**xmc))
           else
                xx = x - 2.57d0
                xkorjaus = stb0 + xx*(stb1 + xx*(stb2 + xx*stb3))
           endif
      else
           xkorjaus = 0.0d0
      endif
C
      return
      end
                     function Dkorj(x)
C   --------------------------------------------------------------------------
C     Derivative of the repulsive pair-potential
C
C       D E_rep = - A*(ro/x)^m*exp(m(aa-(x/rctb)^nc))*m/x*(1+nc(x/rctb)^nc)
C

       implicit real*8(a-h,o-z)

      include 'common_files.inc'

      if(x.lt.2.60d0) then
           if(x.lt.2.57d0) then
                Dkorj=-Am*((do/x)**xmtb)*xmtb/x*(1.0+xmc*(x/dc)**xmc)
     .               *exp(xmtb*(aa1 - (x/dc)**xmc))
           else
                xx = x - 2.57d0
                Dkorj = stb1 + xx*(2.0d0*stb2 + xx*3.0d0*stb3)
           endif
      else
           Dkorj = 0.0d0
      endif

c      Dkorj = -Dkorj
      return
      end

