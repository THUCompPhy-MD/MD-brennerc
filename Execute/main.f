C##############################################################C
C                                                              C
C                                                              C
C  THESE PROGRAMS PERFORM MOLECULAR DYNAMICS SIMULATIONS WITH  C
C  BOND-ORDER POTENTIALS FOR HYDROCARBON, SILICON AND          C
C  GERMANIUM; TIGHT BINDING FOR CARBON; AND LENNARD-JONES      C
C  WITH INPUTTED PARAMETERS. THE LATTER POTENTIALS ARE         C 
C  TRUNCATED AT SMALL DISTANCES FOR C-C, H-H, AND C-H PAIRS    C 
C                                                              C
C   Units: mass = AMU's, length = Angstroms, time = fs         C
C          energy = eV's                                       C
C                                                              C
C   DOCUMENTATION CAN BE FOUND IN:                             C
C      /MD/Documentation                                       C 
C                                                              C
C##############################################################C
c
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'
c 
c************************* 
c set up and input data  * 
c*************************
c 
c open input/output files 
      include 'open.inc' 
c initialize 
      call setin 
C read input data 
      call read_data 
C setup potential parameters  
      call setpp 
c setup predictor-corrector coefficients
      call setpc 
C setup Langevin parameters  
      call setgle 
C initialize random number generator 
      call setran 
c write out data headers 
      call write_data1 
C
c********************
c begin calculation *
c********************

      if(kflag.eq.6) then

c**********************
C Energy minimization *
c**********************
           kvc = 1
           call minimize
           write(6,*) 'minimum energy= ',tote 
      else 
c**********************
C Dynamic Simulation  *
c**********************
            call setmd 
            DO LSTEP=1,KVC

c predictor 
                 call cpred 
c calculate energy and forces 
                 CALL MODEL 
c apply thermostats 
                 call thermos 
c corrector 
                 call ccorr

c write out position file to be post converted to xmol format
                 if(mod(LSTEP,nxmol).eq.0) call xmol
C
                 IF(KFLAG.EQ.5) CALL BERE
c
                 if(mod(LSTEP,maxkb).eq.0) then
c generate and write data 
                 call write_data2 
                 endif 
c volume scaling 
c               call vscale 
           ENDDO 
c
      ENDIF

C
C WRITE OUT POSITIONS FOR RESTART
C
C *****IMPORTANT***********************
C                                     *
C !!! INPUT FILE IS WRITTEN OVER !!!  *
C                                     * 
C**************************************
C
      call write_data3 
C
      include 'close.inc'
      STOP
      END
C
C add included subroutines
C
      include 'subroutines.inc' 
