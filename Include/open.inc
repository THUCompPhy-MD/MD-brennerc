c
c standard input/output files 
      open(13,file='input.d',status='old') 
      open(11,file='coord.d',status='old') 
      open(9,file='output.d',status='unknown') 
c 
c data files for REBO potential 
C 
      open(14,file='../Spline/inter3d_iv_new.d',
     &status='old')
      open(15,file='../Spline/inter2d_iv.d',
     &status='old')
      open(16,file='../Spline/inter3dtors.d',
     &status='old')
      open(17,file='../Spline/inter3d_h.d',
     &status='old')
      open(18,file='../Spline/inter3d_ch.d',
     &status='old')
c 
c optional input/output files
c 
      open(55,file='max_ke.d',status='unknown') 
      open(85,file='pair_energy.d',status='unknown')
      open(50,file='overwrite.d',status='unknown')
      open(51,file='load.d',status='unkown')
      open(1,file='xmol.d',status='unknown') 
C
C output files for tight binding potential 
c 
       open(85,file='eigenvectors.d',status='unknown')
       open(86,file='eigenenergies.d',status='unknown')
       open(87,file='dos.d',status='unknown')
       open(88,file='ldos.d',status='unknown')

