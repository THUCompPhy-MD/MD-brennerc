c source code for these subroutines is available in
c routine calling potentials/forces - application specific 
      include 'model.f' 
c bond order potential subroutines 
      include 'caguts.f'
      include 'bcuint.f'
      include 'mtable.f'
      include 'param.f'
      include 'pibond.f'
      include 'radic.f'
      include 'sili_germ.f'
      include 'tor.f'
c tight binding 
      include 'E_corr.f' 
      include 'dsmooth.f'
      include 'force_tb.f'
      include 'laske_H.f'
      include 'neighbor_tb.f'
      include 'smooth.f'
      include 'tight_bind.f'
      include 'tilatiheys.f'
      include 'xkorjaus.f'
      include 'eispack.f'
      include 'dos.f'
c Lennard-Jones routines 
      include 'ljcont.f'
      include 'ljcset.f'
      include 'ljguts.f'
      include 'ljparam.f'
c general MD routines 
      include 'read_write.f' 
      include 'setin.f' 
      include 'overwr.f'
      include 'pred_corr.f'
      include 'reflect.f'
      include 'thermostats.f'
      include 'random_numbers.f'
c generate data 
      include 'xmol.f'
c energy minimizer
      include 'minimize.f'
c volume scaling 
c      include '../Subroutines/General/vscale.f'
c specialized routines
c
c    indenting
c 
c      include '../Subroutines/Special/move_tip.f'
c      include '../Subroutines/Special/load.f'
c
c    adding energetic particles
c 
c      include '../Subroutines/Special/bombard.f' 
c      include '../Subroutines/Special/reset_itr.f' 
c      include '../Subroutines/Special/max_ke.f' 
