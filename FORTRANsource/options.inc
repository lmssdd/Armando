!!$---------------------------------------------------------
!!$     Including file for parameters and constants used 
!!$     in the entire SPH software packages.
!!$---------------------------------------------------------

	
!!$     maxn    : Maximum number of particles
!!$     max_interation : Maximum number of interaction pairs
      integer maxn,max_interaction,max_grid_data
      integer dim
      integer max_fluka
	
!!$     Control parameters for output 
!!$     int_stat = .true. : Print statistics about SPH particle interactions.
!!$                        including virtual particle information.
!!$     moni_particle: The particle number for information monitoring.
      logical int_stat
      integer moni_particle
	
      double precision pi
      parameter (pi = 3.14159265358979323846)
      
      double precision hdr
      parameter (hdr = 1.2d0)

      parameter (maxn = 200000, max_interaction = maxn*256, &
      max_grid_data = 2*maxn)

      parameter (dim = 3)
      parameter (max_fluka = 12800000)

      logical debug
      parameter (debug = .true.)
