c---------------------------------------------------------
c     Including file for parameters and constants used 
c     in the entire SPH software packages.
c---------------------------------------------------------

	
c     maxn    : Maximum number of particles
c     max_interation : Maximum number of interaction pairs
      integer maxn,max_interaction,max_grid_data
      integer dim
      integer max_fluka
	
c     Control parameters for output 
c     int_stat = .true. : Print statistics about SPH particle interactions.
c                        including virtual particle information.
c     moni_particle: The particle number for information monitoring.
      logical int_stat
      integer moni_particle
	
      double precision pi
      parameter (pi = 3.14159265358979323846)

      parameter (maxn = 600000, max_interaction = maxn*256, 
     &           max_grid_data = 2*maxn)

      parameter (dim = 3)
	  parameter (max_fluka = 12800000)

      logical debug
      parameter (debug = .false.)
