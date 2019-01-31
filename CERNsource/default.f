	subroutine default(filedir, pa_sph, nnps, sle, skf, 
     &	density_method, average_velocity, virtual_part, visc, 
     &	visc_artificial, heat_artificial, heat_external, self_gravity, 
     &	refresh_step, save_step)

c----------------------------------------------------------------------
c     This subroutine is used to define the default values of the 
c     common options
c	
c     filedir -- Directory where the files are located             [out]
c     pa_sph -- SPH algorithm for particle approximation           [out]
c     nnps -- Nearest neighbor particle searching method           [out]
c     sle -- Smoothing length evolution                            [out]
c     skf -- Smoothing kernel function                             [out]
c     density_method -- Density calculation method                 [out]
c	average_velocity                                             [out]
c	virtual_part                                                 [out]
c	visc                                                         [out]
c	visc_artificial                                              [out]
c	heat_artificial                                              [out]
c	heat_external                                                [out]
c	self_gravity                                                 [out]
c	refresh_step                                                 [out]
c	save_step                                                    [out]
	
	implicit none
	include 'options.inc'

	character (len = *) :: filedir
	integer pa_sph, nnps, sle, skf, density_method
	logical average_velocity, virtual_part, visc, 
     &	visc_artificial, heat_artificial, heat_external, self_gravity
	integer refresh_step, save_step
	
	filedir = "./"

c     SPH algorithm for particle approximation (pa_sph)
c     pa_sph = 0 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
c              1 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
      pa_sph = 1

c     Nearest neighbor particle searching (nnps) method
c     nnps = 0 : Simplest and direct searching
c            1 : Sorting grid linked list
c            2 : Tree algorithm
      nnps = 1

c     Smoothing length evolution (sle) algorithm
c     sle = 0 : Keep unchanged,
c           1 : h = fac * (m/rho)^(1/dim)
c           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
c           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
      sle = 0

c     Smoothing kernel function
c     skf = 0, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 1, Gauss kernel   (Gingold and Monaghan 1981)
c         = 2, Quintic kernel (Morris 1997)
      skf = 0

c     Density calculation method
c     density_method = 0, continuity equation
c                    = 1, summation density
c                    = 2, normalized summation density
      density_method = 0
	
c     Switches for different senarios
	
c     average_velocity = .TRUE. : Monaghan treatment on average velocity,
c                        .FALSE.: No average treatment.
c     virtual_part = .TRUE. : Use virtual particle,
c                    .FALSE.: No use of virtual particle.
c     visc = .true. : Consider viscosity,
c            .false.: No viscosity.
c     visc_artificial = .true. : Consider artificial viscosity,
c                       .false.: No considering of artificial viscosity.
c     heat_artificial = .true. : Consider artificial heating,
c                       .false.: No considering of artificial heating.
c     heat_external = .true. : Consider external heating,
c                     .false.: No considering of external heating.
c     self_gravity = .true. : Considering self_gravity,
c                    .false.: No considering of self_gravity
	
      average_velocity  = .true.
      virtual_part  = .true.
      visc  = .true.
      visc_artificial  = .false.
      heat_artificial  = .false.
      heat_external  = .false. ! lucam
      self_gravity  = .false.

      refresh_step = 1

c     Control parameters for output 
c     save_step: Save Timestep (On Disc)
      save_step = 100
	
	end
