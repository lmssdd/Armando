subroutine default(filedir, pa_sph, sle, skf, &
     density_method, average_velocity, virtual_part, ghost_part, visc, &
     visc_artificial, heat_artificial, heat_external, self_gravity, &
     refresh_step, save_step)

!!$------------------------------------------------------------------
!!$ This subroutine is used to define the default values of the 
!!$ common options
!!$ 	
!!$ filedir -- Directory where the files are located             [out]
!!$ pa_sph -- SPH algorithm for particle approximation           [out]
!!$ sle -- Smoothing length evolution                            [out]
!!$ skf -- Smoothing kernel function                             [out]
!!$ density_method -- Density calculation method                 [out]
!!$ average_velocity                                             [out]
!!$ virtual_part                                                 [out]
!!$ visc                                                         [out]
!!$ visc_artificial                                              [out]
!!$ heat_artificial                                              [out]
!!$ heat_external                                                [out]
!!$ self_gravity                                                 [out]
!!$ refresh_step                                                 [out]
!!$ save_step                                                    [out]

  implicit none
  include 'options.inc'

  character (len = *) :: filedir
  integer pa_sph, sle, skf, density_method
  logical average_velocity, virtual_part, ghost_part, visc, &
       visc_artificial, heat_artificial, heat_external, self_gravity
  integer refresh_step, save_step

  filedir = "./"

!!$ SPH algorithm for particle approximation (pa_sph)
!!$ pa_sph = 0 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
!!$          1 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
  pa_sph = 1

!!$ Smoothing length evolution (sle) algorithm
!!$ sle = 0 : Keep unchanged,
!!$       1 : h = fac * (m/rho)^(1/dim)
!!$       2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
!!$       3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
  sle = 0

!!$ Smoothing kernel function
!!$ skf = 0, cubic spline kernel by W4 - Spline (Monaghan 1985)
!!$     = 1, Gauss kernel   (Gingold and Monaghan 1981)
!!$     = 2, Quintic kernel (Morris 1997)
  skf = 0

!!$ Density calculation method
!!$ density_method = 0, continuity equation
!!$                = 1, summation density
!!$                = 2, normalized summation density
  density_method = 0

!!$ Switches for different senarios
!!$	
!!$ average_velocity = .TRUE. : Monaghan treatment on average velocity,
!!$                    .FALSE.: No average treatment.
!!$ virtual_part = .TRUE. : Use virtual particle,
!!$                .FALSE.: No use of virtual particle.
!!$ ghost_part = .TRUE. : Use ghost particle,
!!$              .FALSE.: No use of ghost particle.
!!$ visc = .true. : Consider viscosity,
!!$        .false.: No viscosity.
!!$ visc_artificial = .true. : Consider artificial viscosity,
!!$                   .false.: No considering of artificial viscosity.
!!$ heat_artificial = .true. : Consider artificial heating,
!!$                   .false.: No considering of artificial heating.
!!$ heat_external = .true. : Consider external heating,
!!$                 .false.: No considering of external heating.
!!$ self_gravity = .true. : Considering self_gravity,
!!$                .false.: No considering of self_gravity

  average_velocity  = .false.
  virtual_part  = .false.
  ghost_part  = .false.
  visc  = .false.
  visc_artificial  = .false.
  heat_artificial  = .false.
  heat_external  = .false. ! lucam
  self_gravity  = .false.

  refresh_step = 1

!!$ Control parameters for output 
!!$ save_step: Save Timestep (On Disc)
  save_step = 100

end subroutine default
