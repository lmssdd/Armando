program Armando

!!$------------------------------------------------------------------
!!$ Armando is a SPH code. the followings are the 
!!$ basic parameters needed in this code or calculated by this code
!!$
!!$ mass-- mass of particles                                      [in]
!!$ np-- total particle number used                               [in]
!!$ dt--- Time step used in the time integration                  [in]
!!$ itype-- types of particles                                    [in]
!!$ x-- coordinates of particles                              [in/out]
!!$ v-- velocities of particles                               [in/out]
!!$ rho-- dnesities of particles                              [in/out]
!!$ p-- pressure  of particles                                [in/out]
!!$ u-- internal energy of particles                          [in/out]
!!$ h-- smoothing lengths of particles                        [in/out]
!!$ c-- sound velocity of particles                              [out]
!!$ s-- entropy of particles                                     [out]
!!$ e-- total energy of particles                                [out]
!!$	
!!$ sffluka-- size scale factor (size_sph = sf*size_fluka)       [in]
!!$ effluka-- energy scale factor (energy_sph = ef*energy_fluka) [in]
!!$ npfluka-- number of particles                                [in]
!!$ dtfluka-- deposition time                                    [in]
!!$ ofluka-- the origin of the binning                           [in]
!!$ dfluka-- cell side of the binning                            [in]
!!$ nfluka-- number of cells for each index                      [in]
!!$ tfluka-- identifier for the fluka grid type                  [in]
!!$ datafluka-- data of the binning                              [in]


  implicit none     
  include 'options.inc'

!!$c     SPH algorithm for particle approximation (pa_sph)
!!$c     pa_sph = 0 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
!!$c              1 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
!!$c
!!$c     Smoothing length evolution (sle) algorithm
!!$c     sle = 0 : Keep unchanged,
!!$c           1 : h = fac * (m/rho)^(1/dim)
!!$c           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
!!$c           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
!!$c	
!!$c     Smoothing kernel function
!!$c     skf = 0, cubic spline kernel by W4 - Spline (Monaghan 1985)
!!$c         = 1, Gauss kernel   (Gingold and Monaghan 1981)
!!$c         = 2, Quintic kernel (Morris 1997)
!!$c	
!!$c     Density calculation method
!!$c     density_method = 0, continuity equation
!!$c			       = 1, summation density
!!$c                    = 2, normalized summation density

  character (len = 120) :: filedir
  integer pa_sph, sle, skf, density_method


!!$c     Switches for different senarios
!!$	
!!$c     average_velocity = .TRUE. : Monaghan treatment on average velocity,
!!$c                        .FALSE.: No average treatment.
!!$c     virtual_part = .TRUE. : Use virtual particle,
!!$c                    .FALSE.: No use of virtual particle.
!!$c     ghost_part = .TRUE. : Use ghost particle,
!!$c                  .FALSE.: No use of ghost particle.
!!$c     visc = .true. : Consider viscosity,
!!$c            .false.: No viscosity.
!!$c     visc_artificial = .true. : Consider artificial viscosity,
!!$c                       .false.: No considering of artificial viscosity.
!!$c     heat_artificial = .true. : Consider artificial heating,
!!$c                       .false.: No considering of artificial heating.
!!$c     heat_external = .true. : Consider external heating,
!!$c                     .false.: No considering of external heating.
!!$c     self_gravity = .true. : Considering self_gravity,
!!$c                    .false.: No considering of self_gravity
  logical average_velocity, virtual_part, ghost_part, visc, &
       heat_artificial, heat_external, visc_artificial, gravity

!!$c     Simulation cases
!!$c     shocktube = .true. : carry out shock tube simulation
!!$c     shearcavity = .true. : carry out shear cavity simulation
  logical shocktube, shearcavity

  integer np, nv, steps, refresh_step, save_step, d, m, i
  double precision prop(16, 40), dt
  double precision trfluka(3), sffluka, effluka, npfluka, dtfluka, &
       ofluka(3), dfluka(3)
  integer nfluka(3), tfluka

  integer, dimension(:), allocatable :: mat
  double precision, dimension(:,:), allocatable :: x, v 
  double precision, dimension(:), allocatable :: h, mass, rho, p, &
       u, c, s, e 

  double precision, dimension(:), allocatable :: datafluka

  double precision s1, s2


  allocate(mat(maxn))
  allocate(x(dim, maxn), v(dim, maxn)) 
  allocate(h(maxn), mass(maxn), rho(maxn), p(maxn))
  allocate(u(maxn), c(maxn), s(maxn), e(maxn))
  allocate(datafluka(max_fluka))

  np = 0
  nv = 0
  
  call default(filedir, pa_sph, sle, skf, density_method, &
       average_velocity, virtual_part, ghost_part, visc, &
       visc_artificial, heat_artificial, heat_external, gravity, &
       refresh_step, save_step)

  shocktube = .false.
  shearcavity = .false.

  if (.not. (shearcavity .or. shocktube .or. debug)) then
     call read_input_file(filedir, pa_sph, sle, skf, &
          density_method, average_velocity, &
          virtual_part, ghost_part, visc, &
          visc_artificial, heat_artificial, &
          heat_external, gravity, &
          prop, dt, steps, refresh_step, save_step, &
          trfluka, sffluka, effluka, &
          npfluka, dtfluka)
  endif

  if (shocktube) then
     dt = 0.005
     steps = 40
     prop(1, 1) = 1.0
     prop(2, 1) = 1.4
     call shock_tube(np, x, v, mat, h, mass, rho, p, u)
     
     pa_sph = 1
     sle = 0
     skf = 0
     density_method = 1

     average_velocity  = .true.
     virtual_part  = .false.
     ghost_part  = .true.
     visc  = .false.
     visc_artificial = .true.
     heat_artificial  = .false.
     heat_external  = .false.
     gravity  = .false.

     refresh_step = 1 !100 /2
     save_step = 1000 !100 /2

  endif

  if (debug) then
!!$     dt = 2 * 4.19463E-7
!!$     steps = 2000 / 2
!!$     prop(1,21) = 1.354000e+004
!!$     prop(2,21) = 1.490000e+003
!!$     prop(3,21) = 1.960000
!!$     prop(4,21) = 2.047000
!!$     prop(5,21) = -1.0E30
!!$     prop(6,21) = 0.0d-3
!!$     call partly_heated_rod(np, x, v, mat, h, mass, rho, p, u)
!!$
!!$     pa_sph = 1
!!$     sle = 0
!!$     skf = 0
!!$     density_method = 0
!!$
!!$     average_velocity  = .true.
!!$     virtual_part  = .false.
!!$     ghost_part  = .false.
!!$     visc  = .false.
!!$     visc_artificial  = .true.
!!$     heat_artificial  = .true.
!!$     heat_external  = .false.
!!$     gravity  = .false.
!!$
!!$     refresh_step = steps !100 /2
!!$     save_step = steps !100 /2


!!$     dt = 0.5 * 1.0e-5
!!$     steps = 5000
!!$     prop(1,21) = 1.0e+003
!!$     prop(2,21) = 1.0e+003
!!$     prop(3,21) = 0.0
!!$     prop(4,21) = 0.0
!!$     prop(5,21) = -1.0E30
!!$     prop(6,21) = 0.0d-3
!!$     call bullet(np, nv, x, v, mat, h, mass, rho, p, u)
!!$
!!$     pa_sph = 0
!!$     sle = 0
!!$     skf = 1
!!$     density_method = 0
!!$
!!$     average_velocity  = .false.
!!$     virtual_part  = .true.
!!$     ghost_part  = .false.
!!$     visc  = .false.
!!$     visc_artificial = .true.
!!$     heat_artificial  = .false.
!!$     heat_external  = .false.
!!$     gravity  = .false.
!!$
!!$     refresh_step = 1 !100 /2
!!$     save_step = 1000 !100 /2



!!$     dt = 2.0E-4
!!$     steps = 10000
!!$     
!!$     prop(1,21) = 1.0e+003
!!$     prop(2,21) = 1.0e+001
!!$     prop(3,21) = 0.0
!!$     prop(4,21) = 0.0
!!$     prop(5,21) = -1.0E30
!!$     prop(6,21) = 0.0d-3
!!$     
!!$     prop(1,41) = 1.0e+003
!!$     prop(2,41) = 1.0e+001
!!$     prop(3,41) = -1.0E30
!!$     prop(4,41) = 0.0d-3
!!$     call dam(np, nv, x, v, mat, h, mass, rho, p, u)
!!$
!!$     sle = 0
!!$     skf = 2
!!$     density_method = 0
!!$
!!$     average_velocity  = .true.
!!$     virtual_part  = .true.
!!$     visc  = .false.
!!$     visc_artificial = .true.
!!$     heat_artificial  = .false.
!!$     heat_external  = .false.
!!$     gravity  = .true.
!!$
!!$     refresh_step = 1
!!$     save_step = 500


!!$     dt = 5.0E-8
!!$     steps = 1600
!!$     prop(1,21) = 1.354000e+004
!!$     prop(2,21) = 1.490000e+003
!!$     prop(3,21) = 1.960000
!!$     prop(4,21) = 2.047000
!!$     prop(5,21) = -1.0E30
!!$     prop(6,21) = 0.0d-3
!!$     call partly_heated_disc(np, nv, x, v, mat, h, mass, rho, p, u)
!!$
!!$     pa_sph = 0
!!$     sle = 0
!!$     skf = 0
!!$     density_method = 0
!!$
!!$     average_velocity  = .true.
!!$     virtual_part  = .true.
!!$     ghost_part  = .true.
!!$     visc  = .false.
!!$     visc_artificial  = .true.
!!$     heat_artificial  = .false.
!!$     heat_external  = .false.
!!$     gravity  = .false.
!!$
!!$     refresh_step = 100
!!$     save_step = 100
!!$
!!$
!!$
!!$     dt = 10.0d-9 !5 * 2.0E-8
!!$     steps = 14 !1000 /5
!!$
!!$     prop(1,31) = 1549.000
!!$     prop(2,31) = 1.317e+010
!!$     prop(3,31) = 1.537e+010
!!$     prop(4,31) = 0.0
!!$     prop(5,31) = 0.300000
!!$     prop(6,31) = 0.250000
!!$     prop(7,31) = 1.250000e+006
!!$     prop(8,31) = 1.317e+010
!!$     prop(9,31) = 1.537e+010
!!$     prop(10,31) = -1.50000e+008
!!$     prop(11,31) = 0.0d-3
!!$
!!$     call fluka_heated_disc(np, nv, x, v, mat, h, mass, rho, p, u)
!!$
!!$     pa_sph = 1
!!$     sle = 0
!!$     skf = 0
!!$     density_method = 0
!!$
!!$     average_velocity  = .true. 
!!$     virtual_part  = .false.
!!$     ghost_part  = .true.
!!$     visc  = .false.
!!$     visc_artificial  = .true.
!!$     heat_artificial  = .false.
!!$     heat_external  = .true.
!!$     gravity  = .false.
!!$
!!$     refresh_step = 1! 50 /5
!!$     save_step = 50 /5
!!$
!!$     sffluka = 1.0d-2
!!$     effluka = 1.602d-10
!!$     npfluka = 3.03d13
!!$     dtfluka = 140d-9
!!$     trfluka(3) = 15.d-2 !10d-2
!!$     !trfluka(1) = 1.8d-2
!!$
!!$
!!$
!!$     dt = 10.0E-8
!!$     steps = 400
!!$     prop(1,21) = 1.354000e+004
!!$     prop(2,21) = 1.490000e+003
!!$     prop(3,21) = 1.960000
!!$     prop(4,21) = 2.047000
!!$     prop(5,21) = -1.0E30
!!$     prop(6,21) = 0.0d-3
!!$     call partly_heated_bar(np, nv, x, v, mat, h, mass, rho, p, u)
!!$
!!$     pa_sph = 1
!!$     sle = 0
!!$     skf = 0
!!$     density_method = 0
!!$
!!$     average_velocity  = .false.
!!$     virtual_part  = .false.
!!$     ghost_part  = .false.
!!$     visc  = .false.
!!$     visc_artificial  = .false.
!!$     heat_artificial  = .false.
!!$     heat_external  = .false.
!!$     gravity  = .false.
!!$
!!$     refresh_step = 1
!!$     save_step = 50
!!$
!!$
!!$     dt = 2.0E-7
!!$     steps = 2000
!!$     prop(1,21) = 1.354000e+004
!!$     prop(2,21) = 1.490000e+003
!!$     prop(3,21) = 1.960000
!!$     prop(4,21) = 2.047000
!!$     prop(5,21) = -1.5E5
!!$     prop(6,21) = 0.0d-3
!!$     call partly_heated_cylinder(np, nv, x, v, mat, h, mass, rho, p, u)
!!$
!!$     pa_sph = 0
!!$     sle = 0
!!$     skf = 0
!!$     density_method = 0
!!$
!!$     average_velocity  = .false.
!!$     virtual_part  = .false.
!!$     ghost_part  = .false.
!!$     visc  = .false.
!!$     visc_artificial  = .true.
!!$     heat_artificial  = .false.
!!$     heat_external  = .false.
!!$     gravity  = .false.
!!$     refresh_step = 1
!!$     save_step = 100

     dt = 100.0E-9
     steps = 10000
     prop(1,21) = 1.354000e+004
     prop(2,21) = 1.490000e+003
     prop(3,21) = 1.960000
     prop(4,21) = 2.047000
     prop(5,21) = -1.5E5
     prop(6,21) = 0.0d-3
     call mercury_pot3(np, nv, x, v, mat, h, mass, rho, p, u)

     sle = 0
     skf = 2
     density_method = 0

     average_velocity  = .false.
     virtual_part  = .true.
     visc  = .false.
     visc_artificial  = .true.
     heat_artificial  = .false.
     heat_external  = .false.
     gravity  = .true.
     
     refresh_step = 50
     save_step = 250

  endif

  if (shearcavity) then
     dt = 5.0e-5
     steps = 3000
     prop(1,11) = 1000.0
     prop(2,11) = 1.0d-1
     prop(3,11) = 0.

     prop(4,11) = 0.
     prop(5,11) = 0.
     prop(6,11) = 0.
     prop(7,11) = prop(2,11)
     prop(8,11) = 0.0
     prop(9,11) = -1.0E30
     prop(10,11) = 1.0d-3
     call shear_cavity(np, nv, x, v, mat, h, mass, rho, p, u)
     density_method = 1
     refresh_step = 10
     save_step = 500
  endif

  if (.not. (shearcavity .or. shocktube .or. debug)) then
     call read_initial_conditions(np, x, v, mat, h, &
          mass, rho, p, u, filedir)
     if (virtual_part .or. ghost_part) then
        call read_virtual_particles(np, nv, x, v, mat, h, &
             mass, rho, p, u, filedir)
     endif
  endif

  if (heat_external) then
     call read_binning_bin(trfluka, sffluka, effluka, npfluka, &
          dtfluka, ofluka, dfluka, nfluka, tfluka, &
          datafluka, filedir)
  endif

  call smooth_boundary(np, nv, x, v, h, mass, rho, skf)
  
  call time_print
  call time_elapsed(s1)

  call time_integration(np, nv, x, v, mat, h, &
       mass, rho, p, u, prop, &
       dt, steps, refresh_step, save_step, &
       dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, &
       sle, skf, &
       average_velocity, virtual_part, visc, &
       visc_artificial, heat_artificial, &
       heat_external, gravity, filedir)

  call time_print
  call time_elapsed(s2)      
  write (*,*)'        Elapsed CPU time = ', s2-s1
  write (15,*)'        Elapsed CPU time = ', s2-s1

  deallocate(mat)
  deallocate(x, v) 
  deallocate(h, mass, rho, p)
  deallocate(u, c, s, e)
  deallocate(datafluka)

end program Armando
