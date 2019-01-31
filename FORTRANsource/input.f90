subroutine read_input_file(filedir, pa_sph, sle, skf, &
     density_method, average_velocity, &
     virtual_part, ghost_part, visc, &
     visc_artificial, heat_artificial, &
     heat_external, self_gravity, &
     prop, dt, steps, refresh_step, &
     save_step, trfluka, sffluka, effluka, &
     npfluka, dtfluka)

!!$------------------------------------------------------------------
!!$ Subroutine to read the input file, and set the options for the
!!$ simulation
!!$	
!!$ filedir -- Directory where the files are located            [out]
!!$ pa_sph -- SPH algorithm for particle approximation          [out]
!!$ sle -- Smoothing kernel function                            [out]
!!$ skf -- Smoothing kernel function                            [out]
!!$ density_method -- Density calculation method                [out]
!!$ average_velocity                                            [out]
!!$ virtual_part                                                [out]
!!$ ghost_part                                                  [out]
!!$ visc                                                        [out]
!!$ visc_artificial                                             [out]
!!$ heat_artificial                                             [out]
!!$ heat_external                                               [out]
!!$ self_gravity                                                [out]
!!$ prop-- material properties                                  [out]
!!$ dt-- time step                                              [out]
!!$ steps-- number of time steps                                [out]
!!$ refresh_step-- number of steps between connectivity refresh [out]
!!$ save_step-- number of steps between subsequent saves        [out]
!!$ trfluka-- translation origin of the fluka grid              [out]
!!$ sffluka-- size scale factor (size_sph = sf*size_fluka)      [out]
!!$ effluka-- energy scale factor (energy_sph = ef*energy_fluka)[out]
!!$ npfluka-- number of particles                               [out]
!!$ dtfluka-- deposition time                                   [out]

  implicit none
  include 'options.inc'

  character (len = *) :: filedir
  integer pa_sph, sle, skf, density_method
  logical average_velocity, visc_artificial, heat_artificial
  logical virtual_part, ghost_part, visc, heat_external, self_gravity
  double precision prop(16,40), dt
  integer steps, refresh_step, save_step
  double precision trfluka(3), sffluka, effluka, npfluka, dtfluka

  integer status, intvalue
  character (len = 120) filename
  character (len = 120) line
  double precision doublevalue
  integer i
  integer im, ip1, ip2

  status = 0
  filename = trim(filedir) // "input.sph"
  open(1,file = filename, IOSTAT=status)
  if (status .ne. 0) then
     write(*,'(A)') 'ERROR OPENING input.sph'
  endif

  do
     read(1,'(A)', IOSTAT=status) line

     if (status .ne. 0) exit

     if (line(1:6) == 'PA_SPH') then 
        read(line(7:),*) intvalue
        if ((intvalue >= 0) .and. (intvalue <2)) pa_sph = intvalue
     endif

     if (line(1:3) == 'SLE') then 
        read(line(4:),*) intvalue
        if ((intvalue >= 0) .and. (intvalue <2)) sle = intvalue
     endif

     if (line(1:3) == 'SKF') then 
        read(line(4:),*) intvalue
        if ((intvalue >= 0) .and. (intvalue <2)) skf = intvalue
     endif

     if (line(1:11) == 'DENS_METHOD') then 
        read(line(12:),*) intvalue
        if ((intvalue >= 0) .and. (intvalue <2)) &
             density_method = intvalue
     endif

     if (line(1:7) == 'AVG_VEL') then 
        read(line(8:),*) intvalue
        if (intvalue == 0) average_velocity = .FALSE.
        if (intvalue == 1) average_velocity = .TRUE.
     endif

     if (line(1:6) == 'V_PART') then 
        read(line(7:),*) intvalue
        if (intvalue == 0) virtual_part = .FALSE.
        if (intvalue == 1) virtual_part = .TRUE.
     endif

     if (line(1:6) == 'G_PART') then 
        read(line(7:),*) intvalue
        if (intvalue == 0) ghost_part = .FALSE.
        if (intvalue == 1) ghost_part = .TRUE.
     endif
     
     if (line(1:4) == 'VISC') then 
        read(line(5:),*) intvalue
        if (intvalue == 0) visc = .FALSE.
        if (intvalue == 1) visc = .TRUE.
     endif

     if (line(1:8) == 'ART_VISC') then 
        read(line(9:),*) intvalue
        if (intvalue == 0) visc_artificial = .FALSE.
        if (intvalue == 1) visc_artificial = .TRUE.
     endif

     if (line(1:8) == 'ART_HEAT') then 
        read(line(9:),*) intvalue
        if (intvalue == 0) heat_artificial = .FALSE.
        if (intvalue == 1) heat_artificial = .TRUE.
     endif

     if (line(1:8) == 'EXT_HEAT') then 
        read(line(9:),*) intvalue
        if (intvalue == 0) heat_external = .FALSE.
        if (intvalue == 1) heat_external = .TRUE.
     endif

     if (line(1:7) == 'GRAVITY') then 
        read(line(8:),*) intvalue
        if (intvalue == 0) self_gravity = .FALSE.
        if (intvalue == 1) self_gravity = .TRUE.
     endif

     if (line(1:5) == 'STEPS') then 
        read(line(6:),*) intvalue
        if (intvalue > 0) steps = intvalue
     endif

     if (line(1:7) == 'REFRESH') then 
        read(line(8:),*) intvalue
        if (intvalue > 0) refresh_step = intvalue
     endif

     if (line(1:8) == 'SAVESTEP') then 
        read(line(9:),*) intvalue
        if (intvalue > 0) save_step = intvalue
     endif

     if (line(1:8) == 'TIMESTEP') then 
        read(line(9:),*) doublevalue
        if (doublevalue > 0.0) dt = doublevalue
     endif

     if (line(1:7) == 'TRFLUKA') then 
        read(line(8:),*) (trfluka(i), i = 1,3)
     endif

     if (line(1:7) == 'SFFLUKA') then 
        read(line(8:),*) doublevalue
        if (doublevalue > 0.0) sffluka = doublevalue
     endif

     if (line(1:7) == 'EFFLUKA') then 
        read(line(8:),*) doublevalue
        if (doublevalue > 0.0) effluka = doublevalue
     endif

     if (line(1:7) == 'NPFLUKA') then 
        read(line(8:),*) doublevalue
        if (doublevalue > 0.0) npfluka = doublevalue
     endif

     if (line(1:7) == 'DTFLUKA') then 
        read(line(8:),*) doublevalue
        if (doublevalue > 0.0) dtfluka = doublevalue
     endif

     if (line(1:3) == 'MAT') then 
        read(line(4:),*) im, ip1, ip2
        if (((im .ge. 1) .and. (im .le. 40)) .and. &
             ((ip1 .ge. 1) .and. (ip1 .le. 16)) .and. &
             ((ip2 .gt. ip1) .and. (ip2 .le. 16))) then
           read(line(4:),*) &
                intvalue, intvalue, intvalue, &
                (prop(i, im), i = ip1, ip2)
        endif
     endif

  enddo

  close(1)

end subroutine read_input_file

subroutine read_initial_conditions(np, x, v, mat, h, &
     mass, rho, p, u, filedir)

!!$-------------------------------------------------------------------
!!$ Subroutine for loading initial status of the real particles
!!$
!!$ np-- total particle number                                   [out]
!!$ x-- coordinates of particles                                 [out]
!!$ v-- velocities of particles                                  [out]
!!$ mat-- particle material                                      [out]
!!$ h-- smoothing lengths of particles                           [out]
!!$ mass-- mass of particles                                     [out]
!!$ rho-- dnesities of particles                                 [out]
!!$ p-- pressure  of particles                                   [out]
!!$ u-- internal energy of particles                             [out]
!!$ filedir -- Directory where the files are located              [in]

  implicit none
  include 'options.inc'

  integer np, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn) ,p(maxn), u(maxn)
  character (len = *) :: filedir

  integer i, d, im
  character (len = 120) filename

!!$ load initial particle information from external disk file

  write(*,*)'  **************************************************'
  write(*,*)'      Loading initial particle configuration...   '

  filename = trim(filedir) // "f_xv.dat"
  open(1,file = filename)

  filename = trim(filedir) // "f_state.dat"
  open(2,file = filename)

  filename = trim(filedir) // "f_other.dat"
  open(3,file = filename)

  read (1,*) np

  do i = 1, np
     read(1,*)im, (x(d, i),d = 1, dim), (v(d, i),d = 1, dim)
     read(2,*)im, mass(i), rho(i), p(i), u(i)
     read(3,*)im, mat(i), h(i)
  enddo

  close(1)
  close(2) 
  close(3) 

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine read_initial_conditions

subroutine read_virtual_particles(np, nv, x, v, mat, h, &
     mass, rho, p, u, filedir)

!!$-------------------------------------------------------------------
!!$ Subroutine for loading the data realtive to the boundary particles
!!$
!!$ nv-- virtual particle number                                 [out]
!!$ x-- coordinates of particles                                 [out]
!!$ v-- velocities of particles                                  [out]
!!$ mat-- particle material                                      [out]
!!$ h-- smoothing lengths of particles                           [out]
!!$ mass-- mass of particles                                     [out]
!!$ rho-- dnesities of particles                                 [out]
!!$ p-- pressure  of particles                                   [out]
!!$ u-- internal energy of particles                             [out]
!!$ filedir -- Directory where the files are located              [in]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn) ,p(maxn), u(maxn)
  character (len = *) :: filedir

  integer i, d, im
  character (len = 120) filename

!!$ load initial particle information from external disk file

  write(*,*)'  **************************************************'
  write(*,*)'      Loading initial particle configuration...   '

  filename = trim(filedir) // "v_xv.dat"
  open(1,file = filename)

  filename = trim(filedir) // "v_state.dat"
  open(2,file = filename)

  filename = trim(filedir) // "v_other.dat"
  open(3,file = filename)

  read (1,*) nv

  do i = np +1, np +nv
     read(1,*)im, (x(d, i),d = 1, dim), (v(d, i),d = 1, dim)
     read(2,*)im, mass(i), rho(i), p(i), u(i)
     read(3,*)im, mat(i), h(i)
  enddo

  close(1)
  close(2) 
  close(3) 

  write(*,*)'      Total number of virtual particles   ', nv
  write(*,*)'  **************************************************'

end subroutine read_virtual_particles

