subroutine ensightout_results(np, x, v, rho, p, u, &
     file_step, filedir) 

!!$------------------------------------------------------------------
!!$ Subroutine for saving particle information on position, velocity, 
!!$ density, pressure and enery to external disk file in the EnSight 
!!$ format
!!$
!!$ np-- total particle number                                    [in]
!!$ x-- coordinates of particles                                  [in]
!!$ v-- velocities of particles                                   [in]
!!$ rho-- dnesities of particles                                  [in]
!!$ p-- pressure  of particles                                    [in]
!!$ u-- internal energy of particles                              [in]
!!$ file_step-- step number for file output                       [in]
!!$ filedir -- Directory where the files are located              [in]

  implicit none     
  include 'options.inc'

  integer np, file_step
  double precision x(dim, maxn), v(dim, maxn), &
       rho(maxn), p(maxn), u(maxn)
  character (len = *) :: filedir

  integer i, d
  character (len = 120) file_name

  write(file_name,'(A,A,I5.5,A)') trim(filedir), &
       'CERN_position_', file_step, '.geo'
  open(10,file = file_name)
  write(10,'(A)') "SPH CERN output in EnSight Gold format"
  write(10,'(A)') "EnSight 8.0.7"
  write(10,'(A)') "node id assign"
  write(10,'(A)') "element id assign"
  write(10,'(A)') "extents"
  write(10,'(A)') " 1.00000e+38-1.00000e+38"
  write(10,'(A)') " 1.00000e+38-1.00000e+38"
  write(10,'(A)') " 1.00000e+38-1.00000e+38"
  write(10,'(A)') "part"
  write(10,'(I10)') 1
  write(10,'(A)') "SPH particles"
  write(10,'(A)') "coordinates"
  write(10,'(I10)') np

  do d = 1, dim
     do i = 1, np
        write(10,'(E15.8E3)') x(d,i)
     enddo
  enddo
  do d = dim +1, 3
     do i = 1, np
        write(10,'(E15.8E3)') 0.0d0
     enddo
  enddo
  close(10)

  write(file_name,'(A,A,I5.5,A)') trim(filedir), &
       'CERN_velocity_', file_step, '.dat'
  open(10,file = file_name)
  write(10,'(A)') "particle velocity in EnSight Gold format"
  write(10,'(A)') "part"
  write(10,'(I10)') 1
  write(10,'(A)') "coordinates"
  do d = 1, dim
     do i = 1, np
        write(10,'(E15.8E3)') v(d,i)
     enddo
  enddo
  do d = dim +1, 3
     do i = 1, np
        write(10,'(E15.8E3)') 0.0d0
     enddo
  enddo
  close(10)

  write(file_name,'(A,A,I5.5,A)') trim(filedir), 'CERN_density_', file_step, &
       '.dat'
  open(10,file = file_name)
  write(10,'(A)') "particle density in EnSight Gold format"
  write(10,'(A)') "part"
  write(10,'(I10)') 1
  write(10,'(A)') "coordinates"
  do i = 1, np
     write(10,'(E15.8E3)') rho(i)
  enddo
  close(10)

  write(file_name,'(A,A,I5.5,A)') trim(filedir), &
       'CERN_pressure_', file_step, '.dat'
  open(10,file = file_name)
  write(10,'(A)') "particle pressure in EnSight Gold format"
  write(10,'(A)') "part"
  write(10,'(I10)') 1
  write(10,'(A)') "coordinates"
  do i = 1, np
     write(10,'(E15.8E3)') p(i)
  enddo
  close(10)

  write(file_name,'(A,A,I5.5,A)') trim(filedir), 'CERN_energy_', file_step, &
       '.dat'
  open(10,file = file_name)
  write(10,'(A)') "particle energy in EnSight Gold format"
  write(10,'(A)') "part"
  write(10,'(I10)') 1
  write(10,'(A)') "coordinates"
  do i = 1, np
     write(10,'(E15.8E3)') u(i)
  enddo
  close(10)

end subroutine ensightout_results


subroutine ensightout_case(dt, steps, save_step, filedir) 

!!$-------------------------------------------------------------------
!!$ Subroutine for saving the master file of the results in the 
!!$ EnSight format to external disk file
!!$
!!$ dt-- time step                                                [in]
!!$ steps-- number of time steps                                  [in]
!!$ save_step-- number of time steps between saves                [in]
!!$ filedir -- Directory where the files are located              [in]

  implicit none     
  include 'options.inc'

  double precision dt, time
  integer steps, save_step
  character (len = *) :: filedir

  integer i, its
  character (len = 120) file_name

  write(file_name,'(A,A)')  trim(filedir), 'CERN_SPH.case'
  open(10,file = file_name)
  write(10,'(A)') "# Ensight formatted case file for SPH CERN"
  write(10,'(A)') ""
  write(10,'(A)') "FORMAT"
  write(10,'(A)') "type: ensight gold"
  write(10,'(A)') ""
  write(10,'(A)') "GEOMETRY"
  write(10,'(A,A,A)') "model:    1          ", "CERN", &
       "_position_*****.geo"
  write(10,'(A)') ""
  write(10,'(A)') "VARIABLE"
  write(10,'(A,A,A)') 'vector per node:    1 velocity ', 'CERN', &
       '_velocity_*****.dat'
  write(10,'(A,A,A)') 'scalar per node:    1 density  ', 'CERN', &
       '_density_*****.dat'
  write(10,'(A,A,A)') 'scalar per node:    1 pressure ', 'CERN', &
       '_pressure_*****.dat'
  write(10,'(A,A,A)') 'scalar per node:    1 energy   ', 'CERN', &
       '_energy_*****.dat'
  write(10,'(A)') ""
  write(10,'(A)') "TIME"
  write(10,'(A,I10)') "time set:", 1
  write(10,'(A,I10)') "number of steps:", (steps / save_step) +1
  write(10,'(A,I10)') "filename start number:", 0
  write(10,'(A,I10)') "filename increment:", 1
  write(10,'(A)') "time values:"

  time = 0.
  write(10,'(E15.8E3)') time
  do its = 1, steps
     time = time + dt
     if (mod(its, save_step) == 0) write(10,'(E15.8E3)') time
  enddo

  close(10)

end subroutine ensightout_case
