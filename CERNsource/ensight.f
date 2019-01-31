      subroutine ensightout_results(np, x, v, rho, p, u, 
     &                              file_step, filedir) 
      
c----------------------------------------------------------------------           
c     Subroutine for saving particle information on position, velocity, 
c     density, pressure and enery to external disk file in the EnSight 
c     format

c     np-- total particle number                                    [in]
c     x-- coordinates of particles                                  [in]
c     v-- velocities of particles                                   [in]
c     rho-- dnesities of particles                                  [in]
c     p-- pressure  of particles                                    [in]
c     u-- internal energy of particles                              [in]
c     file_step-- step number for file output                       [in]
c     filedir -- Directory where the files are located              [in]

      implicit none     
      include 'options.inc'
      
      integer np, file_step
      double precision x(dim, maxn), v(dim, maxn),
     &       rho(maxn), p(maxn), u(maxn)
	character (len = *) :: filedir

      integer i, d
	character (len = 120) file_name
      
	write(file_name,1504) trim(filedir), 
     &                      'CERN_position_', file_step, '.geo'
      open(10,file = file_name)
	write(10,1501) "SPH CERN output in EnSight Gold format"
	write(10,1501) "EnSight 8.0.7"
	write(10,1501) "node id assign"
	write(10,1501) "element id assign"
	write(10,1501) "extents"
	write(10,1501) " 1.00000e+38-1.00000e+38"
	write(10,1501) " 1.00000e+38-1.00000e+38"
	write(10,1501) " 1.00000e+38-1.00000e+38"
	write(10,1501) "part"
	write(10,1502) 1
	write(10,1501) "SPH particles"
	write(10,1501) "coordinates"
	write(10,1502) np
	
	do d = 1, dim
	  do i = 1, np
  	    write(10,1503) x(d,i)
	  enddo
	enddo
	do d = dim +1, 3
	  do i = 1, np
  	    write(10,1503) 0.0d0
	  enddo
	enddo
	close(10)
	
	write(file_name,1504) trim(filedir), 
     &                      'CERN_velocity_', file_step, '.dat'
      open(10,file = file_name)
	write(10,1501) "particle velocity in EnSight Gold format"
	write(10,1501) "part"
	write(10,1502) 1
	write(10,1501) "coordinates"
	do d = 1, dim
        do i = 1, np
	    write(10,1503) v(d,i)
	  enddo
	enddo
	do d = dim +1, 3
	  do i = 1, np
  	    write(10,1503) 0.0d0
	  enddo
	enddo
	close(10)
	
	write(file_name,1504) trim(filedir), 'CERN_density_', file_step, '.dat'
      open(10,file = file_name)
	write(10,1501) "particle density in EnSight Gold format"
	write(10,1501) "part"
	write(10,1502) 1
	write(10,1501) "coordinates"
	do i = 1, np
	  write(10,1503) rho(i)
	enddo
	close(10)
	
	write(file_name,1504) trim(filedir), 
     &                      'CERN_pressure_', file_step, '.dat'
      open(10,file = file_name)
	write(10,1501) "particle pressure in EnSight Gold format"
	write(10,1501) "part"
	write(10,1502) 1
	write(10,1501) "coordinates"
	do i = 1, np
	  write(10,1503) p(i)
	enddo
	close(10)
	
	write(file_name,1504) trim(filedir), 'CERN_energy_', file_step, '.dat'
      open(10,file = file_name)
	write(10,1501) "particle energy in EnSight Gold format"
	write(10,1501) "part"
	write(10,1502) 1
	write(10,1501) "coordinates"
	do i = 1, np
	  write(10,1503) u(i)
	enddo
	close(10)
	
1501	format(A)
1502  format(I10)
1503  format(E15.8E3)
1504	format(A,A,I5.5,A)
	
      end
	
	
      subroutine ensightout_case(dt, steps, save_step, filedir) 
      
c----------------------------------------------------------------------           
c     Subroutine for saving the master file of the results in the 
c     EnSight format to external disk file

c     dt-- time step                                                [in]
c     steps-- number of time steps                                  [in]
c     save_step-- number of time steps between saves                [in]
c     filedir -- Directory where the files are located              [in]

      implicit none     
      include 'options.inc'
      
	double precision dt, time
      integer steps, save_step
	character (len = *) :: filedir

      integer i, its
	character (len = 120) file_name
      
	write(file_name,1505)  trim(filedir), 'CERN_SPH.case'
      open(10,file = file_name)
	write(10,1502) "# Ensight formatted case file for SPH CERN"
	write(10,1502) ""
	write(10,1502) "FORMAT"
	write(10,1502) "type: ensight gold"
	write(10,1502) ""
	write(10,1502) "GEOMETRY"
	write(10,1501) "model:    1          ", "CERN", 
     &              "_position_*****.geo"
	write(10,1502) ""
	write(10,1502) "VARIABLE"
	write(10,1501) 'vector per node:    1 velocity ', 'CERN', 
     &              '_velocity_*****.dat'
	write(10,1501) 'scalar per node:    1 density  ', 'CERN', 
     &              '_density_*****.dat'
	write(10,1501) 'scalar per node:    1 pressure ', 'CERN', 
     &              '_pressure_*****.dat'
	write(10,1501) 'scalar per node:    1 energy   ', 'CERN', 
     &              '_energy_*****.dat'
	write(10,1502) ""
	write(10,1502) "TIME"
	write(10,1503) "time set:", 1
	write(10,1503) "number of steps:", (steps / save_step) +1
	write(10,1503) "filename start number:", 0
	write(10,1503) "filename increment:", 1
	write(10,1502) "time values:"
	
	time = 0.
	write(10,1504) time
	do its = 1, steps
        time = time + dt
	  if (mod(its, save_step) == 0) write(10,1504) time
	enddo
	
	close(10)
	
1501	format(A,A,A)
1502  format(A)
1503  format(A,I10)
1504  format(E15.8E3)
1505  format(A,A)
		
      end
