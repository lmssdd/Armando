      subroutine time_integration(np, nv, x, v, mat, h, 
     &           mass, rho, p, u, c, s, e, prop, 
     &           dt, steps, refresh_step, save_step,
     &           dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, 
     &           pa_sph, nnps, sle, skf, density_method, 
     &	         average_velocity, virtual_part, visc, 
     &           visc_artificial, heat_artificial, 
     &           heat_external, self_gravity, filedir)
     
c----------------------------------------------------------------------
c    Subroutine for leapfrog time integration. The particle position, 
c    velocity and energy is updated on the basis of the right hand side 
c    of the momentum and energy equation as calculated by the 
c    single_step subroutine
c
c     np-- particle number                                        [in]
c     nv-- virtual particle number                                [in]
c     x-- coordinates of particles                            [in/out]
c     v-- velocities of particles                             [in/out]
c     mat-- material type                                         [in]
c     h-- smoothing lengths of particles                      [in/out]
c     mass-- mass of particles                                    [in]
c     rho-- densities of particles                            [in/out]
c     p-- pressure  of particles                              [in/out]
c     u-- internal energy of particles                        [in/out]
c     c-- sound velocity of particles                            [out]
c     s-- entropy of particles, not used here                    [out]
c     e-- total energy of particles                              [out]
c     prop-- material properties                                  [in]
c     dt-- timestep                                               [in]
c     steps-- maximum timesteps                                   [in]
c     refresh_step-- timesteps between connectivity refreshing    [in]
c     save_step-- timesteps between saves                         [in]
c	
c     dtfluka-- deposition time for the binning                    [in]
c     ofluka-- the origin of the binning                           [in]
c     dfluka-- cell side of the binning                            [in]
c     nfluka-- number of cells for each index                      [in]
c     tfluka-- identifier for the fluka grid type                  [in]
c     datafluka-- data of the binning                              [in]
c	
c     pa_sph -- SPH algorithm for particle approximation           [in]
c     nnps -- Nearest neighbor particle searching method           [in]
c     sle -- Smoothing kernel function                             [in]
c     skf -- Smoothing kernel function                             [in]
c     density_method -- Density calculation method                 [in]
c	average_velocity                                             [in]
c	virtual_part                                                 [in]
c	visc                                                         [in]
c	visc_artificial                                              [in]
c	heat_artificial                                              [in]
c	heat_external                                                [in]
c	self_gravity                                                 [in]
c     filedir -- Directory where the files are located             [in]

	

      implicit none     
      include 'options.inc'
      
      integer np, nv, mat(maxn), steps, refresh_step, save_step
      double precision x(dim, maxn), v(dim, maxn), h(maxn), 
     &       mass(maxn), rho(maxn), p(maxn), u(maxn), 
     &       c(maxn), s(maxn), e(maxn), prop(16,40), dt
	
	double precision dtfluka
	double precision ofluka(3), dfluka(3), datafluka(max_fluka)
	integer nfluka(3), tfluka

	integer pa_sph, nnps, sle, skf, density_method
	logical average_velocity, virtual_part, visc, 
     &	visc_artificial, heat_artificial, heat_external, self_gravity


	integer niac
      integer, dimension(:), allocatable :: pair_i, pair_j
      integer, dimension(:), allocatable ::  index_ghost

      double precision, dimension(:,:), allocatable :: v_min, dx, dv, av
      double precision, dimension(:), allocatable ::  rho_min, u_min,
     &                  du, drho, ds, t, tdsdt
	character (len = *) :: filedir

      double precision  time

      integer i, j, k, its, d, ng, file_step

     
      allocate(pair_i(max_interaction), pair_j(max_interaction))
      allocate(index_ghost(maxn))
      
      allocate(v_min(dim, maxn), dx(dim, maxn), dv(dim, maxn))
      allocate(av(dim, maxn))
	allocate(rho_min(maxn), u_min(maxn), du(maxn), drho(maxn))
	allocate(ds(maxn), t(maxn), tdsdt(maxn))

	file_step = 0 ! lucam
	time = 0.0
	
c----------------------------------------------------------------------
c---  Interaction parameters, calculating neighboring particles
	ng = 0
      if (nnps .eq. 0) then 
        call direct_find(np +nv +ng, x, h, niac, pair_i, pair_j, skf)
      else if (nnps .eq. 1) then
	  call link_list(np +nv +ng, x, h, niac, pair_i, pair_j, skf)
	endif
	
c---	Calculating ghost particles
	if (virtual_part) then
	  call ghost_particles(np, nv, ng, x, v, h, mass, rho, p, 
     &                       niac, pair_i, pair_j, index_ghost)
	endif

      if (nnps .eq. 0) then 
        call direct_find(np +nv +ng, x, h, niac, pair_i, pair_j, skf)
      else if (nnps .eq. 1) then
	  call link_list(np +nv +ng, x, h, niac, pair_i, pair_j, skf)
	endif
	
	
c     Regularize the density
c	call nor_density(np, h, mass, rho, niac, pair_i, pair_j, w,
c     &                 skf)

	call ensightout_case(dt, steps, save_step, filedir) ! lucam
	call output(np, nv, x, v, mat, h, mass, rho, p, u, c, filedir)
	call ensightout_results(np+nv, x, v, rho, p, u, file_step, filedir) ! lucam

c----------------------------------------------------------------------
      do its = 1, steps

c     If not first time step, then update thermal energy, density and 
c     velocity half a time step

        if (its .ne. 1) then
          do i = 1, np
            u_min(i) = u(i)
            u(i) = u(i) + (dt/2.) * du(i)
            
		  if(u(i).lt.0)  u(i) = 0.
            
            if (density_method .eq. 0) then    
              rho_min(i) = rho(i)
              rho(i) = rho(i) + (dt/2.) * drho(i)
            endif 
           
            do d = 1, dim
              v_min(d, i) = v(d, i)
              v(d, i) = v(d, i) + (dt/2.) * dv(d, i)
            enddo
          enddo

        endif
	  
c---	Calculating ghost particles
	  ng = 0
	  if (virtual_part) then
	    call ghost_particles(np, nv, ng, x, v, h, mass, rho, p, 
     &                         niac, pair_i, pair_j, index_ghost)
	  endif
	  
c---  Interaction parameters, calculating neighboring particles
	  if (mod(its, refresh_step) .eq. 0) then
          if (nnps .eq. 0) then 
            call direct_find(np +nv +ng, x, h, 
     &                       niac, pair_i, pair_j, skf)
          else if (nnps .eq. 1) then
            call link_list(np +nv +ng, x, h, niac, pair_i, pair_j, skf)
          endif
        endif
        
c---  Definition of variables out of the function vector:
        
	  call single_step(np, nv, ng, x, v, mat, h, 
     &       mass, rho, p, u, c, s, e, t, 
     &       niac, pair_i, pair_j, index_ghost, 
     &       dx, dv, drho, du, ds, tdsdt, av, prop, its, dt, 
     &       dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, 
     &       pa_sph, nnps, sle, skf, density_method, 
     &	   average_velocity, virtual_part, visc, 
     &       visc_artificial, heat_artificial, 
     &       heat_external, self_gravity)

        
        if (its .ne. 1) then
		
          do i = 1, np
            u(i) = u_min(i) + dt * du(i)
            if(u(i).lt.0)  u(i) = 0.          
            
            if (density_method .eq. 0) then
              rho(i) = rho_min(i) + dt * drho(i)
            endif
            
            do d = 1, dim
              v(d, i) = v_min(d, i) + dt * dv(d, i) + av(d, i)
              x(d, i) = x(d, i) + dt * v(d, i)                  
            enddo
          enddo
		
        else
          
          do i = 1, np
            u(i) = u(i) + (dt/2.) * du(i)
            if(u(i).lt.0)  u(i) = 0.
	      
            if (density_method .eq. 0) then
              rho(i) = rho(i) + (dt/2.) * drho(i)
            endif
            
            do d = 1, dim
              v(d, i) = v(d, i) + (dt/2.) * dv(d, i) ! + av(d, i)
              x(d, i) = x(d, i) + dt * v(d, i)
            enddo 
          enddo 
          
        endif
	  
        time = time + dt

	  if (mod(its, save_step) == 0) then
          write(*,*)'______________________________________________'
          write(*,*)'  current number of time step =',
     &              its,'     current time=', real(time +dt)
          write(*,*)'______________________________________________'
      
	    call output(np, nv, x, v, mat, h, mass, rho, p, u, c, filedir)
		
		file_step = file_step +1 ! lucam
		call ensightout_results(np, x, v, rho, p, u, 
     &                            file_step, filedir) ! lucam
	  endif
	  
      enddo
	
	call output(np, nv, x, v, mat, h, mass, rho, p, u, c, filedir)
	
	deallocate(pair_i, pair_j)
      deallocate(index_ghost)
      
      deallocate(v_min, dx, dv)
      deallocate(av)
	deallocate(rho_min, u_min, du, drho)
	deallocate(ds, t, tdsdt)

      end
