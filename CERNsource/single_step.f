	subroutine single_step(np, nv, ng, x, v, mat, h, 
     &           mass, rho, p, u, c, s, e, t, 
     &           niac, pair_i, pair_j, index_ghost, 
     &           dx, dv, drho, du, ds, tdsdt, av, prop, its, dt, 
     &           dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, 
     &           pa_sph, nnps, sle, skf, density_method, 
     &	         average_velocity, virtual_part, visc, 
     &           visc_artificial, heat_artificial, 
     &           heat_external, self_gravity)
	
c----------------------------------------------------------------------
c   Subroutine to determine the right hand side of a differential 
c   equation in a single step for performing time integration 

c   In this routine and its subroutines the SPH algorithms are performed.
c     np-- Number of particles                                     [in]
c     nv-- Number of virtual particles                             [in]
c     ng-- Number of ghost particles                               [in]
c     x-- article position                                         [in]
c     v-- particle velocity                                        [in]
c     mat-- particle material                                      [in]
c     h-- smoothing length                                         [in]
c     mass-- particle mass                                         [in]
c     rho-- density                                            [in/out]
c     p-- pressure                                                [out]
c     u-- particle internal energy                                 [in]
c     s-- particle entropy                                         [in]
c     e-- particle internal energy                                 [in]
c     t-- particle temperature                                     [in]
c     dx-- dx = v = dx/dt                                         [out]
c     dv-- dv = dv/dt                                             [out]
c     drho-- drho/dt                                              [out]
c     du-- du/dt                                                  [out]
c     ds-- ds/dt                                                  [out]     
c     tdsdt-- viscous entropy production                          [out]
c     av-- average velocity                                       [out]
c     prop-- material properties                                   [in]
c     its-- current timestep number                                [in]
c     dt- timestep                                                 [in]
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
c     sle -- Smoothing length evolution                            [in]
c     skf -- Smoothing kernel function                             [in]
c     density_method -- Density calculation method                 [in]
c	average_velocity                                             [in]
c	virtual_part                                                 [in]
c	visc                                                         [in]
c	visc_artificial                                              [in]
c	heat_artificial                                              [in]
c	heat_external                                                [in]
c	self_gravity                                                 [in]

      implicit none
      include 'options.inc'
	
      integer np, nv, ng, its, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), 
     &       mass(maxn), rho(maxn), p(maxn), u(maxn), 
     &       c(maxn), s(maxn), e(maxn), t(maxn)
      double precision dx(dim, maxn), dv(dim, maxn), 
     &       drho(maxn), du(maxn), ds(maxn), tdsdt(maxn), av(dim, maxn), 
     &       prop(16,40), dt
	
      integer niac, pair_i(max_interaction), pair_j(max_interaction)
      integer index_ghost(maxn)

	double precision dtfluka
	double precision ofluka(3), dfluka(3), datafluka(max_fluka)
	integer nfluka(3), tfluka

	integer pa_sph, nnps, sle, skf, density_method
	logical average_velocity, virtual_part, visc, 
     &	visc_artificial, heat_artificial, heat_external, self_gravity
	
	double precision indvdt(dim,maxn),exdvdt(dim,maxn), 
     &       avdvdt(dim,maxn), indudt(maxn), avdudt(maxn), 
     &       ahdudt(maxn), exdudt(maxn), eta(maxn)
	
      integer i, d
	double precision gamma, pshift, rho0, eta0, pmin, 
     &                 a1, a2, a3, b0, b1, t1, t2, 
     &                 c0, gamma0, s0, h0, es
	
      do i = 1, np
        indudt(i) = 0.
        ahdudt(i) = 0.
        avdudt(i) = 0.
        exdudt(i) = 0.
        
	  do d = 1, dim
          indvdt(d,i) = 0.
          avdvdt(d,i) = 0.
          exdvdt(d,i) = 0.
          av(d,i) = 0.
        enddo 
      enddo
	

c---  Density approximation or change rate
      if (density_method .eq. 0) then
        call con_density(np, x, v, h, mass, drho, 
     &                   niac, pair_i, pair_j, skf)
      else if (density_method .eq. 1) then
        call sum_density(np, x, h, mass, rho, 
     &                   niac, pair_i, pair_j, skf)
      else if (density_method .eq. 2) then
        call nor_density(np, x, h, mass, rho, 
     &                   niac, pair_i, pair_j, skf)
      endif
	
c     Equation of state
      do i = 1, np
        if ((abs(mat(i)) .gt. 0) .and. (abs(mat(i)) .le. 10)) then
c     Ideal gas EOS
          rho0 = prop(1, abs(mat(i)))
          gamma = prop(2, abs(mat(i)))
          pshift = prop(3, abs(mat(i)))
          eta0 = prop(4, abs(mat(i)))
	    
		call eos_gas(rho(i), u(i), p(i), gamma)
	    
		p(i) = p(i) +pshift
          eta(i)= eta0

        else if ((abs(mat(i)) .gt. 10) .and. (abs(mat(i)) .le. 20)) then
c     Mie-Gruneise polyynomial EOS
          rho0 = prop(1, abs(mat(i)))
          a1 = prop(2, abs(mat(i)))
          a2 = prop(3, abs(mat(i)))
          a3 = prop(4, abs(mat(i)))
          b0 = prop(5, abs(mat(i)))
          b1 = prop(6, abs(mat(i)))
          t1 = prop(7, abs(mat(i)))
          t2 = prop(8, abs(mat(i)))
          pmin = prop(9, abs(mat(i)))
          eta0 = prop(10, abs(mat(i)))
	    
		call eos_poly(rho(i), u(i), p(i), 
     &                  rho0, a1, a2, a3, b0, b1, t1, t2, pmin)

c		call p_art_water(rho(i), p(i), c(i))
	    
          eta(i)= eta0

        else if ((abs(mat(i)) .gt. 20) .and. (abs(mat(i)) .le. 30)) then
c     Mie-Gruneise shock EOS
          rho0 = prop(1, abs(mat(i)))
          c0 = prop(2, abs(mat(i)))
          gamma0 = prop(3, abs(mat(i)))
          s0 = prop(4, abs(mat(i)))
          pmin = prop(5, abs(mat(i)))
          eta0 = prop(6, abs(mat(i)))
	    
		call eos_shock(rho(i), u(i), p(i), 
     &                   rho0, c0, gamma0, s0, pmin)

          eta(i)= eta0

        else if ((abs(mat(i)) .gt. 30) .and. (abs(mat(i)) .le. 40)) then
c     Mie-Gruneise Puff EOS
	    rho0 = prop(1, abs(mat(i)))
		a1 = prop(2, abs(mat(i)))
		a2 = prop(3, abs(mat(i)))
		a3 = prop(4, abs(mat(i)))
		gamma0 = prop(5, abs(mat(i)))
		h0 = prop(6, abs(mat(i)))
		es = prop(7, abs(mat(i)))
		t1 = prop(8, abs(mat(i)))
		t2 = prop(9, abs(mat(i)))
          pmin = prop(10, abs(mat(i)))
          eta0 = prop(11, abs(mat(i)))
	    
		call eos_puff(rho(i), u(i), p(i), 
     &               rho0, a1, a2, a3, gamma0, h0, es, t1, t2, pmin)
	    
          eta(i)= eta0
        endif
      enddo
	
c---	Update ghost particles
	if (virtual_part) then
	  do i = 1, np
	    if (index_ghost(i) .ne. 0) then
	      rho(index_ghost(i)) = rho(i)
	      p(index_ghost(i)) = p(i)
	      u(index_ghost(i)) = u(i)
	    endif
	  enddo
	endif

c---  Internal forces:
	call int_force(np, x, v, mat, h, 
     &     mass, rho, p, eta, indvdt, indudt, tdsdt, 
     &     niac, pair_i, pair_j, dt, skf, pa_sph, visc)

c---  Artificial viscosity:
      if (visc_artificial) call art_visc(np, x, v, h, 
     &    mass, rho, c, avdvdt, avdudt,
     &    niac, pair_i, pair_j, skf)
      
c---  Calculate the artificial heating on the particles
      if (heat_artificial) call art_heat(np, x, v, h, 
     &    mass, rho, u, c, ahdudt,
     &    niac, pair_i, pair_j, skf)
	
c     Calculating average velocity of each partile for avoiding penetration
      if (average_velocity) call av_vel(np, x, v, h, mass, rho, av,
     &                                  niac, pair_i, pair_j, skf)

c---  Calculate the external heating on the particles
	if ((heat_external) .and. (its*dt .le. dtfluka)) then
	  call fluka_heat(np, x, rho, exdudt,
     &        ofluka, dfluka, nfluka, tfluka, datafluka)
	endif
	
c     Calculating the neighboring particles and undating H
      if (sle .ne. 0) call h_upgrade(np, x, v, h, mass, rho, 
     &                     niac, pair_i, pair_j, dt, sle, skf)

c---  Convert velocity, force, and energy to f and dfdt  
      do i = 1, np
        do d = 1, dim
          dv(d,i) = indvdt(d,i) + exdvdt(d,i) + avdvdt(d,i)
        enddo
        du(i) = indudt(i) + avdudt(i) + ahdudt(i) + exdudt(i)
      enddo

      end
