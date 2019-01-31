	subroutine smooth_particles(np, nv, x, v, h, mass, rho, 
     &                            virtual_part)
c----------------------------------------------------------------------
c   Subroutine to optimize the mass of the particles close to the 
c   boundary to have a uniform density in the computational domain

c     np-- Number of particles                                     [in]
c     nv-- Number of virtual particles                             [in]
c     x-- article position                                         [in]
c     v-- particle velocity                                        [in]
c     h-- smoothing length                                         [in]
c     mass-- particle mass                                     [in/out]
c     rho-- density                                                [in]
c	virtual_part                                                 [in]
c	
	
      implicit none
      include 'options.inc'
	
      integer np, nv
      double precision x(dim, maxn), v(dim, maxn), h(maxn), 
     &       mass(maxn), rho(maxn)
	logical virtual_part

c	double precision, dimension(:,:), allocatable :: dx
	double precision, dimension(:), allocatable :: phi
	integer, dimension(:), allocatable :: pair_i, pair_j, index_ghost
      integer ng, niac, i, j, k, d, its, nits
c	double precision dxm, dxmax, hv(dim), dd
	double precision dr, drmax, drold
	double precision dd1, dd2

	allocate(phi(maxn))
	allocate(pair_i(max_interaction), pair_j(max_interaction))
	allocate(index_ghost(maxn))


	call link_list(np +nv +ng, x, h, niac, pair_i, pair_j, 0)
	
c	nits = 10
c	do its = 1, nits
c
c	  ng = 0
c	  if (virtual_part) then
c	    call ghost_particles(np, nv, ng, x, v, h, mass, rho, phi, 
c     &                       niac, pair_i, pair_j, index_ghost)
c	  endif
c	  
c	  call link_list(np +nv +ng, x, h, niac, pair_i, pair_j,
c     &                 w, dwdx, 0)
c	  
c	  do k = 1, niac
c	    i = pair_i(k)
c	    j = pair_j(k)
c
c	    if ((i .le. np) .and. 
c     &      ((j .le. np) .or. (j .gt. (np +nv)))) then
c	        do d = 1, dim
c	          hv(d) = x(d, j) - x(d, i)
c	        enddo
c	        
c	        dd = 0.
c	        do d = 1, dim
c	          dd = dd + hv(d) * hv(d)
c	        enddo
c	        dd = dsqrt(dd)
c
c	      do d = 1, dim
c	        dx(d, i) = dx(d, i) - hv(d) / dd**2
c	        dx(d, j) = dx(d, j) + hv(d) / dd**2
c	      enddo
c	    endif
c	  enddo
c
c	  do i = 1, np
c	    if (index_ghost(i) .ne. 0) then
c	      j = index_ghost(i)
c
c	      do d = 1, dim
c	        hv(d) = x(d, j) - x(d, i)
c	      enddo
c
c	      dd = 0.
c	      do d = 1, dim
c	        dd = dd + hv(d) * hv(d)
c	      enddo
c	      dd = dsqrt(dd)
c
c	      do d = 1, dim
c	        hv(d) = hv(d) / dd
c	      enddo
c
c	      dxm = 0.
c	      do d = 1, dim
c	        dxm = dxm + dx(d, i) * hv(d)
c	      enddo
c	      
c	      do d = 1, dim
c	        dx(d, i) = dxm * hv(d)
c	      enddo
c
c	    else
c	      do d = 1, dim
c	        dx(d, i) = 0.0
c	      enddo
c	    endif
c
c	  enddo
c
c
c	  dxmax = 0.
c	  do i = 1, np
c	    if (index_ghost(i) .ne. 0) then
c	      dxm = 0.
c	      do d = 1, dim
c	        dxm = dxm + dx(d, i)**2
c	      enddo
c	      dxm = dsqrt(dxm)
c	      
c	      if (dxm .gt. dxmax) dxmax = dxm
c	    endif
c	  enddo
c
c	  do i = 1, np
c	    if (index_ghost(i) .ne. 0) then
c	      do d = 1, dim
c	        dx(d, i) = h(i) / dxmax / nits * dx(d, i)
c	      enddo
c	    endif
c	  enddo
c
c	  do i = 1, np
c	    if (index_ghost(i) .ne. 0) then
c	      j = index_ghost(i)
c
c	      do d = 1, dim
c	        hv(d) = x(d, j) - x(d, i)
c	      enddo
c
c	      dd = 0.
c	      do d = 1, dim
c	        dd = dd + hv(d) * hv(d)
c	      enddo
c	      dd = dsqrt(dd)
c	      
c	      if (dd .lt. 0.5*h(i)) then
c	        do d = 1, dim
c	          dx(d, i) = 0.5 * (1. -0.5*h(i)/dd) * hv(d)
c	        enddo
c	      endif
c
c	      do d = 1, dim
cc	        x(d, i) = x(d, i) + dx(d, i)
c	      enddo
c
c	    endif
c	  enddo
c
c	enddo
      
	write(*,'(A)') 'Optimizing densities' 
	
	nits = np

	do i = 1, np
	  index_ghost(i) = 0
	enddo
	
      do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
	  if (i <= np) then
	    if ((j > np) .and. (j <= (np + nv))) then
	      if (index_ghost(i) .eq. 0) then
	        index_ghost(i) = j
	      else
	        dd1 = 0.
	        dd2 = 0.
	        do d = 1, dim
	          dd1 = dd1 + (x(d, index_ghost(i)) - x(d, i))**2
	          dd2 = dd2 + (x(d, j) - x(d, i))**2
	        enddo
			if (dd2 .lt. dd1) index_ghost(i) = j
	      endif
	    endif
	  endif
	enddo
	
	ng = 0
	do i = 1, np
	  if (index_ghost(i) .ne. 0) then
	    ng = ng +1
	    j = index_ghost(i)

	    dd1 = 0.
	    do d = 1, dim
	      dd1 = dd1 + ((x(d, i) - x(d, j)) * v(d, j))
	    enddo
	    
	    do d = 1, dim
	      x(d, np +nv +ng) = x(d, i) -2.0d0 * v(d, j) * dd1
	    enddo

	    h(np +nv +ng) = h(i)
	    index_ghost(i) = np +nv +ng
	  endif
	enddo
	
	call link_list(np +nv +ng, x, h, niac, pair_i, pair_j, 0)
	
	drold = 1.0d30

	do its = 1, nits
	  do i = 1, np
	    if (index_ghost(i) .ne. 0) then
	      rho(index_ghost(i)) = rho(i)
	      mass(index_ghost(i)) = mass(i)
	    endif
	  enddo

	  do i = 1, np +nv +ng
	    phi(i) = rho(i)
	  enddo

	  call sum_density(np, x, h, mass, phi, 
     &                   niac, pair_i, pair_j, 0)
        
	  drmax = 0.0
	  do i = 1, np
	    if (index_ghost(i) .ne. 0) then
	      mass(i) = mass(i) * rho(i) / phi(i)
	      dr = dabs(rho(i) - phi(i)) / (rho(i))
	      if (dr > drmax) drmax = dr
	    endif
	  enddo
	  if (drmax .gt. drold) exit
	  drold = drmax

	enddo

	write(*,'(A, I)') 'Optimization cycles :', its 

	deallocate(phi)
	deallocate(pair_i, pair_j)
	deallocate(index_ghost)
	
	end
