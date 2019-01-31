      subroutine ghost_particles(np, nv, ng, x, v, h, mass, rho, p, 
     &                           niac, pair_i, pair_j, index_ghost)

c----------------------------------------------------------------------
c     Subroutine to calculate the number of ghost particles and to 
c     define an array that linking real and ghost particles

c     np-- number of particles                                     [in]
c     nv-- number of virtual boundary particles                    [in]
c     ng-- number of ghost particles                               [in]
c     x-- article position                                         [in]
c     v-- particle velocity                                        [in]
c     h-- smoothing length                                         [in]
c     mass-- particle mass                                         [in]
c     rho-- density                                                [in]
c     p-- pressure                                                 [in]
c     niac-- number of interaction pairs                           [in]
c     pair_i-- list of first partner of interaction pair           [in]
c     pair_j-- list of second partner of interaction pair          [in]

      implicit none
      include 'options.inc'
      
      integer np, nv, ng
      integer niac, pair_i(max_interaction), pair_j(max_interaction)
      integer index_ghost(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn),
     &                 mass(maxn), rho(maxn), p(maxn)
	
	integer i, j, k, d
	double precision dd1, dd2
	
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
	    j = index_ghost(i)
		
	    ng = ng +1
	    
	    dd1 = 0.
	    dd2 = 0.
	    do d = 1, dim
	      dd1 = dd1 + ((x(d, i) - x(d, j)) * v(d, j))
	      dd2 = dd2 + (v(d, i) * v(d, j))
	    enddo
	    
	    do d = 1, dim
	      x(d, np +nv +ng) = x(d, i) -2.0d0 * v(d, j) * dd1
	      v(d, np +nv +ng) = v(d, i) -2.0d0 * v(d, j) * dd2
	    enddo
	    
	    mass(np +nv +ng) = mass(i)
	    h(np +nv +ng) = h(i)
	    rho(np +nv +ng) = rho(i)
	    p(np +nv +ng) = p(i)
          
	    index_ghost(i) = np +nv +ng
          
	  endif
	enddo

	end
	