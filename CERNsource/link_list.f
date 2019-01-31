	subroutine link_list(nt, x, h, niac, pair_i, pair_j, skf)
c	
c----------------------------------------------------------------------
c   Subroutine to calculate the interaction pairs between the particles
c   the pairs are determined by by using a sorting grid linked list  
c	
c     nt-- total number of particles                               [in]
c     x-- coordinates of all particles                             [in]
c     h-- smoothing length                                         [in]
c     niac-- number of interaction pairs                          [out]
c     pair_i-- list of first partner of interaction pair          [out]
c     pair_j-- list of second partner of interaction pair         [out]
c     skf -- Smoothing kernel function                             [in]

      implicit none
      include 'options.inc'
      
      integer nt, skf, 
     &        niac, pair_i(max_interaction), pair_j(max_interaction)
      double precision x(dim,maxn), h(maxn)
      
	integer scale_k, i, j, k, d, g, counter, ig, jg, kg, grid_nt, 
     &        grid_n(dim), grid_ijk(dim), grid_index(max_grid_data), 
     &        grid_part(maxn), grid_min_ijk(dim), grid_max_ijk(dim)
      double precision dxiac(dim), dr2iac, mh
      double precision grid_min(dim), grid_max(dim), tdwdx(dim)     


      if (skf .eq. 0) then 
        scale_k = 2 
      else if (skf .eq. 1) then 
        scale_k = 3 
      else if (skf .eq. 2) then 
         scale_k = 3 
      endif 
	
	do d = 1, dim
	  grid_min(d) = +1.0d30
	  grid_max(d) = -1.0d30
	enddo
	mh = 0.0d0
	
	do i = 1, nt
	  do d = 1, dim
	    if (x(d,i) < grid_min(d)) grid_min(d) = x(d,i)
	    if (x(d,i) > grid_max(d)) grid_max(d) = x(d,i)
	  enddo
	  
	  if (h(i) > mh) mh = h(i)
	  grid_part(i) = 0
	enddo
	
	do d = 1, dim
	  grid_min(d) = grid_min(d) -mh
	  grid_max(d) = grid_max(d) +mh
	  grid_n(d) = int((grid_max(d) - grid_min(d)) / (scale_k * mh)) +1
	enddo

	grid_nt = 0
	if (dim .eq. 2) grid_nt = grid_n(1) * grid_n(2)
	if (dim .eq. 3) grid_nt = grid_n(1) * grid_n(2) * grid_n(3)
	
	if (grid_nt > max_grid_data) then
        print *,
     &  ' >>> ERROR <<< : Not enough memory for grid' 
	  stop
	endif
	
c	Set to zero the index vector of the grid data structure
	do i = 1, grid_nt
	  grid_index(i) = 0
	enddo
	
c	The index vector is used to store the number of particles in each grid cell
	do i = 1, nt
	  do d = 1, dim
	    grid_ijk(d) = int((x(d,i) - grid_min(d)) / (scale_k *mh))
	  enddo
	  
	  if (dim .eq. 2) then
	    g = grid_n(1)*grid_ijk(2) +grid_ijk(1) +1
	  else if (dim .eq. 3) then
	    g = grid_n(1)*grid_n(2)*grid_ijk(3) 
     &       +grid_n(1)*grid_ijk(2) +grid_ijk(1) +1
	  endif

	  grid_index(g) = grid_index(g) +1
	enddo
	
c	The index vector points at the beginning of the particle list 
c	in the grid data structure
	counter = 1
	do i = 1, grid_nt
	  g = grid_index(i)
	  grid_index(i) = counter
	  counter = counter + g
	enddo
	
c	if (counter <> nt) then
c        print *,
c     &  ' >>> ERROR <<< : Particles out of grid' 
c	  stop
c	endif
	
c	The data vector for particles is filled
	do i = 1, nt
	  do d = 1, dim
	    grid_ijk(d) = int((x(d,i) - grid_min(d)) / (scale_k *mh))
	  enddo
	  
	  if (dim .eq. 2) then
	    g = grid_n(1)*grid_ijk(2) +grid_ijk(1) +1
	  else if (dim .eq. 3) then
	    g = grid_n(1)*grid_n(2)*grid_ijk(3) 
     &       +grid_n(1)*grid_ijk(2) +grid_ijk(1) +1
	  endif

	  do k = grid_index(g), grid_index(g +1) -1
	    if (grid_part(k) == 0) then
	      grid_part(k) = i
	      exit
	    endif
	  enddo
	enddo
	
	
c     Search grid:
      
      niac = 0
      do i = 1, nt -1
	  do d = 1, dim
	    grid_ijk(d) = int((x(d,i) - grid_min(d)) / (scale_k *mh))
	    
		grid_min_ijk(d) = grid_ijk(d) -1
		grid_max_ijk(d) = grid_ijk(d) +1
	    
		if (grid_min_ijk(d) < 0) then
		  grid_min_ijk(d) = 0
	    endif
		if (grid_max_ijk(d) >= grid_n(d)) then
		  grid_max_ijk(d) = grid_n(d) -1
	    endif
	  enddo
	  
	  if (dim .eq. 2) then
	    do jg = grid_min_ijk(2), grid_max_ijk(2)
	      do ig = grid_min_ijk(1), grid_max_ijk(1)
	        g = grid_n(1)*jg +ig +1
	        
			do k = grid_index(g), grid_index(g +1) -1
	          j = grid_part(k)
	          if (j.gt.i) then
	            dr2iac = 0.0
                  do d = 1, 2
                    dxiac(d) = x(d,i) - x(d,j)
                    dr2iac = dr2iac + dxiac(d)*dxiac(d)
                  enddo
                  mh = (h(i) +h(j)) / 2.
                  if (dr2iac < (scale_k * mh)**2) then
                    if (niac < max_interaction) then    
	                
c     Neighboring pair list, and total interaction number and
c     the interaction number for each particle 
                      
                      niac = niac + 1
                      pair_i(niac) = i
                      pair_j(niac) = j
                      
	              else
	                print *,
     &                ' >>> ERROR <<< : Too many interactions' 
	                stop
	              endif
	            endif
		      endif
	        enddo

	      enddo
	    enddo
	  
	  else if (dim .eq. 3) then
	    
	    do kg = grid_min_ijk(3), grid_max_ijk(3)
	      do jg = grid_min_ijk(2), grid_max_ijk(2)
	        do ig = grid_min_ijk(1), grid_max_ijk(1)
	          g = grid_n(1)*grid_n(2)*kg +grid_n(1)*jg +ig +1
	          
	          do k = grid_index(g), grid_index(g +1) -1
	            j = grid_part(k)
	            if (j.gt.i) then
	              dr2iac = 0.0
                    do d = 1, 3
                      dxiac(d) = x(d,i) - x(d,j)
                      dr2iac = dr2iac + dxiac(d)*dxiac(d)
                    enddo
                    mh = (h(i) +h(j)) / 2.
                    if (dr2iac < (scale_k * mh)**2) then
                      if (niac < max_interaction) then    
	                  
c     Neighboring pair list, and total interaction number and
c     the interaction number for each particle 
                      
                        niac = niac + 1
                        pair_i(niac) = i
                        pair_j(niac) = j
                        
	                else
	                  print *,
     &                  ' >>> ERROR <<< : Too many interactions' 
	                  stop
	                endif
	              endif
		        endif
	          enddo
	          
	        enddo
	      enddo
	    enddo
	  endif

	enddo

	end

