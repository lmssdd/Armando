      subroutine sum_density(np, x, h, mass, rho, niac, pair_i, pair_j, 
     &                       skf)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH summation algorithm.

c     np-- number of particles                                     [in]
c     x-- article position                                         [in]
c     h-- smoothing length                                         [in]
c     mass-- particle masses                                       [in]
c     rho-- density                                               [out]
c     niac-- number of interaction pairs                           [in]
c     pair_i-- list of first partner of interaction pair           [in]
c     pair_j-- list of second partner of interaction pair          [in]
c     skf-- smoothing kernel function                              [in]
	
      implicit none
      include 'options.inc'
      
      integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
      integer skf
      double precision x(dim, maxn), h(maxn), mass(maxn), rho(maxn)
      
      integer i, j, k, d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision selfdens, r, hv(dim)
	

      do d = 1, dim
        hv(d) = 0.0d0
      enddo
	
c     Self density of each particle
	
      r=0.
      
c     Calculate the rho integration over the space
	
      do i = 1, np
	  call kernel(r, hv, h(i), selfdens, hv, skf)   
        rho(i) = selfdens * mass(i)
      enddo
	
c     Calculate SPH sum for rho:
      do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        
        dr2iac = 0.0
        do d = 1, dim
          dxiac(d) = x(d,i) - x(d,j)
          dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
	  if (i <= np) rho(i) = rho(i) + mass(j)*w
	  if (j <= np) rho(j) = rho(j) + mass(i)*w
      enddo
      
      end
      
	
      subroutine nor_density(np, x, h, mass, rho, niac, pair_i, pair_j,
     &           skf)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH summation algorithm
c     with density normalization to avoid the boundary deficiency
c     problem (Randles and Libersky 1996).

c     np-- number of particles                                     [in]
c     x-- article position                                         [in]
c     h-- smoothing length                                         [in]
c     mass-- particle masses                                       [in]
c     rho-- density                                            [in/out]
c     niac-- number of interaction pairs                           [in]
c     pair_i-- list of first partner of interaction pair           [in]
c     pair_j-- list of second partner of interaction pair          [in]
c     skf-- smoothing kernel function                              [in]
	
      implicit none
      include 'options.inc'
      
      integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
      integer skf
      
      double precision x(dim, maxn), h(maxn), mass(maxn), rho(maxn)
      
      integer i, j, k, d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision selfdens, r, hv(dim), wi(maxn)
	
c     wi(maxn)---integration of the kernel itself
        
      do d = 1, dim
        hv(d) = 0.d0
      enddo
	
      
c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.

c     Firstly calculate the integration of the kernel over the space
	
      do i = 1, np
	  call kernel(r, hv, h(i), selfdens, hv, skf)   
        wi(i) = selfdens * mass(i) / rho(i)
      enddo

c     Secondly calculate the rho integration over the space
	
      do i = 1, np
	  call kernel(r, hv, h(i), selfdens, hv, skf)   
        rho(i) = selfdens * mass(i)
      enddo
	
      do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        
        dr2iac = 0.0
        do d = 1, dim
          dxiac(d) = x(d,i) - x(d,j)
          dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
	  if (i <= np) wi(i) = wi(i) + mass(j) / rho(j) * w
	  if (j <= np) wi(j) = wi(j) + mass(i) / rho(i) * w
      enddo
	
c     Calculate SPH sum for rho:
      do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        
        dr2iac = 0.0
        do d = 1, dim
          dxiac(d) = x(d,i) - x(d,j)
          dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
	  if (i <= np) rho(i) = rho(i) + mass(j) * w
	  if (j <= np) rho(j) = rho(j) + mass(i) * w
      enddo
	
c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
     
      do i = 1, np
        rho(i) = rho(i) / wi(i)
      enddo
      
      end
	
	
	
      subroutine con_density(np, x, v, h, mass, drho, 
     &                       niac, pair_i, pair_j, skf)

c----------------------------------------------------------------------
c     Subroutine to calculate the density variation based on the 
c     continuity equation.

c     np-- number of particles                                     [in]
c     x-- article position                                         [in]
c     v-- particle velocities                                      [in]
c     h-- smoothing length                                         [in]
c     mass-- particle masses                                       [in]
c     drho-- density change rate of each particle                 [out]
c     niac-- number of interaction pairs                           [in]
c     pair_i-- list of first partner of interaction pair           [in]
c     pair_j-- list of second partner of interaction pair          [in]
c     skf-- smoothing kernel function                              [in]
	
      implicit none
      include 'options.inc'
      
      integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
      integer skf
      
      double precision x(dim, maxn), v(dim, maxn), h(maxn),
     &                 mass(maxn), drho(maxn)
     
      integer i, j, k, d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision vcc, dv(dim)
      
      do i = 1, np
        drho(i) = 0.
      enddo
     
      do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)

        dr2iac = 0.0
        do d = 1, dim
          dxiac(d) = x(d,i) - x(d,j)
          dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
	  vcc = 0.0d0
	  do d = 1, dim
	    vcc = vcc + (v(d,i) - v(d,j)) * dwdx(d)
	  enddo

	  if (i <= np) drho(i) = drho(i) + mass(j) * vcc
	  if (j <= np) drho(j) = drho(j) + mass(i) * vcc
      enddo
	
      end

