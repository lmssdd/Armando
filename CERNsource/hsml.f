      subroutine h_upgrade(nt, x, v, h, mass, rho, 
     &           niac, pair_i, pair_j, dt, sle, skf)

c-----------------------------------------------------------------------
c     Subroutine to evolve smoothing length

c     nt-- number of particles                                      [in]
c     x-- article position                                          [in]
c     v-- velocities of all particles                               [in]
c     h-- smoothing length                                      [in/out]
c     mass-- particle masses                                        [in]
c     rho-- density                                                 [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dt-- timestep                                                 [in]
c     sle-- smoothing length evolution                              [in]
c     skf-- smoothing kernel function                               [in]
   
      implicit none
      include 'options.inc'
      
      integer nt, niac, pair_i(max_interaction), 
     &        pair_j(max_interaction), sle, skf
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &        rho(maxn), dt
      integer i,j,k,d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision dv, vdwdx, dh, divv(maxn)

      if (sle.eq.0 ) then     
c---  Keep smoothing length unchanged. 
        return
      
      else if(sle.eq.1) then
        do i = 1, nt          
          h(i) = 2.0d0 * (mass(i)/rho(i))**(1./dim)
        enddo
       
      else if (sle.eq.2) then
c---  dh/dt = (-1/dim)*(h/rho)*(drho/dt).
        do i = 1, nt
          divv(i) = 0.0d0
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
          
          vdwdx = 0.0d0
          do d = 1, dim
            dv = v(d,j) - v(d,i) 
            vdwdx = vdwdx + dv * dwdx(d)
          enddo
          divv(i) = divv(i) + mass(j)/rho(j)* vdwdx
          divv(j) = divv(j) + mass(i)/rho(i)* vdwdx
        enddo  
        
        do i = 1, nt
          dh = (h(i) / dim) * divv(i)
          h(i) = h(i) + dt * dh
          if (h(i) .le. 0.0d0) h(i) = h(i) - dt*dh
        enddo
        
      endif 
       
      end
      