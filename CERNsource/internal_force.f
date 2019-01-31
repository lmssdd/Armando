	subroutine int_force(nt, x, v, mat, h, 
     &           mass, rho, p, eta, dvdt, dudt, tdsdt, 
     &           niac, pair_i, pair_j, dt, skf, pa_sph, visc)
	
c----------------------------------------------------------------------
c   Subroutine to calculate the internal forces on the right hand side 
c   of the Navier-Stokes equations, i.e. the pressure gradient and the
c   gradient of the viscous stress tensor, used by the time integration. 
c   Moreover the entropy production due to viscous dissipation, tds/dt, 
c   and the change of internal energy per mass, du/dt, are calculated. 

c     nt-- total number of particles                               [in]
c     x-- article position                                         [in]
c     v-- particle velocity                                        [in]
c     mat-- particle material                                      [in]
c     h-- smoothing length                                         [in]
c     mass-- particle mass                                         [in]
c     rho-- density                                                [in]
c     p-- pressure                                                 [in]
c     u-- particle internal energy                                 [in]
c     c-- particle sound velocity                                  [in]
c     eta-- dynamic viscosity                                     [in]
c     dvdt-- avcceleration                                        [out]
c     dudt-- internal energy change rate                          [out]
c     tdsdt-- viscous entropy production                          [out]
c     niac-- number of interaction pairs                           [in]
c     pair_i-- list of first partner of interaction pair           [in]
c     pair_j-- list of second partner of interaction pair          [in]
c     dt- timestep                                                 [in]
c     skf-- smoothing kernel function                              [in]
c	pa_sph                                                       [in]
c	visc                                                         [in]

      implicit none
      include 'options.inc'
	
      integer nt, mat(maxn),
     &        niac, pair_i(max_interaction), pair_j(max_interaction)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), 
     &       mass(maxn), rho(maxn), p(maxn), u(maxn), 
     &       c(maxn), eta(maxn), dt
      double precision dvdt(dim, maxn), dudt(maxn), tdsdt(maxn)
	integer skf, pa_sph
	logical visc
	
      integer i, j, k, d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision  txx(maxn), tyy(maxn), tzz(maxn), 
     &                  txy(maxn), txz(maxn), tyz(maxn), 
     &                  vcc(maxn), dv(dim),
     &                  hxx, hyy, hzz, hxy, hxz, hyz, hv, he, hvcc
	
	
c     Initialization of shear tensor, velocity divergence, 
c     viscous energy, internal energy, acceleration 

      do i = 1, nt
        txx(i) = 0.0d0
        tyy(i) = 0.0d0
        tzz(i) = 0.0d0
        txy(i) = 0.0d0
        txz(i) = 0.0d0
        tyz(i) = 0.0d0
        
	  vcc(i) = 0.0d0
        tdsdt(i) = 0.0d0
        dudt(i) = 0.0d0
        do d = 1 ,dim
          dvdt(d,i) = 0.0d0
        enddo
      enddo
	
c     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

      if (visc) then 
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
          
          do d = 1, dim
            dv(d) = v(d,j) -v(d,i)
          enddo
		
          if (dim .eq. 1) then 
            hxx = 2.0d00 * dv(1)*dwdx(1)        
          else if (dim .eq. 2) then           
            hxx = 2.0d0*dv(1)*dwdx(1) -dv(2)*dwdx(2)
            hyy = 2.0d0*dv(2)*dwdx(2) -dv(1)*dwdx(1)
            hxy = dv(1)*dwdx(2) +dv(2)*dwdx(1)
          else if (dim .eq. 3) then
            hxx = 2.0d0*dv(1)*dwdx(1) -dv(2)*dwdx(2)
     &                                -dv(3)*dwdx(3)
            hyy = 2.0d0*dv(2)*dwdx(2) -dv(1)*dwdx(1)
     &                                -dv(3)*dwdx(3)
            hzz = 2.0d0*dv(3)*dwdx(3) -dv(1)*dwdx(1)
     &                                -dv(2)*dwdx(2)
            hyz = dv(2)*dwdx(3) +dv(3)*dwdx(2)
            hxz = dv(1)*dwdx(3) +dv(3)*dwdx(1)
            hxy = dv(1)*dwdx(2) +dv(2)*dwdx(1)
          endif
	    
          hxx = 2.0d0/3.0d0 * hxx
          hyy = 2.0d0/3.0d0 * hyy
          hzz = 2.0d0/3.0d0 * hzz
	    
          if (dim .eq. 1) then 
            txx(i) = txx(i) +mass(j)*hxx/rho(j)
	      
            txx(j) = txx(j) +mass(i)*hxx/rho(i)
          else if (dim .eq. 2) then
            txx(i) = txx(i) +mass(j)*hxx/rho(j)
            tyy(i) = tyy(i) +mass(j)*hyy/rho(j)
            txy(i) = txy(i) +mass(j)*hxy/rho(j)
	      
            txx(j) = txx(j) +mass(i)*hxx/rho(i)
            tyy(j) = tyy(j) +mass(i)*hyy/rho(i)
            txy(j) = txy(j) +mass(i)*hxy/rho(i)
          else if (dim .eq. 3) then
            txx(i) = txx(i) + mass(j)*hxx/rho(j)
            tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
            tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
            tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
            txz(i) = txz(i) + mass(j)*hxz/rho(j)
            txy(i) = txy(i) + mass(j)*hxy/rho(j)
	      
            txx(j) = txx(j) + mass(i)*hxx/rho(i)   
            tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
            tzz(j) = tzz(j) + mass(i)*hzz/rho(i)
            tyz(j) = tyz(j) + mass(i)*hyz/rho(i)
            txz(j) = txz(j) + mass(i)*hxz/rho(i)
            txy(j) = txy(j) + mass(i)*hxy/rho(i)
          endif
	    
c     Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz:
	    
          hvcc = 0.0d0
          do d = 1, dim
            hvcc = hvcc + dv(d)*dwdx(d)
          enddo
          vcc(i) = vcc(i) +mass(j)*hvcc/rho(j)
          vcc(j) = vcc(j) +mass(i)*hvcc/rho(i)
        enddo
      endif
	
      do i = 1, nt
      
c     Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab

        if (visc) then
          if (dim .eq. 1) then
	      tdsdt(i) = txx(i)*txx(i)                             
          else if (dim .eq. 2) then
            tdsdt(i) = txx(i)*txx(i) + tyy(i)*tyy(i)
     &               + 2.0d00*txy(i)*txy(i)
          else if (dim .eq. 3) then
            tdsdt(i) = txx(i)*txx(i) + tyy(i)*tyy(i) + tzz(i)*tzz(i)
     &               + 2.0d0*tyz(i)*tyz(i) + 2.0d0*txz(i)*txz(i) 
     &               + 2.0d0*txy(i)*txy(i)
          endif
          tdsdt(i) = 0.5e0*eta(i)/rho(i)*tdsdt(i)
        endif
      enddo
	
	
c      Calculate SPH sum for pressure force -p,a/rho
c      and viscous force (eta Tab),b/rho
c      and the internal energy change de/dt due to -p/rho vc,c
	
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
        
c     For SPH algorithm 1
	  
        if(pa_sph .eq. 0) then
	    he = 0.0d0
		
c     DIMENSION = 1
c	
	    if (dim == 1) then
c     X component
c     Pressure part
            hv = -(p(i) + p(j))*dwdx(1)
            he = he + (v(1,j) - v(1,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1)
	      endif
	      
	      hv = hv / (rho(i)*rho(j))
	      dvdt(1,i) = dvdt(1,i) + mass(j)*hv
	      dvdt(1,j) = dvdt(1,j) - mass(i)*hv
	      
c	
c     DIMENSION = 2
c	
	    else if (dim == 2) then
c     X component
c     Pressure part
            hv = -(p(i) + p(j))*dwdx(1)
            he = he + (v(1,j) - v(1,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1)
     &                + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2)
	      endif
	      
	      hv = hv / (rho(i)*rho(j))
	      dvdt(1,i) = dvdt(1,i) + mass(j)*hv
	      dvdt(1,j) = dvdt(1,j) - mass(i)*hv
	      
c     Y component
c     Pressure part
            hv = -(p(i) + p(j))*dwdx(2)
            he = he + (v(2,j) - v(2,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(1)
     &                + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2)
	      endif
	      
	      hv = hv / (rho(i)*rho(j))
	      dvdt(2,i) = dvdt(2,i) + mass(j)*hv
	      dvdt(2,j) = dvdt(2,j) - mass(i)*hv
	      
c	
c     DIMENSION = 3
c	
	    else if (dim == 3) then
c     X component
c     Pressure part
            hv = -(p(i) + p(j))*dwdx(1)
            he = he + (v(1,j) - v(1,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1)
     &                + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2)
     &                + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3)
	      endif
	      
	      hv = hv / (rho(i)*rho(j))
	      dvdt(1,i) = dvdt(1,i) + mass(j)*hv
	      dvdt(1,j) = dvdt(1,j) - mass(i)*hv
	      
c     Y component
c     Pressure part
            hv = -(p(i) + p(j))*dwdx(2)
            he = he + (v(2,j) - v(2,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(1)
     &                + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2)
     &                + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3)
	      endif
	      
	      hv = hv / (rho(i)*rho(j))
	      dvdt(2,i) = dvdt(2,i) + mass(j)*hv
	      dvdt(2,j) = dvdt(2,j) - mass(i)*hv
	      
c     Z component
c     Pressure part
            hv = -(p(i) + p(j))*dwdx(3)
            he = he + (v(3,j) - v(3,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(1)
     &                + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(2)
     &                + (eta(i)*tzz(i) + eta(j)*tzz(j))*dwdx(3)
	      endif
	      
	      hv = hv / (rho(i)*rho(j))
	      dvdt(3,i) = dvdt(3,i) + mass(j)*hv
	      dvdt(3,j) = dvdt(3,j) - mass(i)*hv
	    
		endif
	  
c     ENERGY component
	    he = he / (rho(i)*rho(j))
	    dudt(i) = dudt(i) + mass(j)*he
	    dudt(j) = dudt(j) + mass(i)*he
	  
	  
c     For SPH algorithm 2
        
        else if (pa_sph .eq. 1) then
	    he = 0.0d0
		
c     DIMENSION = 1
c	
	    if (dim == 1) then
c     X component
c     Pressure part
            hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(1)
            he = he + (v(1,j) - v(1,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txx(i)/rho(i)**2 
     &                +  eta(j)*txx(j)/rho(j)**2)*dwdx(1)
	      endif
	      
	      dvdt(1,i) = dvdt(1,i) + mass(j)*hv
	      dvdt(1,j) = dvdt(1,j) - mass(i)*hv
	      
c	
c     DIMENSION = 2
c	
	    else if (dim == 2) then
c     X component
c     Pressure part
            hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(1)
            he = he + (v(1,j) - v(1,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txx(i)/rho(i)**2 
     &                 + eta(j)*txx(j)/rho(j)**2)*dwdx(1)
     &                + (eta(i)*txy(i)/rho(i)**2 
     &                 + eta(j)*txy(j)/rho(j)**2)*dwdx(2)
	      endif
	      
	      dvdt(1,i) = dvdt(1,i) + mass(j)*hv
	      dvdt(1,j) = dvdt(1,j) - mass(i)*hv
	      
c     Y component
c     Pressure part
            hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(2)
            he = he + (v(2,j) - v(2,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txy(i)/rho(i)**2 
     &                 + eta(j)*txy(j)/rho(j)**2)*dwdx(1)
     &                + (eta(i)*tyy(i)/rho(i)**2 
     &                 + eta(j)*tyy(j)/rho(j)**2)*dwdx(2)
	      endif
	      
	      dvdt(2,i) = dvdt(2,i) + mass(j)*hv
	      dvdt(2,j) = dvdt(2,j) - mass(i)*hv
	      
c	
c     DIMENSION = 3
c	
	    else if (dim == 3) then
c     X component
c     Pressure part
            hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(1)
            he = he + (v(1,j) - v(1,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txx(i)/rho(i)**2 
     &                 + eta(j)*txx(j)/rho(j)**2)*dwdx(1)
     &                + (eta(i)*txy(i)/rho(i)**2 
     &                 + eta(j)*txy(j)/rho(j)**2)*dwdx(2)
     &                + (eta(i)*txz(i)/rho(i)**2 
     &                 + eta(j)*txz(j)/rho(j)**2)*dwdx(3)
	      endif
	      
	      dvdt(1,i) = dvdt(1,i) + mass(j)*hv
	      dvdt(1,j) = dvdt(1,j) - mass(i)*hv
	      
c     Y component
c     Pressure part
            hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(2)
            he = he + (v(2,j) - v(2,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txy(i)/rho(i)**2 
     &                 + eta(j)*txy(j)/rho(j)**2)*dwdx(1)
     &                + (eta(i)*tyy(i)/rho(i)**2 
     &                 + eta(j)*tyy(j)/rho(j)**2)*dwdx(2)
     &                + (eta(i)*tyz(i)/rho(i)**2 
     &                 + eta(j)*tyz(j)/rho(j)**2)*dwdx(3)
	      endif
	      
	      dvdt(2,i) = dvdt(2,i) + mass(j)*hv
	      dvdt(2,j) = dvdt(2,j) - mass(i)*hv
	      
c     Z component
c     Pressure part
            hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(3)
            he = he + (v(3,j) - v(3,i))*hv
	      
c     Viscous force
            if (visc) then
              hv = hv + (eta(i)*txz(i)/rho(i)**2 
     &                 + eta(j)*txz(j)/rho(j)**2)*dwdx(1)
     &                + (eta(i)*tyz(i)/rho(i)**2 
     &                 + eta(j)*tyz(j)/rho(j)**2)*dwdx(2)
     &                + (eta(i)*tzz(i)/rho(i)**2 
     &                 + eta(j)*tzz(j)/rho(j)**2)*dwdx(3)
	      endif
	      
	      dvdt(3,i) = dvdt(3,i) + mass(j)*hv
	      dvdt(3,j) = dvdt(3,j) - mass(i)*hv
	    
		endif
	  
c     ENERGY component
	    dudt(i) = dudt(i) + mass(j)*he
	    dudt(j) = dudt(j) + mass(i)*he
	  endif
	enddo
	
c     Change of specific internal energy de/dt = T ds/dt - p/rho vc,c:

      do i = 1, nt
         dudt(i) = tdsdt(i) + 0.5d0*dudt(i)
      enddo
	
      end
