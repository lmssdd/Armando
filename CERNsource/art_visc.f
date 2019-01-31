      subroutine art_visc(nt, x, v, h, mass, rho, c, dvdt, dudt,
     &           niac, pair_i, pair_j, skf)

c----------------------------------------------------------------------
c     Subroutine to calculate the artificial viscosity (Monaghan, 1992) 

c     nt-- Number of particles (including virtual particles        [in]
c     x-- article position                                         [in]
c     v-- particle velocity                                        [in]
c     h-- smoothing length                                         [in]
c     mass-- particle mass                                         [in]
c     rho-- density                                                [in]
c     c-- sound speed                                              [in]
c     dvdt-- = dv/dt                                              [out]
c     dudt-- = du/dt                                              [out]
c     niac-- number of interaction pairs                           [in]
c     pair_i-- list of first partner of interaction pair           [in]
c     pair_j-- list of second partner of interaction pair          [in]
c     skf -- Smoothing kernel function                             [in]
	
	
      implicit none
      include 'options.inc'
      
      integer nt, niac, pair_i(max_interaction), pair_j(max_interaction)
      integer skf
      double precision x(dim, maxn), v(dim, maxn), h(maxn), 
     &       mass(maxn), rho(maxn), c(maxn), dvdt(dim, maxn), 
     &       dudt(maxn)
	
      integer i,j,k,d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision dv, alpha, beta, eta, piv, phi,
     &       vr, vdwdx, mc, mrho
     
c     Parameter for the artificial viscosity:
c     Shear viscosity
      parameter(alpha = 1.0d0)
	
c     Bulk viscosity
      parameter(beta = 1.0d0)
      
c     Parameter to avoid singularities
      parameter(eta = 0.1d0 )
	
c      do i = 1, nt
c        do d = 1, dim
c          dvdt(d,i) = 0.0d0
c        enddo
c        dudt(i) = 0.0d00
c      enddo   
	
c     Calculate SPH sum for artificial viscosity
      
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
        
        vr = 0.0d0
	  vdwdx = 0.0d0
        do d = 1, dim
          dv = v(d,i) - v(d,j)
          vr = vr + dv * dxiac(d)
	    vdwdx = vdwdx + dv*dwdx(d)
        enddo
	  
c     Artificial viscous force only if v_ij * r_ij < 0
	  
        if (vr < 0.0d0) then
	    
c     Calculate phi_ij = h v_ij * r_ij / ( r_ij^2 + h^2 eta^2 )
            
          phi = (mh * vr) / (dr2iac + (mh*eta)**2)
          
c     Calculate PI_ij = (-alpha c_ij phi_ij + beta phi_ij^2) / rho_ij
	    
          mc = 0.5d0 * (c(i) + c(j))
          mrho = 0.5d0 * (rho(i) + rho(j))
          piv = (beta*phi - alpha*mc) * phi/mrho
	    
c     Calculate SPH sum for artificial viscous force 
c	and specific internal energy:
	    
          do d = 1, dim
            dvdt(d, i) = dvdt(d, i) - mass(j) * piv * dwdx(d)
            dvdt(d, j) = dvdt(d, j) + mass(i) * piv * dwdx(d)
          enddo
          
          dudt(i) = dudt(i) + 0.5d0 * mass(j) * piv * vdwdx
          dudt(j) = dudt(j) + 0.5d0 * mass(i) * piv * vdwdx
        endif
      enddo
	
      end
