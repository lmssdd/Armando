      subroutine art_heat(nt, x, v, h, mass, rho, u, c, dudt,
     &           niac, pair_i, pair_j, skf)

c----------------------------------------------------------------------
c     Subroutine to calculate the artificial heat(Noh, 1978) 

c     nt-- Number of particles (including virtual particles        [in]
c     x-- article position                                         [in]
c     v-- particle velocity                                        [in]
c     h-- smoothing length                                         [in]
c     mass-- particle mass                                         [in]
c     rho-- density                                                [in]
c     u-- internal energy                                          [in]
c     c-- sound speed                                              [in]
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
     &       mass(maxn), rho(maxn), u(maxn), c(maxn), dvdt(dim, maxn),
     &       dudt(maxn)
	
      integer i,j,k,d
      double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
      double precision divv(maxn), dx, dv, g1, g2, rdwdx, vdwdx,
     &       psi_i, psi_j, psi_ij, psi, eta, mrho
     
	
c---  Parameter for the artificial heat conduction:
	
      g1 = 0.1d0
      g2 = 1.0d0
	eta = 0.1d0
	
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
	    vdwdx = vdwdx + dv*dwdx(d)
	  enddo
	  
	  divv(i) = divv(i) + mass(j) / rho(i) * vdwdx
	  divv(j) = divv(j) + mass(i) / rho(j) * vdwdx
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
        
        mrho = 0.5d0 * (rho(i) + rho(j))
	  
        rdwdx = 0.0d0
        do d = 1, dim
          rdwdx  = rdwdx + dxiac(d)*dwdx(d)
        enddo
	  
	  psi_i = g1*h(i)*c(i) + g2*h(i)*h(i)*(abs(divv(i)) -divv(i))
	  psi_j = g1*h(j)*c(j) + g2*h(j)*h(j)*(abs(divv(j)) -divv(j)) 
	  psi_ij = 0.5d0 * (psi_i + psi_j)
	  
	  psi = (2.0d0*psi_ij) / (mrho*(dr2iac + (mh*eta)**2)) * rdwdx
	  
	  dudt(i) = dudt(i) + mass(j) * psi * (u(i) - u(j))
	  dudt(j) = dudt(j) + mass(i) * psi * (u(j) - u(i))
	enddo
	
      end
