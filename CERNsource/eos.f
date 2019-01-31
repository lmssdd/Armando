	subroutine eos_puff(rho, e, p, 
     &               rho0, a1, a2, a3, gamma0, h0, es, t1, t2, pmin)
      
c----------------------------------------------------------------------
c     Mie-Gruneisen Puff EOS

c     rho    : Density                                              [in]
c     e      : Energy                                               [in]
c     p      : Pressure                                            [out]
c     rho0   : Material density                                     [in]
c     a      : Region I coefficients                                [in]
c     gamma0 : Gruneisen coefficient                                [in]
c     h0     : Expansion coefficient                                [in]
c     es     : Sublimation energy                                   [in]
c     t      : Region II coefficients                               [in]
c     pmin   : Minimum pressure                                     [in]
	
      implicit none
      double precision rho, e, p
      double precision rho0, a1, a2, a3, gamma0, h0, es, t1, t2, pmin
	double precision mu, eta, gamma

	eta = rho / rho0
	mu = rho / rho0 -1.0
	gamma = gamma0 * rho0 / rho

	if (mu .ge. 0.0) then
	  p = (a1*mu +a2*mu**2 +a3*mu**3)*(1 -0.5*gamma*mu) 
     &    + gamma*rho*e
	else if ((mu .lt. 0.0) .and. (e .lt. es)) then
	  p = (t1*mu +t2*mu**2)*(1 -0.5*gamma*mu) 
     &    + gamma*rho*e
	else if ((mu .lt. 0.0) .and. (e .lt. es)) then
	  p = rho*(h0 +(gamma0 -h0)*eta**0.5) 
     &    * (e -es*(1.0 -dexp((a1*(eta -1.0))/(rho0*gamma*es*eta**2))))
	endif
	
	if (p < pmin) p = 0.0d0
	
      end
	

	subroutine eos_shock(rho, e, p, 
     &                     rho0, c0, gamma0, s0, pmin)
      
c----------------------------------------------------------------------
c     Mie-Gruneisen Shock EOS

c     rho    : Density                                              [in]
c     e      : Energy                                               [in]
c     p      : Pressure                                            [out]
c     rho0   : Material density                                     [in]
c     c0     : Material velocity                                    [in]
c     gamma0 : Gruneisen coefficient                                [in]
c     s0     : Gruneisen parameter                                  [in]
c     pmin   : Minimum pressure                                     [in]

      implicit none
      double precision rho, e, p
      double precision rho0, c0, gamma0, s0, pmin
	double precision mu, ph, eh, gamma

	mu = rho / rho0 -1.0
	gamma = gamma0 * rho0 / rho
	
	ph = (rho0 * c0 * c0 * mu * (1.0 +mu)) / (1.0 - (s0 -1.0) * mu)**2
	eh = 0.5 * ph / rho0 * (mu / (1.0 +mu))
	
	p = ph + gamma * rho * (e -eh)
	
	if (p < pmin) p = 0.0d0
	
      end
	

	subroutine eos_poly(rho, e, p, 
     &                    rho0, a1, a2, a3, b0, b1, t1, t2, pmin)
      
c----------------------------------------------------------------------
c     Mie-Gruneisen polynomial EOS

c     p = a1 mu + a2 mu2 + a3 mu3 + (b0 + b1 mu) rho0 e   in compression
c     p = t1 mu + t2 mu2 + b0 rho0 e                      in tension
c     rho    : Density                                              [in]
c     e      : Energy                                               [in]
c     p      : Pressure                                            [out]
c     rho0   : Material density                                     [in]
c     a      : Compression coefficients                             [in]
c     t      : Tension coefficients                                 [in]
c     b      : Energy coefficients                                  [in]
c     pmin   : Minimum pressure                                     [in]

      implicit none
      double precision rho, e, p
      double precision rho0, a1, a2, a3, b0, b1, t1, t2, pmin
	double precision mu

	mu = rho / rho0 -1.0

	if (mu < 0) then
	  p = (t1*mu + t2*mu**2)
     &    + (b0 * rho0*e)
	else
	  p = (a1*mu + a2*mu**2 + a3*mu**3)
     &    + ((b0 + b1*mu) * rho0*e)
	endif

	if (p < pmin) p = 0.0d0
	
      end

	subroutine eos_gas(rho, e, p, 
     &                   gamma)
      
c----------------------------------------------------------------------
c     Ideal gas EOS

c     rho    : Density                                              [in]
c     e      : Energy                                               [in]
c     p      : Pressure                                            [out]
c     gamma  : Adiabatic exponent gamma = cp/cv                     [in]

      implicit none
      double precision rho, e, p
      double precision gamma
	
	p = (gamma -1.0) * rho * e

      end
