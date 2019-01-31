subroutine eos(np, mat, rho, u, p, c, eta, prop)
!!$
!!$----------------------------------------------------------------------
!!$ Equation Of State computation
!!$
!!$ np     : Particle number                                      [in]
!!$ mat    : Material number                                      [in]
!!$ rho    : Density                                              [in]
!!$ u      : Energy                                               [in]
!!$ p      : Pressure                                            [out]
!!$ c      : Speed of sound                                      [out]
!!$ prop   : Material properties                                  [in]

  implicit none
  include 'options.inc'
  
  integer np, mat(maxn)
  double precision prop(16,40)
  double precision rho(maxn), u(maxn), p(maxn), c(maxn), eta(maxn)
  double precision gamma, pshift, rho0, eta0, pmin, &
       a1, a2, a3, b0, b1, t1, t2, &
       c0, gamma0, s0, h0, es
  integer i

  do i = 1, np
     if ((abs(mat(i)) .gt. 0) .and. (abs(mat(i)) .le. 10)) then
!!$ Ideal gas EOS
        rho0 = prop(1, abs(mat(i)))
        gamma = prop(2, abs(mat(i)))
        pshift = prop(3, abs(mat(i)))
        eta0 = prop(4, abs(mat(i)))

        call eos_gas(rho(i), u(i), p(i), c(i), gamma)

        p(i) = p(i) +pshift
        eta(i)= eta0

     else if ((abs(mat(i)) .gt. 10) .and. (abs(mat(i)) .le. 20)) then
!!$ Mie-Gruneise polyynomial EOS
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

        call eos_poly(rho(i), u(i), p(i), c(i), &
             rho0, a1, a2, a3, b0, b1, t1, t2, pmin)

!!$		call p_art_water(rho(i), p(i), c(i))

        eta(i)= eta0

     else if ((abs(mat(i)) .gt. 20) .and. (abs(mat(i)) .le. 30)) then
!!$ Mie-Gruneise shock EOS
        rho0 = prop(1, abs(mat(i)))
        c0 = prop(2, abs(mat(i)))
        gamma0 = prop(3, abs(mat(i)))
        s0 = prop(4, abs(mat(i)))
        pmin = prop(5, abs(mat(i)))
        eta0 = prop(6, abs(mat(i)))
        
        call eos_shock(rho(i), u(i), p(i), c(i), &
             rho0, c0, gamma0, s0, pmin)

        eta(i)= eta0

     else if ((abs(mat(i)) .gt. 30) .and. (abs(mat(i)) .le. 40)) then
!!$ Mie-Gruneise Puff EOS
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

        call eos_puff(rho(i), u(i), p(i), c(i), &
             rho0, a1, a2, a3, gamma0, h0, es, t1, t2, pmin)

        eta(i)= eta0
        
     else if ((abs(mat(i)) .gt. 40) .and. (abs(mat(i)) .le. 50)) then
!!$ Tait EOS
        rho0 = prop(1, abs(mat(i)))
        c0 = prop(2, abs(mat(i)))
        pmin = prop(4, abs(mat(i)))
        eta0 = prop(5, abs(mat(i)))

        call eos_tait(rho(i), u(i), p(i), c(i), &
             rho0, c0, pmin)

        eta(i)= eta0
     endif
  enddo
  
end subroutine eos


subroutine eos_puff(rho, e, p, c, &
     rho0, a1, a2, a3, gamma0, h0, es, t1, t2, pmin)
!!$
!!$----------------------------------------------------------------------
!!$ Mie-Gruneisen Puff EOS
!!$
!!$ rho    : Density                                              [in]
!!$ e      : Energy                                               [in]
!!$ p      : Pressure                                            [out]
!!$ rho0   : Material density                                     [in]
!!$ a      : Region I coefficients                                [in]
!!$ gamma0 : Gruneisen coefficient                                [in]
!!$ h0     : Expansion coefficient                                [in]
!!$ es     : Sublimation energy                                   [in]
!!$ t      : Region II coefficients                               [in]
!!$ pmin   : Minimum pressure                                     [in]

  implicit none
  double precision rho, e, p, c
  double precision rho0, a1, a2, a3, gamma0, h0, es, t1, t2, pmin
  double precision mu, eta, gamma

  eta = rho / rho0
  mu = rho / rho0 -1.0
  gamma = gamma0 * rho0 / rho

  if (mu .ge. 0.0) then
     p = (a1*mu +a2*mu**2 +a3*mu**3)*(1 -0.5*gamma*mu) &
          + gamma*rho*e
  else if ((mu .lt. 0.0) .and. (e .lt. es)) then
     p = (t1*mu +t2*mu**2)*(1 -0.5*gamma*mu) &
          + gamma*rho*e
  else if ((mu .lt. 0.0) .and. (e .lt. es)) then
     p = rho*(h0 +(gamma0 -h0)*eta**0.5) &
          * (e -es*(1.0 -dexp((a1*(eta -1.0))/(rho0*gamma*es*eta**2))))
  endif

  c = dsqrt(a1 / rho)
  
  if (p < pmin) then
     p = 0.0d0
     rho = rho0
  endif
  
end subroutine eos_puff


subroutine eos_shock(rho, e, p, c, &
     rho0, c0, gamma0, s0, pmin)

!!$----------------------------------------------------------------------
!!$ Mie-Gruneisen Shock EOS
!!$
!!$ rho    : Density                                              [in]
!!$ e      : Energy                                               [in]
!!$ p      : Pressure                                            [out]
!!$ rho0   : Material density                                     [in]
!!$ c0     : Material velocity                                    [in]
!!$ gamma0 : Gruneisen coefficient                                [in]
!!$ s0     : Gruneisen parameter                                  [in]
!!$ pmin   : Minimum pressure                                     [in]

  implicit none
  double precision rho, e, p, c
  double precision rho0, c0, gamma0, s0, pmin
  double precision mu, ph, eh, gamma

  mu = rho / rho0 -1.0
  gamma = gamma0 * rho0 / rho

  ph = (rho0 * c0 * c0 * mu * (1.0 +mu)) &
       / (1.0 - (s0 -1.0) * mu)**2
  eh = 0.5 * ph / rho0 * (mu / (1.0 +mu))

  p = ph + gamma * rho * (e -eh)
  
  c = c0
  
  if (p < pmin)then
     p = 0.0d0
     rho = rho0
  endif
  
end subroutine eos_shock


subroutine eos_poly(rho, e, p, c, &
     rho0, a1, a2, a3, b0, b1, t1, t2, pmin)

!!$-------------------------------------------------------------------
!!$ Mie-Gruneisen polynomial EOS
!!$
!!$ p = a1 mu + a2 mu2 + a3 mu3 + (b0 + b1 mu) rho0 e   in compression
!!$ p = t1 mu + t2 mu2 + b0 rho0 e                      in tension
!!$ rho    : Density                                              [in]
!!$ e      : Energy                                               [in]
!!$ p      : Pressure                                            [out]
!!$ rho0   : Material density                                     [in]
!!$ a      : Compression coefficients                             [in]
!!$ t      : Tension coefficients                                 [in]
!!$ b      : Energy coefficients                                  [in]
!!$ pmin   : Minimum pressure                                     [in]

  implicit none
  double precision rho, e, p, c
  double precision rho0, a1, a2, a3, b0, b1, t1, t2, pmin
  double precision mu

  mu = rho / rho0 -1.0

  if (mu < 0) then
     p = (t1*mu + t2*mu**2) &
          + (b0 * rho0*e)
  else
     p = (a1*mu + a2*mu**2 + a3*mu**3) &
          + ((b0 + b1*mu) * rho0*e)
  endif

  c = dsqrt(a1 / rho)
  
  if (p < pmin) then
     p = 0.0d0
     rho = rho0
  endif
  
end subroutine eos_poly

subroutine eos_gas(rho, e, p, c, gamma)

!!$-------------------------------------------------------------------
!!$ Ideal gas EOS
!!$
!!$ rho    : Density                                              [in]
!!$ e      : Energy                                               [in]
!!$ p      : Pressure                                            [out]
!!$ gamma  : Adiabatic exponent gamma = cp/cv                     [in]

  implicit none
  double precision rho, e, p, c
  double precision gamma

  p = (gamma -1.0) * rho * e
  c = dsqrt((gamma -1.0) * e)

end subroutine eos_gas


subroutine eos_tait(rho, e, p, c, &
     rho0, c0, pmin)

!!$----------------------------------------------------------------------
!!$ Tait EOS for water
!!$
!!$ rho    : Density                                              [in]
!!$ e      : Energy                                               [in]
!!$ p      : Pressure                                            [out]
!!$ rho0   : Material density                                     [in]
!!$ c0     : Material velocity                                    [in]
!!$ gamma0 : Gruneisen coefficient                                [in]
!!$ pmin   : Minimum pressure                                     [in]

  implicit none
  double precision rho, e, p, c
  double precision rho0, c0, pmin
  double precision mu, ph, eh, gamma
  
  p = rho0 * c0 * c0 / 7.0d0 * ((rho / rho0)**7 -1.0d0)
  
  c = c0
  
  if (p < pmin) then
     p = 0.0d0
     rho = rho0
  endif
  
end subroutine eos_tait


