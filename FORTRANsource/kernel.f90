subroutine kernel(r, dx, h, w, dwdx, skf)   

!!$--------------------------------------------------------------------
!!$ Subroutine to calculate the smoothing kernel wij and its 
!!$ derivatives dwdxij.
!!$ if skf = 0, cubic spline kernel by W4 - Spline (Monaghan 1985)
!!$        = 1, Gauss kernel   (Gingold and Monaghan 1981) 
!!$        = 2, Quintic kernel (Morris 1997)
!!$
!!$ r-- distance between particles i and j                       [in]
!!$ dx-- x-, y- and z-distance between i and j                   [in]
!!$ h-- smoothing length                                         [in]
!!$ w-- kernel for all interaction pairs                        [out]
!!$ dwdx-- derivative of kernel with respect to x, y and z      [out]
!!$ skf -- smoothing kernel function                             [in]

  implicit none
  include 'options.inc'

  double precision r, dx(dim), h, w, dwdx(dim)
  integer skf, i, j, d
  double precision q, dw, factor

  q = r / h
  w = 0.0d0
  do d = 1, dim
     dwdx(d) = 0.0d0
  enddo

  if (skf .eq. 0) then
     if (dim .eq. 1) then
        factor = 1.0d0 / h
     else if (dim .eq. 2) then
        factor = 15.0d0 / (7.0d0 * pi * h**2)
     else if (dim .eq. 3) then
        factor = 3.0d0 / (2.0d0 * pi * h**3)
     endif

     if (q >= 0.0d0 .and. q <= 1.0d0) then
        w = factor * (2./3. - q*q + q**3 / 2.)
        do d = 1, dim
           dwdx(d) = factor * (-2.+3./2.*q) / h**2 * dx(d)
        enddo
     else if (q > 1.0d0 .and. q <= 2.0d0) then
        w = factor * 1.0d0 / 6.0d0 * (2.0d0 -q)**3 
        do d = 1, dim
           dwdx(d) =-factor * 1.0d0 / 6.0d0 * 3.*(2.-q)**2 / h &
                * (dx(d) / r)        
        enddo
     else
        w=0.0d0
        do d = 1, dim
           dwdx(d) = 0.0d0
        enddo
     endif

  else if (skf .eq. 1) then
     factor = 1.0d00 / (h**dim * pi**(dim/2.))
     if(q >= 0.0d0 .and. q <= 3.0d0) then
        w = factor * exp(-q*q)
        do d = 1, dim
           dwdx(d) = w * ( -2.* dx(d)/ (h**2))
        enddo
     else
        w = 0.0d0
        do d = 1, dim
           dwdx(d) = 0.0d0
        enddo
     endif

  else if (skf .eq. 2) then

     if (dim .eq. 1) then
        factor = 1.0d0 / (120.0d0 * h**1) !1.0d0 / (120.0d0 * pi * h**1)
     else if (dim .eq. 2) then
        factor = 7.0d0 / (478.0d0 * pi * h**2)
     elseif (dim.eq.3) then
        factor = 1.0d0 / (120.0d0 * pi * h**3)
     endif

     if(q >= 0.0d0 .and. q <= 1.0d0) then
        w = factor * ( (3 -q)**5 - 6*(2 -q)**5 + 15*(1 -q)**5 )
        do d = 1, dim
           dwdx(d) = factor * ( (-120 + 120*q - 50*q**2) &
                / h**2 * dx(d))
        enddo
     else if(q >= 1.0d0 .and. q <= 2.0d0) then
        w = factor * ( (3 -q)**5 - 6*(2 -q)**5 )
        do d = 1, dim
           dwdx(d) = factor * (-5*(3 -q)**4 + 30*(2 -q)**4) &
                / h * (dx(d)/r) 
        enddo
     else if(q >= 2.0d0 .and. q <= 3.0d0) then
        w = factor * (3 -q)**5
        do d = 1, dim
           dwdx(d) = factor * (-5*(3 -q)**4) / h * (dx(d)/r) 
        enddo
     else
        w = 0.0d0
        do d = 1, dim
           dwdx(d) = 0.0d0
        enddo
     endif
  endif

end subroutine kernel
