subroutine sum_density(np, x, h, mass, rho, niac, pair_i, pair_j, &
     skf)

!!$----------------------------------------------------------------
!!$ Subroutine to calculate the density with SPH summation algorithm.
!!$
!!$ np-- number of particles                                     [in]
!!$ x-- article position                                         [in]
!!$ h-- smoothing length                                         [in]
!!$ mass-- particle masses                                       [in]
!!$ rho-- density                                               [out]
!!$ niac-- number of interaction pairs                           [in]
!!$ pair_i-- list of first partner of interaction pair           [in]
!!$ pair_j-- list of second partner of interaction pair          [in]
!!$ skf-- smoothing kernel function                              [in]

  implicit none
  include 'options.inc'

  integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
  integer skf
  double precision x(dim, maxn), h(maxn), mass(maxn), rho(maxn)

  integer i, j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
  
  
  do d = 1, dim
     dxiac(d) = 0.0d0
  enddo
  driac = 0.

!!$ Calculate the rho integration over the space

  do i = 1, np
     call kernel(driac, dxiac, h(i), w, dwdx, skf)   
     rho(i) = mass(i)*w
  enddo
  
!!$ Calculate SPH sum for rho:

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

end subroutine sum_density


subroutine nor_density(np, x, h, mass, rho, niac, pair_i, pair_j, &
     skf)

!!$------------------------------------------------------------------
!!$ Subroutine to calculate the density with SPH summation algorithm
!!$ with density normalization to avoid the boundary deficiency
!!$ problem (Randles and Libersky 1996).
!!$
!!$ np-- number of particles                                     [in]
!!$ x-- article position                                         [in]
!!$ h-- smoothing length                                         [in]
!!$ mass-- particle masses                                       [in]
!!$ rho-- density                                            [in/out]
!!$ niac-- number of interaction pairs                           [in]
!!$ pair_i-- list of first partner of interaction pair           [in]
!!$ pair_j-- list of second partner of interaction pair          [in]
!!$ skf-- smoothing kernel function                              [in]

  implicit none
  include 'options.inc'

  integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
  integer skf

  double precision x(dim, maxn), h(maxn), mass(maxn), rho(maxn)

  integer i, j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
  double precision rn(maxn), wn(maxn)


!!$ Calculate SPH sum for rho:
  
  do d = 1, dim
     dxiac(d) = 0.0d0
  enddo
  driac = 0.

  do i = 1, np
     call kernel(driac, dxiac, h(i), w, dwdx, skf)   
     rn(i) = mass(i)*w
     wn(i) = mass(i)/rho(i)*w
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
     
     if (i <= np) rn(i) = rn(i) + mass(j)*w
     if (j <= np) rn(j) = rn(j) + mass(i)*w
     
     if (i <= np) wn(i) = wn(i) + mass(j)/rho(j)*w
     if (j <= np) wn(j) = wn(j) + mass(i)/rho(i)*w
  enddo

!!$ Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
  
  do i = 1, np
     rho(i) = rn(i) / wn(i)
  enddo

end subroutine nor_density



subroutine con_density(np, x, v, h, mass, drho, &
     niac, pair_i, pair_j, skf)

!!$------------------------------------------------------------------
!!$ Subroutine to calculate the density variation based on the 
!!$ continuity equation.
!!$
!!$ np-- number of particles                                     [in]
!!$ x-- article position                                         [in]
!!$ v-- particle velocities                                      [in]
!!$ h-- smoothing length                                         [in]
!!$ mass-- particle masses                                       [in]
!!$ drho-- density change rate of each particle                 [out]
!!$ niac-- number of interaction pairs                           [in]
!!$ pair_i-- list of first partner of interaction pair           [in]
!!$ pair_j-- list of second partner of interaction pair          [in]
!!$ skf-- smoothing kernel function                              [in]

  implicit none
  include 'options.inc'

  integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
  integer skf

  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), drho(maxn)

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

end subroutine con_density

subroutine unity(np, x, h, mass, rho, wn, niac, pair_i, pair_j, &
     skf)

!!$------------------------------------------------------------------
!!$ Subroutine to calculate the unity with SPH summation algorithm
!!$ with density normalization to avoid the boundary deficiency
!!$ problem (Randles and Libersky 1996).
!!$
!!$ np-- number of particles                                     [in]
!!$ x-- article position                                         [in]
!!$ h-- smoothing length                                         [in]
!!$ mass-- particle masses                                       [in]
!!$ rho-- density                                            [in/out]
!!$ niac-- number of interaction pairs                           [in]
!!$ pair_i-- list of first partner of interaction pair           [in]
!!$ pair_j-- list of second partner of interaction pair          [in]
!!$ skf-- smoothing kernel function                              [in]

  implicit none
  include 'options.inc'

  integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
  integer skf

  double precision x(dim, maxn), h(maxn), mass(maxn), rho(maxn), wn(maxn)

  integer i, j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh


!!$ Calculate SPH sum for rho:
  
  do d = 1, dim
     dxiac(d) = 0.0d0
  enddo
  driac = 0.

  do i = 1, np
     call kernel(driac, dxiac, h(i), w, dwdx, skf)   
     wn(i) = mass(i)/rho(i)*w
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .le. np)) then
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        wn(i) = wn(i) + mass(j)/rho(j)*w
        wn(j) = wn(j) + mass(i)/rho(i)*w
     endif
  enddo

end subroutine unity

