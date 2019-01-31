subroutine av_vel(np, x, v, h, mass, rho, av, &
     niac, pair_i, pair_j, skf)

!!$------------------------------------------------------------------
!!$ Subroutine to calculate the average velocity to correct velocity
!!$ for preventing.penetration (monaghan, 1992)
!!$
!!$ np-- number of particles                                     [in]
!!$ x-- article position                                         [in]
!!$ v-- particle velocity                                        [in]
!!$ h-- smoothing length                                         [in]
!!$ mass-- particle mass                                         [in]
!!$ rho-- density                                                [in]
!!$ av-- average velocity                                       [out]
!!$ niac-- number of interaction pairs                           [in]
!!$ pair_i-- list of first partner of interaction pair           [in]
!!$ pair_j-- list of second partner of interaction pair          [in]
!!$ skf -- Smoothing kernel function                             [in]

  implicit none
  include 'options.inc'

  integer np, niac, pair_i(max_interaction), pair_j(max_interaction)
  integer skf
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn), u(maxn), av(dim, maxn)

  integer i,j,k,d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
  double precision dv, eps, mrho

!!$ epsilon --- a small constants chosen by experence, may lead to instability.
!!$ for example, for the 1 dimensional shock tube problem, the E <= 0.3

  eps = 0.3d0

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

     mrho = (rho(i) + rho(j)) / 2.0d0
     do d = 1, dim
        dv = v(d,i) - v(d,j)
        av(d,i) = av(d,i) - eps * mass(j) / mrho * dv * w
        av(d,j) = av(d,j) + eps * mass(i) / mrho * dv * w
     enddo
  enddo

end subroutine av_vel
