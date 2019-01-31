subroutine art_visc(nt, x, v, h, mass, rho, c, dvdt, dudt, & 
     niac, pair_i, pair_j, skf)

!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
!!$
!!$ nt-- Number of particles (including virtual particles        [in]
!!$ x-- article position                                         [in]
!!$ v-- particle velocity                                        [in]
!!$ h-- smoothing length                                         [in]
!!$ mass-- particle mass                                         [in]
!!$ rho-- density                                                [in]
!!$ c-- sound speed                                              [in]
!!$ dvdt-- = dv/dt                                              [out]
!!$ dudt-- = du/dt                                              [out]
!!$ niac-- number of interaction pairs                           [in]
!!$ pair_i-- list of first partner of interaction pair           [in]
!!$ pair_j-- list of second partner of interaction pair          [in]
!!$ skf -- Smoothing kernel function                             [in]


  implicit none
  include 'options.inc'

  integer nt, niac, pair_i(max_interaction), pair_j(max_interaction)
  integer skf
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn), c(maxn), dvdt(dim, maxn), &
       dudt(maxn)

  integer i,j,k,d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
  double precision dv(dim), alpha, beta, eta, piv, phi, &
       vr, mc, mrho

!!$ Parameter for the artificial viscosity:
!!$ Shear viscosity
  parameter(alpha = 1.0)

!!$ Bulk viscosity
  parameter(beta = 1.0)

!!$ Parameter to avoid singularities
  parameter(eta = 0.1d0 )

!!$      do i = 1, nt
!!$        do d = 1, dim
!!$          dvdt(d,i) = 0.0d0
!!$        enddo
!!$        dudt(i) = 0.0d00
!!$      enddo   

!!$ Calculate SPH sum for artificial viscosity

  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     dr2iac = 0.0
     vr = 0.0d0
     do d = 1, dim
        dxiac(d) = x(d,i) - x(d,j)
        dr2iac = dr2iac + dxiac(d)*dxiac(d)
        dv(d) = v(d,i) - v(d,j)
        vr = vr + dv(d) * dxiac(d)
     enddo
     mh = 0.5*(h(i) +h(j))
     driac = sqrt(dr2iac)
     call kernel(driac, dxiac, mh, w, dwdx, skf)

!!$ Artificial viscous force only if v_ij * r_ij < 0

     if (vr < 0.0d0) then

!!$ Calculate phi_ij = h v_ij * r_ij / ( r_ij^2 + h^2 eta^2 )

        phi = (mh * vr) / (dr2iac + mh*mh*eta*eta)

!!$ Calculate PI_ij = (-alpha c_ij phi_ij + beta phi_ij^2) / rho_ij

        mc = 0.5d0 * (c(i) + c(j))
        mrho = 0.5d0 * (rho(i) + rho(j))
        piv = (beta*phi - alpha*mc) * phi/mrho

!!$ Calculate SPH sum for artificial viscous force 
!!$ and specific internal energy:
        
        do d = 1, dim
           dvdt(d, i) = dvdt(d, i) - mass(j) * piv * dwdx(d)
           dvdt(d, j) = dvdt(d, j) + mass(i) * piv * dwdx(d)

           dudt(i) = dudt(i) + 0.5d0 * mass(j) * dv(d) * piv * dwdx(d)
           dudt(j) = dudt(j) + 0.5d0 * mass(i) * dv(d) * piv * dwdx(d)
        enddo
     endif
  enddo

end subroutine art_visc
