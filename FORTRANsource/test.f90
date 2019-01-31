subroutine test(np, nv, x, v, mat, h, &
     mass, rho, p, eta, dvdt, dudt, tdsdt, &
     niac, pair_i, pair_j, dt, skf, pa_sph, visc)

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn), &
       niac, pair_i(max_interaction), pair_j(max_interaction)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn), p(maxn), u(maxn), &
       c(maxn), eta(maxn), dt
  double precision dvdt(dim, maxn), dudt(maxn), tdsdt(maxn)
  integer skf, pa_sph
  logical visc

  integer i, j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, mh
  double precision  hv,wm

!!$ Initialization of shear tensor, velocity divergence, 
!!$ viscous energy, internal energy, acceleration 

  do i = 1, np
     do d = 1 ,dim
        dvdt(d,i) = 0.0d0
     enddo
  enddo

!!$ Calculate SPH sum for pressure force -p,a/rho
!!$ and viscous force (eta Tab),b/rho
!!$ and the internal energy change de/dt due to -p/rho vc,c

  do i = 1, np
     dr2iac = 0.0
     do d = 1, dim
        dxiac(d) = x(d,i) - (-0.51)
        dr2iac = dr2iac + dxiac(d)*dxiac(d)
     enddo
     driac = sqrt(dr2iac)
     call kernel(0.5*h(i), dxiac, 0.5*h(i), wm, dwdx, skf)
     call kernel(driac, dxiac, 0.5*h(i), w, dwdx, skf)
     
     if (dim == 1) then
        hv = -1.0d3 *4*w**3/wm**4 * dwdx(1)
        
        dvdt(1,i) = dvdt(1,i) + mass(i)/rho(i)*hv
     endif
  enddo
  
end Subroutine test
