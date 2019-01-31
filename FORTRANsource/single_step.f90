subroutine particle_update(np, nv, x, v, mat, h, &
     mass, rho, p, u, &
     niac, pair_i, pair_j, &
     dvdt, drhodt, dudt, av, &
     prop, its, dt, &
     dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, &
     sle, skf, density_method, &
     average_velocity, virtual_part, ghost_part, visc, &
     visc_artificial, heat_artificial, &
     heat_external, self_gravity)


  implicit none
  include 'options.inc'

  integer np, nv, ng, its, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn), p(maxn), u(maxn)
  double precision dvdt(dim, maxn), drhodt(maxn), dudt(maxn), &
       av(dim, maxn)
  double precision prop(16,40), dt
  
  integer niac, pair_i(max_interaction), pair_j(max_interaction)
  
  double precision dtfluka
  double precision ofluka(3), dfluka(3), datafluka(max_fluka)
  integer nfluka(3), tfluka
  
  double precision c(maxn), eta(maxn)
  double precision dx(dim, maxn), exdvdt(dim, maxn), &
       divv(maxn), epsxx(maxn), epsyy(maxn), epszz(maxn), &
       epsxy(maxn), epsyz(maxn), epszx(maxn)
  integer sle, skf, density_method
  logical average_velocity, virtual_part, ghost_part, visc, &
       visc_artificial, heat_artificial, heat_external, self_gravity
  
  integer i, d
  
  do i = 1, np
     dudt(i) = 0.0d0
     drhodt(i) = 0.0d0
     
     do d = 1, dim
        dvdt(d,i) = 0.0d0
        av(d,i) = 0.0d0
     enddo
  enddo
  
!!$ Density approximation
  if (density_method .eq. 0) then
     call con_density(np, x, v, h, mass, drhodt, &
          niac, pair_i, pair_j, skf)
  elseif(density_method .eq. 1) then
     call sum_density(np, x, h, mass, rho, &
          niac, pair_i, pair_j, skf)
  elseif (density_method .eq. 2) then
     call nor_density(np, x, h, mass, rho, &
          niac, pair_i, pair_j, skf)
  endif
  
!!$ Equation of state
  call eos(np, mat, rho, u, p, c, eta, prop)
  
!!$ Interparticle forces
  call pressure_action(np, x, v, h, mass, rho, p, &
       dvdt, dudt, niac, pair_i, pair_j, skf)
  
  if (visc) then
     call viscous_action(np, x, v, h, mass, rho, eta, &
          dvdt, dudt, niac, pair_i, pair_j, skf)
  endif
  
  if (visc_artificial) then
     call shock_correction(np, x, v, h, mass, rho, c, &
          dvdt, dudt, niac, pair_i, pair_j, skf)
  endif
  
  if (heat_artificial) then
     call heat_correction(np, x, v, h, mass, rho, u, c, &
          dudt, niac, pair_i, pair_j, skf)
  endif
  
!!$ External forces:
  if (virtual_part) then
     do i = 1, np
        do d = 1, dim
           dx(d, i) = 0.0d0
        enddo
     enddo
     
     call contact_displacement(np, nv, x, v, h, dx, &
          niac, pair_i, pair_j)
     
     call smooth_vector(np, x, h, mass, rho, dx, exdvdt, &
          niac, pair_i, pair_j, skf)
     
     do i = 1, np
        do d = 1, dim
           exdvdt(d, i) = exdvdt(d, i) / (2.0d0*dt*dt)
        enddo
     enddo
     
!!$     call virtual_force(np, nv, x, v, h, mass, rho, exdvdt, &
!!$          niac, pair_i, pair_j, dt, skf)
     
     do i = 1, np
        do d = 1, dim
           dvdt(d, i) = dvdt(d, i) + exdvdt(d, i)
        enddo
     enddo
  end if
  
!!$ Calculate the external heating on the particles
  if ((heat_external) .and. (its*dt .le. dtfluka)) then
     call fluka_heat(np, x, rho, dudt, &
          ofluka, dfluka, nfluka, tfluka, datafluka)
  endif

!!$ Calculating average velocity of each partile for avoiding penetration
  if (average_velocity) call av_vel(np, x, v, h, mass, rho, av, &
       niac, pair_i, pair_j, skf)

!!$ Calculating the neighboring particles and undating H
  if (sle .ne. 0) call h_upgrade(np, x, v, h, mass, rho, &
       niac, pair_i, pair_j, dt, sle, skf)

!!$ Convert velocity, force, and energy to f and dfdt  
!!$  do i = 1, np
!!$     do d = 1, dim
!!$        dvdt(d,i) = indvdt(d,i) + exdvdt(d,i) + avdvdt(d,i)
!!$     enddo
!!$     du(i) = indudt(i) + avdudt(i) + ahdudt(i) + exdudt(i)
!!$  enddo
  
end subroutine particle_update

