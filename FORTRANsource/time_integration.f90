subroutine time_integration(np, nv, x, v, mat, h, &
     mass, rho, p, u, prop, &
     dt, steps, refresh_step, save_step, &
     dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, &
     sle, skf, &
     average_velocity, virtual_part, visc, &
     visc_artificial, heat_artificial, &
     heat_external, gravity, filedir)



  implicit none     
  include 'options.inc'

  integer np, nv, mat(maxn), steps, refresh_step, save_step
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn), p(maxn), u(maxn), &
       prop(16,40), dt

  double precision dtfluka
  double precision ofluka(3), dfluka(3), datafluka(max_fluka)
  integer nfluka(3), tfluka

  integer sle, skf
  logical average_velocity, virtual_part, visc, &
       visc_artificial, heat_artificial, heat_external, gravity


  integer niac
  integer pair_i(max_interaction), pair_j(max_interaction)

  double precision x0(dim, maxn), v0(dim, maxn), &
       dxdt(dim, maxn), dvdt(dim, maxn), dx(dim, maxn)
  double precision rho0(maxn), u0(maxn), drhodt(maxn), dudt(maxn)
  double precision c(maxn), eta(maxn)
  
  character (len = *) :: filedir

  double precision  time, ut, tt, up, mvt(dim)

  integer i, j, k, its, d, ng, file_step

  
  file_step = 0
  time = 0.0

!!$-----------------------------------------------------------------
!!$ Interaction parameters, calculating neighboring particles
  call link_list(np +nv, x, h, niac, pair_i, pair_j)
  
  call ensightout_case(dt, steps, save_step, filedir)
  call output(np, nv, x, v, mat, h, mass, rho, p, u, filedir)
  call ensightout_results(np +nv, x, v, rho, p, u, file_step, filedir)

  open(21, file = 'X.dat')
  open(22, file = 'V.dat')
  open(23, file = 'P.dat')

!!$------------------------------------------------------------------
  
  do its = 1, steps
     
!!$ Interaction parameters, calculating neighboring particles
     if (mod(its, refresh_step) .eq. 0) then
        call link_list(np +nv, x, h, niac, pair_i, pair_j)
!!$        call direct_find(np +nv, x, h, niac, pair_i, pair_j)
     endif
     
!!$     do i = 1, niac
!!$        if (pair_j(i) .eq. np+1) write(41,*) pair_i(i)
!!$     enddo
     
!!$ Half step predictor
     do i = 1, np
        do d = 1, dim
           x0(d, i) = x(d, i)
           v0(d, i) = v(d, i)
        enddo
        rho0(i) = rho(i)
        u0(i) = u(i)
     enddo
     
     do i = 1, np
        do d = 1, dim
           v(d, i) = v(d, i) + dt/2.0d0 * dvdt(d, i)
           x(d, i) = x(d, i) + dt/2.0d0 * dxdt(d, i)
        enddo
        rho(i) = rho(i) + dt/2.0d0 * drhodt(i)
        u(i) = u(i) + dt/2.0d0 * dudt(i)
     enddo
     
     do i = 1, np
        do d = 1, dim
           dxdt(d, i) = 0.0d0
           dvdt(d, i) = 0.0d0
        enddo
        drhodt(i) = 0.0d0
        dudt(i) = 0.0d0
     enddo
     
!!$ Particle quantities change rate:
     call con_density(np, x, v, h, mass, drhodt, &
          niac, pair_i, pair_j, skf)
     
     call eos(np, mat, rho, u, p, c, eta, prop)
     
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
     
!!$     call tensile_correction(np, x, h, mass, rho, p, &
!!$          dvdt, niac, pair_i, pair_j, skf)
     
     if (virtual_part) then
        do i = 1, np
           do d = 1, dim
              dx(d, i) = 0.0d0
           enddo
        enddo
        
        call contact_displacement(np, nv, x, v, h, dx, &
             niac, pair_i, pair_j)
        
        call smooth_vector(np, x, h, mass, rho, dx, &
             niac, pair_i, pair_j, skf)
        
        do i = 1, np
           do d = 1, dim
              dvdt(d, i) = dvdt(d, i) + 0.5d0*dx(d, i)/(dt*dt)
           enddo
        enddo
     end if
     
     if ((gravity) .and. (dim .ge. 2)) then
        do i = 1, np
           dvdt(2, i) = dvdt(2, i) - 9.806
        enddo
     endif
     
!!$ Calculate the external heating on the particles
     if ((heat_external) .and. (its*dt .le. dtfluka)) then
        call fluka_heat(np, x, rho, dudt, &
             ofluka, dfluka, nfluka, tfluka, datafluka)
     endif
     
     do i = 1, np
        do d = 1, dim
           v(d, i) = v0(d, i) + dt*dvdt(d, i)
           dxdt(d, i) = v(d, i)
        enddo
        rho(i) = rho0(i) + dt*drhodt(i)
        u(i) = u0(i) + dt*dudt(i)
     enddo
     
!!$ Calculating average velocity of each partile for avoiding penetration
     if (average_velocity) call av_vel(np, x, v, h, mass, rho, dxdt, &
          niac, pair_i, pair_j, skf)
     
     do i = 1, np
        do d = 1, dim
           x(d, i) = x0(d, i) + dt*dxdt(d, i)
        enddo
     enddo
     
!!$ Calculating the neighboring particles and undating H
     if (sle .ne. 0) call h_upgrade(np, x, v, h, mass, rho, &
          niac, pair_i, pair_j, dt, sle, skf)
     
     time = time + dt
     
     ut = 0.0
     tt = 0.0
     up = 0.0
     do d = 1, dim
        mvt(d) = 0.0
     enddo
     
     do i = 1, np
        ut = ut + mass(i)*u(i)
        do d = 1, dim
           tt = tt + 0.5*mass(i)*v(d,i)*v(d,i)
           mvt(d) = mvt(d) + mass(i) * v(d,i)
        enddo
     enddo
     
!!$     DEBUG
     write(31,*) time, x(1,1), v(1,1), p(1), p(2), p(3), p(4), p(5)
     write(32,*) time, x(1,np +2), v(1,np +2)
     write(24,*) time, ut, tt, ut + tt
!!$     
     if (mod(its, save_step) == 0) then
        write(*,*)'______________________________________________'
        write(*,*)'  current number of time step =', &
             its,'     current time=', real(time +dt)
        write(*,*)'______________________________________________'

        call output(np, nv, x, v, mat, h, mass, rho, p, u, filedir)

        file_step = file_step +1
        call ensightout_results(np +nv, x, v, rho, p, u, &
             file_step, filedir)
     endif

  enddo

  close(21)
  close(22)
  close(23)
  close(24)
  
  call output(np, nv, x, v, mat, h, mass, rho, p, u, filedir)

end subroutine time_integration
