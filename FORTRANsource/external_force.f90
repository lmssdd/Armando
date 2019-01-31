subroutine virtual_force(np, nv, x, v, h, mass, rho, dvdt, &
     niac, pair_i, pair_j, dt, skf)

!!$  ------------------------------------------------------------------
!!$  Subroutine to calculate the external forces on the right hand side 
!!$  of the Navier-Stokes equations
!!$
!!$  np-- number of real particles                                [in]
!!$  nv-- number of virtual particles                             [in]
!!$  x-- article position                                         [in]
!!$  v-- particle velocity                                        [in]
!!$  h-- smoothing length                                         [in]
!!$  dvdt-- avcceleration                                        [out]
!!$  niac-- number of interaction pairs                           [in]
!!$  pair_i-- list of first partner of interaction pair           [in]
!!$  pair_j-- list of second partner of interaction pair          [in]
!!$  dt- timestep                                                 [in]

  implicit none
  include 'options.inc'

  integer np, nv, &
       niac, pair_i(max_interaction), pair_j(max_interaction), &
       skf
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn), dt
  double precision dvdt(dim, maxn)
  
  integer i, j, k, d
  double precision pni, pnj, pn2, mh
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac
  double precision pv(dim,maxn)
  
  do i = 1, np +nv
     do d = 1 ,dim
        dvdt(d,i) = 0.0d0
        
        pv(d,i) = 0.0d0
     enddo
  enddo
  
  
!!$  Calculate penetration on the boundary
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((j .gt. np) .and. (j .le. np +nv) .and. (i .le. np)) then
        
        pni = 0.0d0
        pn2 = 0.0d0
        do d = 1, dim
           pni = pni + (x(d,i) +dt*v(d,i) -x(d,j)) * v(d,j)
           pn2 = pn2 + pv(d,i)*pv(d,i)
        enddo
        pni = 0.5*h(i) -pni
        
        if ((pni .gt. 0.0) .and. (pni*pni .gt. pn2)) then
           do d = 1, dim
              pv(d,i) = pni * v(d,j)
           enddo
        endif
     endif
  enddo
  
  
!!$  Smooth out contact forces
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        pni = 0.0d0
        pnj = 0.0d0
        do d = 1, dim
           pni = pni + pv(d,i)*pv(d,i)
           pnj = pnj + pv(d,j)*pv(d,j)
        enddo
        
        if ((pni .gt. 1.0d3* epsilon(pni)) &
             .or. (pnj .gt. 1.0d3* epsilon(pnj))) then
           
           dr2iac = 0.0
           do d = 1, dim
              dxiac(d) = x(d,i) - x(d,j)
              dr2iac = dr2iac + dxiac(d)*dxiac(d)
           enddo
           mh = 0.5d0 * (h(i) + h(j))
           driac = sqrt(dr2iac)
           call kernel(driac, dxiac, mh, w, dwdx, skf)
           
           if (pni .gt. 1.0d3* epsilon(pni)) then
              do d = 1, dim
                 dvdt(d,j) = dvdt(d,j) &
                      + w * mass(j)/rho(j) * pv(d,i)/(2.0*dt*dt)
              enddo
           endif
           
           if (pnj .gt. 1.0d3* epsilon(pnj)) then
              do d = 1, dim
                 dvdt(d,i) = dvdt(d,i) &
                      + w * mass(i)/rho(i) * pv(d,j)/(2.0*dt*dt)
              enddo
           endif
           
        endif
     endif
  enddo
  
  do i = 1, np
     pni = 0.0d0
     do d = 1, dim
        pni = pni + pv(d,i)*pv(d,i)
     enddo
     
     if (pni .gt. 1.0d3* epsilon(pni)) then
        do d = 1, dim
           dxiac(d) = 0.0d0
        enddo
        driac = 0.0d0
        call kernel(driac, dxiac, h(i), w, dwdx, skf)
        
        do d = 1, dim
           dvdt(d,i) = dvdt(d,i) &
                + w * mass(i)/rho(i) * pv(d,i)/(2.0*dt*dt)
        enddo
     endif
  enddo
  
end subroutine virtual_force



subroutine contact_displacement_old(np, nv, x, v, h, dx, &
     niac, pair_i, pair_j)

!!$  ------------------------------------------------------------------
!!$  Subroutine to calculate the contact displacemetn of the particles
!!$
!!$  np-- number of real particles                                [in]
!!$  nv-- number of virtual particles                             [in]
!!$  x-- article position                                         [in]
!!$  v-- particle velocity                                        [in]
!!$  h-- smoothing length                                         [in]
!!$  dx-- displacement                                        [in/out]
!!$  niac-- number of interaction pairs                           [in]
!!$  pair_i-- list of first partner of interaction pair           [in]
!!$  pair_j-- list of second partner of interaction pair          [in]
  
  
  implicit none
  include 'options.inc'

  integer np, nv, &
       niac, pair_i(max_interaction), pair_j(max_interaction)
  
  double precision x(dim, maxn), v(dim, maxn), h(maxn), dx(dim, maxn)
  
  integer i, j, k, d
  double precision pni
  
  
!!$  Calculate penetration on the boundary
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((j .gt. np) .and. (j .le. np +nv) .and. (i .le. np)) then
        
        pni = 0.0d0
        do d = 1, dim
           pni = pni + (x(d,i) -x(d,j)) * v(d,j)
        enddo
        pni = h(i) / hdr -pni
        
        if (pni .gt. 0.0) then
           do d = 1, dim
              dx(d,i) = dx(d,i) + pni * v(d,j)
           enddo
        endif
     endif
  enddo
  
end subroutine contact_displacement_old



subroutine contact_displacement(np, nv, x, v, h, dx, &
     niac, pair_i, pair_j)

!!$  ------------------------------------------------------------------
!!$  Subroutine to calculate the contact displacemetn of the particles
!!$
!!$  np-- number of real particles                                [in]
!!$  nv-- number of virtual particles                             [in]
!!$  x-- article position                                         [in]
!!$  v-- particle velocity                                        [in]
!!$  h-- smoothing length                                         [in]
!!$  dx-- displacement                                        [in/out]
!!$  niac-- number of interaction pairs                           [in]
!!$  pair_i-- list of first partner of interaction pair           [in]
!!$  pair_j-- list of second partner of interaction pair          [in]
  
  
  implicit none
  include 'options.inc'

  integer np, nv, &
       niac, pair_i(max_interaction), pair_j(max_interaction)
  
  double precision x(dim, maxn), v(dim, maxn), h(maxn), dx(dim, maxn)
  
  integer i, j, k, d
  double precision pni, pnm
  
  
!!$  Calculate penetration on the boundary
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((j .gt. np) .and. (j .le. np +nv) .and. (i .le. np)) then
        
        pni = 0.0d0
        pnm = 0.0d0
        do d = 1, dim
           pni = pni + (x(d,i) -x(d,j)) * v(d,j)
           pnm = pnm + dx(d,i) * v(d,j)
        enddo
        pni = h(i) / hdr -pni
        
        if ((pni .gt. 0.0) .and. (pni .gt. pnm)) then
           do d = 1, dim
              dx(d,i) = dx(d,i) + (pni -pnm) * v(d,j)
           enddo
        endif
     endif
  enddo
  
end subroutine contact_displacement


subroutine smooth_vector(np, x, h, mass, rho, vec, &
     niac, pair_i, pair_j, skf)

!!$  ------------------------------------------------------------------
!!$  Smooth out a vector on the particles
!!$
!!$  np-- number of real particles                                [in]
!!$  x-- article position                                         [in]
!!$  h-- smoothing length                                         [in]
!!$  vec-- vector                                                 [in]
!!$  niac-- number of interaction pairs                           [in]
!!$  pair_i-- list of first partner of interaction pair           [in]
!!$  pair_j-- list of second partner of interaction pair          [in]

  implicit none
  include 'options.inc'

  integer np, &
       niac, pair_i(max_interaction), pair_j(max_interaction), &
       skf
  double precision x(dim, maxn), h(maxn), mass(maxn), rho(maxn) 
  double precision vec(dim, maxn)
  
  integer i, j, k, d
  double precision vout(dim, maxn)
  double precision vmi, vmj, mh
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac
  
  
  do i = 1, np
     do d = 1, dim
        vout(d, i) = 0.0d0
     enddo
  enddo
  
  do i = 1, np
     vmi = 0.0d0
     do d = 1, dim
        vmi = vmi + vec(d,i)*vec(d,i)
     enddo
     
     do d = 1, dim
        dxiac(d) = 0.0d0
     enddo
     driac = 0.0d0
     call kernel(driac, dxiac, h(i), w, dwdx, skf)
     
     do d = 1, dim
        vout(d,i) = w * mass(i)/rho(i) * vec(d,i)
     enddo
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        vmi = 0.0d0
        vmj = 0.0d0
        do d = 1, dim
           vmi = vmi + vec(d,i)*vec(d,i)
           vmj = vmj + vec(d,j)*vec(d,j)
        enddo
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5d0 * (h(i) + h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        do d = 1, dim
           vout(d,j) = vout(d,j) &
                + w * mass(j)/rho(j) * vec(d,i)
        enddo
        
        do d = 1, dim
           vout(d,i) = vout(d,i) &
                + w * mass(i)/rho(i) * vec(d,j)
        enddo
        
     endif
  enddo
  
  do i = 1, np
     do d = 1, dim
        vec(d, i) = vout(d, i)
     enddo
  enddo
  
end subroutine smooth_vector



subroutine ghost_force(np, nv, x, v, mat, h, &
     mass, rho, p, eta, dvdt, dudt, tdsdt, &
     niac, pair_i, pair_j, dt, skf, pa_sph, visc)

!!$  ------------------------------------------------------------------
!!$  Subroutine to calculate the external forces on the right hand side 
!!$  of the Navier-Stokes equations
!!$
!!$  np-- number of real particles                                [in]
!!$  nv-- number of virtual particles                             [in]
!!$  x-- article position                                         [in]
!!$  v-- particle velocity                                        [in]
!!$  mat-- particle material                                      [in]
!!$  h-- smoothing length                                         [in]
!!$  mass-- particle mass                                         [in]
!!$  rho-- density                                                [in]
!!$  p-- pressure                                                 [in]
!!$  u-- particle internal energy                                 [in]
!!$  c-- particle sound velocity                                  [in]
!!$  eta-- dynamic viscosity                                     [in]
!!$  dvdt-- avcceleration                                        [out]
!!$  dudt-- internal energy change rate                          [out]
!!$  tdsdt-- viscous entropy production                          [out]
!!$  niac-- number of interaction pairs                           [in]
!!$  pair_i-- list of first partner of interaction pair           [in]
!!$  pair_j-- list of second partner of interaction pair          [in]
!!$  dt- timestep                                                 [in]
!!$  skf-- smoothing kernel function                              [in]
!!$  pa_sph                                                       [in]
!!$  visc                                                         [in]

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
  double precision  hv

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

  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .gt. np +nv)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)

!!$ For SPH algorithm 1

        if(pa_sph .eq. 0) then

!!$ DIMENSION = 1

           if (dim == 1) then
!!$ X component
              hv = -(p(i) + p(j))*dwdx(1)
              hv = hv / (rho(i)*rho(j))

              dvdt(1,i) = dvdt(1,i) + mass(j)*hv
              dvdt(1,j) = dvdt(1,j) - mass(i)*hv

!!$ DIMENSION = 2

           else if (dim == 2) then
!!$ X component
              hv = -(p(i) + p(j))*dwdx(1)
              hv = hv / (rho(i)*rho(j))

              dvdt(1,i) = dvdt(1,i) + mass(j)*hv
              dvdt(1,j) = dvdt(1,j) - mass(i)*hv

!!$ Y component
              hv = -(p(i) + p(j))*dwdx(2)
              hv = hv / (rho(i)*rho(j))

              dvdt(2,i) = dvdt(2,i) + mass(j)*hv
              dvdt(2,j) = dvdt(2,j) - mass(i)*hv


!!$ DIMENSION = 3

           else if (dim == 3) then
!!$ X component
              hv = -(p(i) + p(j))*dwdx(1)
              hv = hv / (rho(i)*rho(j))

              dvdt(1,i) = dvdt(1,i) + mass(j)*hv
              dvdt(1,j) = dvdt(1,j) - mass(i)*hv

!!$ Y component
              hv = -(p(i) + p(j))*dwdx(2)
              hv = hv / (rho(i)*rho(j))

              dvdt(2,i) = dvdt(2,i) + mass(j)*hv
              dvdt(2,j) = dvdt(2,j) - mass(i)*hv

!!$ Z component
              hv = -(p(i) + p(j))*dwdx(3)
              hv = hv / (rho(i)*rho(j))

              dvdt(3,i) = dvdt(3,i) + mass(j)*hv
              dvdt(3,j) = dvdt(3,j) - mass(i)*hv

           endif

!!$ For SPH algorithm 2

        else if (pa_sph .eq. 1) then

!!$ DIMENSION = 1

           if (dim == 1) then
!!$ X component
              hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(1)

              dvdt(1,i) = dvdt(1,i) + mass(j)*hv
              dvdt(1,j) = dvdt(1,j) - mass(i)*hv


!!$ DIMENSION = 2

           else if (dim == 2) then
!!$ X component
              hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(1)

              dvdt(1,i) = dvdt(1,i) + mass(j)*hv
              dvdt(1,j) = dvdt(1,j) - mass(i)*hv

!!$ Y component
              hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(2)

              dvdt(2,i) = dvdt(2,i) + mass(j)*hv
              dvdt(2,j) = dvdt(2,j) - mass(i)*hv


!!$ DIMENSION = 3

           else if (dim == 3) then
!!$ X component
              hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(1)

              dvdt(1,i) = dvdt(1,i) + mass(j)*hv
              dvdt(1,j) = dvdt(1,j) - mass(i)*hv

!!$ Y component
              hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(2)

              dvdt(2,i) = dvdt(2,i) + mass(j)*hv
              dvdt(2,j) = dvdt(2,j) - mass(i)*hv

!!$ Z component
              hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(3)

              dvdt(3,i) = dvdt(3,i) + mass(j)*hv
              dvdt(3,j) = dvdt(3,j) - mass(i)*hv

           endif

        endif

     endif

  enddo

end subroutine ghost_force
