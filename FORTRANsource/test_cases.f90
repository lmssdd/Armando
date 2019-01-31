subroutine partly_heated_rod(np, x, v, mat, h, mass, rho, p, u)

!!$-----------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d
  double precision dx, d0, e0, d1


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  np = 800
  dx = 1.0/np

  do i = 1, np
     mass(i) = 13540.0/np
     h(i) = 2.0*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0

     do d = 1, dim
        x(d,i) = 0. 
        v(d,i) = 0.
     enddo
  enddo

  do i = 1, np
     x(1,i) = -0.5 + dx * (i -1)
  enddo

  d0 = 12.0d-2
  d1 = 13.0d-2
  e0 = 13540.0

  do i = 1, np
     if (dabs(x(1,i)) < d0) u(i) = e0 / rho(i)
     if ((dabs(x(1,i)) < d1) .and. (dabs(x(1,i)) >= d0)) then
        u(i) = (e0 / rho(i)) * (d1 - dabs(x(1,i))) / (d1 -d0)
     endif
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine partly_heated_rod

subroutine bullet(np, nv, x, v, mat, h, &
     mass, rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d
  double precision dx, d0, e0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  np = 100
  nv = 1
  dx = 1.0/np

  do i = 1, np
     mass(i) = 1000.0/np
     h(i) = 2.0*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 1000.0
     p(i) = 0.0
  enddo

  do i = 1, np
     x(1,i) = -0.5 + dx * (i -1) +0.5*dx
     v(1,i) = -1.0
  enddo
  
  x(1, np +1) = x(1,1) -1.5*dx
  v(1, np +1) = 1.0
  mass(np +1) = 0.0
  h(np +1) = 0.0
  mat(np +1) = 0

!!$  x(1, np +2) = +0.5 !+0.5*dx
!!$  v(1, np +2) = -1.0 !+0.5*dx

!!$  d0 = 5.0d-2
!!$  e0 = 13540.0
!!$
!!$  do i = 1, np
!!$     if (dabs(x(1,i)) < d0) u(i) = e0 / rho(i)
!!$  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine bullet


subroutine drop(np, nv, x, v, mat, h, &
     mass, rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 2d drop test
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, nx, ny
  double precision dx, dy, r, xp, yp 


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 20
  ny = 20
  dx = 0.01/nx
  dy = 0.01/ny

  np = 0
  do j = 1, ny
     do i = 1, nx
        xp = -0.5d-2 +dx * (i -1)
        yp = -0.5d-2 +dy * (j -1)
        r = dsqrt(xp*xp + yp*yp)
        if (r < 0.5d-2) then
           np = np +1
           x(1,np) = xp
           x(2,np) = yp
        endif
     enddo
  enddo
  
  do i = 1, np
     mass(i) = 1000.0*dx*dy
     h(i) = 2.0*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 1000.0
     p(i) = 0.0
     v(1,i) = 1.0d0
  enddo
  
  nv = 0
  
!!$  do i = 1, 2*ny
!!$     xp = -1.0d-2
!!$     yp = -1.0d-2 + dy * (i -1)
!!$     nv = nv +1
!!$     x(1,np +nv) = xp
!!$     x(2,np +nv) = yp
!!$     v(1,np +nv) = 1.0d0
!!$     v(2,np +nv) = 0.0d0
!!$     mass(np +nv) = 0.0d0
!!$     h(np +nv) = 2.0d0*dy
!!$     mat(np +nv) = 0
!!$  enddo

  do i = 1, 2*ny
     xp = 1.0d-2
     yp = -1.0d-2 + dy * (i -1)
     nv = nv +1
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = -1.0d0
     v(2,np +nv) = 0.0d0
     mass(np +nv) = 0.0d0
     h(np +nv) = 2.0d0*dy
     mat(np +nv) = 0
  enddo

  do i = 1, 2*nx
     xp = -1.0d-2 + dx * (i -1)
     yp = -1.0d-2
     nv = nv +1
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = 0.0d0
     v(2,np +nv) = 1.0d0
     mass(np +nv) = 0.0d0
     h(np +nv) = 2.0d0*dy
     mat(np +nv) = 0
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine drop


subroutine dam(np, nv, x, v, mat, h, &
     mass, rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 2d drop test
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, nx, ny
  double precision dx, dy, r, xp, yp 


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 40
  ny = 80
  dx = 1.0/nx
  dy = 2.0/ny

  np = 0
  do j = 1, ny
     do i = 1, nx
        xp = dx * (i -1)
        yp = dy * (j -1)
        np = np +1
        x(1,np) = xp
        x(2,np) = yp
     enddo
  enddo
  
  do i = 1, np
     mass(i) = 1000.0*dx*dy
     h(i) = hdr*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 1000.0
     p(i) = 0.0e5
     v(1,i) = 0.0d0
     v(2,i) = 0.0d0
  enddo
  
  nv = 0
  
  do i = 1, 2*ny +1
     xp = -dx
     yp = dy * (i -1)
     nv = nv +1
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = 1.0d0
     v(2,np +nv) = 0.0d0
     mass(np +nv) = 0.0d0
     rho (np +nv) = 1000.0d0
     h(np +nv) = hdr*dy
     mat(np +nv) = 0
  enddo

  do i = 1, 2*ny +1
     xp = 4.0d0 +dx
     yp = dy * (i -1)
     nv = nv +1
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = -1.0d0
     v(2,np +nv) = 0.0d0
     mass(np +nv) = 0.0d0
     rho (np +nv) = 1000.0d0
     h(np +nv) = hdr*dy
     mat(np +nv) = 0
  enddo
  
  do i = 1, 4*nx +1
     xp = dx * (i -1)
     yp = -dy
     nv = nv +1
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = 0.0d0
     v(2,np +nv) = 1.0d0
     mass(np +nv) = 0.0d0
     rho (np +nv) = 1000.0d0
     h(np +nv) = hdr*dy
     mat(np +nv) = 0
  enddo

  do i = 1, 4*nx +1
     xp = dx * (i -1)
     yp = 4.0 +dy
     nv = nv +1
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = 0.0d0
     v(2,np +nv) = -1.0d0
     mass(np +nv) = 0.0d0
     rho (np +nv) = 1000.0d0
     h(np +nv) = hdr*dy
     mat(np +nv) = 0
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine dam


subroutine partly_heated_rod_b(np, nv, x, v, mat, h, &
     mass, rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d
  double precision dx, d0, e0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  np = 800
  nv = 2
  dx = 1.0/np

  do i = 1, np + nv
     mass(i) = 13540.0/np
     h(i) = 2.0*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0

     do d = 1, dim
        x(d,i) = 0. 
        v(d,i) = 0.
     enddo
  enddo

  do i = 1, np
     x(1,i) = -0.5 + dx * (i -1) +0.5*dx
  enddo

  x(1, np +1) = -0.5 !-0.5*dx
  x(1, np +2) = +0.5 !+0.5*dx

  d0 = 5.0d-2
  e0 = 13540.0

  do i = 1, np
     if (dabs(x(1,i)) < d0) u(i) = e0 / rho(i)
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine partly_heated_rod_b


subroutine partly_heated_disc(np, nv, x, v, mat, h, mass, &
     rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, nx, ny
  double precision dx, dy, xp, yp, r, r0, e, e0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 200
  ny = 200
  dx = 0.1/nx
  dy = 0.1/ny
  nv = 3*nx

  np = 0
  do j = 1, ny
     do i = 1, nx
        xp = -0.5d-1 +dx * (i -1)
        yp = -0.5d-1 +dy * (j -1)
        r = dsqrt(xp*xp + yp*yp)
        if (r < 5.0d-2) then
           np = np +1
           x(1,np) = xp
           x(2,np) = yp
        endif
     enddo
  enddo

  do i = 1, nv
     x(1, np +i) = 0.05 * cos(2.0d0*pi/nv * i)
     x(2, np +i) = 0.05 * sin(2.0d0*pi/nv * i)

     v(1, np +i) = cos(2.0d0*pi/nv * i)
     v(2, np +i) = sin(2.0d0*pi/nv * i)
  enddo

  do i = 1, np
     mass(i) = 1.0d-2*13540.0/(nx*ny)
     h(i) = 2.0*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  do i = np +1, np +nv
     mass(i) = 0.0
     h(i) = 2.0*dx
     mat(i) = -21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  r0 = 2.0d-3
  e0 = 13540.0

  do i = 1, np
     !	  r = dsqrt(x(1,i)**2 +(x(2,i)-5.0*r0)**2)
     r = dsqrt(x(1,i)**2 +x(2,i)**2)
     e = e0 * exp(-(r/r0)**2)
     e = e / rho(i)

     u(i) = e
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine partly_heated_disc


subroutine partly_heated_bar(np, nv, x, v, mat, h, mass, &
     rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1 d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ nv-- virtual particle number                                [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, nx, ny
  double precision dx, dy, xp, yp, r, r0, e, e0, vx, vy, r1, r2


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 400
  ny = 400
  dx = 0.12/nx
  dy = 0.12/ny

  vx = cos(60./180. * pi)
  vy = sin(60./180. * pi)

  np = 0
  do j = 1, ny
     do i = 1, nx
        xp = -0.6d-1 +dx * (i -0.5)
        yp = -0.6d-1 +dy * (j -0.5) 
        r1 = xp*vx + yp*vy
        r2 = -xp*vy + yp*vx
        if ((dabs(r1) < 0.5d-2) .and. (dabs(r2) < 5.0d-2)) then
           np = np +1
           x(1,np) = xp
           x(2,np) = yp
        endif
     enddo
  enddo

  do i = 1, np
     mass(i) = 0.12*0.12*13540.0/(nx*ny)
     h(i) = 2.01*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  nv = 0

  do i = 1, nx
     nv = nv +1
     xp = 0.5d-2 * vx - dx * (i -0.5 -0.5*nx) * vy
     yp = 0.5d-2 * vy + dx * (i -0.5 -0.5*nx) * vx
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = vx
     v(2,np +nv) = vy
  enddo

  do i = 1, nx
     nv = nv +1
     xp = -0.5d-2 * vx - dx * (i -0.5 -0.5*nx) * vy
     yp = -0.5d-2 * vy + dx * (i -0.5 -0.5*nx) * vx
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = vx
     v(2,np +nv) = vy
  enddo

  do i = 1, nx
     nv = nv +1
     xp = -5.0d-2 * vy + dx * (i -0.5 -0.5*nx) * vx
     yp = 5.0d-2 * vx + dx * (i -0.5 -0.5*nx) * vy
     x(1,np +nv) = xp
     x(2,np +nv) = yp
     v(1,np +nv) = -vy
     v(2,np +nv) = vx
  enddo

  do i = np +1, np +nv
     mass(i) = 0.0
     h(i) = 2.01*dx
     mat(i) = -21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  r0 = 2.0d-3
  e0 = 13540.0
  do i = 1, np
     r = dabs(-x(1,i)*vy + x(2,i)*vx)
     e = e0 * exp(-(r/r0)**2)
     e = e / rho(i)

     u(i) = e
  enddo

!!$	j = 12345
!!$      do i = 1, np
!!$	  x(2,i) = x(2,i) +0.1*dy*(ran(j) -0.5)
!!$      enddo


  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine partly_heated_bar


subroutine partly_heated_cylinder(np, nv, x, v, mat, h, mass, &
     rho, p, u)

!!$ ------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ nv-- virtual particle number                                [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, k, nxy, nz
  double precision dxy, dz, xp, yp, zp, rp, e, e0, r, r0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nxy = 20
  nz = 160
  dxy = 0.02/nxy
  dz = 0.16/nz

  np = 0
  do k = 1, nz
     do j = 1, nxy
        do i = 1, nxy
           xp = -0.01 +dxy * (i -0.5)
           yp = -0.01 +dxy * (j -0.5) 
           zp = -0.08 +dz * (k -0.5)

           rp = dsqrt(xp*xp + yp*yp)

           if (rp .le. 0.01) then
              np = np +1
              x(1,np) = xp
              x(2,np) = yp
              x(3,np) = zp
           endif
        enddo
     enddo
  enddo

  do i = 1, np
     mass(i) = 0.02*0.02*0.16* 13540.0/(nxy*nxy*nz)
     h(i) = 2.01*dxy
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  r0 = 0.01
  e0 = 139.0*200.0
  do i = 1, np
     r = dsqrt(x(3,i)**2 + x(2,i)**2)
     e = e0 * exp(-(r/r0)**2)

     u(i) = e
  enddo


  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine partly_heated_cylinder


subroutine fluka_heated_disc(np, nv, x, v, mat, h, mass, &
     rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ nv-- virtual particle number                                [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, nx
  double precision dx, xp, yp, r, r0, e, e0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 200
  nv = 3 * nx
  dx = 0.06/nx
  np = 0
  do j = 1, nx
     do i = 1, nx
        xp = -0.03 +dx * (i -1)
        yp = -0.03 +dx * (j -1)
        r = dsqrt(xp*xp + yp*yp)
        if ((r < (3.0d-2 -0.5*dx)) .and. (xp < 0.02)) then
           np = np +1
           x(1,np) = xp
           x(2,np) = yp
        endif
     enddo
  enddo

  do i = 1, nv
     x(1, np +i) = 0.03 * cos(2.0d0*pi/nv * i)
     x(2, np +i) = 0.03 * sin(2.0d0*pi/nv * i)

     v(1, np +i) = cos(2.0d0*pi/nv * i)
     v(2, np +i) = sin(2.0d0*pi/nv * i)
  enddo

  do i = 1, np
     mass(i) = 0.06*0.06*1549.0/(nx*nx)
     h(i) = 3.0*dx
     mat(i) = 31
     u(i) = 0.0
     rho(i) = 1549.0
     p(i) = 0.0
  enddo

  do i = np +1, np +nv
     mass(i) = 0.0
     h(i) = 3.0*dx
     mat(i) = -31
     u(i) = 0.0
     rho(i) = 1549.0
     p(i) = 0.0
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine fluka_heated_disc


subroutine mercury_pot2(np, nv, x, v, mat, h, mass, &
     rho, p, u)

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, nx
  double precision dx, xp, yp, r, r0, e, e0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 50
  dx = 0.012/nx
  
  np = 0
  do j = 1, nx +1
     do i = 1, nx +1
        xp = -0.006 +dx * (i -1)
        yp = -0.006 +dx * (j -1)
        r = dsqrt(xp*xp + yp*yp)
        if (.not. ((r > (0.006 - 0.5*dx)) .and. (yp < 0.00))) then
           np = np +1
           x(1,np) = xp
           x(2,np) = yp
        endif
     enddo
  enddo
  
  nv = 0
  do i = 1, nx
     nv = nv +1
     x(1, np +nv) = -0.006 * cos(1.0d0*pi/nx * i)
     x(2, np +nv) = -0.006 * sin(1.0d0*pi/nx * i)

     v(1, np +nv) = cos(1.0d0*pi/nx * i)
     v(2, np +nv) = sin(1.0d0*pi/nx * i)
  enddo

  do i = 1, nx/2
     nv = nv +1
     x(1, np +nv) = -0.006
     x(2, np +nv) = dx * (i -1)
     v(1, np +nv) = 1.0d0
     v(2, np +nv) = 0.0d0
  enddo

  do i = 1, nx/2
     nv = nv +1
     x(1, np +nv) = 0.006
     x(2, np +nv) = dx * (i -1)
     v(1, np +nv) = -1.0d0
     v(2, np +nv) = 0.0d0
  enddo

  do i = 1, np
     mass(i) = 0.012*0.012*13540.0/(nx*nx)
     h(i) = hdr*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 13540.0*9.81*(0.006 -x(2,i))
  enddo
  
  do i = np +1, np +nv
     mass(i) = 0.0
     h(i) = hdr*dx
     mat(i) = -21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  r0 = 0.00315
  e0 = 139.0*160.0
  do i = 1, np
     r = dsqrt(x(1,i)**2 + (x(2,i))**2)
     e = e0 * exp(-(r/r0)**2)

     u(i) = e
  enddo
  
  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine mercury_pot2


subroutine mercury_pot3(np, nv, x, v, mat, h, mass, &
     rho, p, u)

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, k, d, nx, ni, nj
  double precision dx, xp, yp, zp, rp, r, r0, e, e0


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  nx = 20
  dx = 0.012/nx
  
  np = 0
  do k = 1, nx +1
     do j = 1, nx +1
        do i = 1, nx +1
           xp = -0.006 +dx * (i -1)
           yp = -0.006 +dx * (j -1)
           zp = -0.006 +dx * (k -1)
           if (yp .ge. 0.0d0) then
              r = dsqrt(xp*xp + zp*zp)
              if (r .lt. 0.006) then
                 np = np +1
                 x(1,np) = xp
                 x(2,np) = yp
                 x(3,np) = zp
              endif
           else
              r = dsqrt(xp*xp + yp*yp + zp*zp)
              if (r .lt. 0.006) then
                 np = np +1
                 x(1,np) = xp
                 x(2,np) = yp
                 x(3,np) = zp
              endif
           endif
        enddo
     enddo
  enddo
  
  nv = 0
  
  nj = nx
  do j = 1, nj
     yp = -0.006*sin(pi/2.0/nj*j)
     rp = dsqrt(0.006**2 - yp*yp)
     ni = int(2*pi*rp/dx)
     
     do i = 1, ni
        nv = nv +1
        xp = -rp * cos(2.0d0*pi/ni * i)
        zp = -rp * sin(2.0d0*pi/ni * i)
        
        r = dsqrt(xp*xp + yp*yp + zp*zp)
        
        x(1, np +nv) = xp
        x(2, np +nv) = yp
        x(3, np +nv) = zp
        
        v(1, np +nv) = -xp/r
        v(2, np +nv) = -yp/r
        v(3, np +nv) = -zp/r
     enddo
  enddo
  

  nj = nx/2
  do j = 1, nj +1
     ni = 3*nx
     do i = 1, ni
        nv = nv +1
        x(1, np +nv) = -0.006 * cos(2.0d0*pi/ni * i)
        x(2, np +nv) = 0.006/nj*(j -1)
        x(3, np +nv) = -0.006 * sin(2.0d0*pi/ni * i)
        
        v(1, np +nv) = cos(2.0d0*pi/ni * i)
        v(2, np +nv) = 0.0d0
        v(3, np +nv) = sin(2.0d0*pi/ni * i)
     enddo
  enddo

  do i = 1, np
     mass(i) = 13540.0 * (0.012/nx)**3
     h(i) = hdr*dx
     mat(i) = 21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 13540.0*9.81*(0.006 -x(2,i))
  enddo
  
  do i = np +1, np +nv
     mass(i) = 0.0
     h(i) = hdr*dx
     mat(i) = -21
     u(i) = 0.0
     rho(i) = 13540.0
     p(i) = 0.0
  enddo

  r0 = 0.00315
  e0 = 139.0*160.0
  do i = 1, np
     r = dsqrt(x(1,i)**2 + (x(2,i))**2)
     e = e0 * exp(-(r/r0)**2)

     u(i) = e
  enddo
  
  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine mercury_pot3

subroutine shock_tube(np, x, v, mat, h, mass, rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 1d noh shock tube problem
!!$ np-- particle number                                        [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d
  double precision space_x


  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  np = 400
  space_x = 0.6/80.

  do i = 1, np
     mass(i) = 0.75/400.
     h(i) = 0.015
     mat(i) = 1
     do d = 1, dim
        x(d,i) = 0. 
        v(d,i) = 0.
     enddo
  enddo

  do i = 1,320
     x(1,i) = -0.6 + space_x/4. * (i -1)
  enddo

  do i = 320+1, np
     x(1,i) = 0. + space_x * (i-320)
  enddo

  do i = 1, np
     if (x(1,i).le.1.e-8) then
        u(i) = 2.5
        rho(i) = 1.
        p(i) = 1.
     endif
     if (x(1,i).gt.1.e-8)  then
        u(i) = 1.795
        rho(i) = 0.25
        p(i) = 0.1795
     endif
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'  **************************************************'

end subroutine shock_tube

subroutine shear_cavity(np, nv, x, v, mat, h, mass, rho, p, u)

!!$------------------------------------------------------------------     
!!$ This subroutine is used to generate initial data for the 
!!$ 2d shear driven cavity probem with Re = 1
!!$ np-- particle number                                        [out]
!!$ nv-- virtual particle number                                [out]
!!$ x-- coordinates of particles                                [out]
!!$ v-- velocities of particles                                 [out]
!!$ mat-- material of particles                                 [out]
!!$ h-- smoothing lengths of particles                          [out]
!!$ mass-- mass of particles                                    [out]
!!$ rho-- dnesities of particles                                [out]
!!$ p-- pressure  of particles                                  [out]
!!$ u-- internal energy of particles                            [out]

  implicit none
  include 'options.inc'

  integer np, nv, mat(maxn)
  double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn), &
       rho(maxn), p(maxn), u(maxn)
  integer i, j, d, m, n, k
  double precision xl, yl, dx, dy, v_inf

!!$ Giving mass and smoothing length as well as other data.

  write(*,*)'  **************************************************'
  write(*,*)'      Generating initial particle configuration   '

  m = 41
  n = 41
  xl = 1.0d-3
  yl = 1.0d-3
  dx = xl / (m -1)
  dy = yl / (n -1)
  v_inf = -1.e-3

  np = 0
  do i = 2, m -1
     do j = 2, n -1
        np = np +1
        x(1, np) = (i -1)*dx 
        x(2, np) = (j -1)*dy 
     enddo
  enddo

  nv = 0
!!$ Monaghan type virtual particle on the Upper side

!!$    do i = 1, m
!!$       nv = nv +1
!!$       x(1, np +nv) = (i -1)*dx
!!$       x(2, np +nv) = yl +dy
!!$       v(1, np +nv) = v_inf
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$    do i = 1, m
!!$       nv = nv +1
!!$       x(1, np +nv) = (i -1)*dx +0.5*dx
!!$       x(2, np +nv) = yl +1.5*dy
!!$       v(1, np +nv) = v_inf
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$ Monaghan type virtual particle on the Lower side

!!$    do i = 1, m
!!$       nv = nv +1
!!$       x(1, np +nv) = (i -1)*dx
!!$       x(2, np +nv) = 0. -dy
!!$       v(1, np +nv) = 0.
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$    do i = 1, m
!!$       nv = nv +1
!!$       x(1, np +nv) = (i -1)*dx -0.5*dx
!!$       x(2, np +nv) = 0. -1.5*dy
!!$       v(1, np +nv) = 0.
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$ Monaghan type virtual particle on the Left side

!!$    do j = -1, n +2
!!$       nv = nv +1
!!$       x(1, np +nv) = 0. -dx
!!$       x(2, np +nv) = (j -1)*dy
!!$       v(1, np +nv) = 0.
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$    do j = -1, n +2
!!$       nv = nv +1
!!$       x(1, np +nv) = 0. -1.5*dx
!!$       x(2, np +nv) = (j -1)*dy +0.5*dy
!!$       v(1, np +nv) = 0.
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$ Monaghan type virtual particle on the Right side

!!$    do j = -1, n +2
!!$       nv = nv +1
!!$       x(1, np +nv) = xl +dx
!!$       x(2, np +nv) = (j -1)*dy
!!$       v(1, np +nv) = 0.
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$    do j = -1, n +2
!!$       nv = nv +1
!!$       x(1, np +nv) = xl +1.5*dx
!!$       x(2, np +nv) = (j -1)*dy -0.5*dy
!!$       v(1, np +nv) = 0.
!!$       v(2, np +nv) = 0.
!!$    enddo

!!$ Monaghan type virtual particle on the Upper side

  do i = 1, 2*(m -1) +1
     nv = nv +1
     x(1, np +nv) = (i -1)*dx/2
     x(2, np +nv) = yl
     v(1, np +nv) = v_inf
     v(2, np +nv) = 0.
  enddo

!!$ Monaghan type virtual particle on the Lower side

  do i = 1, 2*(m -1) +1
     nv = nv +1
     x(1, np +nv) = (i -1)*dx/2
     x(2, np +nv) = 0.
     v(1, np +nv) = 0.
     v(2, np +nv) = 0.
  enddo

!!$ Monaghan type virtual particle on the Left side

  do i = 1, 2*(n -1)
     nv = nv +1
     x(1, np +nv) = 0.
     x(2, np +nv) = i*dy/2
     v(1, np +nv) = 0.
     v(2, np +nv) = 0.
  enddo

!!$ Monaghan type virtual particle on the Right side

  do i = 1, 2*(n -1)
     nv = nv +1
     x(1, np +nv) = xl
     x(2, np +nv) = i*dy/2
     v(1, np +nv) = 0.
     v(2, np +nv) = 0.
  enddo

  do i = 1, np
     v(1, i) = 0.
     v(2, i) = 0.
     mat(i) = 11
     h(i) = dx
     rho (i) = 1000.
     mass(i) = dx*dy*rho(i)
     p(i)= 0.
     u(i)= 0. ! 357.1
  enddo

  do i = np +1, np +nv
     mat(i) = -11
     h(i) = dx
     rho (i) = 1000.
     mass(i) = 0.75*dx*dy*rho(i)
     p(i)= 0.
     u(i)= 0. ! 357.1
  enddo

  write(*,*)'      Total number of particles   ', np
  write(*,*)'      Total number of virtual particles   ', nv
  write(*,*)'  **************************************************'

end subroutine shear_cavity
