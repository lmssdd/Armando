      subroutine partly_heated_rod(np, x, v, mat, h, mass, rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
	
      end
      

      subroutine partly_heated_rod_b(np, nv, x, v, mat, h, 
     &                               mass, rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
	
      end
      

      subroutine partly_heated_disc(np, nv, x, v, mat, h, mass, 
     &                              rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
c	  r = dsqrt(x(1,i)**2 +(x(2,i)-5.0*r0)**2)
	  r = dsqrt(x(1,i)**2 +x(2,i)**2)
	  e = e0 * exp(-(r/r0)**2)
	  e = e / rho(i)

	  u(i) = e
      enddo
      
	write(*,*)'      Total number of particles   ', np
      write(*,*)'  **************************************************'
	
      end


      subroutine partly_heated_bar(np, nv, x, v, mat, h, mass, 
     &                             rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     nv-- virtual particle number                                [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
	
c	j = 12345
c      do i = 1, np
c	  x(2,i) = x(2,i) +0.1*dy*(ran(j) -0.5)
c      enddo


	write(*,*)'      Total number of particles   ', np
      write(*,*)'  **************************************************'
	
      end
      

      subroutine partly_heated_cylinder(np, nv, x, v, mat, h, mass, 
     &                             rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     nv-- virtual particle number                                [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
	
      end
      

      subroutine fluka_heated_disc(np, nv, x, v, mat, h, mass, 
     &                             rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     nv-- virtual particle number                                [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
	
      end
      
      subroutine shock_tube(np, x, v, mat, h, mass, rho, p, u)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     np-- particle number                                        [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]

      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
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
	
      end
      
      subroutine shear_cavity(np, nv, x, v, mat, h, mass, rho, p, u)
	
c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2 d shear driven cavity probem with Re = 1
c     np-- particle number                                        [out]
c     nv-- virtual particle number                                [out]
c     x-- coordinates of particles                                [out]
c     v-- velocities of particles                                 [out]
c     mat-- material of particles                                 [out]
c     h-- smoothing lengths of particles                          [out]
c     mass-- mass of particles                                    [out]
c     rho-- dnesities of particles                                [out]
c     p-- pressure  of particles                                  [out]
c     u-- internal energy of particles                            [out]
	
      implicit none
      include 'options.inc'
      
      integer np, nv, mat(maxn)
      double precision x(dim, maxn), v(dim, maxn), h(maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn)
      integer i, j, d, m, n, k
      double precision xl, yl, dx, dy, v_inf
	
c     Giving mass and smoothing length as well as other data.

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
c     Monaghan type virtual particle on the Upper side

c	do i = 1, m
c	  nv = nv +1
c	  x(1, np +nv) = (i -1)*dx
c	  x(2, np +nv) = yl +dy
c	  v(1, np +nv) = v_inf
c	  v(2, np +nv) = 0.
c	enddo
	
c	do i = 1, m
c	  nv = nv +1
c	  x(1, np +nv) = (i -1)*dx +0.5*dx
c	  x(2, np +nv) = yl +1.5*dy
c	  v(1, np +nv) = v_inf
c	  v(2, np +nv) = 0.
c	enddo
	
c     Monaghan type virtual particle on the Lower side

c	do i = 1, m
c	  nv = nv +1
c	  x(1, np +nv) = (i -1)*dx
c	  x(2, np +nv) = 0. -dy
c	  v(1, np +nv) = 0.
c	  v(2, np +nv) = 0.
c	enddo

c	do i = 1, m
c	  nv = nv +1
c	  x(1, np +nv) = (i -1)*dx -0.5*dx
c	  x(2, np +nv) = 0. -1.5*dy
c	  v(1, np +nv) = 0.
c	  v(2, np +nv) = 0.
c	enddo

c     Monaghan type virtual particle on the Left side

c	do j = -1, n +2
c	  nv = nv +1
c	  x(1, np +nv) = 0. -dx
c	  x(2, np +nv) = (j -1)*dy
c	  v(1, np +nv) = 0.
c	  v(2, np +nv) = 0.
c	enddo
	
c	do j = -1, n +2
c	  nv = nv +1
c	  x(1, np +nv) = 0. -1.5*dx
c	  x(2, np +nv) = (j -1)*dy +0.5*dy
c	  v(1, np +nv) = 0.
c	  v(2, np +nv) = 0.
c	enddo
	
c     Monaghan type virtual particle on the Right side

c	do j = -1, n +2
c	  nv = nv +1
c	  x(1, np +nv) = xl +dx
c	  x(2, np +nv) = (j -1)*dy
c	  v(1, np +nv) = 0.
c	  v(2, np +nv) = 0.
c	enddo
	
c	do j = -1, n +2
c	  nv = nv +1
c	  x(1, np +nv) = xl +1.5*dx
c	  x(2, np +nv) = (j -1)*dy -0.5*dy
c	  v(1, np +nv) = 0.
c	  v(2, np +nv) = 0.
c	enddo

c     Monaghan type virtual particle on the Upper side
	
	do i = 1, 2*(m -1) +1
	  nv = nv +1
	  x(1, np +nv) = (i -1)*dx/2
	  x(2, np +nv) = yl
	  v(1, np +nv) = v_inf
	  v(2, np +nv) = 0.
	enddo
	
c     Monaghan type virtual particle on the Lower side

	do i = 1, 2*(m -1) +1
	  nv = nv +1
	  x(1, np +nv) = (i -1)*dx/2
	  x(2, np +nv) = 0.
	  v(1, np +nv) = 0.
	  v(2, np +nv) = 0.
	enddo

c     Monaghan type virtual particle on the Left side

	do i = 1, 2*(n -1)
	  nv = nv +1
	  x(1, np +nv) = 0.
	  x(2, np +nv) = i*dy/2
	  v(1, np +nv) = 0.
	  v(2, np +nv) = 0.
	enddo
	
c     Monaghan type virtual particle on the Right side

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
	
      end	 
