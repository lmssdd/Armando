      subroutine output(np, nv, x, v, mat, h, mass, rho, p, u, c, 
     &                  filedir) 
      
c----------------------------------------------------------------------           
c     Subroutine for saving particle information to external disk file

c     np-- total particle number                                    [in]
c     x-- coordinates of particles                                  [in]
c     v-- velocities of particles                                   [in]
c     mat-- material type of particles                              [in]
c     h-- smoothing lengths of particles                            [in]
c     mass-- mass of particles                                      [in]
c     rho-- dnesities of particles                                  [in]
c     p-- pressure  of particles                                    [in]
c     u-- internal energy of particles                              [in]
c     c-- sound velocity of particles                               [in]
c     filedir -- Directory where the files are located             [in]

      implicit none     
      include 'options.inc'
      
      integer mat(maxn), np, nv
      double precision x(dim, maxn), v(dim, maxn), h(maxn), 
     &       mass(maxn), rho(maxn),p(maxn), u(maxn), c(maxn)
	character (len = *) :: filedir

	character (len = 120) :: filename
      integer i, d     
      
	filename = trim(filedir) // "f_xv.dat"
      open(1,file = filename)
	
	filename = trim(filedir) // "f_state.dat"
      open(2,file = filename)
	
	filename = trim(filedir) // "f_other.dat"
      open(3,file = filename) 
	
      write(1,*) np
      do i = 1, np
        write(1,1001) i, (x(d, i), d = 1, dim), (v(d, i), d = 1, dim)
        write(2,1002) i, mass(i), rho(i), p(i), u(i)
        write(3,1003) i, mat(i), h(i)
      enddo 
      
      close(1)
      close(2)
      close(3)
      
	filename = trim(filedir) // "v_xv.dat"
      open(1,file = filename)
	
	filename = trim(filedir) // "v_state.dat"
      open(2,file = filename)
	
	filename = trim(filedir) // "v_other.dat"
      open(3,file = filename) 
	
      write(1,*) nv
      do i = np +1, np +nv
        write(1,1001) i, (x(d, i), d = 1, dim), (v(d, i), d = 1, dim)
        write(2,1002) i, mass(i), rho(i), p(i), u(i)
        write(3,1003) i, mat(i), h(i)
      enddo 
      
      close(1)
      close(2)
      close(3)
      
1001  format(1x, I6, 6(2x, e14.8))
1002  format(1x, I6, 7(2x, e14.8))
1003  format(1x, I6, 2x, I4, 2x, e14.8)
      
      end
