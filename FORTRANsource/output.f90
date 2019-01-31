subroutine output(np, nv, x, v, mat, h, mass, rho, p, u, &
     filedir) 

!!$------------------------------------------------------------------           
!!$ Subroutine for saving particle information to external disk file
!!$
!!$ np-- total particle number                                    [in]
!!$ x-- coordinates of particles                                  [in]
!!$ v-- velocities of particles                                   [in]
!!$ mat-- material type of particles                              [in]
!!$ h-- smoothing lengths of particles                            [in]
!!$ mass-- mass of particles                                      [in]
!!$ rho-- dnesities of particles                                  [in]
!!$ p-- pressure  of particles                                    [in]
!!$ u-- internal energy of particles                              [in]
!!$ c-- sound velocity of particles                               [in]
!!$ filedir -- Directory where the files are located             [in]

  implicit none     
  include 'options.inc'

  integer mat(maxn), np, nv
  double precision x(dim, maxn), v(dim, maxn), h(maxn), &
       mass(maxn), rho(maxn),p(maxn), u(maxn)
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
     write(1,'(1x, I6, 6(2x, e14.8))') i, (x(d, i), d = 1, dim), &
          (v(d, i), d = 1, dim)
     write(2,'(1x, I6, 4(2x, e14.8))') i, mass(i), rho(i), p(i), u(i)
     write(3,'(1x, I6, 2x, I4, 2x, e14.8)') i, mat(i), h(i)
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
     write(1,'(1x, I6, 6(2x, e14.8))') i, (x(d, i), d = 1, dim), &
          (v(d, i), d = 1, dim)
     write(2,'(1x, I6, 4(2x, e14.8))') i, mass(i), rho(i), p(i), u(i)
     write(3,'(1x, I6, 2x, I4, 2x, e14.8)') i, mat(i), h(i)
  enddo

  close(1)
  close(2)
  close(3)

end subroutine output
