subroutine time_elapsed(s)
!!$  use dfport
!!$  implicit none
!!$
!!$  integer, parameter :: output = 6
!!$  real(8) :: s
!!$
!!$  s = rtc()
!!$
  implicit none
  double precision :: s
  integer :: count, rate

  call system_clock(count, rate)
  s = count / rate

end subroutine time_elapsed
