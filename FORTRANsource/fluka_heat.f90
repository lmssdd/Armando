subroutine read_binning_bin(trfluka, sffluka, effluka, npfluka, &
     dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, filedir)

!!$-------------------------------------------------------------------
!!$ Subroutine to read and store the FLUKA results from a binning 
!!$ binary file
!!$
!!$ trfluka-- translation of the fluka grid origin               [in]
!!$ sffluka-- size scale factor (size_sph = sf*size_fluka)       [in]
!!$ effluka-- energy scale factor (energy_sph = ef*energy_fluka) [in]
!!$ npfluka-- number of particles                                [in]
!!$ dtfluka-- deposition time                                    [in]
!!$ ofluka-- the origin of the binning                          [out]
!!$ dfluka-- cell side of the binning                           [out]
!!$ nfluka-- number of cells for each index                     [out]
!!$ tfluka-- identifier for the fluka grid type                 [out]
!!$ datafluka-- data of the binning                             [out]
!!$ filedir -- Directory where the files are located             [in]

  implicit none
  include 'options.inc'

  double precision trfluka(3), sffluka, effluka, npfluka, dtfluka, &
       ofluka(3), dfluka(3), datafluka(max_fluka)
  integer nfluka(3), tfluka
  character (len = *) :: filedir

  character (len = 120) filename
  character runtitle*80, runhead*32, namebin*10
  real weight, minbin(3), maxbin(3), dbin(3), b1bin, b2bin, tcbin
  real realfluka(max_fluka)
  integer ncase, mbin, typebin, distbin, nbin(3), lntzer
  integer tb, i

  filename = trim(filedir) // "fluka_binning"
  open (1, file = filename, form="UNFORMATTED")
!!$*----------- read 1st record -------------------------------------------
  read(1) runtitle, runhead, weight, ncase
!!$*----------- loop on binning detector data in the present file ---------
!!$* ---------------- read and write 2nd record --------------------

  read(1) mbin, namebin, typebin, distbin, &
       minbin(1), maxbin(1), nbin(1), dbin(1), &
       minbin(2), maxbin(2), nbin(2), dbin(2), &
       minbin(3), maxbin(3), nbin(3), dbin(3), &
       lntzer, b1bin, b2bin, tcbin

  tb = mod(typebin,10)

  if (tb .eq. 2) then
!!$* Region binning
  else if (tb .eq. 8) then
!!$* Region/Lattice/User binning
  else if (((tb .eq. 1) .or. (tb .eq. 7)) .and. nbin(2) .lt. 2) then
!!$* R-Z binning
  else if ((tb .eq. 1) .or. (tb .eq. 7)) then
!!$* R-Phi-Z binning
     read(1) (realfluka(i), i = 1, (nbin(1) * nbin(2) * nbin(3)))
     tfluka = tb
     ofluka(1) = sffluka * minbin(1)
     ofluka(2) = minbin(2)
     ofluka(3) = sffluka * minbin(3)
     dfluka(1) = sffluka * dbin(1)
     dfluka(2) = dbin(2)
     dfluka(3) = sffluka * dbin(3)
     nfluka(1) = nbin(1)
     nfluka(2) = nbin(2)
     nfluka(3) = nbin(3)
     do i = 1, 3
        ofluka(i) = ofluka(i) + trfluka(i)
     enddo
     do i = 1, nbin(1) * nbin(2) * nbin(3)
        datafluka(i) = ((npfluka * effluka) / (dtfluka * sffluka**3)) &
             * datafluka(i)
     enddo
  else if (tb .eq. 0) then
!!$* Cartesian binning
     read(1) (realfluka(i), i = 1, (nbin(1) * nbin(2) * nbin(3)))
     tfluka = tb
     ofluka(1) = sffluka * minbin(1)
     ofluka(2) = sffluka * minbin(2)
     ofluka(3) = sffluka * minbin(3)
     dfluka(1) = sffluka * dbin(1)
     dfluka(2) = sffluka * dbin(2)
     dfluka(3) = sffluka * dbin(3)
     nfluka(1) = nbin(1)
     nfluka(2) = nbin(2)
     nfluka(3) = nbin(3)
     do i = 1, 3
        ofluka(i) = ofluka(i) + trfluka(i)
     enddo
     do i = 1, nbin(1) * nbin(2) * nbin(3)
        datafluka(i) = ((npfluka * effluka) / (dtfluka * sffluka**3)) &
             * realfluka(i)
     enddo
  endif
  close(1)
end subroutine read_binning_bin

subroutine read_binning_txt(trfluka, sffluka, effluka, npfluka, &
     dtfluka, ofluka, dfluka, nfluka, tfluka, datafluka, filedir)

!!$------------------------------------------------------------------
!!$ Subroutine to read and store the FLUKA results from a binning 
!!$ text file
!!$
!!$ trfluka-- translation of the fluka grid origin               [in]
!!$ sffluka-- size scale factor (size_sph = sf*size_fluka)       [in]
!!$ effluka-- energy scale factor (energy_sph = ef*energy_fluka) [in]
!!$ npfluka-- number of particles                                [in]
!!$ dtfluka-- deposition time                                    [in]
!!$ ofluka-- the origin of the binning                          [out]
!!$ dfluka-- cell side of the binning                           [out]
!!$ nfluka-- number of cells for each index                     [out]
!!$ tfluka-- identifier for the fluka grid type                 [out]
!!$ datafluka-- data of the binning                             [out]
!!$ filedir -- Directory where the files are located             [in]

  implicit none
  include 'options.inc'

  double precision trfluka(3), sffluka, effluka, npfluka, dtfluka, &
       ofluka(3), dfluka(3), datafluka(max_fluka)
  integer nfluka(3), tfluka
  character (len = *) :: filedir

  character (len = 120) filename
  real minbin(3), maxbin(3), dbin(3)
  real realfluka(max_fluka)
  integer typebin, nbin(3)
  integer tb, i

  filename = trim(filedir) // "fluka_binning.FORM"
  open (1, file = filename, form="FORMATTED")

  read(1,*) ! title
  read(1,*) ! subtitle
  read(1,*) ! weight and ncase
  read(1,*) ! bin number
  read(1,*) ! name of the bin
  read(1,*) typebin

  read(1,*) minbin(1), maxbin(1), nbin(1), dbin(1)
  read(1,*) minbin(2), maxbin(2), nbin(2), dbin(2)
  read(1,*) minbin(3), maxbin(3), nbin(3), dbin(3)

  read(1,*)

  tb = mod(typebin,10)

  if (tb .eq. 2) then
!!$* Region binning
  else if (tb .eq. 8) then
!!$* Region/Lattice/User binning
  else if (((tb .eq. 1) .or. (tb .eq. 7)) .and. nbin(2) .lt. 2) then
!!$* R-Z binning
  else if ((tb .eq. 1) .or. (tb .eq. 7)) then
!!$* R-Phi-Z binning
     read(1,*) (realfluka(i), i = 1, (nbin(1) * nbin(2) * nbin(3)))
     tfluka = tb
     ofluka(1) = sffluka * minbin(1)
     ofluka(2) = minbin(2)
     ofluka(3) = sffluka * minbin(3)
     dfluka(1) = sffluka * dbin(1)
     dfluka(2) = dbin(2)
     dfluka(3) = sffluka * dbin(3)
     nfluka(1) = nbin(1)
     nfluka(2) = nbin(2)
     nfluka(3) = nbin(3)
     do i = 1, 3
        ofluka(i) = ofluka(i) + trfluka(i)
     enddo
     do i = 1, nbin(1) * nbin(2) * nbin(3)
        datafluka(i) = ((npfluka * effluka) / (dtfluka * sffluka**3)) &
             * datafluka(i)
     enddo
  else if (tb .eq. 0) then
!!$* Cartesian binning
     read(1,*) (realfluka(i), i = 1, (nbin(1) * nbin(2) * nbin(3)))
     tfluka = tb
     ofluka(1) = sffluka * minbin(1)
     ofluka(2) = sffluka * minbin(2)
     ofluka(3) = sffluka * minbin(3)
     dfluka(1) = sffluka * dbin(1)
     dfluka(2) = sffluka * dbin(2)
     dfluka(3) = sffluka * dbin(3)
     nfluka(1) = nbin(1)
     nfluka(2) = nbin(2)
     nfluka(3) = nbin(3)
     do i = 1, 3
        ofluka(i) = ofluka(i) + trfluka(i)
     enddo
     do i = 1, nbin(1) * nbin(2) * nbin(3)
        datafluka(i) = ((npfluka * effluka) / (dtfluka * sffluka**3)) &
             * realfluka(i)
     enddo
  endif
  close(1)
end subroutine read_binning_txt

subroutine fluka_heat(np, x, rho, dudt, &
     ofluka, dfluka, nfluka, tfluka, datafluka)


!!$------------------------------------------------------------------
!!$ Subroutine to calculate the power deposition on the particles 
!!$ from the FLUKA results
!!$	
!!$ np-- umber of particles                                      [in]
!!$ x-- particle position                                        [in]
!!$ rho-- density                                                [in]
!!$ dedt   : produced external heat, adding to energy Eq.       [out]
!!$ ofluka-- the origin of the binning                           [in]
!!$ dfluka-- cell side of the binning                            [in]
!!$ nfluka-- number of cells for each index                      [in]
!!$ tfluka-- identifier for the fluka grid type                  [in]
!!$ datafluka-- data of the binning                              [in]

  implicit none
  include 'options.inc'

  integer np
  double precision x(dim, maxn), rho(maxn), dudt(maxn), &
       ofluka(3), dfluka(3), datafluka(max_fluka)
  integer nfluka(3), tfluka

  integer ip, ic, d, imin, imax, jmin, jmax, kmin, kmax, ir(3)
  double precision r, phi, z, xr(3)
  double precision f1, f2, f3, f4, f5, f6, f7, f8, &
       v1, v2, v3, v4, v5, v6, v7, v8
  integer i, j, k


  if (dim .eq. 3) then
     do ip = 1, np
        if (tfluka .eq. 0) then
           if ((x(1,ip) >= ofluka(1)) .and. &
               (x(2,ip) >= ofluka(2)) .and. &
               (x(3,ip) >= ofluka(3)) .and. &
               (x(1,ip) <= (ofluka(1) + nfluka(1)*dfluka(1))) .and. &
               (x(2,ip) <= (ofluka(2) + nfluka(2)*dfluka(2))) .and. &
               (x(3,ip) <= (ofluka(3) + nfluka(3)*dfluka(3)))) then
              
              do d = 1, dim
                 ir(d) = int((x(d, ip) - ofluka(d)) / dfluka(d))
              enddo

              do d = 1, dim
                 xr(d) = (x(d, ip) - ofluka(d)) / dfluka(d) - ir(d)
              enddo

              if (xr(1) .ge. 0.5) then
                 imin = ir(1)
                 imax = ir(1) +1
                 xr(1) = 2.0 * xr(1) -2.0
              else
                 imin = ir(1) -1
                 imax = ir(1)
                 xr(1) = 2.0 * xr(1)
              endif

              if (xr(2) .ge. 0.5) then
                 jmin = ir(2)
                 jmax = ir(2) +1
                 xr(2) = 2.0 * xr(2) -2.0
              else
                 jmin = ir(2) -1
                 jmax = ir(2)
                 xr(2) = 2.0 * xr(2)
              endif

              if (xr(3) .ge. 0.5) then
                 kmin = ir(3)
                 kmax = ir(3) +1
                 xr(3) = 2.0 * xr(3) -2.0
              else
                 kmin = ir(3) -1
                 kmax = ir(3)
                 xr(3) = 2.0 * xr(3)
              endif

              if (imin .lt. 0) imin = 0
              if (jmin .lt. 0) jmin = 0
              if (kmin .lt. 0) kmin = 0
              if (imax .ge. nfluka(1)) imax = nfluka(1) -1
              if (jmax .ge. nfluka(2)) jmax = nfluka(2) -1
              if (kmax .ge. nfluka(3)) kmax = nfluka(3) -1

              f1 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v1 = datafluka(ic)

              f2 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v2 = datafluka(ic)

              f3 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v3 = datafluka(ic)

              f4 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v4 = datafluka(ic)

              f5 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v5 = datafluka(ic)

              f6 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v6 = datafluka(ic)

              f7 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v7 = datafluka(ic)

              f8 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v8 = datafluka(ic)

              dudt(ip) = dudt(ip) + (f1*v1 +f2*v2 +f3*v3 +f4*v4 &
                   +f5*v5 +f6*v6 +f7*v7 +f8*v8) / rho(ip)
           endif

        else if ((tfluka .eq. 1) .or. (tfluka .eq. 7)) then
           r = dsqrt(x(1,ip)**2 +x(2,ip)**2)
           phi = 180.0d0 / pi * datan2(x(2,ip),x(1,ip))
           z = x(3, ip)

           if ((r >= ofluka(1)) .and. &
                (phi >= ofluka(2)) .and. &
                (z >= ofluka(3)) .and. &
                (r <= (ofluka(1) + nfluka(1)*dfluka(1))) .and. &
                (phi <= (ofluka(2) + nfluka(2)*dfluka(2))) .and. &
                (z <= (ofluka(3) + nfluka(3)*dfluka(3)))) then

              ir(1) = int((r - ofluka(1)) / dfluka(1))
              ir(2) = int((phi - ofluka(2)) / dfluka(2))
              ir(3) = int((z - ofluka(3)) / dfluka(3))

              xr(1) = (r - ofluka(1)) / dfluka(1) - ir(1)
              xr(2) = (phi - ofluka(2)) / dfluka(2) - ir(2)
              xr(3) = (z - ofluka(3)) / dfluka(2) - ir(3)

              if (xr(1) .ge. 0.5) then
                 imin = ir(1)
                 imax = ir(1) +1
                 xr(1) = 2.0 * xr(1) -2.0
              else
                 imin = ir(1) -1
                 imax = ir(1)
                 xr(1) = 2.0 * xr(1)
              endif

              if (xr(2) .ge. 0.5) then
                 jmin = ir(2)
                 jmax = ir(2) +1
                 xr(2) = 2.0 * xr(2) -2.0
              else
                 jmin = ir(2) -1
                 jmax = ir(2)
                 xr(2) = 2.0 * xr(2)
              endif

              if (xr(3) .ge. 0.5) then
                 kmin = ir(3)
                 kmax = ir(3) +1
                 xr(3) = 2.0 * xr(3) -2.0
              else
                 kmin = ir(3) -1
                 kmax = ir(3)
                 xr(3) = 2.0 * xr(3)
              endif

              if (imin .lt. 0) imin = 0
              if (jmin .lt. 0) jmin = 0
              if (kmin .lt. 0) kmin = 0
              if (imax .ge. nfluka(1)) imax = nfluka(1) -1
              if (jmax .ge. nfluka(2)) jmax = nfluka(2) -1
              if (kmax .ge. nfluka(3)) kmax = nfluka(3) -1

              if (xr(1) .ge. 0.0) then
                 imin = ir(1)
                 imax = ir(1) +1
              else
                 imin = ir(1) -1
                 imax = ir(1)
              endif

              if (xr(2) .ge. 0.0) then
                 jmin = ir(2)
                 jmax = ir(2) +1
              else
                 jmin = ir(2) -1
                 jmax = ir(2)
              endif

              if (xr(3) .ge. 0.0) then
                 kmin = ir(3)
                 kmax = ir(3) +1
              else
                 kmin = ir(3) -1
                 kmax = ir(3)
              endif

              f1 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v1 = datafluka(ic)

              f2 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v2 = datafluka(ic)

              f3 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v3 = datafluka(ic)

              f4 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v4 = datafluka(ic)

              f5 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v5 = datafluka(ic)

              f6 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v6 = datafluka(ic)

              f7 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v7 = datafluka(ic)

              f8 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v8 = datafluka(ic)

              dudt(ip) = dudt(ip) + (f1*v1 +f2*v2 +f3*v3 +f4*v4 &
                   +f5*v5 +f6*v6 +f7*v7 +f8*v8) / rho(ip)

!!$	        i = int((r -ofluka(1))/dfluka(1))
!!$	        j = int((phi -ofluka(2))/dfluka(2))
!!$	        k = int((x(3,ip) -ofluka(3))/dfluka(3))
!!$	        ic = i +j*nfluka(1) +k*nfluka(1)*nfluka(2) +1
!!$	        
!!$	        dedt(ip) = datafluka(ic) / rho(ip)
           endif
        endif
     enddo

  else if (dim .eq. 2) then

     do ip = 1, np
        if (tfluka .eq. 0) then
           if ((x(1,ip) >= ofluka(1)) .and. &
                (x(2,ip) >= ofluka(2)) .and. &
                (x(1,ip) <= (ofluka(1) + nfluka(1)*dfluka(1))) .and. &
                (x(2,ip) <= (ofluka(2) + nfluka(2)*dfluka(2)))) then

              do d = 1, dim
                 ir(d) = int((x(d, ip) - ofluka(d)) / dfluka(d))
              enddo
              ir(3) = int((-ofluka(3)) / dfluka(3))

              do d = 1, dim
                 xr(d) = (x(d, ip) - ofluka(d)) / dfluka(d) - ir(d)
              enddo
              xr(3) = (-ofluka(3)) / dfluka(3) - ir(3)

              if (xr(1) .ge. 0.5) then
                 imin = ir(1)
                 imax = ir(1) +1
                 xr(1) = 2.0 * xr(1) -2.0
              else
                 imin = ir(1) -1
                 imax = ir(1)
                 xr(1) = 2.0 * xr(1)
              endif

              if (xr(2) .ge. 0.5) then
                 jmin = ir(2)
                 jmax = ir(2) +1
                 xr(2) = 2.0 * xr(2) -2.0
              else
                 jmin = ir(2) -1
                 jmax = ir(2)
                 xr(2) = 2.0 * xr(2)
              endif

              if (xr(3) .ge. 0.5) then
                 kmin = ir(3)
                 kmax = ir(3) +1
                 xr(3) = 2.0 * xr(3) -2.0
              else
                 kmin = ir(3) -1
                 kmax = ir(3)
                 xr(3) = 2.0 * xr(3)
              endif

              if (imin .lt. 0) imin = 0
              if (jmin .lt. 0) jmin = 0
              if (kmin .lt. 0) kmin = 0
              if (imax .ge. nfluka(1)) imax = nfluka(1) -1
              if (jmax .ge. nfluka(2)) jmax = nfluka(2) -1
              if (kmax .ge. nfluka(3)) kmax = nfluka(3) -1

              f1 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v1 = datafluka(ic)

              f2 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v2 = datafluka(ic)

              f3 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v3 = datafluka(ic)

              f4 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v4 = datafluka(ic)

              f5 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v5 = datafluka(ic)

              f6 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v6 = datafluka(ic)

              f7 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v7 = datafluka(ic)

              f8 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v8 = datafluka(ic)

              dudt(ip) = dudt(ip) + (f1*v1 +f2*v2 +f3*v3 +f4*v4 &
                   +f5*v5 +f6*v6 +f7*v7 +f8*v8) / rho(ip)
           endif

        else if ((tfluka .eq. 1) .or. (tfluka .eq. 7)) then
           r = dsqrt(x(1,ip)**2 +x(2,ip)**2)
           phi = 180.0d0 / pi * datan2(x(2,ip),x(1,ip))

           if ((r >= ofluka(1)) .and. &
                (phi >= ofluka(2)) .and. &
                (r <= (ofluka(1) + nfluka(1)*dfluka(1))) .and. &
                (phi <= (ofluka(2) + nfluka(2)*dfluka(2)))) then

              ir(1) = int((r - ofluka(1)) / dfluka(1))
              ir(2) = int((phi - ofluka(2)) / dfluka(2))
              ir(3) = int((-ofluka(3)) / dfluka(3))

              xr(1) = (r - ofluka(1)) / dfluka(1) - ir(1)
              xr(2) = (phi - ofluka(2)) / dfluka(2) - ir(2)
              xr(3) = (-ofluka(3)) / dfluka(2) - ir(3)

              if (xr(1) .ge. 0.5) then
                 imin = ir(1)
                 imax = ir(1) +1
                 xr(1) = 2.0 * xr(1) -2.0
              else
                 imin = ir(1) -1
                 imax = ir(1)
                 xr(1) = 2.0 * xr(1)
              endif

              if (xr(2) .ge. 0.5) then
                 jmin = ir(2)
                 jmax = ir(2) +1
                 xr(2) = 2.0 * xr(2) -2.0
              else
                 jmin = ir(2) -1
                 jmax = ir(2)
                 xr(2) = 2.0 * xr(2)
              endif

              if (xr(3) .ge. 0.5) then
                 kmin = ir(3)
                 kmax = ir(3) +1
                 xr(3) = 2.0 * xr(3) -2.0
              else
                 kmin = ir(3) -1
                 kmax = ir(3)
                 xr(3) = 2.0 * xr(3)
              endif

              if (imin .lt. 0) imin = 0
              if (jmin .lt. 0) jmin = 0
              if (kmin .lt. 0) kmin = 0
              if (imax .ge. nfluka(1)) imax = nfluka(1) -1
              if (jmax .ge. nfluka(2)) jmax = nfluka(2) -1
              if (kmax .ge. nfluka(3)) kmax = nfluka(3) -1

              f1 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v1 = datafluka(ic)

              f2 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmin*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v2 = datafluka(ic)

              f3 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imax +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v3 = datafluka(ic)

              f4 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * -1.0)
              ic = imin +jmax*nfluka(1) +kmin*nfluka(1)*nfluka(2) +1
              v4 = datafluka(ic)

              f5 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v5 = datafluka(ic)

              f6 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * -1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmin*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v6 = datafluka(ic)

              f7 = 0.125 * (1.0 + xr(1) * +1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imax +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v7 = datafluka(ic)

              f8 = 0.125 * (1.0 + xr(1) * -1.0) &
                   * (1.0 + xr(2) * +1.0) &
                   * (1.0 + xr(3) * +1.0)
              ic = imin +jmax*nfluka(1) +kmax*nfluka(1)*nfluka(2) +1
              v8 = datafluka(ic)

              dudt(ip) = dudt(ip) + (f1*v1 +f2*v2 +f3*v3 +f4*v4 &
                   +f5*v5 +f6*v6 +f7*v7 +f8*v8) / rho(ip)

           endif
        endif
     enddo
  endif

end subroutine fluka_heat

