subroutine pressure_action_new(np, x, v, h, mass, rho, p, &
     dvdt, dudt, niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the shear deformation tensor 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ p    :  pressure                                             [in]
!!$ dvdt :  particle acceleration                               [out]
!!$ dudt :  particle energy increase                            [out]
  
  
  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn), p(maxn)
  double precision dvdt(dim, maxn), dudt(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i,j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, dv(dim)
  double precision mh, hvi, hvj
  
  
!!$ Calculate pressure force and the internal energy change

  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .le. np)) then

        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)

        do d = 1, dim
           dv(d) = v(d,j) -v(d,i)
        enddo
        
        hvi = -(p(j)-p(i)) / (rho(i) * rho(j))
        hvj = -(p(i)-p(j)) / (rho(i) * rho(j))
        
        do d = 1, dim
           dvdt(d,i) = dvdt(d,i) + mass(j)*hvi*dwdx(d)
           dvdt(d,j) = dvdt(d,j) - mass(i)*hvj*dwdx(d)
           
           dudt(i) = dudt(i) + 0.5d0*dv(d)*mass(j)*hvi*dwdx(d)
           dudt(j) = dudt(j) + 0.5d0*dv(d)*mass(i)*hvj*dwdx(d)
        enddo
     endif
     
  enddo

end Subroutine pressure_action_new


subroutine pressure_action(np, x, v, h, mass, rho, p, &
     dvdt, dudt, niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the shear deformation tensor 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ p    :  pressure                                             [in]
!!$ dvdt :  particle acceleration                               [out]
!!$ dudt :  particle energy increase                            [out]
  
  
  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn), p(maxn)
  double precision dvdt(dim, maxn), dudt(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i,j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, dv(dim)
  double precision mh, hv
  
  
!!$ Calculate pressure force and the internal energy change

  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .le. np)) then

        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)

        do d = 1, dim
           dv(d) = v(d,j) -v(d,i)
        enddo
        
        hv = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)
        do d = 1, dim
           dvdt(d,i) = dvdt(d,i) + mass(j)*hv*dwdx(d)
           dvdt(d,j) = dvdt(d,j) - mass(i)*hv*dwdx(d)
           
           dudt(i) = dudt(i) + 0.5d0*dv(d)*mass(j)*hv*dwdx(d)
           dudt(j) = dudt(j) + 0.5d0*dv(d)*mass(i)*hv*dwdx(d)
        enddo
     endif
     
  enddo

end Subroutine pressure_action


subroutine shear_deformation(np, x, v, h, mass, rho, &
     epsxx, epsyy, epszz, epsxy, epsyz, epszx, &
     niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the shear deformation tensor 
!!$
!!$ np   : number of particles                                   [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ eps  : shear deformations                                   [out]

  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn)
  double precision epsxx(maxn), epsyy(maxn), epszz(maxn), &
       epsxy(maxn), epsyz(maxn), epszx(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i,j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, dv(dim)
  double precision mrho, mh
  double precision hxx, hyy, hzz, hyz, hzx, hxy
  
  do i = 1, np
     if (dim .eq. 1) then
        epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
     elseif (dim .eq. 2) then
        epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
        epsyy(i) = epsyy(i) +mass(j)*hyy/rho(j)
        epsxy(i) = epsxy(i) +mass(j)*hxy/rho(j)
     elseif (dim .eq. 3) then
        epsxx(i) = epsxx(i) + mass(j)*hxx/rho(j)
        epsyy(i) = epsyy(i) + mass(j)*hyy/rho(j)
        epszz(i) = epszz(i) + mass(j)*hzz/rho(j)
        epsxy(i) = epsxy(i) + mass(j)*hxy/rho(j)
        epsyz(i) = epsyz(i) + mass(j)*hyz/rho(j)
        epszx(i) = epszx(i) + mass(j)*hzx/rho(j)
     endif
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        do d = 1, dim
           dv(d) = v(d,j) -v(d,i)
        enddo
        
        if (dim .eq. 1) then
           hxx = 2.0d00 * dv(1)*dwdx(1)
           
           hxx = 2.0d0/3.0d0 * hxx
           
           epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
           
           epsxx(j) = epsxx(j) +mass(i)*hxx/rho(i)
           
        elseif (dim .eq. 2) then
           hxx = 2.0d0*dv(1)*dwdx(1) -dv(2)*dwdx(2)
           hyy = 2.0d0*dv(2)*dwdx(2) -dv(1)*dwdx(1)
           hxy = dv(1)*dwdx(2) +dv(2)*dwdx(1)
           
           hxx = 2.0d0/3.0d0 * hxx
           hyy = 2.0d0/3.0d0 * hyy
           
           epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
           epsyy(i) = epsyy(i) +mass(j)*hyy/rho(j)
           epsxy(i) = epsxy(i) +mass(j)*hxy/rho(j)
           
           epsxx(j) = epsxx(j) +mass(i)*hxx/rho(i)
           epsyy(j) = epsyy(j) +mass(i)*hyy/rho(i)
           epsxy(j) = epsxy(j) +mass(i)*hxy/rho(i)
        
        elseif (dim .eq. 3) then
           hxx = 2.0d0*dv(1)*dwdx(1) -dv(2)*dwdx(2) &
                -dv(3)*dwdx(3)
           hyy = 2.0d0*dv(2)*dwdx(2) -dv(1)*dwdx(1) &
                -dv(3)*dwdx(3)
           hzz = 2.0d0*dv(3)*dwdx(3) -dv(1)*dwdx(1) &
                -dv(2)*dwdx(2)
           hxy = dv(1)*dwdx(2) +dv(2)*dwdx(1)
           hyz = dv(2)*dwdx(3) +dv(3)*dwdx(2)
           hzx = dv(1)*dwdx(3) +dv(3)*dwdx(1)
        
           hxx = 2.0d0/3.0d0 * hxx
           hyy = 2.0d0/3.0d0 * hyy
           hzz = 2.0d0/3.0d0 * hzz
           
           epsxx(i) = epsxx(i) + mass(j)*hxx/rho(j)
           epsyy(i) = epsyy(i) + mass(j)*hyy/rho(j)
           epszz(i) = epszz(i) + mass(j)*hzz/rho(j)
           epsxy(i) = epsxy(i) + mass(j)*hxy/rho(j)
           epsyz(i) = epsyz(i) + mass(j)*hyz/rho(j)
           epszx(i) = epszx(i) + mass(j)*hzx/rho(j)
           
           epsxx(j) = epsxx(j) + mass(i)*hxx/rho(i)
           epsyy(j) = epsyy(j) + mass(i)*hyy/rho(i)
           epszz(j) = epszz(j) + mass(i)*hzz/rho(i)
           epsxy(j) = epsxy(j) + mass(i)*hxy/rho(i)
           epsyz(j) = epsyz(j) + mass(i)*hyz/rho(i)
           epszx(j) = epszx(j) + mass(i)*hzx/rho(i)
           
        endif
     endif
     
  enddo
  
end subroutine shear_deformation


subroutine viscous_action(np, x, v, h, mass, rho, eta, &
     dvdt, dudt, niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the shear deformation tensor 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ eta  :  viscosity                                            [in]
!!$ dvdt :  particle acceleration                               [out]
!!$ dudt :  particle energy increase                            [out]
  
  
  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn), eta(maxn)
  double precision epsxx(maxn), epsyy(maxn), epszz(maxn), &
       epsxy(maxn), epsyz(maxn), epszx(maxn)
  double precision dvdt(dim, maxn), dudt(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i,j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac
  double precision mh, hv, tdsdt, dv(dim)
  double precision hxx, hyy, hzz, hyz, hzx, hxy
  
  
!!$ Calculate acceleration due to viscous effects
  
  do i = 1, np
     if (dim .eq. 1) then
        epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
     elseif (dim .eq. 2) then
        epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
        epsyy(i) = epsyy(i) +mass(j)*hyy/rho(j)
        epsxy(i) = epsxy(i) +mass(j)*hxy/rho(j)
     elseif (dim .eq. 3) then
        epsxx(i) = epsxx(i) + mass(j)*hxx/rho(j)
        epsyy(i) = epsyy(i) + mass(j)*hyy/rho(j)
        epszz(i) = epszz(i) + mass(j)*hzz/rho(j)
        epsxy(i) = epsxy(i) + mass(j)*hxy/rho(j)
        epsyz(i) = epsyz(i) + mass(j)*hyz/rho(j)
        epszx(i) = epszx(i) + mass(j)*hzx/rho(j)
     endif
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        do d = 1, dim
           dv(d) = v(d,j) -v(d,i)
        enddo
        
        if (dim .eq. 1) then
           hxx = 2.0d00 * dv(1)*dwdx(1)
           
           hxx = 2.0d0/3.0d0 * hxx
           
           epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
           
           epsxx(j) = epsxx(j) +mass(i)*hxx/rho(i)
           
        elseif (dim .eq. 2) then
           hxx = 2.0d0*dv(1)*dwdx(1) -dv(2)*dwdx(2)
           hyy = 2.0d0*dv(2)*dwdx(2) -dv(1)*dwdx(1)
           hxy = dv(1)*dwdx(2) +dv(2)*dwdx(1)
           
           hxx = 2.0d0/3.0d0 * hxx
           hyy = 2.0d0/3.0d0 * hyy
           
           epsxx(i) = epsxx(i) +mass(j)*hxx/rho(j)
           epsyy(i) = epsyy(i) +mass(j)*hyy/rho(j)
           epsxy(i) = epsxy(i) +mass(j)*hxy/rho(j)
           
           epsxx(j) = epsxx(j) +mass(i)*hxx/rho(i)
           epsyy(j) = epsyy(j) +mass(i)*hyy/rho(i)
           epsxy(j) = epsxy(j) +mass(i)*hxy/rho(i)
        
        elseif (dim .eq. 3) then
           hxx = 2.0d0*dv(1)*dwdx(1) -dv(2)*dwdx(2) &
                -dv(3)*dwdx(3)
           hyy = 2.0d0*dv(2)*dwdx(2) -dv(1)*dwdx(1) &
                -dv(3)*dwdx(3)
           hzz = 2.0d0*dv(3)*dwdx(3) -dv(1)*dwdx(1) &
                -dv(2)*dwdx(2)
           hxy = dv(1)*dwdx(2) +dv(2)*dwdx(1)
           hyz = dv(2)*dwdx(3) +dv(3)*dwdx(2)
           hzx = dv(1)*dwdx(3) +dv(3)*dwdx(1)
        
           hxx = 2.0d0/3.0d0 * hxx
           hyy = 2.0d0/3.0d0 * hyy
           hzz = 2.0d0/3.0d0 * hzz
           
           epsxx(i) = epsxx(i) + mass(j)*hxx/rho(j)
           epsyy(i) = epsyy(i) + mass(j)*hyy/rho(j)
           epszz(i) = epszz(i) + mass(j)*hzz/rho(j)
           epsxy(i) = epsxy(i) + mass(j)*hxy/rho(j)
           epsyz(i) = epsyz(i) + mass(j)*hyz/rho(j)
           epszx(i) = epszx(i) + mass(j)*hzx/rho(j)
           
           epsxx(j) = epsxx(j) + mass(i)*hxx/rho(i)
           epsyy(j) = epsyy(j) + mass(i)*hyy/rho(i)
           epszz(j) = epszz(j) + mass(i)*hzz/rho(i)
           epsxy(j) = epsxy(j) + mass(i)*hxy/rho(i)
           epsyz(j) = epsyz(j) + mass(i)*hyz/rho(i)
           epszx(j) = epszx(j) + mass(i)*hzx/rho(i)
           
        endif
     endif
     
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        if (dim .eq. 1) then
           hv = (eta(i)*epsxx(i)/rho(i)**2 &
                +  eta(j)*epsxx(j)/rho(j)**2)*dwdx(1)
           
           dvdt(1,i) = dvdt(1,i) + mass(j)*hv
           dvdt(1,j) = dvdt(1,j) - mass(i)*hv
        
        elseif (dim .eq. 2) then
           hv = (eta(i)*epsxx(i)/rho(i)**2 &
              +  eta(j)*epsxx(j)/rho(j)**2)*dwdx(1) &
              + (eta(i)*epsxy(i)/rho(i)**2 &
              +  eta(j)*epsxy(j)/rho(j)**2)*dwdx(2)
           
           dvdt(1,i) = dvdt(1,i) + mass(j)*hv
           dvdt(1,j) = dvdt(1,j) - mass(i)*hv
           
           hv = (eta(i)*epsxy(i)/rho(i)**2 &
              +  eta(j)*epsxy(j)/rho(j)**2)*dwdx(1) &
              + (eta(i)*epsyy(i)/rho(i)**2 &
              +  eta(j)*epsyy(j)/rho(j)**2)*dwdx(2)
           
           dvdt(2,i) = dvdt(2,i) + mass(j)*hv
           dvdt(2,j) = dvdt(2,j) - mass(i)*hv
           
        elseif (dim .eq. 3) then
           hv = (eta(i)*epsxx(i)/rho(i)**2 &
              +  eta(j)*epsxx(j)/rho(j)**2)*dwdx(1) &
              + (eta(i)*epsxy(i)/rho(i)**2 &
              +  eta(j)*epsxy(j)/rho(j)**2)*dwdx(2) &
              + (eta(i)*epszx(i)/rho(i)**2 &
              +  eta(j)*epszx(j)/rho(j)**2)*dwdx(3)
           
           dvdt(1,i) = dvdt(1,i) + mass(j)*hv
           dvdt(1,j) = dvdt(1,j) - mass(i)*hv
           
           hv = (eta(i)*epsxy(i)/rho(i)**2 &
              +  eta(j)*epsxy(j)/rho(j)**2)*dwdx(1) &
              + (eta(i)*epsyy(i)/rho(i)**2 &
              +  eta(j)*epsyy(j)/rho(j)**2)*dwdx(2) &
              + (eta(i)*epsyz(i)/rho(i)**2 &
              +  eta(j)*epsyz(j)/rho(j)**2)*dwdx(3)
           
           dvdt(2,i) = dvdt(2,i) + mass(j)*hv
           dvdt(2,j) = dvdt(2,j) - mass(i)*hv
           
           hv = (eta(i)*epszx(i)/rho(i)**2 &
              +  eta(j)*epszx(j)/rho(j)**2)*dwdx(1) &
              + (eta(i)*epsyz(i)/rho(i)**2 &
              +  eta(j)*epsyz(j)/rho(j)**2)*dwdx(2) &
              + (eta(i)*epszz(i)/rho(i)**2 &
              +  eta(j)*epszz(j)/rho(j)**2)*dwdx(3)
           
           dvdt(3,i) = dvdt(3,i) + mass(j)*hv
           dvdt(3,j) = dvdt(3,j) - mass(i)*hv
        endif
     endif
  enddo

!!$ Change of specific internal energy
  do i = 1, np
     if (dim .eq. 1) then
        tdsdt = epsxx(i)*epsxx(i)
        dudt(i) = dudt(i) + 0.5e0 * eta(i)/rho(i) * tdsdt
     elseif (dim .eq. 2) then
        tdsdt = epsxx(i)*epsxx(i) + epsyy(i)*epsyy(i) &
              + 2.0d00*epsxy(i)*epsxy(i)
        dudt(i) = dudt(i) + 0.5e0 * eta(i)/rho(i) * tdsdt
     elseif (dim .eq. 3) then
        tdsdt = epsxx(i)*epsxx(i) + epsyy(i)*epsyy(i) &
              + epszz(i)*epszz(i) &
              + 2.0d0*epsyz(i)*epsyz(i) + 2.0d0*epszx(i)*epszx(i) &
              + 2.0d0*epsxy(i)*epsxy(i)
        dudt(i) = dudt(i) + 0.5e0 * eta(i)/rho(i) * tdsdt
     endif
  enddo
  
end subroutine viscous_action


subroutine velocity_divergence(np, x, v, h, mass, rho, &
     divv, &
     niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the shear deformation tensor 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ divv :  divergence of the velocity field                    [out]

  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn), divv(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i, j, k, d
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac
  double precision mh
  

!!$ Calculate the velocity divervence
  
  do i = 1, np
     divv(i) =0.0d0
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        do d = 1, dim
           divv(i) = divv(i) &
                   + mass(j)/rho(i)*(v(d,j) -v(d,i))*dwdx(d)
           divv(j) = divv(j) &
                   + mass(i)/rho(j)*(v(d,j) -v(d,i))*dwdx(d)
        enddo
        
     endif
     
  enddo
  
end subroutine velocity_divergence



subroutine shock_correction(np, x, v, h, mass, rho, c, &
     dvdt, dudt, niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the artificial viscous correction 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ p    :  pressure                                             [in]
!!$ c    :  sound speed                                          [in]
!!$ dvdt :  particle acceleration                               [out]
!!$ dudt :  particle energy increase                            [out]
  
  
  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn), c(maxn)
  double precision dvdt(dim, maxn), dudt(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i, j, k, d
  double precision mh, mrho, mc, phi, vr, pi_ij
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, dv(dim)
  double precision rdwdx, psi_i, psi_j, psi_ij, psi
  double precision alpha, beta, eta
  
  
!!$ Parameter for the artificial viscosity:
  parameter(alpha = 1.0d0) !!$ Shear viscosity
  parameter(beta = 1.0d0) !!$ Bulk viscosity
  parameter(eta = 0.1d0) !!$ Parameter to avoid singularities
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        vr = 0.0d0
        do d = 1, dim
           dv(d) = v(d,i) -v(d,j)
           vr = vr + dv(d) * dxiac(d)
        enddo
        
!!$ Artificial viscous force
        if (vr < 0.0d0) then
           phi = (mh * vr) / (dr2iac + mh*mh*eta*eta)
           
           mc = 0.5d0 * (c(i) + c(j))
           mrho = 0.5d0 * (rho(i) + rho(j))
           pi_ij = (beta*phi - alpha*mc) * phi/mrho
           
           do d = 1, dim
              dvdt(d, i) = dvdt(d, i) - mass(j) * pi_ij * dwdx(d)
              dvdt(d, j) = dvdt(d, j) + mass(i) * pi_ij * dwdx(d)
              
              dudt(i) = dudt(i) + 0.5d0 * mass(j) * dv(d) * pi_ij * dwdx(d)
              dudt(j) = dudt(j) + 0.5d0 * mass(i) * dv(d) * pi_ij * dwdx(d)
           enddo
        endif
        
     endif
  enddo
  
end Subroutine shock_correction


subroutine heat_correction(np, x, v, h, mass, rho, u, c, &
     dudt, niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the artificial heat correction 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ v    :  particle velocity                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ u    :  energy                                               [in]
!!$ c    :  sound speed                                          [in]
!!$ dudt :  particle energy increase                            [out]
  
  
  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), v(dim,maxn), h(maxn), &
       mass(maxn), rho(maxn), u(maxn), c(maxn)
  double precision dudt(maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i, j, k, d
  double precision divv(maxn)
  double precision mh, mrho
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, dv(dim)
  double precision rdwdx, psi_i, psi_j, psi_ij, psi
  double precision g1, g2, g3
  
  
!!$ Parameter for the artificial heat conductivity:
  parameter(g1 = 0.1d0)
  parameter(g2 = 1.0d0)
  parameter(g3 = 0.1d0)
  
!!$ Calculate the velocity divervence
  
  do i = 1, np
     divv(i) =0.0d0
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)
     
     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        do d = 1, dim
           divv(i) = divv(i) &
                   + mass(j)/rho(i)*(v(d,j) -v(d,i))*dwdx(d)
           divv(j) = divv(j) &
                   + mass(i)/rho(j)*(v(d,j) -v(d,i))*dwdx(d)
        enddo
        
     endif
     
  enddo
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .le. np)) then
        
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = 0.5*(h(i) +h(j))
        driac = sqrt(dr2iac)
        call kernel(driac, dxiac, mh, w, dwdx, skf)
        
        mrho = 0.5d0 * (rho(i) + rho(j))
        
        rdwdx = 0.0d0
        do d = 1, dim
           rdwdx  = rdwdx + dxiac(d)*dwdx(d)
        enddo
        
        psi_i = g1*h(i)*c(i) + g2*h(i)*h(i)*(abs(divv(i)) -divv(i))
        psi_j = g1*h(j)*c(j) + g2*h(j)*h(j)*(abs(divv(j)) -divv(j)) 
        psi_ij = 0.5d0 * (psi_i + psi_j)
        
        psi = (2.0d0*psi_ij) / (mrho*(dr2iac + (mh*g3)**2)) * rdwdx
        
        dudt(i) = dudt(i) + mass(j) * psi * (u(i) - u(j))
        dudt(j) = dudt(j) + mass(i) * psi * (u(j) - u(i))
        
     endif
  enddo
  
end Subroutine heat_correction


subroutine tensile_correction(np, x, h, mass, rho, p, &
     dvdt, niac, pair_i, pair_j, skf)
  
!!$----------------------------------------------------------------------
!!$ Subroutine to calculate the tensile instability correction 
!!$
!!$ np   :  number of particles                                  [in]
!!$ x    :  particle position                                    [in]
!!$ h    :  smoothing length                                     [in]
!!$ mass :  particle mass                                        [in]
!!$ rho  :  density                                              [in]
!!$ p    :  pressure                                             [in]
!!$ dvdt :  particle acceleration                               [out]
  
  
  implicit none
  include 'options.inc'
  
  integer np
  double precision x(dim,maxn), h(maxn), mass(maxn), rho(maxn), p(maxn)
  double precision dvdt(dim, maxn)
  integer niac, pair_i(max_interaction), pair_j(max_interaction), skf
  
  integer i, j, k, d
  double precision mh, fij, Rij, epsti, w0
  double precision w, dwdx(dim), dxiac(dim), driac, dr2iac, dv(dim)
  integer nti
  
  
!!$ Parameter for the tensile instability:
  parameter(nti = 4)
  parameter(epsti = 0.2)
  
  do k = 1, niac
     i = pair_i(k)
     j = pair_j(k)

     if ((i .le. np) .and. (j .le. np)) then
        
        Rij = 0.01*(dabs(p(j)) / (rho(j)*rho(j)) + dabs(p(i)) / (rho(i)*rho(i)))
        if (p(i) .le. 0.0d0) then
           Rij = Rij + epsti * dabs(p(i)) / (rho(i)*rho(i))
        endif
        
        if (p(j) .le. 0.0d0) then
           Rij = Rij + epsti * dabs(p(j)) / (rho(j)*rho(j))
        endif
        
        if (Rij .gt. 0.0d0) then
           
           dr2iac = 0.0
           do d = 1, dim
              dxiac(d) = x(d,i) - x(d,j)
              dr2iac = dr2iac + dxiac(d)*dxiac(d)
           enddo
           mh = 0.5*(h(i) +h(j))
           driac = sqrt(dr2iac)
           
           call kernel(mh / hdr, dxiac, mh, w0, dwdx, skf)
           
           call kernel(driac, dxiac, mh, w, dwdx, skf)
           
           fij = - (w / w0)**nti * Rij
           
           do d = 1, dim
              dvdt(d, i) = dvdt(d, i) + mass(j) * fij * dwdx(d)
              dvdt(d, j) = dvdt(d, j) - mass(i) * fij * dwdx(d)
           enddo
        endif
        
     endif
  enddo
  
end subroutine tensile_correction


