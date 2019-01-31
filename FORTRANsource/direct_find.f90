subroutine direct_find(nt, x, h, niac, pair_i, pair_j)

!!$-------------------------------------------------------------------
!!$ Subroutine to calculate the interaction pairs between the particles
!!$ the pairs are determined by directly comparing the particle distance 
!!$ with the corresponding smoothing length.
!!$
!!$   nt-- total number of particles                               [in]
!!$   h-- smoothing length                                         [in]
!!$   x-- coordinates of all particles                             [in]
!!$   niac-- number of interaction pairs                          [out]
!!$   pair_i-- list of first partner of interaction pair          [out]
!!$   pair_j-- list of second partner of interaction pair         [out]

  implicit none
  include 'options.inc'

  integer nt, niac, pair_i(max_interaction), pair_j(max_interaction)
  double precision x(dim,maxn), h(maxn)

  integer scale_k, i, j, d
  double precision dxiac(dim), dr2iac, mh
  
  scale_k = 3
  
  niac = 0
  do i = 1, nt -1 
     do j = i +1, nt 
        dr2iac = 0.0
        do d = 1, dim
           dxiac(d) = x(d,i) - x(d,j)
           dr2iac = dr2iac + dxiac(d)*dxiac(d)
        enddo
        mh = (h(i) +h(j)) / 2.
        if (sqrt(dr2iac) < scale_k * mh) then
           if (niac < max_interaction) then    

!!$ Neighboring pair list, and total interaction number and
!!$ the interaction number for each particle 

              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
           else
              print *, &
                   ' >>> ERROR <<< : Too many interactions' 
              stop
           endif
        endif
     enddo
  enddo

end subroutine direct_find
