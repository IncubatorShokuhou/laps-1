!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis

subroutine model_pblhgt(theta, thsfc, psig, zsig, topo, nx, ny, nz, pblhgt)

   !  subroutine to estimate height agl in meters of pbl from native
   !  coordinate model data.  adapted from the laps routine for
   !  terrain-following model coordinates.

   implicit none

   integer, intent(in)   :: nx, ny, nz
   real, intent(in)      :: theta(nx, ny, nz)
   real, intent(in)      :: thsfc(nx, ny)
   real, intent(in)      :: psig(nx, ny, nz)
   real, intent(in)      :: zsig(nx, ny, nz)
   real, intent(in)      :: topo(nx, ny)
   real, intent(out)     :: pblhgt(nx, ny)

   integer  :: i, j, k, ktop
   real     :: thresh_k, topwgt, botwgt
   logical  :: found_pbl_top

   print *, 'generating pbl height using theta and surface temp.'
   loop_j: do j = 1, ny
      loop_i: do i = 1, nx

         ! compute threshold value that theta needs to exceed
         ! to be above pbl.  we use surface theta plus an
         ! additional 3k for slop.

         thresh_k = thsfc(i, j) + 3.0

         ! now begin at the bottom and work our way up until
         ! we find the first level with a theta exceeding the
         ! threshold

         found_pbl_top = .false.
         loop_k: do k = 1, nz

            if (theta(i, j, k) .ge. thresh_k) then
               ktop = k
               found_pbl_top = .true.
               exit loop_k
            end if
         end do loop_k

         ! if we did not find a good pbl, set pbl to first level
         ! and print out some diagnostics
         if (.not. found_pbl_top) then
            print *, 'pbl height not found at i/j = ', i, j
            print *, 'surface theta = ', thsfc(i, j)
            print *, 'theta in the column:'
            print *, 'pressure height  theta'
            print *, '-------- ------- --------'
            diag_loop: do k = 1, nz
               print '(f8.0,f8.0,f8.2)', psig(i, j, k), zsig(i, j, k), theta(i, j, k)
            end do diag_loop
            ktop = 1
            pblhgt(i, j) = zsig(i, j, 1) - topo(i, j)

         else

            ! we found the top k-level bounding the pbl so interpolate
            ! to the actual level

            if (ktop .eq. 1) then
               pblhgt(i, j) = zsig(i, j, 1) - topo(i, j)
            else
               ! interpolate to get height at thresh_k
               botwgt = ((theta(i, j, ktop) - thresh_k)/ &
                         (theta(i, j, ktop) - theta(i, j, ktop - 1)))
               topwgt = 1.0 - botwgt
               pblhgt(i, j) = botwgt*zsig(i, j, ktop - 1) + &
                              topwgt*zsig(i, j, ktop) - topo(i, j)
            end if
         end if
      end do loop_i
   end do loop_j
   return
end subroutine model_pblhgt
